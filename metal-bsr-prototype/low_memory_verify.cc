#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"
#include "chromabase.h"
#include "production_bridge.hpp"
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <sys/resource.h>

static_assert(BASE_PRECISION == 64, "Native outer operator must be FP64");
namespace {
  using namespace QDP;
  using namespace Chroma;
  using Z = chroma_bridge::Complex;
  using Gauge = multi1d<LatticeColorMatrixD>;
  using GaugeF = multi1d<LatticeColorMatrixF>;
  constexpr uint32_t seed = 20261107;
  using Clock = std::chrono::steady_clock;
  double elapsed(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now() - t).count();
  }
  void save(const std::string &path, const LatticeFermionD &x,
            const std::vector<int> &local) {
    if (std::filesystem::exists(path))
      throw std::runtime_error("solution already exists");
    std::ofstream out(path, std::ios::binary);
    std::array<Z, 12> row;
    for (int site : local) {
      for (int c = 0; c < 3; ++c)
        for (int s = 0; s < 4; ++s) {
          auto z = x.elem(site).elem(s).elem(c);
          row[c * 4 + s] = {z.real(), z.imag()};
        }
      chroma_bridge::write(out, row.data(), 12);
    }
    out.close();
    if (!out)
      throw std::runtime_error("solution write failed");
  }
  void load(const std::string &path, LatticeFermionD &x,
            const std::vector<int> &local) {
    std::ifstream in(path, std::ios::binary);
    std::array<Z, 12> row;
    for (int site : local) {
      chroma_bridge::read(in, row.data(), 12);
      for (int c = 0; c < 3; ++c)
        for (int s = 0; s < 4; ++s) {
          auto &z = x.elem(site).elem(s).elem(c);
          auto a = row[c * 4 + s];
          z.real() = a.re;
          z.imag() = a.im;
        }
    }
    if (in.peek() != EOF)
      throw std::runtime_error("wrong spinor size");
  }
  void source(LatticeFermionD &b, const std::vector<int> &local, uint32_t tag) {
    double scale = 1 / std::sqrt(double(local.size()) * 12);
    for (size_t site = 0; site < local.size(); ++site)
      for (int c = 0; c < 3; ++c)
        for (int s = 0; s < 4; ++s) {
          auto z = chroma_bridge::source(site * 12 + c * 4 + s, tag);
          auto &v = b.elem(local[site]).elem(s).elem(c);
          v.real() = scale * z.re;
          v.imag() = scale * z.im;
        }
  }
} // namespace

int main(int argc, char **argv) {
  if (argc != 7) {
    std::cerr << "Usage: low_memory_verify ACTION_XML BRIDGE_PREFIX "
                 "OUTPUT_PREFIX METHOD RHS_COUNT CSV\n";
    return 2;
  }
  std::string xml = argv[1], prefix = argv[2], output = argv[3],
              method = argv[4], csvpath = argv[6];
  int count = std::stoi(argv[5]);
  int qargc = 1;
  QDP_initialize(&qargc, &argv);
  int status = 0;
  try {
    if (method != "metal-mg" && method != "metal-bicgstab")
      throw std::runtime_error("method");
    XMLReader metadata(xml);
    multi1d<int> dims, boundary;
    read(metadata, "/DBMetaData/lattSize", dims);
    read(metadata,
         "/DBMetaData/Params/Propagator/FermionAction/FermState/FermionBC/"
         "boundary",
         boundary);
    CloverFermActParams params(metadata,
                               "/DBMetaData/Params/Propagator/FermionAction");
    Layout::setLattSize(dims);
    Layout::create();
    std::ifstream input(prefix + ".bridge", std::ios::binary);
    auto h = production::read_header(input);
    for (int i = 0; i < 4; ++i)
      if (dims[i] != int(h.shape[i]))
        throw std::runtime_error("shape mismatch");
    if (toDouble(params.Mass) != h.mass ||
        toDouble(params.anisoParam.xi_0) != h.xi0 ||
        toDouble(params.anisoParam.nu) != h.nu ||
        toDouble(params.clovCoeffR) != h.clover_r ||
        toDouble(params.clovCoeffT) != h.clover_t)
      throw std::runtime_error("action mismatch");
    std::vector<int> local(h.volume());
    multi1d<int> coord(4);
    for (uint32_t site = 0; site < h.volume(); ++site) {
      uint32_t q = site;
      for (int mu = 0; mu < 4; ++mu) {
        coord[mu] = q % h.shape[mu];
        q /= h.shape[mu];
      }
      local[site] = Layout::linearSiteIndex(coord);
    }
    Gauge links(4);
    std::array<Z, 36> row;
    for (uint32_t site = 0; site < h.volume(); ++site) {
      chroma_bridge::read(input, row.data(), 36);
      uint32_t stride = 1;
      for (int mu = 0; mu < 4; ++mu) {
        int sign =
            (site / stride % h.shape[mu] == h.shape[mu] - 1) ? boundary[mu] : 1;
        stride *= h.shape[mu];
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) {
            auto z = row[mu * 9 + i * 3 + j];
            auto &v = links[mu].elem(local[site]).elem().elem(i, j);
            v.real() = sign * z.re;
            v.imag() = sign * z.im;
          }
      }
    }
    std::unique_ptr<UnprecCloverLinOp> op;
    {
      Handle<FermState<LatticeFermionD, Gauge, Gauge>> state(
          new PeriodicFermState<LatticeFermionD, Gauge, Gauge>(links));
      op.reset(new UnprecCloverLinOp(state, params));
    }
    links.resize(0);
    LatticeFermionD b, x, ax, r;
    source(x, local, chroma_bridge::probe_rhs);
    x *= std::sqrt(double(h.volume()) * 12);
    (*op)(ax, x, PLUS);
    double error = 0, norm = 0;
    std::array<Z, 12> values;
    for (uint32_t site = 0; site < h.volume(); ++site) {
      chroma_bridge::read(input, values.data(), 12);
      for (int c = 0; c < 3; ++c)
        for (int s = 0; s < 4; ++s) {
          auto v = ax.elem(local[site]).elem(s).elem(c);
          auto z = values[c * 4 + s];
          error += (v.real() - z.re) * (v.real() - z.re) +
                   (v.imag() - z.im) * (v.imag() - z.im);
          norm += z.re * z.re + z.im * z.im;
        }
    }
    double match = std::sqrt(error / norm);
    std::cout << std::setprecision(17)
              << "NATIVE original_UnprecCloverLinOp probe_relative=" << match
              << std::endl;
    if (!std::isfinite(match) || match > 1e-13)
      throw std::runtime_error("original native probe mismatch");
    std::ofstream csv(csvpath);
    csv << std::setprecision(17) << "rhs,method,full_relative\n";
    for (int rhs = 0; rhs < count; ++rhs) {
      source(b, local, seed + rhs);
      load(output + "." + method + ".rhs" + std::to_string(rhs) + ".bin", x,
           local);
      (*op)(ax, x, PLUS);
      r = b - ax;
      double rho = std::sqrt(toDouble(norm2(r)) / toDouble(norm2(b)));
      std::cout << "CHECK rhs=" << rhs << " method=" << method
                << " original_full_FP64_relative=" << rho << std::endl;
      if (!std::isfinite(rho) || rho > 1e-12)
        throw std::runtime_error("original Chroma full FP64 residual failed");
      csv << rhs << ',' << method << ',' << rho << '\n';
      csv.flush();
    }
    rusage u{};
    getrusage(RUSAGE_SELF, &u);
    std::cout << "PASS original Chroma FP64 verification peak_RSS_GiB="
              << double(u.ru_maxrss) / (1ull << 30) << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "ERROR " << e.what() << std::endl;
    status = 1;
  }
  QDP_finalize();
  return status;
}
