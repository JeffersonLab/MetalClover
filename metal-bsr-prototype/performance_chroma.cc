#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/reliable_bicgstab.h"
#include "actions/ferm/linop/eoprec_clover_dumb_linop_w.h"
#include "actions/ferm/linop/eoprec_clover_linop_w.h"
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
  if (argc != 6) {
    std::cerr << "Usage: performance_chroma ACTION_XML BRIDGE_PREFIX "
                 "OUTPUT_PREFIX cpu|verify RHS_COUNT\n";
    return 2;
  }
  std::string xml = argv[1], bridge = argv[2], output = argv[3], mode = argv[4];
  int count = std::stoi(argv[5]);
  int qargc = 1;
  QDP_initialize(&qargc, &argv);
  int status = 0;
  try {
    if (count < 1 || count > 192 || (mode != "cpu" && mode != "verify"))
      throw std::runtime_error("invalid benchmark mode/count");
    auto setup = Clock::now();
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
    if (Layout::numNodes() != 1)
      throw std::runtime_error("single rank required");
    std::ifstream input(bridge + ".bridge", std::ios::binary);
    auto h = production::read_header(input);
    for (int mu = 0; mu < 4; ++mu)
      if (dims[mu] != int(h.shape[mu]))
        throw std::runtime_error("geometry mismatch");
    if (toDouble(params.Mass) != h.mass ||
        toDouble(params.clovCoeffR) != h.clover_r ||
        toDouble(params.clovCoeffT) != h.clover_t ||
        toDouble(params.anisoParam.xi_0) != h.xi0 ||
        toDouble(params.anisoParam.nu) != h.nu)
      throw std::runtime_error("action/bridge mismatch");
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
    std::array<Z, 36> raw;
    // The bridge was made by native QIO, reunitarization, and native FP64
    // stout. It contains physical links with no anisotropy or fermion BC.
    // Apply BC once here; PeriodicFermState then preserves those links.
    for (uint32_t site = 0; site < h.volume(); ++site) {
      chroma_bridge::read(input, raw.data(), 36);
      uint32_t stride = 1;
      for (int mu = 0; mu < 4; ++mu) {
        int sign =
            (site / stride % h.shape[mu] == h.shape[mu] - 1) ? boundary[mu] : 1;
        stride *= h.shape[mu];
        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c) {
            auto a = raw[mu * 9 + r * 3 + c];
            auto &z = links[mu].elem(local[site]).elem().elem(r, c);
            z.real() = sign * a.re;
            z.imag() = sign * a.im;
          }
      }
    }
    std::unique_ptr<EvenOddPrecCloverLinOp> op;
    {
      Handle<FermState<LatticeFermionD, Gauge, Gauge>> state(
          new PeriodicFermState<LatticeFermionD, Gauge, Gauge>(links));
      op.reset(new EvenOddPrecCloverLinOp(state, params));
    }
    std::unique_ptr<EvenOddPrecDumbCloverFLinOp> single;
    if (mode == "cpu") {
      GaugeF low(4);
      for (int mu = 0; mu < 4; ++mu)
        low[mu] = links[mu];
      Handle<FermState<LatticeFermionF, GaugeF, GaugeF>> state(
          new PeriodicFermState<LatticeFermionF, GaugeF, GaugeF>(low));
      single.reset(new EvenOddPrecDumbCloverFLinOp(state, params));
    }
    links.resize(0);
    LatticeFermionD b, x, ax, r;
    // Independent original Chroma full operator, compared to the original
    // UnprecCloverLinOp probe exported separately from the LIME pipeline.
    source(x, local, chroma_bridge::probe_rhs);
    x *= std::sqrt(double(h.volume()) * 12);
    op->unprecLinOp(ax, x, PLUS);
    double error = 0, norm = 0;
    std::array<Z, 12> row;
    for (uint32_t site = 0; site < h.volume(); ++site) {
      chroma_bridge::read(input, row.data(), 12);
      for (int c = 0; c < 3; ++c)
        for (int s = 0; s < 4; ++s) {
          auto z = ax.elem(local[site]).elem(s).elem(c);
          auto a = row[c * 4 + s];
          error += (z.real() - a.re) * (z.real() - a.re) +
                   (z.imag() - a.im) * (z.imag() - a.im);
          norm += a.re * a.re + a.im * a.im;
        }
    }
    double match = std::sqrt(error / norm);
    if (!std::isfinite(match) || match > 1e-13)
      throw std::runtime_error("native EO/full operator mismatch");
    std::cout << std::setprecision(17) << "NATIVE full_operator_match=" << match
              << " setup_seconds=" << elapsed(setup) << std::endl;
    std::ofstream csv(
        output + (mode == "cpu" ? ".chroma-solves.csv" : ".chroma-checks.csv"));
    csv << std::setprecision(17);
    if (mode == "cpu")
      csv << "rhs,iterations,attempts,seconds,full_relative,save_seconds\n";
    else
      csv << "rhs,method,full_relative,relative_solution_difference\n";
    for (int rhs = 0; rhs < count; ++rhs) {
      source(b, local, seed + rhs);
      double bn = toDouble(norm2(b));
      auto residual = [&]() {
        op->unprecLinOp(ax, x, PLUS);
        r = b - ax;
        return std::sqrt(toDouble(norm2(r)) / bn);
      };
      if (mode == "cpu") {
        x = zero;
        r = b;
        int iterations = 0, attempts = 0;
        double rho = 1;
        auto start = Clock::now();
        // Original Chroma reliable mixed-precision BiCGStab, with the
        // original EO clover blocks for reduction and reconstruction.
        for (; attempts < 3 && rho > 1e-12; ++attempts) {
          LatticeFermionD reduced, dx, tmp, work;
          op->evenEvenInvLinOp(tmp, r, PLUS);
          op->oddEvenLinOp(work, tmp, PLUS);
          reduced[rb[1]] = r - work;
          dx = zero;
          auto info = InvBiCGStabReliable(*op, *single, reduced, dx,
                                          Real(5e-14), Real(.1), 20000, PLUS);
          iterations += info.n_count + 1;
          op->evenOddLinOp(tmp, dx, PLUS);
          work[rb[0]] = r - tmp;
          op->evenEvenInvLinOp(tmp, work, PLUS);
          dx[rb[0]] = tmp;
          x += dx;
          rho = residual();
          std::cout << "CHROMA rhs=" << rhs << " attempt=" << attempts + 1
                    << " iterations=" << info.n_count + 1
                    << " full_relative=" << rho << std::endl;
        }
        double seconds = elapsed(start);
        if (!std::isfinite(rho) || rho > 1e-12)
          throw std::runtime_error("native full FP64 residual failed");
        auto io = Clock::now();
        save(output + ".chroma-bicgstab.rhs" + std::to_string(rhs) + ".bin", x,
             local);
        csv << rhs << ',' << iterations << ',' << attempts << ',' << seconds
            << ',' << rho << ',' << elapsed(io) << '\n';
        csv.flush();
        std::cout << "CPU SOLVED rhs=" << rhs << " seconds=" << seconds
                  << std::endl;
      } else {
        LatticeFermionD reference;
        load(output + ".chroma-bicgstab.rhs" + std::to_string(rhs) + ".bin",
             reference, local);
        double xn = toDouble(norm2(reference));
        for (auto method : {"chroma-bicgstab", "metal-bicgstab", "metal-mg"}) {
          load(output + "." + method + ".rhs" + std::to_string(rhs) + ".bin", x,
               local);
          double rho = residual();
          r = x - reference;
          double diff = std::sqrt(toDouble(norm2(r)) / xn);
          std::cout << "CHECK rhs=" << rhs << " method=" << method
                    << " full_relative=" << rho
                    << " solution_difference=" << diff << std::endl;
          if (!std::isfinite(rho) || rho > 1e-12 || !std::isfinite(diff))
            throw std::runtime_error("native verification failed");
          csv << rhs << ',' << method << ',' << rho << ',' << diff << '\n';
          csv.flush();
        }
      }
    }
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "PASS mode=" << mode
              << " peak_rss_GiB=" << double(usage.ru_maxrss) / (1ull << 30)
              << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    status = 1;
  }
  QDP_finalize();
  return status;
}
