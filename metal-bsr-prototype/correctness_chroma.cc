#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"
#include "actions/ferm/invert/invbicgstab.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"
#include "chromabase.h"
#include "production_bridge.hpp"

#include "util/gauge/reunit.h"
#include <chrono>
#include <iomanip>

using namespace MetalCloverEnv;

static_assert(BASE_PRECISION == 64, "Independent reference must be FP64");
namespace {
  using namespace QDP;
  using namespace Chroma;
  using Z = chroma_bridge::Complex;
  using Gauge = multi1d<LatticeColorMatrix>;
  constexpr uint32_t seed = 20261107;
  void save(const std::string &path, const LatticeFermion &x,
            const std::vector<int> &local) {
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
  }
  void load(const std::string &path, LatticeFermion &x,
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
} // namespace
int main(int argc, char **argv) {
  if (argc != 6) {
    std::cerr << "Usage: correctness_chroma GAUGE ACTION_XML BRIDGE_PREFIX "
                 "OUTPUT_PREFIX reference|verify|state-only\n";
    return 2;
  }
  std::string gauge = argv[1], xml = argv[2], bridge = argv[3],
              output = argv[4], mode = argv[5];
  int qargc = 1;
  QDP_initialize(&qargc, &argv);
  int status = 0;
  try {
    XMLReader metadata(xml);
    multi1d<int> dims, boundary;
    read(metadata, "/DBMetaData/lattSize", dims);
    read(metadata,
         "/DBMetaData/Params/Propagator/FermionAction/FermState/FermionBC/"
         "boundary",
         boundary);
    CloverFermActParams params(metadata,
                               "/DBMetaData/Params/Propagator/FermionAction");
    StoutFermStateParams stout(
        metadata, "/DBMetaData/Params/Propagator/FermionAction/FermState");
    Layout::setLattSize(dims);
    Layout::create();
    if (Layout::numNodes() != 1)
      throw std::runtime_error("single rank required");
    std::ifstream interchange(bridge + ".bridge", std::ios::binary);
    auto h = production::read_header(interchange);
    for (int mu = 0; mu < 4; ++mu)
      if (dims[mu] != int(h.shape[mu]))
        throw std::runtime_error("geometry mismatch");
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
    XMLReader file_xml, record_xml;
    QDPFileReader reader(file_xml, gauge, QDPIO_SERIAL);
    read(reader, record_xml, links);
    if (reader.bad())
      throw std::runtime_error("QIO read/checksum failure");
    reader.close();
    for (int mu = 0; mu < 4; ++mu)
      reunit(links[mu]);
    Handle<FermBC<LatticeFermion, Gauge, Gauge>> bc(
        new SimpleFermBC<LatticeFermion, Gauge, Gauge>(boundary));
    Handle<FermState<LatticeFermion, Gauge, Gauge>> state(
        new StoutFermState<LatticeFermion, Gauge, Gauge>(bc, stout, links));
    links.resize(0);
    double difference = 0, norm = 0;
    std::array<Z, 36> raw;
    for (uint32_t site = 0; site < h.volume(); ++site) {
      chroma_bridge::read(interchange, raw.data(), 36);
      for (int mu = 0; mu < 4; ++mu)
        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c) {
            auto z = state->getLinks()[mu].elem(local[site]).elem().elem(r, c);
            auto a = raw[mu * 9 + r * 3 + c];
            uint32_t stride = 1;
            for (int k = 0; k < mu; ++k)
              stride *= h.shape[k];
            int sign = (site / stride % h.shape[mu] == h.shape[mu] - 1)
                           ? boundary[mu]
                           : 1;
            difference += (z.real() - sign * a.re) * (z.real() - sign * a.re) +
                          (z.imag() - sign * a.im) * (z.imag() - sign * a.im);
            norm += a.re * a.re + a.im * a.im;
          }
    }
    double link_error = std::sqrt(difference / norm);
    std::cout << std::setprecision(17)
              << "NATIVE StoutFermState link_match=" << link_error << std::endl;
    if (!std::isfinite(link_error) || link_error > 1e-14)
      throw std::runtime_error("independent stout/BC pipeline differs");
    if (mode == "state-only") {
      std::cout << "STATE PASS independent LIME/QIO, reunitarization, native "
                   "StoutFermState and boundary conditions"
                << std::endl;
    } else {
      UnprecCloverLinOp op(state, params);
      LatticeFermion b, x, ax, r, reference;
      double scale = 1 / std::sqrt(double(h.volume()) * 12);
      for (uint32_t site = 0; site < h.volume(); ++site)
        for (int c = 0; c < 3; ++c)
          for (int s = 0; s < 4; ++s) {
            auto z = chroma_bridge::source(size_t(site) * 12 + c * 4 + s, seed);
            auto &v = b.elem(local[site]).elem(s).elem(c);
            v.real() = scale * z.re;
            v.imag() = scale * z.im;
            auto p = chroma_bridge::source(size_t(site) * 12 + c * 4 + s,
                                           chroma_bridge::probe_rhs);
            auto &a = x.elem(local[site]).elem(s).elem(c);
            a.real() = p.re;
            a.imag() = p.im;
          }
      op(ax, x, PLUS);
      double error = 0, reference_norm = 0;
      std::array<Z, 12> row;
      for (uint32_t site = 0; site < h.volume(); ++site) {
        chroma_bridge::read(interchange, row.data(), 12);
        for (int c = 0; c < 3; ++c)
          for (int s = 0; s < 4; ++s) {
            auto z = ax.elem(local[site]).elem(s).elem(c);
            auto a = row[c * 4 + s];
            error += (z.real() - a.re) * (z.real() - a.re) +
                     (z.imag() - a.im) * (z.imag() - a.im);
            reference_norm += a.re * a.re + a.im * a.im;
          }
      }
      double operator_error = std::sqrt(error / reference_norm);
      std::cout << "NATIVE original UnprecCloverLinOp match=" << operator_error
                << std::endl;
      if (!std::isfinite(operator_error) || operator_error > 1e-13)
        throw std::runtime_error("original Chroma full operator mismatch");
      const double bn = toDouble(norm2(b));
      auto residual = [&]() {
        op(ax, x, PLUS);
        r = b - ax;
        return std::sqrt(toDouble(norm2(r)) / bn);
      };
      if (mode == "reference") {
        save(bridge + ".rhs.bin", b, local);
        x = zero;
        auto start = std::chrono::steady_clock::now();
        int iterations = 0;
        double rho = 1;
        for (int attempt = 0; attempt < 3 && rho > 1e-12; ++attempt) {
          auto result = InvBiCGStab(op, b, x, Real(5e-14), 20000, PLUS);
          iterations += result.n_count;
          rho = residual();
          std::cout << "CHROMA BICGSTAB attempt=" << attempt
                    << " iterations=" << result.n_count
                    << " full_relative=" << rho << std::endl;
        }
        if (!std::isfinite(rho) || rho > 1e-12)
          throw std::runtime_error(
              "native FP64 BiCGStab did not reach full target");
        double seconds = std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - start)
                             .count();
        save(output + ".chroma-bicgstab.bin", x, local);
        std::ofstream out(output + ".chroma-reference.csv");
        out << std::setprecision(17)
            << "iterations,seconds,full_relative,stout_link_match,original_"
               "operator_match\n"
            << iterations << ',' << seconds << ',' << rho << ',' << link_error
            << ',' << operator_error << '\n';
        std::cout << "REFERENCE PASS seconds=" << seconds << std::endl;
      } else if (mode == "verify") {
        load(output + ".chroma-bicgstab.bin", reference, local);
        std::ofstream out(output + ".chroma-checks.csv");
        out << std::setprecision(17)
            << "method,full_relative,relative_solution_difference\n";
        double reference_size = toDouble(norm2(reference));
        for (auto method : {"chroma-bicgstab", "metal-bicgstab", "metal-mg"}) {
          load(output + "." + method + ".bin", x, local);
          double rho = residual();
          r = x - reference;
          double diff = std::sqrt(toDouble(norm2(r)) / reference_size);
          std::cout << "CHROMA CHECK method=" << method
                    << " full_relative=" << rho
                    << " relative_solution_difference=" << diff << std::endl;
          if (!std::isfinite(rho) || rho > 1e-12 || !std::isfinite(diff))
            throw std::runtime_error("independent full residual failed");
          out << method << ',' << rho << ',' << diff << '\n';
          out.flush();
        }
        std::cout << "VERIFY PASS all three full native FP64 residuals <= 1e-12"
                  << std::endl;
      } else
        throw std::runtime_error("unknown mode");
    }
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    status = 1;
  }
  QDP_finalize();
  return status;
}
