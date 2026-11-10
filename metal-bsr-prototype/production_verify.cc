#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/linop/clover_term_qdp_w.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "chromabase.h"
#include "production_bridge.hpp"
#include <iomanip>

static_assert(BASE_PRECISION == 64,
              "Verification requires FP64 action parameters");
int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Usage: production_verify BRIDGE_PREFIX SOLUTION_PREFIX\n";
    return 2;
  }
  std::string prefix = argv[1], solutions = argv[2];
  int qargc = 1;
  QDP::QDP_initialize(&qargc, &argv);
  int status = 0;
  try {
    using namespace QDP;
    using namespace Chroma;
    std::ifstream in(prefix + ".bridge", std::ios::binary);
    auto h = production::read_header(in);
    multi1d<int> nrow(4), boundary(4);
    boundary = 1;
    boundary[3] = -1;
    for (int mu = 0; mu < 4; ++mu)
      nrow[mu] = h.shape[mu];
    Layout::setLattSize(nrow);
    Layout::create();
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
    using Gauge = multi1d<LatticeColorMatrixD>;
    Gauge links(4);
    std::array<chroma_bridge::Complex, 36> raw;
    for (uint32_t site = 0; site < h.volume(); ++site) {
      chroma_bridge::read(in, raw.data(), raw.size());
      for (int mu = 0; mu < 4; ++mu)
        for (int row = 0; row < 3; ++row)
          for (int col = 0; col < 3; ++col) {
            auto z = raw[mu * 9 + row * 3 + col];
            auto &v = links[mu].elem(local[site]).elem().elem(row, col);
            v.real() = z.re;
            v.imag() = z.im;
          }
    }
    in.close();
    CloverFermActParams params;
    params.Mass = h.mass;
    params.clovCoeffR = h.clover_r;
    params.clovCoeffT = h.clover_t;
    params.u0 = 1;
    params.anisoParam.anisoP = true;
    params.anisoParam.t_dir = 3;
    params.anisoParam.xi_0 = h.xi0;
    params.anisoParam.nu = h.nu;
    Handle<FermState<LatticeFermionD, Gauge, Gauge>> state(
        new SimpleFermState<LatticeFermionD, Gauge, Gauge>(boundary, links));
    links.resize(0);
    QDPCloverTermD clover;
    clover.create(state, params);
    QDPWilsonDslashD dslash(state, params.anisoParam);
    std::ifstream eig(prefix + ".eigen0", std::ios::binary);
    std::array<char, 8> tag;
    chroma_bridge::read(eig, tag.data(), 8);
    if (tag != production::eigen_magic)
      throw std::runtime_error("bad eigenvector tag");
    uint32_t first, count, spatial;
    for (auto *n : {&first, &count, &spatial})
      chroma_bridge::read(eig, n, 1);
    if (first != h.first_time || count != h.time_count ||
        spatial != h.shape[0] * h.shape[1] * h.shape[2])
      throw std::runtime_error("eigenvector shape mismatch");
    eig.seekg(std::streamoff(h.source_time - first) * spatial * 6 *
                  sizeof(float),
              std::ios::cur);
    std::vector<float> source(size_t(spatial) * 6);
    chroma_bridge::read(eig, source.data(), source.size());
    LatticeFermionD x, rhs, full, hop;
    std::array<chroma_bridge::Complex, 12> values;
    std::ofstream csv(solutions + ".chroma-residuals.csv");
    csv << "spin_source,full_relative\n" << std::setprecision(17);
    for (uint32_t spin = 0; spin < 4; ++spin) {
      std::ifstream solution(solutions + ".spin" + std::to_string(spin) +
                                 ".bin",
                             std::ios::binary);
      for (uint32_t site = 0; site < h.volume(); ++site) {
        chroma_bridge::read(solution, values.data(), values.size());
        for (int color = 0; color < 3; ++color)
          for (int s = 0; s < 4; ++s) {
            auto z = values[color * 4 + s];
            auto &v = x.elem(local[site]).elem(s).elem(color);
            v.real() = z.re;
            v.imag() = z.im;
          }
      }
      if (solution.peek() != EOF)
        throw std::runtime_error("solution payload length mismatch");
      rhs = zero;
      for (uint32_t site = 0; site < spatial; ++site)
        for (int color = 0; color < 3; ++color) {
          auto &v = rhs.elem(local[size_t(h.source_time) * spatial + site])
                        .elem(spin)
                        .elem(color);
          v.real() = source[size_t(site) * 6 + color * 2];
          v.imag() = source[size_t(site) * 6 + color * 2 + 1];
        }
      clover(full, x, PLUS);
      for (int cb = 0; cb < 2; ++cb)
        dslash.apply(hop, x, PLUS, cb);
      full -= Double(.5) * hop;
      full -= rhs;
      double relative = std::sqrt(toDouble(norm2(full)) / toDouble(norm2(rhs)));
      if (!std::isfinite(relative) || relative > 1e-10)
        throw std::runtime_error("independent Chroma residual failed");
      csv << spin << ',' << relative << '\n';
      csv.flush();
      std::cout << std::setprecision(17) << "CHROMA full residual spin=" << spin
                << " relative=" << relative << std::endl;
    }
    std::cout << "CHROMA VERIFY PASS" << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    status = 1;
  }
  QDP::QDP_finalize();
  return status;
}
