#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/linop/clover_term_qdp_w.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "chromabase.h"
#include "production_bridge.hpp"
#include "qdp_disk_map_slice.h"
#include "qdp_map_obj_disk.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include "util/gauge/reunit.h"
#include "util/gauge/stout_utils.h"
#include <filesystem>
#include <iomanip>
#include <sstream>

static_assert(BASE_PRECISION == 64,
              "Production action parameters must remain FP64");

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "Usage: production_export GAUGE_LIME PERAM_XML EIGEN_MOD "
                 "OUTPUT_PREFIX\n";
    return 2;
  }
  std::string gauge_path = argv[1], xml_path = argv[2], eigen_path = argv[3],
              prefix = argv[4];
  if (std::filesystem::exists(prefix + ".bridge") ||
      std::filesystem::exists(prefix + ".eigen0"))
    throw std::runtime_error("production export already exists");
  int qargc = 1;
  QDP::QDP_initialize(&qargc, &argv);
  int status = 0;
  try {
    using namespace QDP;
    using namespace Chroma;
    XMLReader metadata(xml_path);
    multi1d<int> nrow;
    read(metadata, "/DBMetaData/lattSize", nrow);
    if (nrow.size() != 4 || Layout::numNodes() != 1)
      throw std::runtime_error("expected single-rank 4D field");
    CloverFermActParams params(metadata,
                               "/DBMetaData/Params/Propagator/FermionAction");
    production::Header h;
    for (int mu = 0; mu < 4; ++mu)
      h.shape[mu] = nrow[mu];
    h.mass = toDouble(params.Mass);
    h.xi0 = toDouble(params.anisoParam.xi_0);
    h.nu = toDouble(params.anisoParam.nu);
    h.clover_r = toDouble(params.clovCoeffR);
    h.clover_t = toDouble(params.clovCoeffT);
    multi1d<int> sources;
    read(metadata, "/DBMetaData/Params/Contractions/t_sources", sources);
    if (sources.size() != 1)
      throw std::runtime_error("expected one perambulator source time");
    h.source_time = h.first_time = sources[0];
    read(metadata, "/DBMetaData/Params/Contractions/Nt_forward", h.time_count);
    h.validate();
    const std::string action = "/DBMetaData/Params/Propagator/FermionAction/";
    std::string fermact;
    read(metadata, action + "FermAct", fermact);
    if (fermact != "CLOVER" || toDouble(params.u0) != 1)
      throw std::runtime_error("expected CLOVER action with u0=1");
    multi1d<int> boundary;
    read(metadata, action + "FermState/FermionBC/boundary", boundary);
    if (boundary.size() != 4 || boundary[0] != 1 || boundary[1] != 1 ||
        boundary[2] != 1 || boundary[3] != -1 || !params.anisoParam.anisoP ||
        params.anisoParam.t_dir != 3)
      throw std::runtime_error("expected the production spatial-periodic "
                               "temporal-antiperiodic action");
    Real rho;
    int steps, orthogonal;
    std::string state_name;
    read(metadata, action + "FermState/rho", rho);
    read(metadata, action + "FermState/n_smear", steps);
    read(metadata, action + "FermState/orthog_dir", orthogonal);
    read(metadata, action + "FermState/Name", state_name);
    if (orthogonal != 3 || steps < 0 || state_name != "STOUT_FERM_STATE")
      throw std::runtime_error("unsupported stout state");
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
    std::cout << std::setprecision(17) << "ACTION mass=" << h.mass
              << " xi0=" << h.xi0 << " nu=" << h.nu << " cR=" << h.clover_r
              << " cT=" << h.clover_t << " rho=" << toDouble(rho)
              << " steps=" << steps << std::endl;
    {
      multi1d<LatticeColorMatrix> links(4);
      XMLReader file_xml, record_xml;
      QDPFileReader reader(file_xml, gauge_path, QDPIO_SERIAL);
      read(reader, record_xml, links);
      if (reader.bad())
        throw std::runtime_error("LIME/QIO read or checksum failed");
      reader.close();
      std::cout << "LIME gauge read and QIO checksum checked" << std::endl;
      // Match gaugeStartup: promote the disk field to FP64, then project
      // every link to SU(3) before the fermion state's stout smearing.
      for (int mu = 0; mu < 4; ++mu)
        reunit(links[mu]);
      std::cout << "Chroma gaugeStartup FP64 reunitarization complete"
                << std::endl;
      multi1d<bool> smear(4);
      smear = true;
      smear[3] = false;
      multi2d<Real> weights(4, 4);
      weights = rho;
      for (int mu = 0; mu < 4; ++mu) {
        weights(mu, mu) = 0;
        weights(3, mu) = weights(mu, 3) = 0;
      }
      for (int step = 0; step < steps; ++step) {
        multi1d<LatticeColorMatrix> next(4);
        Stouting::smear_links(links, next, smear, weights);
        links = next;
        std::cout << "STOUT step=" << step + 1 << " complete" << std::endl;
      }
      std::ofstream out(prefix + ".bridge", std::ios::binary);
      if (std::getenv("METALCLOVER_ALIGNED_BRIDGE")) {
        std::ostringstream header(std::ios::binary);
        production::write_header(header, h);
        auto bytes = header.str();
        bytes[7] = '3';
        out.write(bytes.data(), bytes.size());
        uint32_t padding = 0;
        chroma_bridge::write(out, &padding, 1);
      } else
        production::write_header(out, h);
      std::array<chroma_bridge::Complex, 36> site_links;
      for (uint32_t site = 0; site < h.volume(); ++site) {
        for (int mu = 0; mu < 4; ++mu)
          for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col) {
              const auto &z = links[mu].elem(local[site]).elem().elem(row, col);
              site_links[mu * 9 + row * 3 + col] = {z.real(), z.imag()};
            }
        chroma_bridge::write(out, site_links.data(), site_links.size());
      }
      using Gauge = multi1d<LatticeColorMatrix>;
      Handle<FermState<LatticeFermion, Gauge, Gauge>> state(
          new SimpleFermState<LatticeFermion, Gauge, Gauge>(boundary, links));
      links.resize(0);
      QDPCloverTermD clover;
      clover.create(state, params);
      QDPWilsonDslashD dslash(state, params.anisoParam);
      LatticeFermion x, full, hop;
      for (uint32_t site = 0; site < h.volume(); ++site)
        for (int color = 0; color < 3; ++color)
          for (int spin = 0; spin < 4; ++spin) {
            auto z = chroma_bridge::source(size_t(site) * 12 + color * 4 + spin,
                                           chroma_bridge::probe_rhs);
            auto &v = x.elem(local[site]).elem(spin).elem(color);
            v.real() = z.re;
            v.imag() = z.im;
          }
      clover(full, x, PLUS);
      for (int cb = 0; cb < 2; ++cb)
        dslash.apply(hop, x, PLUS, cb);
      full -= Double(.5) * hop;
      std::array<chroma_bridge::Complex, 12> values;
      for (uint32_t site = 0; site < h.volume(); ++site) {
        for (int color = 0; color < 3; ++color)
          for (int spin = 0; spin < 4; ++spin) {
            const auto &z = full.elem(local[site]).elem(spin).elem(color);
            values[color * 4 + spin] = {z.real(), z.imag()};
          }
        chroma_bridge::write(out, values.data(), values.size());
      }
      out.close();
      if (!out)
        throw std::runtime_error("bridge close failed");
      std::cout << "FP64 full-operator reference exported" << std::endl;
    }
    if (eigen_path == "-") {
      std::cout << "EXPORT reference-only PASS" << std::endl;
      QDP_finalize();
      return 0;
    }
    MapObjectDisk<KeyTimeSliceColorVec_t, TimeSliceIO<LatticeColorVectorF>>
        eigen;
    eigen.open(eigen_path, std::ios_base::in);
    std::ofstream out(prefix + ".eigen0", std::ios::binary);
    chroma_bridge::write(out, production::eigen_magic.data(), 8);
    chroma_bridge::write(out, &h.first_time, 1);
    chroma_bridge::write(out, &h.time_count, 1);
    const uint32_t spatial = h.shape[0] * h.shape[1] * h.shape[2];
    chroma_bridge::write(out, &spatial, 1);
    LatticeColorVectorF field = zero;
    std::vector<float> values(size_t(spatial) * 6);
    for (uint32_t time = h.first_time; time < h.first_time + h.time_count;
         ++time) {
      TimeSliceIO<LatticeColorVectorF> slice(field, time);
      if (eigen.get(KeyTimeSliceColorVec_t(time, 0), slice))
        throw std::runtime_error("missing first eigenvector");
      double norm = 0;
      for (uint32_t site = 0; site < spatial; ++site)
        for (int color = 0; color < 3; ++color) {
          const auto &z = field.elem(local[size_t(time) * spatial + site])
                              .elem()
                              .elem(color);
          size_t q = size_t(site) * 6 + color * 2;
          values[q] = z.real();
          values[q + 1] = z.imag();
          norm += double(values[q]) * values[q] +
                  double(values[q + 1]) * values[q + 1];
        }
      if (std::abs(norm - 1) > 1e-5)
        throw std::runtime_error("eigenvector normalization mismatch");
      chroma_bridge::write(out, values.data(), values.size());
    }
    eigen.close();
    out.close();
    if (!out)
      throw std::runtime_error("eigenvector close failed");
    std::cout << "EXPORT PASS source=" << h.source_time
              << " sink_count=" << h.time_count << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    status = 1;
  }
  QDP::QDP_finalize();
  return status;
}
