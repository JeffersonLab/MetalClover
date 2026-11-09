#include "chromabase.h"
#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/linop/clover_term_qdp_w.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "qdp_iogauge.h"

#include "chroma_bridge.hpp"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

namespace {
using namespace QDP;
using namespace Chroma;
using Gauge = multi1d<LatticeColorMatrixD>;
using BridgeComplex = chroma_bridge::Complex;

uint32_t nersc_extent(const char* path) {
    std::ifstream file(path, std::ios::binary);
    std::string line;
    std::array<uint32_t, 4> dims{};
    if (!std::getline(file, line) || line != "BEGIN_HEADER")
        throw std::runtime_error("expected a NERSC gauge archive");
    for (int n = 0; n < 64 && std::getline(file, line); ++n) {
        if (line == "END_HEADER") {
            chroma_bridge::volume(dims[0]);
            for (uint32_t d : dims)
                if (d != dims[0]) throw std::runtime_error("expected four equal NERSC dimensions");
            return dims[0];
        }
        std::istringstream entry(line);
        std::string name, equals;
        entry >> name >> equals;
        if (name.size() == 11 && name.substr(0, 10) == "DIMENSION_" &&
            name.back() >= '1' && name.back() <= '4' && equals == "=")
            entry >> dims[name.back() - '1'];
    }
    throw std::runtime_error("incomplete NERSC header");
}

void process_archive(const char* gauge_path, const char* payload, double requested_mass,
                     bool verify, const char* rhs_path) {
    const uint32_t extent = nersc_extent(gauge_path);
    const uint32_t sites = chroma_bridge::volume(extent);
    if (Layout::numNodes() != 1 || Nd != 4 || Nc != 3 || Ns != 4)
        throw std::runtime_error("expected single-rank four-dimensional SU(3) Chroma");
    multi1d<int> dims(4), boundary(4);
    dims = extent;
    boundary = 1;
    boundary[3] = -1;
    Layout::setLattSize(dims);
    Layout::create();

    multi1d<LatticeColorMatrix> archive(4);
    ArchivGauge_t header;
    readArchiv(header, archive, gauge_path);
    Gauge links(4);
    for (int mu = 0; mu < 4; ++mu) links[mu] = archive[mu];

    CloverFermActParams params;
    params.Mass = requested_mass;
    params.clovCoeffR = params.clovCoeffT = params.u0 = 1;
    params.anisoParam.anisoP = false;
    const double mass = toDouble(params.Mass);
    if (!std::isfinite(mass) || !(mass > -.5 && mass < 1))
        throw std::runtime_error("mass outside compact solver range");

    using Fermion = LatticeFermionD;
    Handle<FermState<Fermion, Gauge, Gauge>> state(
        new SimpleFermState<Fermion, Gauge, Gauge>(boundary, links));
    QDPCloverTermD clover;
    clover.create(state, params);
    QDPWilsonDslashD dslash(state, params.anisoParam);

    std::vector<int> local(sites);
    std::vector<bool> odd(sites);
    multi1d<int> coord(4);
    for (uint32_t site = 0; site < sites; ++site) {
        uint32_t q = site;
        int parity = 0;
        for (int mu = 0; mu < 4; ++mu) {
            coord[mu] = q % extent;
            q /= extent;
            parity += coord[mu];
        }
        odd[site] = parity % 2;
        local[site] = Layout::linearSiteIndex(coord);
    }

    auto fill = [&](Fermion& field, const BridgeComplex* values, uint32_t rhs) {
        for (uint32_t site = 0; site < sites; ++site)
            for (int color = 0; color < 3; ++color)
                for (int spin = 0; spin < 4; ++spin) {
                    size_t q = size_t(site) * 12 + color * 4 + spin;
                    auto z = values ? values[q] : chroma_bridge::source(q, rhs);
                    auto& v = field.elem(local[site]).elem(spin).elem(color);
                    v.real() = z.re;
                    v.imag() = z.im;
                }
    };
    Fermion x, full, hop;
    auto apply_full = [&] {
        clover(full, x, PLUS);
        for (int cb = 0; cb < 2; ++cb) dslash.apply(hop, x, PLUS, cb);
        full -= Double(0.5) * hop;
    };

    if (verify) {
        std::ifstream solutions(payload, std::ios::binary | std::ios::ate);
        const std::streamoff bytes_per_rhs = std::streamoff(sites) * 12 * sizeof(BridgeComplex);
        const std::streamoff bytes = solutions.tellg();
        if (!solutions || bytes <= 0 || bytes % bytes_per_rhs ||
            bytes / bytes_per_rhs > chroma_bridge::probe_rhs)
            throw std::runtime_error("invalid Chroma solutions file size");
        solutions.seekg(0);
        std::ifstream source_file;
        if (rhs_path) {
            source_file.open(rhs_path, std::ios::binary | std::ios::ate);
            if (!source_file || source_file.tellg() != bytes)
                throw std::runtime_error("RHS and solutions must have the same size");
            source_file.seekg(0);
        }
        std::vector<BridgeComplex> buffer(size_t(sites) * 12);
        Fermion rhs;
        double worst = 0;
        for (uint32_t r = 0; r < bytes / bytes_per_rhs; ++r) {
            chroma_bridge::read(solutions, buffer.data(), buffer.size());
            fill(x, buffer.data(), r);
            apply_full();
            if (rhs_path) chroma_bridge::read(source_file, buffer.data(), buffer.size());
            fill(rhs, rhs_path ? buffer.data() : nullptr, r);
            full -= rhs;
            double norm = std::sqrt(toDouble(norm2(rhs)));
            double relative = std::sqrt(toDouble(norm2(full))) / (norm > 0 ? norm : 1);
            if (!std::isfinite(relative) || relative > 1e-10)
                throw std::runtime_error("Chroma FP64 residual failed on RHS " + std::to_string(r));
            worst = std::max(worst, relative);
        }
        std::cout << "CHROMA SOLUTIONS PASS rhs=" << bytes / bytes_per_rhs
                  << " worst_full_relative=" << worst << std::endl;
        return;
    }

    std::vector<BridgeComplex> packed_links(size_t(sites) * 36);
    for (uint32_t site = 0; site < sites; ++site) {
        for (int mu = 0; mu < 4; ++mu)
            for (int row = 0; row < 3; ++row)
                for (int col = 0; col < 3; ++col) {
                    const auto& z = links[mu].elem(local[site]).elem().elem(row, col);
                    packed_links[(size_t(site) * 4 + mu) * 9 + row * 3 + col] =
                        {z.real(), z.imag()};
                }
    }

    // Same decomposition as UnprecCloverLinOp and unscaled EvenOddPrecCloverLinOp.
    fill(x, nullptr, chroma_bridge::probe_rhs);
    apply_full();
    QDPCloverTermD inverse;
    inverse.create(state, params, clover);
    inverse.choles(0);
    Fermion schur, even, inverse_even;
    clover.apply(schur, x, PLUS, 1);
    dslash.apply(even, x, PLUS, 0);
    inverse.apply(inverse_even, even, PLUS, 0);
    dslash.apply(hop, inverse_even, PLUS, 1);
    schur[rb[1]] -= Double(0.25) * hop;

    auto pack = [&](const Fermion& field, bool odd_only) {
        std::vector<BridgeComplex> values;
        values.reserve(size_t(odd_only ? sites / 2 : sites) * 12);
        for (uint32_t site = 0; site < sites; ++site) {
            if (odd_only && !odd[site]) continue;
            for (int color = 0; color < 3; ++color)
                for (int spin = 0; spin < 4; ++spin) {
                    const auto& z = field.elem(local[site]).elem(spin).elem(color);
                    values.push_back({z.real(), z.imag()});
                }
        }
        return values;
    };
    std::ofstream out(payload, std::ios::binary);
    chroma_bridge::write(out, chroma_bridge::magic.data(), chroma_bridge::magic.size());
    chroma_bridge::write(out, &extent, 1);
    chroma_bridge::write(out, &mass, 1);
    chroma_bridge::write(out, packed_links.data(), packed_links.size());
    for (bool odd_only : {false, true}) {
        const auto values = pack(odd_only ? schur : full, odd_only);
        chroma_bridge::write(out, values.data(), values.size());
    }
    out.close();
    if (!out) throw std::runtime_error("could not close Chroma bridge output");
    std::cout << std::setprecision(12) << "CHROMA EXPORT L=" << extent
              << " mass=" << mass << " c_sw=1 boundary=(1,1,1,-1)"
              << " archive_plaquette=" << toDouble(header.w_plaq)
              << " reference=FP64" << std::endl;
}
} // namespace

int main(int argc, char** argv) {
    bool verify = argc > 1 && std::string(argv[1]) == "--check-solutions";
    if ((!verify && argc != 4) || (verify && argc != 5 && argc != 6)) {
        std::cerr << "Usage: chroma_export NERSC_CONFIG OUTPUT_BRIDGE MASS\n"
                  << "       chroma_export --check-solutions NERSC_CONFIG SOLUTIONS MASS [RHS_FILE]\n";
        return 2;
    }
    const char* gauge = argv[1 + verify];
    const char* payload = argv[2 + verify];
    const std::string mass = argv[3 + verify];
    const char* rhs_path = verify && argc == 6 ? argv[5] : nullptr;
    // QDP's own command-line parser is independent of the exporter arguments.
    int qdp_argc = 1;
    QDP::QDP_initialize(&qdp_argc, &argv);
    int status = 0;
    try {
        size_t parsed = 0;
        double requested_mass = std::stod(mass, &parsed);
        if (parsed != mass.size()) throw std::runtime_error("invalid mass");
        process_archive(gauge, payload, requested_mass, verify, rhs_path);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        status = 1;
    }
    QDP::QDP_finalize();
    return status;
}
