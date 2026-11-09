#include "chromabase.h"
#include "util/ferm/key_prop_matelem.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include "qdp_map_obj_disk.h"
#include "qdp_disk_map_slice.h"
#include <sstream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <set>
#include <tuple>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 3 && argc != 4) {
        std::cerr << "Usage: anisotropic_inspect PERAMBULATOR OUTPUT_PREFIX [EIGEN_MOD]\n";
        return 2;
    }
    const std::string input = argv[1], prefix = argv[2];
    const std::string eigen_path = argc == 4 ? argv[3] : "";
    int qdp_argc = 1;
    QDP::QDP_initialize(&qdp_argc, &argv);
    int status = 0;
    try {
        using namespace Chroma;
        using Key = SerialDBKey<KeyPropElementalOperator_t>;
        using Value = SerialDBData<ValPropElementalOperator_t>;
        QDP::BinaryStoreDB<Key, Value> db;
        db.open(input, O_RDONLY, 0);
        std::string metadata;
        db.getUserdata(metadata);
        std::ofstream(prefix + "-peram-metadata.xml") << metadata;
        std::vector<Key> keys;
        db.keys(keys);
        std::sort(keys.begin(), keys.end(), [](const Key& a, const Key& b) {
            const auto& x = a.key();
            const auto& y = b.key();
            return std::tie(x.t_source, x.t_slice, x.spin_src, x.spin_snk, x.mass_label) <
                   std::tie(y.t_source, y.t_slice, y.spin_src, y.spin_snk, y.mass_label);
        });
        std::set<int> sources, sinks;
        std::set<std::string> masses;
        std::ofstream csv(prefix + "-peram-00.csv");
        csv << "t_source,t_sink,spin_source,spin_sink,mass_label,real,imag\n"
            << std::setprecision(17);
        for (const auto& key : keys) {
            const auto& k = key.key();
            sources.insert(k.t_source);
            sinks.insert(k.t_slice);
            masses.insert(k.mass_label);
            Value value;
            if (db.get(key, value)) throw std::runtime_error("cannot read perambulator value");
            const auto& matrix = value.data().mat;
            if (matrix.size1() != 162 || matrix.size2() != 162)
                throw std::runtime_error("unexpected perambulator dimensions");
            auto z = matrix(0, 0);
            csv << k.t_source << ',' << k.t_slice << ',' << k.spin_src << ',' << k.spin_snk
                << ',' << k.mass_label << ',' << toDouble(real(z)) << ',' << toDouble(imag(z)) << '\n';
        }
        std::cout << "PERAM keys=" << keys.size() << " metadata_bytes=" << metadata.size();
        std::cout << " source_times=";
        for (auto t : sources) std::cout << t << ' ';
        std::cout << " sink_times=";
        for (auto t : sinks) std::cout << t << ' ';
        std::cout << " mass_labels=";
        for (const auto& m : masses) std::cout << m << ' ';
        std::cout << std::endl;
        db.close();
        if (!eigen_path.empty()) {
            std::istringstream text(metadata);
            XMLReader xml(text);
            multi1d<int> dimensions;
            read(xml, "/DBMetaData/lattSize", dimensions);
            QDP::Layout::setLattSize(dimensions);
            QDP::Layout::create();
            QDP::MapObjectDisk<KeyTimeSliceColorVec_t, TimeSliceIO<LatticeColorVectorF>> eigen;
            eigen.open(eigen_path, std::ios_base::in);
            std::vector<KeyTimeSliceColorVec_t> eigen_keys;
            eigen.keys(eigen_keys);
            if (eigen_keys.size() != 162 * 256)
                throw std::runtime_error("unexpected number of stored eigenvectors");
            std::cout << "EIGEN keys=" << eigen_keys.size() << std::endl;
            for (int t : {*sources.begin(), *sinks.rbegin()}) {
                LatticeColorVectorF field = zero;
                TimeSliceIO<LatticeColorVectorF> slice(field, t);
                if (eigen.get(KeyTimeSliceColorVec_t(t, 0), slice))
                    throw std::runtime_error("cannot read first eigenvector");
                double norm = toDouble(norm2(field));
                if (!std::isfinite(norm) || std::abs(norm - 1.) > 1e-5)
                    throw std::runtime_error("unexpected first-eigenvector norm");
                std::cout << std::setprecision(17) << "EIGEN vec=0 t=" << t
                          << " norm_squared=" << norm << std::endl;
            }
            eigen.close();
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        status = 1;
    }
    QDP::QDP_finalize();
    return status;
}
