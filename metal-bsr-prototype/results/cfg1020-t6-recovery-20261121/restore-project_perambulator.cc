// Project verified same-time point sinks into the supplied LapH basis and
// write the standard Chroma propElemOp key/value database.
#include <Accelerate/Accelerate.h>

#include "chromabase.h"
#include "qdp_disk_map_slice.h"
#include "qdp_map_obj_disk.h"
#include "util/ferm/key_prop_matelem.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace Chroma;
namespace {
  using ProjectComplex = std::complex<double>;
  static_assert(sizeof(ProjectComplex) == 16);
  static_assert(std::endian::native == std::endian::little);
  std::string text_file(const std::string &path) {
    std::ifstream f(path, std::ios::binary);
    if (!f)
      throw std::runtime_error("cannot read " + path);
    return {std::istreambuf_iterator<char>(f), {}};
  }
  std::string fragment(std::string s) {
    auto first = s.find_first_not_of(" \t\r\n");
    if (first != std::string::npos && s.compare(first, 5, "<?xml") == 0) {
      auto end = s.find("?>", first);
      if (end == std::string::npos)
        throw std::runtime_error("invalid XML declaration");
      s.erase(first, end + 2 - first);
    }
    return s;
  }
  std::string fnv(const std::string &s) {
    uint64_t h = 14695981039346656037ull;
    for (unsigned char c : s)
      h = (h ^ c) * 1099511628211ull;
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16) << h;
    return out.str();
  }
  std::string header(int n, const multi1d<int> &dims) {
    std::ostringstream out;
    out << "{'descr': '<c16', 'fortran_order': False, 'shape': (" << n
        << ", 4, " << dims[2] << ", " << dims[1] << ", " << dims[0]
        << ", 3, 4), }";
    auto s = out.str();
    s.append((64 - ((10 + s.size() + 1) % 64)) % 64, ' ');
    s += '\n';
    std::string h("\x93NUMPY\x01\x00", 8);
    h += char(s.size() & 255);
    h += char((s.size() >> 8) & 255);
    return h + s;
  }
  struct PointFile {
    int fd = -1;
    void *base = MAP_FAILED;
    size_t bytes = 0;
    const ProjectComplex *data = nullptr;
    PointFile(const std::string &path, int n, const multi1d<int> &dims) {
      const auto h = header(n, dims);
      bytes =
          h.size() + uint64_t(n) * 4 * dims[0] * dims[1] * dims[2] * 12 * 16;
      fd = open(path.c_str(), O_RDONLY);
      struct stat st{};
      if (fd < 0 || fstat(fd, &st) || uint64_t(st.st_size) != bytes)
        throw std::runtime_error("point-file size/open mismatch");
      base = mmap(nullptr, bytes, PROT_READ, MAP_PRIVATE, fd, 0);
      if (base == MAP_FAILED || std::memcmp(base, h.data(), h.size()))
        throw std::runtime_error("point-file mapping/header mismatch");
      data = reinterpret_cast<const ProjectComplex *>(
          static_cast<const char *>(base) + h.size());
    }
    ~PointFile() {
      if (base != MAP_FAILED)
        munmap(base, bytes);
      if (fd >= 0)
        close(fd);
    }
  };
  using EigenDB = QDP::MapObjectDisk<KeyTimeSliceColorVec_t,
                                     TimeSliceIO<LatticeColorVectorF>>;
  using Key = KeyPropElementalOperator_t;
  using Value = ValPropElementalOperator_t;
  using Store = Datalib::TensorStorage::StorageRecords;
  ProjectComplex value(const ComplexD &z) {
    return {toDouble(real(z)), toDouble(imag(z))};
  }
  Key key(int source, int sink, int time_slice) {
    KeyPropElementalOperator_t k;
    k.t_source = k.t_slice = time_slice;
    k.spin_src = source;
    k.spin_snk = sink;
    k.mass_label = "light";
    return k;
  }
  void get(Store &db, int source, int sink, Value &v, int n, int time_slice) {
    std::string encoded;
    if (db.get_record(encodeStorageRecord(key(source, sink, time_slice)),
                      encoded))
      throw std::runtime_error("standard Chroma record readback failed");
    QDP::BinaryBufferReader reader(encoded);
    read(reader, v);
    if (reader.fail() || v.mat.size1() != n || v.mat.size2() != n)
      throw std::runtime_error("standard Chroma value decoding failed");
  }
} // namespace

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "Usage: project_perambulator POINT_PREFIX EIG_MOD OUTPUT_SDB "
                 "NUM_VECS\n";
    return 2;
  }
  const std::string prefix = argv[1], eig = argv[2], output = argv[3];
  const int n = std::stoi(argv[4]);
  if (n < 1 || n > 5184)
    throw std::runtime_error("invalid vector count");
  for (const auto &suffix :
       {"", ".partial", ".projection.csv", ".projection.xml"})
    if (std::filesystem::exists(output + suffix))
      throw std::runtime_error(
          "refusing to overwrite existing projection output");
  int qargc = 1;
  QDP_initialize(&qargc, &argv);
  int status = 0;
  try {
    static_assert(BASE_PRECISION == 64 && Ns == 4 && Nc == 3);
    StopWatch total;
    total.start();
    XMLReader meta(prefix + ".meta.xml");
    multi1d<int> dims;
    int input_n, first, source_time, sink_time, precision;
    read(meta, "/MetalCloverUnsmeared/lattSize", dims);
    read(meta, "/MetalCloverUnsmeared/num_vectors", input_n);
    read(meta, "/MetalCloverUnsmeared/first_vector", first);
    read(meta, "/MetalCloverUnsmeared/source_time", source_time);
    read(meta, "/MetalCloverUnsmeared/sink_time", sink_time);
    read(meta, "/MetalCloverUnsmeared/input_eigen_precision_bits", precision);
    if (dims.size() != 4 || dims[0] != 12 || dims[1] != 12 || dims[2] != 12 ||
        dims[3] != 96 || input_n < n || first != 0 || source_time < 0 ||
        source_time >= dims[3] || sink_time != source_time || precision != 32)
      throw std::runtime_error("unsupported point-file metadata");
    Layout::setLattSize(dims);
    Layout::create();
    if (Layout::numNodes() != 1)
      throw std::runtime_error("single rank required");
    std::ifstream state(prefix + ".state");
    std::string magic, signature;
    uint64_t completed, csv_bytes;
    double residual;
    if (!(state >> magic >> signature >> completed >> csv_bytes >> residual) ||
        magic != "METAL_UNSMEARED_V1" || completed != uint64_t(input_n) * 4 ||
        !(residual >= 0 && residual <= 1e-12))
      throw std::runtime_error("point collection is incomplete/unaccepted");
    PointFile points(prefix + ".npy", input_n, dims);
    const int k = dims[0] * dims[1] * dims[2] * 3;
    std::vector<int> sites(k / 3);
    multi1d<int> coord(4);
    coord[3] = sink_time;
    for (size_t i = 0; i < sites.size(); ++i) {
      size_t x = i;
      for (int d = 0; d < 3; ++d) {
        coord[d] = x % dims[d];
        x /= dims[d];
      }
      sites[i] = Layout::linearSiteIndex(coord);
    }
    const auto point = [&](int vector, int src, int p, int snk) {
      return points.data[(uint64_t(vector) * 4 + src) * k * 4 + p * 4 + snk];
    };
    EigenDB eigen;
    eigen.open(eig, std::ios_base::in);
    std::string eigen_xml;
    eigen.getUserdata(eigen_xml);
    std::istringstream eigen_text(eigen_xml);
    XMLReader eigen_meta(eigen_text);
    int available, eigen_time;
    read(eigen_meta, "/MODMetaData/num_vecs", available);
    read(eigen_meta, "/MODMetaData/t_slice", eigen_time);
    if (available < n || eigen_time != sink_time)
      throw std::runtime_error("eigenvector metadata mismatch");
    std::vector<ProjectComplex> basis(uint64_t(k) * n), rhs(uint64_t(k) * n),
        projected(uint64_t(n) * n);
    LatticeColorVectorF vf = zero;
    for (int v = 0; v < n; ++v) {
      TimeSliceIO<LatticeColorVectorF> slice(vf, sink_time);
      if (eigen.get(KeyTimeSliceColorVec_t(sink_time, v), slice))
        throw std::runtime_error("cannot read eigenvector");
      for (size_t i = 0; i < sites.size(); ++i)
        for (int c = 0; c < 3; ++c) {
          auto z = vf.elem(sites[i]).elem().elem(c);
          basis[uint64_t(v) * k + i * 3 + c] = {double(z.real()),
                                                double(z.imag())};
        }
    }
    XMLBufferWriter metadata;
    push(metadata, "DBMetaData");
    write(metadata, "id", std::string("propElemOp"));
    write(metadata, "lattSize", dims);
    write(metadata, "decay_dir", 3);
    push(metadata, "Params");
    XMLReader prop(
        meta, "/MetalCloverUnsmeared/ActionXML/DBMetaData/Params/Propagator");
    write(metadata, "Propagator", prop);
    push(metadata, "Contractions");
    write(metadata, "num_vecs", n);
    multi1d<int> times(1), phase(3);
    times[0] = source_time;
    phase = 0;
    write(metadata, "t_sources", times);
    write(metadata, "Nt_forward", 1);
    write(metadata, "Nt_backward", 0);
    write(metadata, "decay_dir", 3);
    write(metadata, "mass_label", std::string("light"));
    write(metadata, "phase", phase);
    write(metadata, "use_superb_format", false);
    pop(metadata);
    pop(metadata);
    push(metadata, "ProjectionProvenance");
    write(metadata, "definition", std::string("V_dagger M_inverse V"));
    write(metadata, "point_prefix", prefix);
    write(metadata, "eigenvector_file", eig);
    write(metadata, "source_order",
          std::string("preserved supplied MOD ordering and phases"));
    write(metadata, "original_full_FP64_tolerance", Double(1e-12));
    write(metadata, "original_max_full_relative", Double(residual));
    write(metadata, "point_metadata",
          fragment(text_file(prefix + ".meta.xml")));
    pop(metadata);
    pop(metadata);
    const auto metadata_text = metadata.str();
    Datalib::TensorStorage::StorageRecords writer;
    writer.record_options().metadata_capacity = metadata_text.size();
    writer.record_options().cache_bytes = 50 * 1024 * 1024;
    writer.record_options().configurations = 1;
    const auto partial = output + ".partial";
    if (writer.open_records(partial, O_RDWR | O_CREAT, 0644) ||
        writer.write_metadata(metadata_text))
      throw std::runtime_error("create standard perambulator failed");
    std::ofstream journal(output + ".projection.csv");
    journal << "spin_source,spin_sink,num_vectors,gemm_seconds,qdp_dot_scaled_"
               "error,encoded_fnv1a64\n"
            << std::setprecision(17);
    std::array<std::array<std::string, 4>, 4> hashes;
    double max_check = 0, gemm_seconds = 0;
    LatticeColorVectorD v = zero, w = zero;
    for (int src = 0; src < 4; ++src)
      for (int snk = 0; snk < 4; ++snk) {
        for (int j = 0; j < n; ++j)
          for (int p = 0; p < k; ++p)
            rhs[uint64_t(j) * k + p] = point(j, src, p, snk);
        const ProjectComplex alpha{1, 0}, beta{0, 0};
        StopWatch multiply;
        multiply.start();
        cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, n, n, k,
                    &alpha, basis.data(), k, rhs.data(), k, &beta,
                    projected.data(), n);
        multiply.stop();
        gemm_seconds += multiply.getTimeInSeconds();
        for (const auto &z : projected)
          if (!std::isfinite(z.real()) || !std::isfinite(z.imag()))
            throw std::runtime_error("nonfinite projected element");
        double block_check = 0;
        for (int check = 0; check < 8; ++check) {
          int i = check == 0 ? 0 : check == 1 ? n - 1 : (check * 1543 + 17) % n;
          int j = check == 0 ? 0 : check == 1 ? n - 1 : (check * 727 + 13) % n;
          TimeSliceIO<LatticeColorVectorF> slice(vf, sink_time);
          if (eigen.get(KeyTimeSliceColorVec_t(sink_time, i), slice))
            throw std::runtime_error("independent eigenvector read failed");
          v = vf;
          w = zero;
          for (size_t site = 0; site < sites.size(); ++site)
            for (int c = 0; c < 3; ++c) {
              auto z = point(j, src, site * 3 + c, snk);
              w.elem(sites[site]).elem().elem(c) =
                  RComplex<double>(z.real(), z.imag());
            }
          ProjectComplex reference = value(innerProduct(v, w));
          double scale = std::sqrt(toDouble(norm2(v)) * toDouble(norm2(w)));
          double error =
              std::abs(reference - projected[uint64_t(j) * n + i]) / scale;
          if (!(scale > 0 && std::isfinite(error) && error <= 1e-12))
            throw std::runtime_error("independent QDP sink contraction failed");
          block_check = std::max(block_check, error);
        }
        max_check = std::max(max_check, block_check);
        ValPropElementalOperator_t matrix;
        matrix.mat.resize(n, n);
        for (int i = 0; i < n; ++i)
          for (int j = 0; j < n; ++j) {
            auto z = projected[uint64_t(j) * n + i];
            matrix.mat(i, j).elem().elem().elem() =
                RComplex<double>(z.real(), z.imag());
          }
        auto encoded = encodeStorageRecord(matrix);
        hashes[src][snk] = fnv(encoded);
        if (writer.put_record(encodeStorageRecord(key(src, snk, source_time)),
                              encoded))
          throw std::runtime_error("write standard perambulator block failed");
        journal << src << ',' << snk << ',' << n << ','
                << multiply.getTimeInSeconds() << ',' << block_check << ','
                << hashes[src][snk] << '\n';
        journal.flush();
        if (!journal)
          throw std::runtime_error("write projection journal failed");
        std::cout << std::setprecision(17) << "PROJECTED source_spin=" << src
                  << " sink_spin=" << snk << " vectors=" << n
                  << " gemm_seconds=" << multiply.getTimeInSeconds()
                  << " qdp_dot_scaled_error=" << block_check << std::endl;
      }
    eigen.close();
    if (writer.close_records())
      throw std::runtime_error("close perambulator failed");
    basis.clear();
    basis.shrink_to_fit();
    rhs.clear();
    rhs.shrink_to_fit();
    projected.clear();
    projected.shrink_to_fit();
    Store db;
    db.record_options().cache_bytes = 50 * 1024 * 1024;
    if (db.open_records(partial, O_RDONLY, 0))
      throw std::runtime_error("legacy Chroma reader open failed");
    auto keys = db.record_keys();
    if (keys.size() != 16)
      throw std::runtime_error("standard perambulator must have sixteen keys");
    std::string stored_meta;
    db.read_metadata(stored_meta);
    if (stored_meta != metadata_text)
      throw std::runtime_error("metadata roundtrip failed");
    for (int src = 0; src < 4; ++src)
      for (int snk = 0; snk < 4; ++snk) {
        Value block;
        get(db, src, snk, block, n, source_time);
        if (fnv(encodeStorageRecord(block)) != hashes[src][snk])
          throw std::runtime_error("complete standard block checksum failed");
      }
    long double difference = 0, norm = 0;
    double max_abs = 0;
    for (int src = 0; src < 4; ++src)
      for (int snk = src; snk < 4; ++snk) {
        Value a, b;
        get(db, src, snk, a, n, source_time);
        if (src != snk)
          get(db, snk, src, b, n, source_time);
        const auto &am = a.mat;
        const auto &bm = src == snk ? am : b.mat;
        const int sign = (src < 2) == (snk < 2) ? 1 : -1;
        for (int i = 0; i < n; ++i)
          for (int j = 0; j < n; ++j) {
            ProjectComplex x = value(am(i, j)), y = value(bm(j, i));
            double error = std::abs(x - double(sign) * std::conj(y));
            difference += (src == snk ? 1 : 2) * (long double)error * error;
            norm += std::norm(x);
            if (src != snk)
              norm += std::norm(y);
            max_abs = std::max(max_abs, error);
          }
      }
    double hermiticity = std::sqrt(difference / norm);
    if (!(std::isfinite(hermiticity) && hermiticity <= 1e-10))
      throw std::runtime_error(
          "full projected gamma5 Hermiticity check failed");
    db.close_records();
    if (rename(partial.c_str(), output.c_str()))
      throw std::runtime_error("publish perambulator failed");
    total.stop();
    XMLBufferWriter report;
    push(report, "PerambulatorProjection");
    write(report, "status", std::string("pass"));
    write(report, "num_vectors", n);
    write(report, "spin_records", 16);
    write(report, "source_time", source_time);
    write(report, "sink_time", sink_time);
    write(report, "max_QDP_dot_scaled_error", Double(max_check));
    write(report, "full_gamma5_relative_Frobenius", Double(hermiticity));
    write(report, "full_gamma5_max_absolute", Double(max_abs));
    write(report, "gemm_seconds", Double(gemm_seconds));
    write(report, "total_seconds", Double(total.getTimeInSeconds()));
    write(report, "metadata", fragment(metadata_text));
    pop(report);
    std::ofstream(output + ".projection.xml") << report.str();
    std::cout << std::setprecision(17)
              << "PASS projected_perambulator vectors=" << n
              << " records=16 full_gamma5_relative=" << hermiticity
              << " qdp_dot_scaled_error=" << max_check
              << " total_seconds=" << total.getTimeInSeconds() << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "ERROR " << e.what() << std::endl;
    status = 1;
  } catch (const std::string &e) {
    std::cerr << "ERROR " << e << std::endl;
    status = 1;
  }
  QDP_finalize();
  return status;
}
