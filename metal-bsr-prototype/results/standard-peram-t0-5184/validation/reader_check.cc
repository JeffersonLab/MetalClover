// Check the independent analysis-side reader, including spin/index ordering.
#define USE_SUPERBBLAS
#define USE_SUPERBBLAS_WITH_GPU_SUPPORT
#include <datalib/superb/prop_superb.h>
#include <fstream>
#include <iomanip>
#include <mpi.h>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int status = 0;
  try {
    if (argc != 4)
      throw std::runtime_error(
          "Usage: reader_check SDB NUM_VECS CSV_OR_KEYS_ONLY");
    const int n = std::stoi(argv[2]);
    Datalib::TensorStorage::StorageProp4 reader(false);
    if (reader.open({argv[1]}, n))
      throw std::runtime_error("analysis reader open failed");
    std::vector<Hadron::KeyProp4ElementalOperator_t> keys;
    reader.keys(keys);
    if (keys.size() != 1 || keys[0].t_source != 0 || keys[0].t_slice != 0 ||
        keys[0].mass_label != "light")
      throw std::runtime_error("analysis reader keys mismatch");
    if (std::string(argv[3]) != "KEYS_ONLY") {
      Hadron::ValProp4ElementalOperator_t val;
      if (reader.get(keys[0], val))
        throw std::runtime_error("analysis reader value failed");
      std::ofstream csv(argv[3]);
      csv << "source_spin,sink_spin,sink_vector,source_vector,real,imag\n"
          << std::setprecision(17);
      for (int src = 0; src < 4; ++src)
        for (int snk = 0; snk < 4; ++snk)
          for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) {
              auto z = val.op(i, j, snk, src);
              csv << src << ',' << snk << ',' << i << ',' << j << ','
                  << z.real() << ',' << z.imag() << '\n';
            }
      if (!csv)
        throw std::runtime_error("write analysis-reader comparison failed");
    }
    std::cout << "PASS analysis_reader t_source=0 t_sink=0 mass=light vectors="
              << n << std::endl;
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    status = 1;
  }
  MPI_Finalize();
  return status;
}
