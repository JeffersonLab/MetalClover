from pathlib import Path
import subprocess
repo=Path('/Users/admin/Downloads/chroma.metalclover')
base='b3a380ea6'
names=['lib/meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.cc','lib/meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.h','lib/util/ferm/mgproton.cc','lib/util/ferm/mgproton.h']
for name in names:
 raw=subprocess.run(['git','show',base+':'+name],cwd=repo,check=True,capture_output=True).stdout
 (repo/name).write_bytes(raw)
# Preserve useful accounting without changing Eloy's spin-major scheduling.
p=repo/'lib/util/ferm/mgproton.cc';s=p.read_text()
a='''      QDPIO::cout << "Time to compute inversions for " << num_spins << " spin sources and "
		  << num_vecs << " colorvecs : " << snarss1.getTimeInSeconds() << " secs"
		  << std::endl;'''
b='''      QDPIO::cout << "Time to compute inversions for " << num_spins
                  << " spin sources and " << num_vecs << " colorvecs ("
                  << size_t(num_spins) * num_vecs
                  << " physical RHS) : " << snarss1.getTimeInSeconds()
                  << " secs" << std::endl;'''
assert a in s;s=s.replace(a,b,1);p.write_text(s)
p=repo/'mainprogs/tests/t_metal_clover.cc';s=p.read_text()
a=s.index('      const SB::DistillationRHSPlan rounded(2, 7);')
b=s.index('        for (int col = 0; col < vectors; ++col)',a)
s=s[:a]+'''      for (int vectors : {1, 3}) {
        const int max_rhs = vectors == 1 ? 1 : 2;
        recording->widths.clear();
        auto input_colorvec = colorvec.kvslice_from_size({}, {{'n', vectors}});
        auto spin_solutions = SB::doInversion(
            chimera, input_colorvec, source_time, 0, dims[3], {0, 1, 2, 3},
            max_rhs, "cxyzXnSst");
        SB::Tensor<5, SB::Complex> projected(
            "NnSst", SB::Coor<5>{vectors, vectors, Ns, Ns, 1}, SB::OnHost,
            SB::OnEveryoneReplicated);
        projected.contract(
            input_colorvec, {{'n', 'N'}}, SB::Conjugate,
            spin_solutions.kvslice_from_size({{'t', source_time}}, {{'t', 1}}),
            {}, SB::NotConjugate);
        const std::vector<size_t> expected_widths = vectors == 1
            ? std::vector<size_t>{1, 1, 1, 1}
            : std::vector<size_t>{2, 1, 2, 1, 2, 1, 2, 1};
        if (recording->widths != expected_widths)
          throw std::runtime_error("Eloy spin-major qprop batching failed");
        std::cout << "ELOY_BATCH solver=" << name << " num_vecs=" << vectors
                  << " source_spins=" << Ns << " physical_rhs=" << Ns * vectors
                  << " max_rhs=" << max_rhs << " qprop_calls=" << recording->widths.size()
                  << " widths=";
        for (size_t i = 0; i < recording->widths.size(); ++i)
          std::cout << (i ? "," : "") << recording->widths[i];
        std::cout << std::endl;
''' + s[b:]
s=s.replace('"four-spin full Chroma FP64 residual failed"','"spin-major full Chroma FP64 residual failed"',1)
s=s.replace('"four-spin perambulator index/contraction failed"','"spin-major perambulator index/contraction failed"',1)
s=s.replace('"DISTILLATION_SPIN solver="','"ELOY_SPIN solver="',1)
s=s.replace('<< " physical_rhs=" << plan.total_rhs','<< " physical_rhs=" << Ns * vectors',1)
s=s.replace('<< " effective_max_rhs=" << plan.max_rhs','<< " max_rhs=" << max_rhs',1)
s=s.replace('"full/odd solves, batched qprop RHS tiles=3+2, four-spin "\n           "distillation, "',
            '"full/odd solves, batched qprop RHS tiles=3+2, Eloy spin-major "\n           "distillation max_rhs=1/2, "',1)
p.write_text(s)
p=repo/'mainprogs/tests/t_metal_clover_batch.cc';s=p.read_text()
a=s.index('        for (int i = 0; i < count; i += tile) {')
b=s.index('        watch.stop();',a)
s=s[:a]+'''        // Match Eloy's production measurement: batch eigenvectors within
        // each source spin. max_rhs=1 remains a strict one-field caller cap.
        for (int spin = 0; spin < Ns; ++spin)
          for (int first = 0; first < count / Ns; first += tile) {
            int end = std::min(count / Ns, first + tile);
            std::vector<std::shared_ptr<T>> out;
            std::vector<std::shared_ptr<const T>> in;
            for (int vector = first; vector < end; ++vector) {
              size_t column = size_t(vector) * Ns + spin;
              out.push_back(solutions[column]);
              in.push_back(sources[column]);
            }
            auto result = (*prop)(out, in);
            if (result.size() != size_t(end - first))
              throw std::runtime_error("array result count mismatch");
            for (int vector = first; vector < end; ++vector)
              results[size_t(vector) * Ns + spin] = result[vector - first];
          }
''' + s[b:]
s=s.replace(' seed_base=20261113 spin_fast=true',' seed_base=20261113 spin_fast_index=true submission_order=spin_major',1)
s=s.replace('"PASS real 12^3x96 four-spin distillation qprop "',
            '"PASS real 12^3x96 Eloy spin-major distillation qprop "',1)
p.write_text(s)
for name in ['mainprogs/tests/t_metal_clover.cc','mainprogs/tests/t_metal_clover_batch.cc']:
 subprocess.run(['/opt/homebrew/bin/clang-format','--style=file','-i',str(repo/name)],check=True)
print('Restored the original measurement and spin-major helper; kept Metal backend and bounded temporary pools.')
