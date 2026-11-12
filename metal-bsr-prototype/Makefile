CXX := xcrun clang++
.DEFAULT_GOAL := all
CXXFLAGS := -std=c++17 -O3 -mcpu=apple-m1 -ffp-contract=fast -Wall -Wextra -fobjc-arc -fblocks -DACCELERATE_NEW_LAPACK
FRAMEWORKS := -framework Metal -framework Foundation -framework Accelerate
CLANG_FORMAT ?= clang-format
FORMAT_SOURCES := $(wildcard *.hpp *.cc *.mm *.metal)
METAL_SOURCES := kernels.metal solve.metal bicgstab.metal evenodd.metal \
	clover.metal multigrid.metal sparse_coarse.metal compact_clover.metal \
	mg_real_tune.metal packed_clover.metal

.PHONY: format format-check metal-check
format:
	$(CLANG_FORMAT) --style=file -i $(FORMAT_SOURCES)

format-check:
	$(CLANG_FORMAT) --style=file --dry-run --Werror $(FORMAT_SOURCES)

# Compile the shader fragments together, in dependency order, without a GPU.
metal-check: build/standalone-shaders.air

build/standalone-shaders.metal: $(METAL_SOURCES)
	mkdir -p build
	cat $(METAL_SOURCES) > $@

build/standalone-shaders.air: build/standalone-shaders.metal
	xcrun -sdk macosx metal -std=metal3.0 -c $< -o $@

build/metal_three_level: three_level.mm $(wildcard *.hpp) $(wildcard *.metal)
	mkdir -p build
	$(CXX) $(CXXFLAGS) three_level.mm -o $@ $(FRAMEWORKS)

build/compact_probe: compact_probe.mm $(wildcard *.hpp)
	mkdir -p build
	$(CXX) $(CXXFLAGS) compact_probe.mm -o $@ $(FRAMEWORKS)

build/critical_probe: critical_probe.mm sparse_rhs.hpp clover_common.hpp clover_checks.hpp evenodd.hpp wilson_common.hpp bicgstab.hpp solve_common.hpp bsr_common.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) critical_probe.mm -o $@ $(FRAMEWORKS)

.PHONY: all test bench solve-test solve-bench clean
.PHONY: three-test
three-test: build/metal_three_level
	mkdir -p results
	VECLIB_MAXIMUM_THREADS=1 ./build/metal_three_level --self-test --csv results/threelevel-self-test-default.csv
all: build/metal_bsr build/metal_solve build/metal_wilson build/metal_evenodd build/metal_clover build/metal_multigrid build/metal_sparse_multigrid build/metal_solver_compare build/metal_three_level

build/metal_solver_compare: solver_compare.mm sparse_rhs.hpp sparse_coarse.hpp sparse_checks.hpp sparse_coarse.metal multigrid.hpp multigrid.metal fgmres.hpp clover_checks.hpp clover_common.hpp clover.metal evenodd.hpp wilson_common.hpp bicgstab.hpp bicgstab.metal solve_common.hpp bsr_common.hpp kernels.metal
	mkdir -p build
	$(CXX) $(CXXFLAGS) solver_compare.mm -o $@ $(FRAMEWORKS)

.PHONY: solver-compare-test solver-compare-bench
solver-compare-test: build/metal_solver_compare
	VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare --self-test

solver-compare-bench: build/metal_solver_compare
	mkdir -p results
	VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare --csv results/solver-compare16.csv
	python3 solver_compare_summary.py

build/metal_sparse_multigrid: sparse_multigrid.mm sparse_rhs.hpp sparse_coarse.hpp sparse_checks.hpp sparse_coarse.metal multigrid.hpp multigrid.metal fgmres.hpp clover_checks.hpp clover_common.hpp clover.metal evenodd.hpp wilson_common.hpp bicgstab.hpp bicgstab.metal solve_common.hpp bsr_common.hpp kernels.metal
	mkdir -p build
	$(CXX) $(CXXFLAGS) sparse_multigrid.mm -o $@ $(FRAMEWORKS)

.PHONY: sparse-mg-test sparse-mg-bench sparse-mg-study
sparse-mg-test: build/metal_sparse_multigrid
	VECLIB_MAXIMUM_THREADS=1 ./build/metal_sparse_multigrid --self-test

sparse-mg-bench: build/metal_sparse_multigrid
	mkdir -p results
	VECLIB_MAXIMUM_THREADS=1 ./build/metal_sparse_multigrid --tile-rhs 96 --coarse-iters 4 --csv results/sparse-mg16.csv

sparse-mg-study: build/metal_sparse_multigrid build/metal_multigrid build/metal_clover
	python3 run_sparse_study.py
	python3 sparse_summary.py

build/metal_multigrid: multigrid.mm multigrid.hpp fgmres.hpp mg_checks.hpp multigrid.metal clover_checks.hpp clover_common.hpp clover.metal evenodd.hpp wilson_common.hpp bicgstab.hpp bicgstab.metal solve_common.hpp bsr_common.hpp kernels.metal
	mkdir -p build
	$(CXX) $(CXXFLAGS) multigrid.mm -o $@ $(FRAMEWORKS)

.PHONY: mg-test mg-bench mg-tuned-bench mg-tune-screen
mg-test: build/metal_multigrid
	VECLIB_MAXIMUM_THREADS=1 ./build/metal_multigrid --self-test

mg-bench: build/metal_multigrid
	mkdir -p results
	VECLIB_MAXIMUM_THREADS=1 ./build/metal_multigrid --csv results/multigrid.csv

mg-tune-screen: build/metal_multigrid
	python3 mg_tune_sweep.py
	python3 mg_tune_sweep.py --stage2
	python3 mg_tune_sweep.py --stage3

mg-tuned-bench: build/metal_multigrid
	python3 mg_tuned_bench.py
	python3 mg_tuned_summary.py

build/metal_clover: clover.mm clover_checks.hpp clover_common.hpp clover.metal evenodd.hpp wilson_common.hpp bicgstab.hpp bicgstab.metal solve_common.hpp bsr_common.hpp kernels.metal
	mkdir -p build
	$(CXX) $(CXXFLAGS) clover.mm -o $@ $(FRAMEWORKS)

build/clover_stage1: clover_stage1.mm clover_checks.hpp clover_common.hpp clover.metal evenodd.hpp wilson_common.hpp bicgstab.hpp bicgstab.metal solve_common.hpp bsr_common.hpp kernels.metal
	mkdir -p build
	$(CXX) $(CXXFLAGS) clover_stage1.mm -o $@ $(FRAMEWORKS)

.PHONY: clover-test clover-bench
clover-test: build/metal_clover
	./build/metal_clover --self-test

clover-bench: build/metal_clover
	mkdir -p results
	./build/metal_clover --csv results/clover.csv

build/metal_bsr: main.mm bsr_common.hpp kernels.metal
	mkdir -p build
	$(CXX) $(CXXFLAGS) main.mm -o $@ $(FRAMEWORKS)

build/metal_solve: solve.mm solve_common.hpp bsr_common.hpp solve.metal kernels.metal
	mkdir -p build
	$(CXX) $(CXXFLAGS) solve.mm -o $@ $(FRAMEWORKS)

build/metal_wilson: wilson.mm wilson_common.hpp bicgstab.hpp bicgstab.metal solve_common.hpp bsr_common.hpp kernels.metal
	mkdir -p build
	$(CXX) $(CXXFLAGS) wilson.mm -o $@ $(FRAMEWORKS)

build/metal_evenodd: evenodd.mm evenodd.hpp evenodd.metal wilson_common.hpp bicgstab.hpp bicgstab.metal solve_common.hpp bsr_common.hpp kernels.metal
	mkdir -p build
	$(CXX) $(CXXFLAGS) evenodd.mm -o $@ $(FRAMEWORKS)

.PHONY: evenodd-test evenodd-bench
evenodd-test: build/metal_evenodd
	./build/metal_evenodd --self-test

evenodd-bench: build/metal_evenodd
	mkdir -p results
	./build/metal_evenodd --csv results/evenodd.csv

.PHONY: wilson-test wilson-bench
wilson-test: build/metal_wilson
	./build/metal_wilson --self-test

wilson-bench: build/metal_wilson
	mkdir -p results
	./build/metal_wilson --csv results/wilson.csv

solve-test: build/metal_solve
	VECLIB_MAXIMUM_THREADS=1 ./build/metal_solve --self-test

solve-bench: build/metal_solve
	mkdir -p results
	VECLIB_MAXIMUM_THREADS=1 ./build/metal_solve --csv results/solve.csv

test: build/metal_bsr
	VECLIB_MAXIMUM_THREADS=1 ./build/metal_bsr --self-test

bench: build/metal_bsr
	mkdir -p results
	VECLIB_MAXIMUM_THREADS=1 ./build/metal_bsr --csv results/benchmark.csv

clean:
	rm -rf build

# Build the QDP++ bridge with the compiler and flags of the installed Chroma.
CHROMA_CONFIG ?= chroma-config
build/chroma_export: chroma_export.cc chroma_bridge.hpp
	mkdir -p build
	"$$("$(CHROMA_CONFIG)" --cxx)" $$("$(CHROMA_CONFIG)" --cxxflags) $(CHROMA_EXTRA_CPPFLAGS) chroma_export.cc -o $@ $$("$(CHROMA_CONFIG)" --ldflags) $$("$(CHROMA_CONFIG)" --libs)

build/chroma_mg: three_level.mm $(wildcard *.hpp) $(wildcard *.metal)
	mkdir -p build
	$(CXX) $(CXXFLAGS) three_level.mm -o $@ $(FRAMEWORKS)

build/residual_profile: residual_profile.mm $(wildcard *.hpp)
	mkdir -p build
	$(CXX) $(CXXFLAGS) -g residual_profile.mm -o $@ $(FRAMEWORKS)

build/anisotropic_inspect: anisotropic_inspect.cc
	mkdir -p build
	"$$("$(CHROMA_CONFIG)" --cxx)" $$("$(CHROMA_CONFIG)" --cxxflags) $(CHROMA_EXTRA_CPPFLAGS) anisotropic_inspect.cc -o $@ $$("$(CHROMA_CONFIG)" --ldflags) $$("$(CHROMA_CONFIG)" --libs)
