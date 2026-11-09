#include "solve_common.hpp"

int main(int argc,char** argv) {
    @autoreleasepool {try {
        std::vector<uint32_t> sizes{4,8},rhs_values{1,12,48,96,192,384};
        double mass=0.25,tol=1e-10,inner_tol=2e-3;int samples=3,workers=8,chunk=8,limit=256;
        bool self=false;std::string csv_path;
        for(int j=1;j<argc;++j){std::string opt=argv[j];auto next=[&](){if(++j>=argc)throw std::runtime_error("missing option value");return std::string(argv[j]);};
            if(opt=="--sizes")sizes=list_arg(next(),16);else if(opt=="--rhs")rhs_values=list_arg(next(),2048);
            else if(opt=="--samples")samples=std::stoi(next());else if(opt=="--threads")workers=std::stoi(next());
            else if(opt=="--chunk")chunk=std::stoi(next());else if(opt=="--mass")mass=std::stod(next());
            else if(opt=="--tol")tol=std::stod(next());else if(opt=="--inner-tol")inner_tol=std::stod(next());
            else if(opt=="--csv")csv_path=next();else if(opt=="--self-test")self=true;
            else if(opt=="--help"){std::cout<<"Run from prototype directory. --sizes 4,8 --rhs 1,12,48,96,192,384 --samples 3 --threads 8 --chunk 8 --mass 0.25 --tol 1e-10 --inner-tol 0.002 --csv results/solve.csv --self-test\n";return 0;}
            else throw std::runtime_error("unknown option "+opt);
        }
        if(!(mass>0)||!std::isfinite(mass)||!(tol>=1e-13&&tol<1e-3)||!(inner_tol>tol&&inner_tol<0.2)||
           samples<1||samples%2==0||samples>15||workers<1||workers>32||chunk<1||chunk>32)
            throw std::runtime_error("invalid numerical or timing options");
        if(self){sizes={1,2,3};rhs_values={1,3,7};samples=1;}
        for(auto l:sizes)for(auto r:rhs_values)if(uint64_t(l)*l*l*l*12*r>24000000)
            throw std::runtime_error("prototype working-set limit exceeded; use smaller lattice/RHS batch");
        std::ifstream k("kernels.metal"),s("solve.metal");if(!k||!s)throw std::runtime_error("run from prototype directory (shader files missing)");
        std::string source((std::istreambuf_iterator<char>(k)),{});source+='\n';source.append((std::istreambuf_iterator<char>(s)),{});
        auto init=Clock::now();MetalRun metal(source);double library_ms=seconds(init)*1000;
        std::cout<<"GPU: "<<[[metal.device name] UTF8String]<<" | macOS "<<[[NSProcessInfo processInfo].operatingSystemVersionString UTF8String]
                 <<"\nCPU workers="<<workers<<" library_ms="<<library_ms<<" mass="<<mass<<" true_tol="<<tol<<" inner_tol="<<inner_tol<<" GPU_chunk="<<chunk<<std::endl;
        std::ofstream csv,trace;
        if(!csv_path.empty()){csv.open(csv_path);trace.open(csv_path+".trace.csv");if(!csv||!trace)throw std::runtime_error("cannot create results");
            csv<<"extent,rows,rhs,backend,trial,mass,tolerance,inner_tolerance,workers,chunk,total_ms,inner_ms,residual_ms,update_ms,setup_ms,library_ms,outer_steps,inner_max_sum,dispatched_steps,max_true_residual,solution_error,ms_per_rhs,rhs_per_second\n";
            trace<<"extent,rhs,backend,trial,outer_step,max_true_residual,inner_min_iterations,inner_max_iterations,dispatched_steps\n";}
        int cases=0;
        for(auto l:sizes)for(auto rhs:rhs_values){
            System a=make_system(l,rhs,mass,workers,self);
            if(self)check_spd_structure(a);
            auto t=Clock::now();CPUCG<float> cf(a,workers);double sf=seconds(t)*1000;
            t=Clock::now();CPUCG<double> cd(a,workers);double sd=seconds(t)*1000;
            t=Clock::now();GPUCG gpu(metal,a,chunk,workers);double sg=seconds(t)*1000;
            t=Clock::now();OuterWorkspace w(a,workers);double common_setup=seconds(t)*1000;
            auto run=[&](int backend){
                if(backend==0)return solve(a,cd,w,workers,tol,tol,limit);
                if(backend==1)return solve(a,cf,w,workers,tol,inner_tol,limit);
                return solve(a,gpu,w,workers,tol,inner_tol,limit);
            };
            const std::string names[]={"cpu_fp64","cpu_mixed","metal_mixed"};double setup[]={sd+common_setup,sf+common_setup,sg+common_setup};
            if(!self)for(int b=0;b<3;++b)run((b+cases)%3); // full-solve warmup, not reported
            std::array<std::vector<SolveResult>,3> results;
            for(int trial=0;trial<samples;++trial)for(int b=0;b<3;++b){int backend=(b+trial+cases)%3;auto r=run(backend);results[backend].push_back(r);
                std::cout<<"L="<<l<<" RHS="<<rhs<<" "<<names[backend]<<" trial="<<trial<<" total_ms="<<r.total_ms
                         <<" outer="<<r.outer<<" inner="<<r.inner_max_sum<<" dispatched="<<r.dispatched<<" true_r="<<r.max_residual
                         <<" x_error="<<r.solution_error<<" PASS"<<std::endl;
                if(csv){csv<<std::setprecision(12)<<l<<','<<a.p.rows<<','<<rhs<<','<<names[backend]<<','<<trial<<','<<mass<<','<<tol<<','
                    <<(backend==0?tol:inner_tol)<<','<<workers<<','<<(backend==2?chunk:1)<<','<<r.total_ms<<','<<r.inner_ms<<','<<r.residual_ms<<','
                    <<r.update_ms<<','<<setup[backend]<<','<<library_ms<<','<<r.outer<<','<<r.inner_max_sum<<','<<r.dispatched<<','<<r.max_residual<<','
                    <<r.solution_error<<','<<r.total_ms/rhs<<','<<rhs*1000/r.total_ms<<'\n';csv.flush();
                    for(auto& tr:r.trace) {
                        trace<<std::setprecision(12)<<l<<','<<rhs<<','<<names[backend]<<','<<trial<<','<<tr.outer<<','
                            <<tr.max_residual<<','<<tr.inner.min_iterations<<','<<tr.inner.max_iterations<<','<<tr.inner.dispatched<<'\n';
                    }
                    trace.flush();}
            }
            ++cases;
        }
        std::cout<<"PASS: "<<cases<<" systems; every RHS has full CPU FP64 true-residual and known-solution checks."<<std::endl;
        return 0;
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}}
}
