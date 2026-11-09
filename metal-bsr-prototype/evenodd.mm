#include "evenodd.hpp"

void validate_evenodd(const WilsonSystem& w,const EvenOdd& e,MetalRun& metal,int workers){
    uint32_t rhs=w.a.p.rhs,width=12*rhs;double inv=1/e.diagonal,inv2=inv*inv;
    std::vector<D> xe(e.eo.p.outputs()),xfull(w.a.p.outputs()),oddonly(xfull.size());
    for(size_t i=0;i<e.even.size();++i)for(uint32_t k=0;k<width;++k){
        xe[i*width+k]=w.a.truth[size_t(e.even[i])*width+k];xfull[size_t(e.even[i])*width+k]=xe[i*width+k];
    }
    auto first=direct_wilson(w,xfull);
    for(size_t i=0;i<e.odd.size();++i)for(uint32_t k=0;k<width;++k)oddonly[size_t(e.odd[i])*width+k]=first[size_t(e.odd[i])*width+k];
    auto second=direct_wilson(w,oddonly);
    std::vector<D> expected(xe.size()),got(xe.size());
    for(size_t i=0;i<e.even.size();++i)for(uint32_t k=0;k<width;++k){size_t j=i*width+k;auto z=second[size_t(e.even[i])*width+k];
        expected[j]={xe[j].re-inv2*z.re,xe[j].im-inv2*z.im};}
    CPUSchur<double> cpu(e,workers);cpu.apply(xe.data(),got.data());double err=0;
    for(size_t j=0;j<got.size();++j)err=std::max(err,std::hypot(got[j].re-expected[j].re,got[j].im-expected[j].im));
    if(err>2e-13)throw std::runtime_error("Schur versus independent direct Wilson mismatch");
    // Reduction of b must equal S*x_even for the known solution.
    std::fill(oddonly.begin(),oddonly.end(),D{});
    for(size_t i=0;i<e.odd.size();++i)for(uint32_t k=0;k<width;++k)oddonly[size_t(e.odd[i])*width+k]=w.a.b[size_t(e.odd[i])*width+k];
    auto hbo=direct_wilson(w,oddonly);err=0;
    for(size_t i=0;i<e.even.size();++i)for(uint32_t k=0;k<width;++k){size_t f=size_t(e.even[i])*width+k,j=i*width+k;
        D ge{(w.a.b[f].re-inv*hbo[f].re)*inv,(w.a.b[f].im-inv*hbo[f].im)*inv};
        err=std::max(err,std::hypot(ge.re-expected[j].re,ge.im-expected[j].im));}
    if(err>2e-13)throw std::runtime_error("Schur RHS reduction mismatch");
    err=0;
    for(size_t i=0;i<e.odd.size();++i)for(uint32_t k=0;k<width;++k){size_t f=size_t(e.odd[i])*width+k;
        D rec{(w.a.b[f].re-first[f].re)*inv,(w.a.b[f].im-first[f].im)*inv};
        err=std::max(err,std::hypot(rec.re-w.a.truth[f].re,rec.im-w.a.truth[f].im));}
    if(err>2e-13)throw std::runtime_error("odd reconstruction mismatch");
    // Test Metal and CPU FP32 on the same rounded input; use direct Wilson as oracle.
    for(size_t i=0;i<e.even.size();++i)for(uint32_t k=0;k<width;++k){size_t j=i*width+k;
        xe[j]={double(float(xe[j].re)),double(float(xe[j].im))};xfull[size_t(e.even[i])*width+k]=xe[j];}
    first=direct_wilson(w,xfull);std::fill(oddonly.begin(),oddonly.end(),D{});
    for(size_t i=0;i<e.odd.size();++i)for(uint32_t k=0;k<width;++k)oddonly[size_t(e.odd[i])*width+k]=first[size_t(e.odd[i])*width+k];
    second=direct_wilson(w,oddonly);
    GPUSchur gpu(metal,e);auto in=gpu.eo.buffers[3],out=gpu.eo.buffers[4];
    auto* x=static_cast<Cx<float>*>([in contents]);auto* y=static_cast<Cx<float>*>([out contents]);
    std::vector<Cx<float>> xf(xe.size()),yf(xe.size());
    for(size_t j=0;j<xe.size();++j)x[j]=xf[j]={float(xe[j].re),float(xe[j].im)};
    auto cmd=[metal.queue commandBuffer];auto enc=[cmd computeCommandEncoder];gpu.encode(enc,in,out);[enc endEncoding];[cmd commit];[cmd waitUntilCompleted];
    if(cmd.status==MTLCommandBufferStatusError)throw std::runtime_error([[cmd.error localizedDescription] UTF8String]);
    CPUSchur<float> cpuf(e,workers);cpuf.apply(xf.data(),yf.data());err=0;
    for(size_t i=0;i<e.even.size();++i)for(uint32_t k=0;k<width;++k){size_t j=i*width+k;auto z=second[size_t(e.even[i])*width+k];
        D ref{xe[j].re-inv2*z.re,xe[j].im-inv2*z.im};
        err=std::max({err,std::hypot(y[j].re-ref.re,y[j].im-ref.im),std::hypot(yf[j].re-ref.re,yf[j].im-ref.im)});}
    if(err>2e-5)throw std::runtime_error("CPU/Metal FP32 Schur mismatch");
    std::cout<<"EO direct-composition/RHS/reconstruction/CPU/Metal checks PASS L="<<w.extent<<" RHS="<<rhs<<" gauge="<<w.strength<<" APBC="<<w.antiperiodic<<std::endl;
}

int main(int argc,char** argv){
    @autoreleasepool{try{
        std::vector<uint32_t> sizes{8},rhs_values{96,192};std::vector<double> masses{.1,.5};
        double strength=.25,tol=1e-10,inner_tol=.002;int workers=8,chunk=4,samples=3,limit=512;
        bool self=false,apbc=true;std::string csv_path,selection="both";
        for(int j=1;j<argc;++j){std::string opt=argv[j];auto next=[&](){if(++j>=argc)throw std::runtime_error("missing value");return std::string(argv[j]);};
            if(opt=="--sizes")sizes=list_arg(next(),12);else if(opt=="--rhs")rhs_values=list_arg(next(),384);
            else if(opt=="--mass"){masses.clear();std::stringstream in(next());std::string s;while(std::getline(in,s,','))masses.push_back(std::stod(s));}
            else if(opt=="--gauge-strength")strength=std::stod(next());else if(opt=="--tol")tol=std::stod(next());
            else if(opt=="--inner-tol")inner_tol=std::stod(next());else if(opt=="--threads")workers=std::stoi(next());
            else if(opt=="--chunk")chunk=std::stoi(next());else if(opt=="--samples")samples=std::stoi(next());
            else if(opt=="--csv")csv_path=next();else if(opt=="--preconditioner")selection=next();
            else if(opt=="--self-test")self=true;else if(opt=="--periodic-time")apbc=false;
            else if(opt=="--help"){std::cout<<"--sizes 8 --rhs 96,192 --mass .1,.5 --samples 3 --threads 8 --chunk 4 --preconditioner both|none|even_odd --gauge-strength .25 --tol 1e-10 --inner-tol .002 --periodic-time --csv results/evenodd.csv --self-test\n";return 0;}
            else throw std::runtime_error("unknown option "+opt);
        }
        if(masses.empty()||!std::all_of(masses.begin(),masses.end(),[](double m){return m>0&&std::isfinite(m);})||!(strength>=0&&strength<=1)||
           !(tol>=1e-13&&tol<1e-3)||!(inner_tol>tol&&inner_tol<.2)||workers<1||workers>32||chunk<1||chunk>32||samples<1||samples>15||samples%2==0||
           (selection!="both"&&selection!="none"&&selection!="even_odd"))throw std::runtime_error("invalid options");
        if(self){sizes={2,4};rhs_values={1,3,7};masses={.25};samples=1;selection="both";}
        for(auto l:sizes){if(l%2)throw std::runtime_error("even/odd requires even periodic lattice extents");
            for(auto r:rhs_values)if(uint64_t(l)*l*l*l*12*r>18000000)throw std::runtime_error("working-set limit exceeded");}
        std::string source;for(const char* name:{"kernels.metal","bicgstab.metal","evenodd.metal"}){
            std::ifstream f(name);if(!f)throw std::runtime_error("run from prototype directory");source.append((std::istreambuf_iterator<char>(f)),{});source+='\n';}
        auto start=Clock::now();MetalRun metal(source);double library_ms=seconds(start)*1000;
        std::cout<<"GPU: "<<[[metal.device name] UTF8String]<<" OS: "<<[[NSProcessInfo processInfo].operatingSystemVersionString UTF8String]
                 <<"\nWilson EO BiCGStab workers="<<workers<<" library_ms="<<library_ms<<" gauge="<<strength<<" APBC="<<apbc<<" tol="<<tol<<" chunk="<<chunk<<std::endl;
        if(self){bool rejected=false;try{auto w=make_wilson(3,1,.25,strength,apbc,workers);EvenOdd invalid(w);}catch(const std::runtime_error&){rejected=true;}
            if(!rejected)throw std::runtime_error("odd periodic extent was not rejected");}
        std::ofstream csv,trace;
        if(!csv_path.empty()){csv.open(csv_path);trace.open(csv_path+".trace.csv");if(!csv||!trace)throw std::runtime_error("cannot create results");
            csv<<"extent,rows,rhs,mass,gauge_strength,antiperiodic,preconditioner,backend,trial,workers,chunk,tolerance,inner_tolerance,total_ms,inner_ms,residual_ms,update_ms,eo_prepare_ms,eo_reconstruct_ms,setup_ms,library_ms,outer_steps,inner_max_sum,dispatched_steps,max_true_residual,solution_error,ms_per_rhs,rhs_per_second\n";
            trace<<"extent,rhs,mass,preconditioner,backend,trial,outer_step,max_true_residual,inner_min_iterations,inner_max_iterations,dispatched_steps\n";}
        int cases=0;
        for(auto mass:masses)for(auto l:sizes)for(auto rhs:rhs_values){
            auto w=make_wilson(l,rhs,mass,strength,apbc,workers,self);auto& a=w.a;
            start=Clock::now();EvenOdd eo(w);double eo_setup=seconds(start)*1000;
            if(self){validate_evenodd(w,eo,metal,workers);auto free=make_wilson(l,rhs,mass,0,!apbc,workers,self);EvenOdd fe(free);validate_evenodd(free,fe,metal,workers);}
            start=Clock::now();CPUBiCGStab<double> plain64(a,workers);double sd=seconds(start)*1000;
            start=Clock::now();CPUBiCGStab<float> plain32(a,workers);double sf=seconds(start)*1000;
            start=Clock::now();GPUBiCGStab plain_gpu(metal,a,chunk,workers);double sg=seconds(start)*1000;
            start=Clock::now();CPUSchur<double> cs64(eo,workers);CPUBiCGStab<double> bi64(eo.eo,workers);
            bi64.apply_override=[&](const D* x,D* y){cs64.apply(x,y);};EOCorrection<CPUBiCGStab<double>> e64(eo,bi64,workers);double ed=eo_setup+seconds(start)*1000;
            start=Clock::now();CPUSchur<float> cs32(eo,workers);CPUBiCGStab<float> bi32(eo.eo,workers);
            bi32.apply_override=[&](const Cx<float>* x,Cx<float>* y){cs32.apply(x,y);};EOCorrection<CPUBiCGStab<float>> e32(eo,bi32,workers);double ef=eo_setup+seconds(start)*1000;
            start=Clock::now();GPUSchur gs(metal,eo);GPUBiCGStab bigpu(metal,eo.eo,chunk,workers);
            bigpu.encode_override=[&](id<MTLComputeCommandEncoder> enc,id<MTLBuffer> x,id<MTLBuffer> y){gs.encode(enc,x,y);};
            EOCorrection<GPUBiCGStab> egpu(eo,bigpu,workers);double es=eo_setup+seconds(start)*1000;
            start=Clock::now();OuterWorkspace outer(a,workers);double common_setup=seconds(start)*1000;
            double setup[]={sd,sf,sg,ed,ef,es};std::string names[]={"cpu_fp64","cpu_mixed","metal_mixed"};
            auto run=[&](int b){
                if(b==0)return solve(a,plain64,outer,workers,tol,tol,limit);
                if(b==1)return solve(a,plain32,outer,workers,tol,inner_tol,limit);
                if(b==2)return solve(a,plain_gpu,outer,workers,tol,inner_tol,limit);
                if(b==3){e64.reset_metrics();return solve(a,e64,outer,workers,tol,tol,limit);}
                if(b==4){e32.reset_metrics();return solve(a,e32,outer,workers,tol,inner_tol,limit);}
                egpu.reset_metrics();return solve(a,egpu,outer,workers,tol,inner_tol,limit);
            };
            std::vector<int> backends;for(int b=0;b<6;++b)if(selection=="both"||(selection=="none"&&b<3)||(selection=="even_odd"&&b>=3))backends.push_back(b);
            if(!self)for(size_t b=0;b<backends.size();++b)run(backends[(b+cases)%backends.size()]);
            for(int trial=0;trial<samples;++trial)for(size_t b=0;b<backends.size();++b){int method=backends[(b+trial+cases)%backends.size()];auto r=run(method);
                bool pre=method>=3;int backend=method%3;double prep=method==3?e64.prepare_ms:method==4?e32.prepare_ms:method==5?egpu.prepare_ms:0;
                double recon=method==3?e64.reconstruct_ms:method==4?e32.reconstruct_ms:method==5?egpu.reconstruct_ms:0;std::string prec=pre?"even_odd":"none";
                std::cout<<"mass="<<mass<<" L="<<l<<" RHS="<<rhs<<" "<<prec<<" "<<names[backend]<<" trial="<<trial<<" ms="<<r.total_ms
                    <<" outer="<<r.outer<<" inner="<<r.inner_max_sum<<" dispatched="<<r.dispatched<<" true_r="<<r.max_residual<<" x_error="<<r.solution_error<<" PASS"<<std::endl;
                if(csv){csv<<std::setprecision(12)<<l<<','<<a.p.rows<<','<<rhs<<','<<mass<<','<<strength<<','<<apbc<<','<<prec<<','<<names[backend]<<','<<trial<<','
                    <<workers<<','<<(backend==2?chunk:1)<<','<<tol<<','<<(backend==0?tol:inner_tol)<<','<<r.total_ms<<','<<r.inner_ms<<','<<r.residual_ms<<','
                    <<r.update_ms<<','<<prep<<','<<recon<<','<<setup[method]+common_setup<<','<<library_ms<<','<<r.outer<<','<<r.inner_max_sum<<','<<r.dispatched<<','
                    <<r.max_residual<<','<<r.solution_error<<','<<r.total_ms/rhs<<','<<1000*rhs/r.total_ms<<'\n';csv.flush();
                    for(auto tr:r.trace){trace<<std::setprecision(12)<<l<<','<<rhs<<','<<mass<<','<<prec<<','<<names[backend]<<','<<trial<<','<<tr.outer<<','<<tr.max_residual
                        <<','<<tr.inner.min_iterations<<','<<tr.inner.max_iterations<<','<<tr.inner.dispatched<<'\n';}trace.flush();}
            }
            ++cases;
        }
        std::cout<<"PASS: "<<cases<<" systems; every RHS checked on the full FP64 Wilson equation after reconstruction."<<std::endl;return 0;
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}}
}
