#include <filesystem>
#include <iomanip>
#include <sstream>
#include <sys/resource.h>
#include "production_bridge.hpp"
#include "three_level.hpp"
#include "compact_solve.hpp"

namespace {
constexpr double target=1e-12;
constexpr int workers=8;
struct CountedSchur:GPUCompactSchur {
    using GPUCompactSchur::GPUCompactSchur;
    uint64_t applications=0;
    void encode(id<MTLComputeCommandEncoder> e,id<MTLBuffer> x,id<MTLBuffer> y){
        ++applications;GPUCompactSchur::encode(e,x,y);
    }
};
template<class T> void release(std::vector<T>& v){std::vector<T>().swap(v);}
void source(CompactClover& c,CompactOuter& w,int rhs){
    double scale=1/std::sqrt(double(c.w.a.p.outputs()));
    chunks(c.w.a.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t q=first;q<last;++q){auto z=chroma_bridge::source(q,20261107+rhs);c.w.a.b[q]={scale*z.re,scale*z.im};}});
    compact_rhs_norms(c.w.a,w,workers);
}
struct Input {
    std::ifstream file;production::Header h;std::streampos probe;
    explicit Input(const std::string& prefix):file(prefix+".bridge",std::ios::binary),h(production::read_header(file)){}
    std::vector<M3> links(){
        std::vector<M3> out(size_t(h.volume())*4);std::array<chroma_bridge::Complex,9> row;
        for(auto& m:out){chroma_bridge::read(file,row.data(),9);for(int k=0;k<9;++k)m[k]={row[k].re,row[k].im};}
        probe=file.tellg();return out;
    }
    double compare(const std::vector<D>& x){
        file.clear();file.seekg(probe);std::vector<D> reference(x.size());chroma_bridge::read(file,reference.data(),reference.size());
        double e=0,n=0;for(size_t j=0;j<x.size();++j){e+=std::norm(CD(x[j].re-reference[j].re,x[j].im-reference[j].im));n+=reference[j].re*reference[j].re+reference[j].im*reference[j].im;}
        double rel=std::sqrt(e/n);if(!std::isfinite(rel)||rel>1e-13)throw std::runtime_error("standalone full operator mismatch");return rel;
    }
};
double full_residual(const CompactClover& c,CompactOuter& w){
    compact_fp64::residual(c,w.x.data(),c.w.a.b.data(),w.residual.data(),w.sums.data(),workers);
    double n=0;for(double v:w.sums)n+=v;return std::sqrt(n)/c.w.a.bnorm[0];
}
void check_schur(CompactClover& c,GPUCompactSchur& op){
    std::vector<Cx<float>> x(c.parity.oe.p.outputs()),reference(x.size());
    for(size_t j=0;j<x.size();++j){auto z=chroma_bridge::source(j,20261108);x[j]={float(z.re),float(z.im)};}
    CPUCloverSchur<float,CompactClover> cpu(c,workers);cpu.apply(x.data(),reference.data());
    auto input=[op.m.device newBufferWithBytes:x.data() length:x.size()*8 options:MTLResourceStorageModeShared];
    auto output=[op.m.device newBufferWithLength:x.size()*8 options:MTLResourceStorageModeShared];
    auto cmd=[op.m.queue commandBuffer];auto enc=[cmd computeCommandEncoder];op.encode(enc,input,output);[enc endEncoding];[cmd commit];[cmd waitUntilCompleted];
    if(cmd.status==MTLCommandBufferStatusError)throw std::runtime_error("Schur probe GPU error");
    auto* values=static_cast<Cx<float>*>([output contents]);double e=0,n=0;
    for(size_t j=0;j<x.size();++j){e+=std::norm(CD(values[j].re-reference[j].re,values[j].im-reference[j].im));n+=reference[j].re*reference[j].re+reference[j].im*reference[j].im;}
    double rel=std::sqrt(e/n);if(!std::isfinite(rel)||rel>3e-6)throw std::runtime_error("CPU/Metal Schur convention mismatch");
    std::cout<<"METAL CPU_FP32 Schur_match="<<rel<<std::endl;
}
template<class Correction> void solve(const CompactClover& c,CompactOuter& w,Correction& correction,
                                      const std::string& method,const std::string& output,
                                      std::ostream& trace,std::ostream& summary,int rhs,CountedSchur& op){
    std::fill(w.x.begin(),w.x.end(),D{});w.residual=c.w.a.b;
    const uint64_t before=op.applications;
    auto start=Clock::now();int total=0,outers=0;double rho=1,residual_s=0;
    for(int outer=1;outer<=12;++outer){
        const uint64_t outer_before=op.applications;
        auto info=correction.run(w.residual,.002,8192);total+=info.max_iterations;auto* dx=correction.correction();
        chunks(c.w.a.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t j=first;j<last;++j){w.x[j].re+=double(dx[j].re);w.x[j].im+=double(dx[j].im);}});
        auto check=Clock::now();rho=full_residual(c,w);residual_s+=seconds(check);outers=outer;
        if(!std::isfinite(rho))throw std::runtime_error("nonfinite full FP64 residual");
        trace<<method<<','<<rhs<<','<<outer<<','<<info.max_iterations<<','<<rho<<','<<op.applications-outer_before<<'\n';trace.flush();
        std::cout<<"OUTER method="<<method<<" rhs="<<rhs<<" step="<<outer<<" K="<<info.max_iterations<<" full_relative="<<rho<<std::endl;
        if(rho<=target)break;
    }
    if(rho>target)throw std::runtime_error("full FP64 1e-12 refinement failed");
    double elapsed=seconds(start);std::string path=output+"."+method+".rhs"+std::to_string(rhs)+".bin";
    if(std::filesystem::exists(path))throw std::runtime_error("solution output already exists");
    auto io=Clock::now();std::ofstream saved(path,std::ios::binary);chroma_bridge::write(saved,w.x.data(),w.x.size());saved.close();if(!saved)throw std::runtime_error("solution write failed");
    summary<<method<<','<<rhs<<','<<outers<<','<<total<<','<<elapsed<<','<<residual_s<<','<<rho<<','<<seconds(io)<<','<<op.applications-before<<'\n';summary.flush();
    std::cout<<"SOLVED method="<<method<<" rhs="<<rhs<<" outer="<<outers<<" K="<<total<<" seconds="<<elapsed<<" full_relative="<<rho<<std::endl;
}
}
int main(int argc,char** argv){@autoreleasepool{
    try{
        if(argc!=5)throw std::runtime_error("Usage: performance_compare BRIDGE_PREFIX OUTPUT_PREFIX RHS_COUNT BX,BY,BZ,BT");
        std::string prefix=argv[1],output=argv[2];int count=std::stoi(argv[3]);
        if(count<1||count>192)throw std::runtime_error("invalid RHS count");
        LatticeShape blocks{};
        {std::stringstream values(argv[4]);for(int mu=0;mu<4;++mu){char comma;if(mu&&(!(values>>comma)||comma!=','))throw std::runtime_error("expected four MG extents");if(!(values>>blocks[mu]))throw std::runtime_error("invalid MG extent");}char tail;if(values>>tail)throw std::runtime_error("MG suffix");}
        if(!std::getenv("VECLIB_MAXIMUM_THREADS")||std::string(std::getenv("VECLIB_MAXIMUM_THREADS"))!="1")throw std::runtime_error("set VECLIB_MAXIMUM_THREADS=1");
        if(std::filesystem::exists(output+".metal-solves.csv"))throw std::runtime_error("choose a fresh output prefix");
        Input in(prefix);auto h=in.h;
        for(int mu=0;mu<4;++mu)if(!blocks[mu]||h.shape[mu]%blocks[mu]||(h.shape[mu]/blocks[mu])%2)throw std::runtime_error("MG blocks must divide the lattice and leave even first-coarse extents");
        CompactClover c(h.shape,h.mass,workers,.25,in.links(),h.xi0,h.nu,h.clover_r,h.clover_t);
        release(c.evenA32);
        c.rhs(1,false);c.w.a.b.resize(c.w.a.p.outputs());c.w.a.bnorm.resize(1);CompactOuter outer(c.w.a,workers);outer.known_solution=false;
        std::cout<<std::setprecision(17);
        std::cout<<"MG_SETTINGS blocks="<<blocks[0]<<','<<blocks[1]<<','<<blocks[2]<<','<<blocks[3]<<std::endl;
        for(size_t q=0;q<outer.x.size();++q){auto z=chroma_bridge::source(q,chroma_bridge::probe_rhs);outer.x[q]={z.re,z.im};}
        apply_cpu(c.w.a,outer.x.data(),outer.residual.data(),workers);
        std::cout<<"STANDALONE original FP64 operator_match="<<in.compare(outer.residual)<<std::endl;
        std::string shader;for(auto name:{"kernels.metal","bicgstab.metal","clover.metal","multigrid.metal","sparse_coarse.metal","compact_clover.metal"}){std::ifstream f(name);if(!f)throw std::runtime_error("missing shader");shader.append((std::istreambuf_iterator<char>(f)),{});shader+='\n';}
        MetalRun metal(shader);std::cout<<"DEVICE "<<[[metal.device name] UTF8String]<<std::endl;
        std::ofstream trace(output+".trace.csv"),summary(output+".metal-solves.csv");
        trace<<std::setprecision(17)<<"method,rhs,outer,krylov,full_relative,fine_schur_applications\n";
        summary<<std::setprecision(17)<<"method,rhs,outer,krylov,seconds,residual_seconds,full_relative,save_seconds,fine_schur_applications\n";
        CountedSchur op(metal,c);
        @autoreleasepool{check_schur(c,op);}
        // Subsequent solves/checks use the GPU copies or the FP64 operator.
        // The FP64 basis construction retains its separate even inverse.
        release(c.links32);release(c.oddA32);release(c.inverseEven32);
        {
            GPUBiCGStab bi(metal,c.parity.oe,1,workers);
            bi.encode_override=[&](id<MTLComputeCommandEncoder> e,id<MTLBuffer> x,id<MTLBuffer> y){op.encode(e,x,y);};
            bi.progress=[](int k,double r){std::cout<<"INNER metal-bicgstab K="<<k<<" relative_FP32="<<r<<std::endl;};
            CompactDeviceCorrection<decltype(bi)> correction(c,op,bi,workers);
            for(int rhs=0;rhs<count;++rhs){source(c,outer,rhs);solve(c,outer,correction,"metal-bicgstab",output,trace,summary,rhs,op);}
        }
        auto start=Clock::now();const uint32_t vectors=12,tile=2;LatticeShape coarse_shape;
        for(int mu=0;mu<4;++mu){if(h.shape[mu]%blocks[mu])throw std::runtime_error("MG block shape mismatch");coarse_shape[mu]=h.shape[mu]/blocks[mu];}
        c.rhs(tile,false);
        Transfer first(c.parity.oe.p.rows,12,vectors,h.shape,blocks,true,workers,&c.parity.odd);
        {CPUCloverSchur<double,CompactClover> schur(c,workers);
         first.build([&](const D* x,D* y){schur.apply(x,y);},[&](const D* x,D* y){local_blocks(c.inverseOdd32,x,y,first.rows,tile,workers);},8,.65,tile);}
        BlockMatrix a1(2*vectors,first.aggregates,workers);assemble_schur(a1,c,first);release(c.inverseEven);
        Transfer second(a1.rows,a1.bs,vectors,coarse_shape,hypercube(2),false,workers);
        second.build([&](const D* x,D* y){a1.apply(x,y,vectors);},[&](const D* x,D* y){a1.apply(x,y,vectors,true);},8);
        BlockMatrix a2(2*vectors,second.aggregates,workers);project_matrix(a2,a1,second);
        first.release64();second.release64();a1.release64();a2.release64();
        c.rhs(1,false);
        std::cout<<"HIERARCHY seconds="<<seconds(start)<<" coarse1="<<a1.n<<" coarse2="<<a2.n<<std::endl;
        GPUCycle pre(first,metal,1,c.inverseOdd32);release(first.basis32);release(c.inverseOdd32);
        pre.product=[&](id<MTLComputeCommandEncoder> e,id<MTLBuffer> x,id<MTLBuffer> y){op.encode(e,x,y);};
        GPUBlockCoarseSolve<BlockMatrix> matrix1(a1,metal,1,1,.1),bottom(a2,metal,1,16,.1);
        GPUBlockOperator op1(matrix1);GPUCycle pre1(second,metal,1,a1.di32);release(second.basis32);
        pre1.product=[&](id<MTLComputeCommandEncoder> e,id<MTLBuffer> x,id<MTLBuffer> y){op1.encode(e,x,y);};
        pre1.coarse=[&](id<MTLComputeCommandEncoder> e,id<MTLBuffer> x,id<MTLBuffer> y){bottom.encode(e,x,y);};pre1.check=[&](){bottom.check();};
        GPUGenericFGMRES<GPUBlockOperator,GPUCycle> coarse(a1.n,1,op1,pre1,8);coarse.allow_inexact=true;
        pre.blocking_coarse=[&](id<MTLBuffer> x,id<MTLBuffer> y){coarse.run_device(x,.1,16);pre.copy(coarse.delta,y,a1.n);};
        release(a1.a32);release(a1.di32);release(a2.a32);release(a2.di32);
        GPUGenericFGMRES<CountedSchur,GPUCycle> mg(first.rows*12,1,op,pre,8);
        mg.blocking_pre=[&](id<MTLBuffer> x,id<MTLBuffer> y){pre.apply_device(x,y);};
        mg.progress=[](int k,double r){std::cout<<"INNER metal-mg K="<<k<<" relative_FP32="<<r<<std::endl;};
        CompactDeviceCorrection<decltype(mg)> correction(c,op,mg,workers);
        std::cout<<"MG_READY setup_seconds="<<seconds(start)<<" allocated_GiB="<<double(metal.device.currentAllocatedSize)/(1ull<<30)<<std::endl;
        const size_t primary_host_bytes=c.w.links.size()*sizeof(M3)+(c.evenA.size()+c.oddA.size()+c.inverseEven.size()+c.w.a.b.size()+outer.x.size()+outer.residual.size())*sizeof(D);
        std::cout<<"MEMORY primary_host_bytes="<<primary_host_bytes<<" metal_allocated_bytes="<<metal.device.currentAllocatedSize<<" primary_combined_GiB="<<double(primary_host_bytes+metal.device.currentAllocatedSize)/(1ull<<30)<<std::endl;
        for(int rhs=0;rhs<count;++rhs){source(c,outer,rhs);solve(c,outer,correction,"metal-mg",output,trace,summary,rhs,op);}
        rusage usage{};getrusage(RUSAGE_SELF,&usage);std::cout<<"PEAK_RSS_GiB="<<double(usage.ru_maxrss)/(1ull<<30)<<std::endl;
        std::cout<<"METAL PASS all full FP64 residuals <= 1e-12; native Chroma verification still required"<<std::endl;
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
    return 0;
}}
