#pragma once
#include "compact_residual.hpp"
#include "sparse_rhs.hpp"

// Three full FP64 fields: b (in the system), x, and r. The old residual is
// overwritten by M64*x after each correction. Manufactured truth is regenerated.
struct CompactOuter {
    std::vector<D> x,residual;std::vector<double> sums,relative;uint32_t offset=0;bool zero=false,known_solution=true;
    std::function<void(int,double,int)> observe;
    CompactOuter(const System& a,int workers):x(a.p.outputs()),residual(x.size()),sums(size_t(workers)*a.p.rhs),relative(a.p.rhs){}
};
inline D compact_truth(size_t q,uint32_t rhs){uint64_t seed=rhs_hash(q)^rhs_hash(uint64_t(rhs)+771337);auto uniform=[](uint64_t v){return 2*double(v>>11)*0x1p-53-1;};return {uniform(rhs_hash(seed)),uniform(rhs_hash(seed+1))};}
inline void compact_rhs_norms(System& a,CompactOuter& outer,int workers){uint32_t rhs=a.p.rhs;std::fill(outer.sums.begin(),outer.sums.end(),0.);
    chunks(a.p.rows*12,workers,[&](int worker,uint32_t first,uint32_t last){auto* p=outer.sums.data()+size_t(worker)*rhs;for(size_t q=first;q<last;++q)for(uint32_t r=0;r<rhs;++r){auto z=a.b[q*rhs+r];p[r]+=z.re*z.re+z.im*z.im;}});a.bnorm.assign(rhs,0);for(uint32_t r=0;r<rhs;++r){for(int w=0;w<workers;++w)a.bnorm[r]+=outer.sums[size_t(w)*rhs+r];a.bnorm[r]=std::sqrt(a.bnorm[r]);}
}
inline void set_compact_rhs(CompactClover& c,CompactOuter& outer,uint32_t offset,int workers,bool zero=false){auto& a=c.w.a;uint32_t rhs=a.p.rhs;outer.offset=offset;outer.zero=zero;
    chunks(a.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t q=first;q<last;++q)for(uint32_t r=0;r<rhs;++r)outer.residual[q*rhs+r]=zero&&r==rhs-1?D{}:compact_truth(q,offset+r);});apply_cpu(a,outer.residual.data(),a.b.data(),workers);compact_rhs_norms(a,outer,workers);
}
template<class Inner> SolveResult compact_solve(const CompactClover& c,Inner& inner,CompactOuter& w,int workers,int limit,bool use_initial_guess=false){const auto& a=c.w.a;uint32_t rhs=a.p.rhs;SolveResult out;auto start=Clock::now();bool initially_converged=false;
    if(!use_initial_guess){std::fill(w.x.begin(),w.x.end(),D{});w.residual=a.b;}
    else{
        auto t=Clock::now();compact_fp64::residual(c,w.x.data(),a.b.data(),w.residual.data(),w.sums.data(),workers);
        for(uint32_t r=0;r<rhs;++r){double n=0;for(int worker=0;worker<workers;++worker)n+=w.sums[size_t(worker)*rhs+r];double rel=std::sqrt(n)/(a.bnorm[r]>0?a.bnorm[r]:1);if(!std::isfinite(rel))throw std::runtime_error("nonfinite initial full residual");w.relative[r]=rel;out.max_residual=std::max(out.max_residual,rel);}
        out.residual_ms+=seconds(t)*1000;initially_converged=out.max_residual<=1e-10;if(w.observe)w.observe(0,out.max_residual,0);
    }
    for(int outer=0;!initially_converged&&outer<12;++outer){auto t=Clock::now();auto info=inner.run(w.residual,.002,limit);out.inner_ms+=seconds(t)*1000;t=Clock::now();auto* dx=inner.correction();chunks(a.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t j=size_t(first)*rhs;j<size_t(last)*rhs;++j){w.x[j].re+=double(dx[j].re);w.x[j].im+=double(dx[j].im);}});out.update_ms+=seconds(t)*1000;
        t=Clock::now();compact_fp64::residual(c,w.x.data(),a.b.data(),w.residual.data(),w.sums.data(),workers);
        out.max_residual=0;for(uint32_t r=0;r<rhs;++r){double n=0;for(int worker=0;worker<workers;++worker)n+=w.sums[size_t(worker)*rhs+r];double rel=std::sqrt(n)/(a.bnorm[r]>0?a.bnorm[r]:1);if(!std::isfinite(rel))throw std::runtime_error("nonfinite compact full residual");w.relative[r]=rel;out.max_residual=std::max(out.max_residual,rel);}
        out.residual_ms+=seconds(t)*1000;out.outer=outer+1;out.inner_max_sum+=info.max_iterations;out.dispatched+=info.dispatched;out.trace.push_back({outer+1,out.max_residual,info});if(w.observe)w.observe(outer+1,out.max_residual,info.max_iterations);if(out.max_residual<=1e-10)break;
        chunks(a.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t q=first;q<last;++q)for(uint32_t r=0;r<rhs;++r)if(w.relative[r]<=1e-10)w.residual[q*rhs+r]={};});if(outer==11)throw std::runtime_error("compact refinement failed");
    }out.total_ms=seconds(start)*1000;
    if(!w.known_solution)return out;
    std::vector<double> errors(workers),refs(workers);chunks(a.p.rows*12,workers,[&](int worker,uint32_t first,uint32_t last){double e=0,n=0;for(size_t q=first;q<last;++q)for(uint32_t r=0;r<rhs;++r){auto exact=w.zero&&r==rhs-1?D{}:compact_truth(q,w.offset+r);auto z=minus_z(w.x[q*rhs+r],exact);e+=z.re*z.re+z.im*z.im;n+=exact.re*exact.re+exact.im*exact.im;}errors[worker]=e;refs[worker]=n;});double e=0,n=0;for(int i=0;i<workers;++i){e+=errors[i];n+=refs[i];}out.solution_error=std::sqrt(e/std::max(n,1e-300));if(!std::isfinite(out.solution_error)||out.solution_error>1e-7)throw std::runtime_error("compact known solution failed");return out;
}

template<class Inner> struct CompactHostCorrection {
    const CompactClover& c;Inner& inner;int workers;std::vector<Cx<float>> re,reduced,temp,aux,delta;
    CompactHostCorrection(const CompactClover& c,Inner& inner,int workers):c(c),inner(inner),workers(workers),re(c.parity.eo.p.outputs()),reduced(re.size()),temp(re.size()),aux(re.size()),delta(c.w.a.p.outputs()){}
    InnerInfo run(const std::vector<D>& residual,double tol,int limit){uint32_t rhs=c.w.a.p.rhs,width=12*rhs;chunks(c.parity.eo.p.rows,workers,[&](int,uint32_t first,uint32_t last){for(size_t i=first;i<last;++i)for(uint32_t q=0;q<width;++q){auto a=residual[size_t(c.parity.even[i])*width+q],b=residual[size_t(c.parity.odd[i])*width+q];re[i*width+q]={float(a.re),float(a.im)};reduced[i*width+q]={float(b.re),float(b.im)};}});
        local_blocks(c.inverseEven32,re.data(),temp.data(),c.parity.eo.p.rows,rhs,workers);apply_cpu(c.parity.oe,temp.data(),aux.data(),workers);chunks(c.parity.eo.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t j=size_t(first)*rhs;j<size_t(last)*rhs;++j)reduced[j]=minus_z(reduced[j],aux[j]);});auto info=inner.run(reduced,tol,limit);auto* odd=inner.correction();apply_cpu(c.parity.eo,odd,temp.data(),workers);chunks(c.parity.eo.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t j=size_t(first)*rhs;j<size_t(last)*rhs;++j)aux[j]=minus_z(re[j],temp[j]);});local_blocks(c.inverseEven32,aux.data(),temp.data(),c.parity.eo.p.rows,rhs,workers);
        chunks(c.parity.eo.p.rows,workers,[&](int,uint32_t first,uint32_t last){for(size_t i=first;i<last;++i)for(uint32_t q=0;q<width;++q){delta[size_t(c.parity.even[i])*width+q]=temp[i*width+q];delta[size_t(c.parity.odd[i])*width+q]=odd[i*width+q];}});return info;
    }
    const Cx<float>* correction()const{return delta.data();}
};

template<class Inner> struct CompactDeviceCorrection {
    const CompactClover& c;GPUCompactSchur& op;Inner& inner;int workers;id<MTLBuffer> input,delta,re,reduced;
    CompactDeviceCorrection(const CompactClover& c,GPUCompactSchur& op,Inner& inner,int workers):c(c),op(op),inner(inner),workers(workers){auto alloc=[&](size_t n){auto b=[op.m.device newBufferWithLength:n options:MTLResourceStorageModeShared];if(!b)throw std::runtime_error("compact correction allocation");return b;};input=alloc(c.w.a.p.outputs()*8);delta=input;re=alloc(c.parity.eo.p.outputs()*8);reduced=alloc(c.parity.eo.p.outputs()*8);}
    void dispatch(id<MTLComputeCommandEncoder> enc,id<MTLComputePipelineState> pipe,std::initializer_list<id<MTLBuffer>> buffers){[enc setComputePipelineState:pipe];size_t i=0;for(auto b:buffers){[enc setBuffer:b offset:0 atIndex:i];++i;}[enc setBytes:&op.params length:sizeof(op.params) atIndex:i];[enc dispatchThreads:MTLSizeMake(op.eo.outputs,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];[enc memoryBarrierWithScope:MTLBarrierScopeBuffers];}
    void finish(id<MTLCommandBuffer> cmd,id<MTLComputeCommandEncoder> enc){[enc endEncoding];[cmd commit];[cmd waitUntilCompleted];if(cmd.status==MTLCommandBufferStatusError)throw std::runtime_error([[cmd.error localizedDescription] UTF8String]);}
    InnerInfo run(const std::vector<D>& residual,double tol,int limit){return run(residual.data(),tol,limit);}
    InnerInfo run(const D* residual,double tol,int limit){uint32_t rhs=c.w.a.p.rhs;auto* in=static_cast<Cx<float>*>([input contents]);chunks(c.w.a.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t j=size_t(first)*rhs;j<size_t(last)*rhs;++j)in[j]={float(residual[j].re),float(residual[j].im)};});return run_prepared(tol,limit);}
    InnerInfo run_prepared(double tol,int limit){
        @autoreleasepool{auto cmd=[op.m.queue commandBuffer];auto enc=[cmd computeCommandEncoder];dispatch(enc,op.gather,{input,op.eo.map,op.oe.map,re,reduced});op.block(enc,op.inv,re,op.inverted);op.hop(enc,op.oe,op.inverted,op.temp);op.difference(enc,reduced,op.temp,reduced);finish(cmd,enc);}
        auto info=inner.run_device(reduced,tol,limit);@autoreleasepool{auto cmd=[op.m.queue commandBuffer];auto enc=[cmd computeCommandEncoder];op.hop(enc,op.eo,inner.delta,op.temp);op.difference(enc,re,op.temp,op.inverted);op.block(enc,op.inv,op.inverted,op.temp);dispatch(enc,op.scatter,{op.temp,inner.delta,input,op.eo.map,op.oe.map});finish(cmd,enc);}return info;
    }
    const Cx<float>* correction()const{return static_cast<Cx<float>*>([delta contents]);}
};
