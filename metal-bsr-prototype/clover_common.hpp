#pragma once
#include "evenodd.hpp"
using M6=std::array<CD,36>;
M3 m3_identity(){M3 a{};a[0]=a[4]=a[8]=1;return a;}
M3 m3_mul(const M3& a,const M3& b){M3 c{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)c[i*3+j]+=a[i*3+k]*b[k*3+j];return c;}
M3 m3_adj(const M3& a){M3 b;for(int i=0;i<3;++i)for(int j=0;j<3;++j)b[i*3+j]=std::conj(a[j*3+i]);return b;}
M3 gauge_loop(const WilsonSystem& w,uint32_t x,std::initializer_list<int> path){
    M3 product=m3_identity();uint32_t at=x;
    for(int step:path){uint32_t mu=uint32_t(std::abs(step)-1);bool f=step>0;M3 u;
        for(int c=0;c<3;++c)for(int k=0;k<3;++k)u[c*3+k]=w.link(at,mu,f,c,k);
        product=m3_mul(product,u);at=w.neighbor(at,mu,f);
    }
    if(at!=x)throw std::runtime_error("nonclosed clover path");return product;
}
std::vector<M3> clover_field(const WilsonSystem& w){
    std::vector<M3> f(w.a.p.rows*6);
    for(uint32_t x=0;x<w.a.p.rows;++x){int pair=0;
        for(int mu=0;mu<4;++mu)for(int nu=mu+1;nu<4;++nu){int m=mu+1,n=nu+1;
            auto p0=gauge_loop(w,x,{m,n,-m,-n}),p1=gauge_loop(w,x,{n,-m,-n,m});
            auto p2=gauge_loop(w,x,{-m,-n,m,n}),p3=gauge_loop(w,x,{-n,m,n,-m});M3 q;
            for(int j=0;j<9;++j)q[j]=p0[j]+p1[j]+p2[j]+p3[j];auto qa=m3_adj(q);
            for(int j=0;j<9;++j)f[x*6+pair][j]=(q[j]-qa[j])*.125;
            ++pair;
        }
    }return f; // mesField convention: anti-Hermitian, without trace subtraction.
}
M6 invert_cholesky(const M6& a,double& min_pivot){
    M6 l{},inv{};
    for(int i=0;i<6;++i)for(int j=0;j<=i;++j){CD z=a[i*6+j];for(int k=0;k<j;++k)z-=l[i*6+k]*std::conj(l[j*6+k]);
        if(i==j){if(!(z.real()>0)||std::abs(z.imag())>1e-12)throw std::runtime_error("clover block is not HPD");min_pivot=std::min(min_pivot,z.real());l[i*6+j]=std::sqrt(z.real());}
        else l[i*6+j]=z/l[j*6+j].real();
    }
    for(int col=0;col<6;++col){CD y[6]{},x[6]{};
        for(int i=0;i<6;++i){CD z=i==col?1:0;for(int k=0;k<i;++k)z-=l[i*6+k]*y[k];y[i]=z/l[i*6+i].real();}
        for(int i=5;i>=0;--i){CD z=y[i];for(int k=i+1;k<6;++k)z-=std::conj(l[k*6+i])*x[k];x[i]=z/l[i*6+i].real();}
        for(int i=0;i<6;++i)inv[i*6+col]=x[i];
    }return inv;
}
template<class T> std::vector<Cx<T>> block_cast(const std::vector<D>& a){std::vector<Cx<T>> b(a.size());for(size_t j=0;j<a.size();++j)b[j]={T(a[j].re),T(a[j].im)};return b;}

struct CloverSystem {
    WilsonSystem w;EvenOdd parity;double coefficient,min_pivot=1e100;
    std::vector<M3> field;
    std::vector<D> fullA,evenA,oddA,inverseEven;
    std::vector<Cx<float>> evenA32,oddA32,inverseEven32;
    CloverSystem(uint32_t l,uint32_t rhs,double mass,double gauge,double csw,bool apbc,int workers,bool zero_rhs=false):
        w(make_wilson(l,rhs,mass,gauge,apbc,workers,zero_rhs,true)),parity(w),coefficient(csw),field(clover_field(w)){
        auto& a=w.a;fullA.resize(a.p.rows*72);a.local64.resize(a.p.rows*72);
        for(uint32_t site=0;site<a.p.rows;++site){
            std::array<CD,144> c{};int pair=0;
            for(int mu=0;mu<4;++mu)for(int nu=mu+1;nu<4;++nu){
                for(int s=0;s<4;++s)for(int t=0;t<4;++t){CD spin=0;for(int k=0;k<4;++k)spin+=w.gamma[mu][s*4+k]*w.gamma[nu][k*4+t];
                    for(int i=0;i<3;++i)for(int j=0;j<3;++j)c[(s*3+i)*12+t*3+j]+=-.5*csw*spin*field[site*6+pair][i*3+j];
                }
                ++pair;
            }
            for(int ch=0;ch<2;++ch)for(int i=0;i<6;++i)for(int j=0;j<6;++j){
                CD z=c[(ch*6+i)*12+ch*6+j];size_t idx=(size_t(site)*2+ch)*36+i*6+j;
                a.local64[idx]={z.real(),z.imag()};if(i==j)z+=mass+4;fullA[idx]={z.real(),z.imag()};
            }
            for(int i=0;i<12;++i)for(int j=0;j<12;++j)if((i/6)!=(j/6)&&std::abs(c[i*12+j])>1e-14)throw std::runtime_error("clover mixes chirality");
        }
        a.local32=block_cast<float>(a.local64);a.solution_error_limit=1e-7;
        for(auto site:parity.even)for(int ch=0;ch<2;++ch){M6 block;
            for(int j=0;j<36;++j){auto z=fullA[(size_t(site)*2+ch)*36+j];block[j]={z.re,z.im};evenA.push_back(z);}
            auto inv=invert_cholesky(block,min_pivot);for(auto z:inv)inverseEven.push_back({z.real(),z.imag()});
        }
        for(auto site:parity.odd)for(int j=0;j<72;++j)oddA.push_back(fullA[size_t(site)*72+j]);
        evenA32=block_cast<float>(evenA);oddA32=block_cast<float>(oddA);inverseEven32=block_cast<float>(inverseEven);
        apply_cpu(a,a.truth.data(),a.b.data(),workers);a.bnorm.assign(rhs,0);
        for(size_t j=0;j<a.p.outputs();++j)a.bnorm[j%rhs]+=a.b[j].re*a.b[j].re+a.b[j].im*a.b[j].im;
        for(auto& n:a.bnorm)n=std::sqrt(n);
    }
    template<class T> const std::vector<Cx<T>>& inv()const{if constexpr(std::is_same_v<T,double>)return inverseEven;else return inverseEven32;}
    template<class T> const std::vector<Cx<T>>& odd()const{if constexpr(std::is_same_v<T,double>)return oddA;else return oddA32;}
};

template<class T,class Clover=CloverSystem> struct CPUCloverSchur {
    const Clover& c;int workers;std::vector<Cx<T>> temp,inverted;
    CPUCloverSchur(const Clover& c,int workers):c(c),workers(workers),temp(c.parity.eo.p.outputs()),inverted(temp.size()){}
    void apply(const Cx<T>* x,Cx<T>* y){
        const auto& e=c.parity;apply_cpu(e.eo,x,temp.data(),workers);
        local_blocks(c.template inv<T>(),temp.data(),inverted.data(),e.eo.p.rows,e.eo.p.rhs,workers);
        apply_cpu(e.oe,inverted.data(),y,workers);
        local_blocks(c.template odd<T>(),x,y,e.oe.p.rows,e.oe.p.rhs,workers,T(-1));
    }
};

struct GPUCloverSchur {
    MetalRun& m;const CloverSystem& c;GPUProblem eo,oe;
    id<MTLBuffer> inv,Aodd,temp,inverted;
    id<MTLComputePipelineState> local,subtract,gather,scatter;
    struct Params{uint32_t rows;float beta;} params;
    GPUCloverSchur(MetalRun& m,const CloverSystem& c):m(m),c(c),eo(m,c.parity.eo.p),oe(m,c.parity.oe.p){
        auto put=[&](const void* data,size_t n){auto b=[m.device newBufferWithLength:n options:MTLResourceStorageModeShared];if(!b)throw std::runtime_error("clover Metal allocation failed");if(data)std::memcpy([b contents],data,n);return b;};
        inv=put(c.inverseEven32.data(),c.inverseEven32.size()*8);Aodd=put(c.oddA32.data(),c.oddA32.size()*8);
        temp=put(nullptr,c.parity.eo.p.outputs()*8);inverted=put(nullptr,c.parity.eo.p.outputs()*8);
        auto values=[MTLFunctionConstantValues new];uint32_t rhs=c.w.a.p.rhs;[values setConstantValue:&rhs type:MTLDataTypeUInt atIndex:0];
        auto pipe=[&](const char* name){NSError* error=nil;auto fn=[m.library newFunctionWithName:[NSString stringWithUTF8String:name] constantValues:values error:&error];
            if(!fn)throw std::runtime_error([[error localizedDescription] UTF8String]);auto p=[m.device newComputePipelineStateWithFunction:fn error:&error];if(!p)throw std::runtime_error([[error localizedDescription] UTF8String]);return p;};
        local=pipe("clover_local");subtract=pipe("clover_subtract");gather=pipe("clover_gather");scatter=pipe("clover_scatter");params={eo.rows,0};
    }
    void hop(id<MTLComputeCommandEncoder> enc,GPUProblem& op,id<MTLBuffer> x,id<MTLBuffer> y){
        [enc setComputePipelineState:op.pipeline];for(size_t j=0;j<op.buffers.size();++j)[enc setBuffer:op.buffers[j] offset:0 atIndex:j];
        [enc setBuffer:x offset:0 atIndex:3];[enc setBuffer:y offset:0 atIndex:4];[enc setBytes:&op.rows length:4 atIndex:7];
        [enc dispatchThreads:MTLSizeMake(op.outputs,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];[enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
    }
    void block(id<MTLComputeCommandEncoder> enc,id<MTLBuffer> a,id<MTLBuffer> x,id<MTLBuffer> y,float beta=0){
        params.beta=beta;[enc setComputePipelineState:local];[enc setBuffer:a offset:0 atIndex:0];[enc setBuffer:x offset:0 atIndex:1];[enc setBuffer:y offset:0 atIndex:2];
        [enc setBytes:&params length:sizeof(params) atIndex:3];[enc dispatchThreads:MTLSizeMake(eo.outputs,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];[enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
    }
    void difference(id<MTLComputeCommandEncoder> enc,id<MTLBuffer> a,id<MTLBuffer> b,id<MTLBuffer> y){
        [enc setComputePipelineState:subtract];[enc setBuffer:a offset:0 atIndex:0];[enc setBuffer:b offset:0 atIndex:1];[enc setBuffer:y offset:0 atIndex:2];
        [enc setBytes:&params length:sizeof(params) atIndex:3];[enc dispatchThreads:MTLSizeMake(eo.outputs,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];[enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
    }
    void encode(id<MTLComputeCommandEncoder> enc,id<MTLBuffer> x,id<MTLBuffer> y){
        hop(enc,eo,x,temp);block(enc,inv,temp,inverted);hop(enc,oe,inverted,y);block(enc,Aodd,x,y,-1);
    }
};

template<class Inner,class Elim,class Clover=CloverSystem> struct HostCloverCorrection {
    const Clover& c;Inner& inner;int workers;
    std::vector<Cx<Elim>> re,ro,reduced,temp,aux,odd_delta,even_delta,delta;
    double prepare_ms=0,reconstruct_ms=0;
    HostCloverCorrection(const Clover& c,Inner& inner,int workers):c(c),inner(inner),workers(workers),re(c.parity.eo.p.outputs()),ro(re.size()),reduced(re.size()),
        temp(re.size()),aux(re.size()),odd_delta(re.size()),even_delta(re.size()),delta(c.w.a.p.outputs()){}
    void reset_metrics(){prepare_ms=reconstruct_ms=0;}
    InnerInfo run(const std::vector<D>& residual,double tol,int limit){
        auto start=Clock::now();const auto& e=c.parity;uint32_t width=12*c.w.a.p.rhs;
        chunks(e.eo.p.rows,workers,[&](int,uint32_t first,uint32_t last){for(size_t i=first;i<last;++i)for(uint32_t k=0;k<width;++k){
            auto a=residual[size_t(e.even[i])*width+k],b=residual[size_t(e.odd[i])*width+k];re[i*width+k]={Elim(a.re),Elim(a.im)};ro[i*width+k]={Elim(b.re),Elim(b.im)};}});
        local_blocks(c.template inv<Elim>(),re.data(),temp.data(),e.eo.p.rows,e.eo.p.rhs,workers);apply_cpu(e.oe,temp.data(),aux.data(),workers);
        chunks(e.eo.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){uint32_t rhs=e.eo.p.rhs;for(size_t j=size_t(first)*rhs;j<size_t(last)*rhs;++j)reduced[j]=minus_z(ro[j],aux[j]);});
        prepare_ms+=seconds(start)*1000;auto info=inner.run(reduced,tol,limit);start=Clock::now();auto* source=inner.correction();
        chunks(e.eo.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){uint32_t rhs=e.eo.p.rhs;for(size_t j=size_t(first)*rhs;j<size_t(last)*rhs;++j)odd_delta[j]={Elim(source[j].re),Elim(source[j].im)};});
        apply_cpu(e.eo,odd_delta.data(),temp.data(),workers);
        chunks(e.eo.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){uint32_t rhs=e.eo.p.rhs;for(size_t j=size_t(first)*rhs;j<size_t(last)*rhs;++j)aux[j]=minus_z(re[j],temp[j]);});
        local_blocks(c.template inv<Elim>(),aux.data(),even_delta.data(),e.eo.p.rows,e.eo.p.rhs,workers);
        chunks(e.eo.p.rows,workers,[&](int,uint32_t first,uint32_t last){for(size_t i=first;i<last;++i)for(uint32_t k=0;k<width;++k){size_t j=i*width+k;
            delta[size_t(e.even[i])*width+k]=even_delta[j];delta[size_t(e.odd[i])*width+k]=odd_delta[j];}});
        reconstruct_ms+=seconds(start)*1000;return info;
    }
    const Cx<Elim>* correction()const{return delta.data();}
};

template<class Inner,class Clover=CloverSystem,class Operator=GPUCloverSchur> struct DeviceCloverCorrectionT {
    const Clover& c;Operator& op;Inner& inner;int workers;
    id<MTLBuffer> input,re,ro,reduced,even_delta,delta,even_map,odd_map;
    double prepare_ms=0,reconstruct_ms=0;
    DeviceCloverCorrectionT(const Clover& c,Operator& op,Inner& inner,int workers):c(c),op(op),inner(inner),workers(workers){
        auto put=[&](const void* data,size_t bytes){auto b=[op.m.device newBufferWithLength:bytes options:MTLResourceStorageModeShared];if(!b)throw std::runtime_error("clover transfer allocation failed");if(data)std::memcpy([b contents],data,bytes);return b;};
        size_t half=c.parity.eo.p.outputs()*8,full=c.w.a.p.outputs()*8;
        input=put(nullptr,full);delta=put(nullptr,full);re=put(nullptr,half);ro=put(nullptr,half);reduced=put(nullptr,half);even_delta=put(nullptr,half);
        even_map=put(c.parity.even.data(),c.parity.even.size()*4);odd_map=put(c.parity.odd.data(),c.parity.odd.size()*4);
    }
    void reset_metrics(){prepare_ms=reconstruct_ms=0;}
    void dispatch(id<MTLComputeCommandEncoder> enc,id<MTLComputePipelineState> pipeline,std::initializer_list<id<MTLBuffer>> buffers){
        [enc setComputePipelineState:pipeline];size_t j=0;for(auto b:buffers){[enc setBuffer:b offset:0 atIndex:j];++j;}
        [enc setBytes:&op.params length:sizeof(op.params) atIndex:j];[enc dispatchThreads:MTLSizeMake(op.eo.outputs,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
    }
    void complete(id<MTLCommandBuffer> cmd,id<MTLComputeCommandEncoder> enc){[enc endEncoding];[cmd commit];[cmd waitUntilCompleted];
        if(cmd.status==MTLCommandBufferStatusError)throw std::runtime_error([[cmd.error localizedDescription] UTF8String]);}
    void prepare(const std::vector<D>& residual){
        uint32_t rhs=c.w.a.p.rhs;auto* ptr=static_cast<Cx<float>*>([input contents]);
        chunks(c.w.a.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t j=size_t(first)*rhs;j<size_t(last)*rhs;++j)ptr[j]={float(residual[j].re),float(residual[j].im)};});
        @autoreleasepool{auto cmd=[op.m.queue commandBuffer];auto enc=[cmd computeCommandEncoder];
            dispatch(enc,op.gather,{input,even_map,odd_map,re,ro});
            op.block(enc,op.inv,re,op.inverted);op.hop(enc,op.oe,op.inverted,op.temp);op.difference(enc,ro,op.temp,reduced);complete(cmd,enc);
        }
    }
    void reconstruct(){
        @autoreleasepool{auto cmd=[op.m.queue commandBuffer];auto enc=[cmd computeCommandEncoder];
            op.hop(enc,op.eo,inner.delta,op.temp);op.difference(enc,re,op.temp,op.inverted);op.block(enc,op.inv,op.inverted,even_delta);
            dispatch(enc,op.scatter,{even_delta,inner.delta,delta,even_map,odd_map});complete(cmd,enc);
        }
    }
    InnerInfo run(const std::vector<D>& residual,double tol,int limit){
        auto start=Clock::now();prepare(residual);prepare_ms+=seconds(start)*1000;
        auto info=inner.run_device(reduced,tol,limit);start=Clock::now();reconstruct();reconstruct_ms+=seconds(start)*1000;return info;
    }
    const Cx<float>* correction()const{return static_cast<Cx<float>*>([delta contents]);}
};
using DeviceCloverCorrection=DeviceCloverCorrectionT<GPUBiCGStab>;
