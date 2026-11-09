#pragma once
#include "sparse_coarse.hpp"

// Same gauge generator, QDP gamma basis and clover normalization as CloverSystem.
// Only forward links, neighbor indices, and local 6x6 blocks are retained.
struct CompactClover {
    WilsonSystem w;
    struct Parity {System eo,oe;std::vector<uint32_t> even,odd,index,parity;} parity;
    std::vector<uint32_t> neighbors;
    std::vector<Cx<float>> links32,evenA32,oddA32,inverseEven32,inverseOdd32;
    std::vector<D> evenA,oddA,inverseEven;
    std::array<uint32_t,16> perm{};std::array<D,16> phase{};
    std::array<double,4> hopping{1,1,1,1};
    double magnetic=1,electric=1;
    double coefficient=1,min_pivot=1e100,setup_ms=0;int workers;
    CompactClover(uint32_t l,double mass,int workers,double strength=.25,std::vector<M3> links={})
        :CompactClover(hypercube(l),mass,workers,strength,std::move(links),1,1,1,1){}
    CompactClover(LatticeShape shape,double mass,int workers,double strength,
                  std::vector<M3> links,double xi0,double nu,double clover_r,double clover_t)
        :magnetic(clover_r/xi0),electric(clover_t),workers(workers){
        auto start=Clock::now();w.extent=shape[0];w.dimensions=shape;w.strength=strength;w.antiperiodic=true;
        hopping={nu/xi0,nu/xi0,nu/xi0,1};
        for(auto& g:w.gamma)for(auto& z:g)z=std::conj(z);
        auto& a=w.a;auto& p=a.p;p.rows=lattice_volume(shape);p.colors=3;p.spins=4;p.rhs=1;a.mass=mass;a.solution_error_limit=1e-7;
        for(uint32_t mu=0;mu<4;++mu)for(uint32_t s=0;s<4;++s)for(uint32_t k=0;k<4;++k)if(std::abs(w.gamma[mu][s*4+k])>0){perm[mu*4+s]=k;auto z=w.gamma[mu][s*4+k];phase[mu*4+s]={z.real(),z.imag()};}
        if(!links.empty()){
            if(links.size()!=size_t(p.rows)*4)throw std::runtime_error("wrong number of imported gauge links");
            w.links=std::move(links);
        }else{
            std::mt19937_64 rng(77347);std::uniform_real_distribution<double> angle(-strength,strength),ph(-std::acos(-1.),std::acos(-1.));
            w.links.resize(size_t(p.rows)*4);
            for(auto& u:w.links){u=m3_identity();for(auto pair:{std::pair<int,int>{0,1},{1,2},{0,2}}){double theta=angle(rng),phi=ph(rng);CD s=std::sin(theta)*std::polar(1.,phi);double c=std::cos(theta);
                for(int k=0;k<3;++k){auto x=u[pair.first*3+k],y=u[pair.second*3+k];u[pair.first*3+k]=c*x+s*y;u[pair.second*3+k]=-std::conj(s)*x+c*y;}}
            double t0=angle(rng),t1=angle(rng);CD d[3]={std::polar(1.,t0),std::polar(1.,t1),std::polar(1.,-t0-t1)};for(int c=0;c<3;++c)for(int k=0;k<3;++k)u[c*3+k]*=d[c];}
        }
        links32.resize(size_t(p.rows)*36);for(size_t i=0;i<w.links.size();++i)for(int j=0;j<9;++j)links32[i*9+j]={float(hopping[i%4]*w.links[i][j].real()),float(hopping[i%4]*w.links[i][j].imag())};
        parity.index.resize(p.rows);parity.parity.resize(p.rows);neighbors.resize(size_t(p.rows)*8);
        for(uint32_t i=0;i<p.rows;++i){uint32_t q=i,pa=0;for(int mu=0;mu<4;++mu){pa+=q%shape[mu];q/=shape[mu];}pa%=2;parity.parity[i]=pa;auto& list=pa?parity.odd:parity.even;parity.index[i]=uint32_t(list.size());list.push_back(i);
            for(uint32_t d=0;d<8;++d)neighbors[size_t(i)*8+d]=w.neighbor(i,d/2,d%2==0);}
        uint32_t half=p.rows/2;evenA.resize(size_t(half)*72);oddA.resize(evenA.size());inverseEven.resize(evenA.size());inverseOdd32.resize(evenA.size());
        std::vector<double> pivots(workers,1e100);
        chunks(p.rows,workers,[&](int worker,uint32_t first,uint32_t last){for(uint32_t site=first;site<last;++site){std::array<CD,144> cl{};
            for(int mu=0;mu<4;++mu)for(int nu=mu+1;nu<4;++nu){int m=mu+1,n=nu+1;auto a=gauge_loop(w,site,{m,n,-m,-n}),b=gauge_loop(w,site,{n,-m,-n,m}),c=gauge_loop(w,site,{-m,-n,m,n}),d=gauge_loop(w,site,{-n,m,n,-m});M3 q;for(int j=0;j<9;++j)q[j]=a[j]+b[j]+c[j]+d[j];auto adj=m3_adj(q);for(int j=0;j<9;++j)q[j]=(q[j]-adj[j])*.125;
                for(int s=0;s<4;++s)for(int t=0;t<4;++t){CD spin=0;for(int k=0;k<4;++k)spin+=w.gamma[mu][s*4+k]*w.gamma[nu][k*4+t];for(int i=0;i<3;++i)for(int j=0;j<3;++j)cl[(s*3+i)*12+t*3+j]+=-.5*(nu==3?electric:magnetic)*spin*q[i*3+j];}}
            bool odd=parity.parity[site];uint32_t ix=parity.index[site];auto& blocks=odd?oddA:evenA;
            for(int ch=0;ch<2;++ch){M6 block;for(int i=0;i<6;++i)for(int j=0;j<6;++j){auto z=cl[(ch*6+i)*12+ch*6+j]+CD(i==j?mass+hopping[0]+hopping[1]+hopping[2]+hopping[3]:0);block[i*6+j]=z;blocks[(size_t(ix)*2+ch)*36+i*6+j]={z.real(),z.imag()};}
                auto inverse=invert_cholesky(block,pivots[worker]);for(int j=0;j<36;++j){auto z=inverse[j];if(odd)inverseOdd32[(size_t(ix)*2+ch)*36+j]={float(z.real()),float(z.imag())};else inverseEven[(size_t(ix)*2+ch)*36+j]={z.real(),z.imag()};}}
        }});min_pivot=*std::min_element(pivots.begin(),pivots.end());
        evenA32=block_cast<float>(evenA);oddA32=block_cast<float>(oddA);inverseEven32=block_cast<float>(inverseEven);
        for(auto* h:{&parity.eo,&parity.oe}){h->mass=mass;h->p.rows=half;h->p.colors=3;h->p.spins=4;h->p.rhs=1;}
        a.apply64=[this](const D* x,D* y,int n){apply(x,y,n,-1,true);};a.apply32=[this](const Cx<float>* x,Cx<float>* y,int n){apply(x,y,n,-1,true);};
        parity.eo.apply64=[this](const D* x,D* y,int n){apply(x,y,n,0,false);};parity.eo.apply32=[this](const Cx<float>* x,Cx<float>* y,int n){apply(x,y,n,0,false);};
        parity.oe.apply64=[this](const D* x,D* y,int n){apply(x,y,n,1,false);};parity.oe.apply32=[this](const Cx<float>* x,Cx<float>* y,int n){apply(x,y,n,1,false);};
        setup_ms=seconds(start)*1000;
    }
    void rhs(uint32_t n,bool allocate=true){w.a.p.rhs=parity.eo.p.rhs=parity.oe.p.rhs=n;if(allocate){w.a.truth.resize(w.a.p.outputs());w.a.b.resize(w.a.p.outputs());w.a.bnorm.resize(n);}}
    template<class T> const std::vector<Cx<T>>& inv()const{if constexpr(std::is_same_v<T,double>)return inverseEven;else return inverseEven32;}
    template<class T> const std::vector<Cx<T>>& odd()const{if constexpr(std::is_same_v<T,double>)return oddA;else return oddA32;}
    template<class T> Cx<T> gauge(uint32_t site,uint32_t d,uint32_t c,uint32_t k)const{
        uint32_t mu=d/2,n=neighbors[size_t(site)*8+d];bool f=d%2==0;Cx<T> z;
        if constexpr(std::is_same_v<T,double>){auto a=w.links[(f?site:n)*4+mu][f?c*3+k:k*3+c];z={a.real(),f?a.imag():-a.imag()};}
        else {z=links32[((size_t(f?site:n)*4+mu)*3+(f?c:k))*3+(f?k:c)];if(!f)z.im=-z.im;}
        T scale=T(-.5*w.bc(site,mu,f));if constexpr(std::is_same_v<T,double>)scale*=hopping[mu];return {scale*z.re,scale*z.im};
    }
    template<class T> void apply(const Cx<T>* x,Cx<T>* y,int workers,int pa,bool onsite)const{
        using V=SIMD<T>;uint32_t rows=pa<0?w.a.p.rows:w.a.p.rows/2,nrhs=w.a.p.rhs;
        chunks(rows,workers,[&](int,uint32_t first,uint32_t last){for(uint32_t row=first;row<last;++row){uint32_t site=pa<0?row:(pa?parity.odd[row]:parity.even[row]);
            for(uint32_t color=0;color<3;++color)for(uint32_t s=0;s<4;++s){uint32_t r=0;size_t out=(size_t(row)*12+color*4+s)*nrhs;
                auto terms=[&](const auto& emit){for(uint32_t d=0;d<8;++d){uint32_t mu=d/2,n=neighbors[size_t(site)*8+d],ni=pa<0?n:parity.index[n],si=perm[mu*4+s];auto ph=phase[mu*4+s];T sign=d%2==0?T(-1):T(1);
                        for(uint32_t k=0;k<3;++k){auto z=gauge<T>(site,d,color,k);emit(z,(size_t(ni)*12+k*4+s)*nrhs);emit(cm(z,Cx<T>{sign*T(ph.re),sign*T(ph.im)}),(size_t(ni)*12+k*4+si)*nrhs);}}
                    if(onsite){uint32_t ch=s/2,q=(s%2)*3+color;const auto& blocks=[&]()->const std::vector<Cx<T>>&{if constexpr(std::is_same_v<T,double>)return parity.parity[site]?oddA:evenA;else return parity.parity[site]?oddA32:evenA32;}();
                        for(uint32_t k=0;k<6;++k)emit(blocks[(size_t(parity.index[site])*2+ch)*36+q*6+k],(size_t(row)*12+(k%3)*4+ch*2+k/3)*nrhs);}};
                for(;r+V::width<=nrhs;r+=V::width){auto re=V::zero(),im=V::zero();terms([&](Cx<T> z,size_t q){auto v=V::load(x+q+r);re=V::ma(re,v.first,z.re);re=V::ms(re,v.second,z.im);im=V::ma(im,v.second,z.re);im=V::ma(im,v.first,z.im);});V::store(y+out+r,re,im);}
                for(;r<nrhs;++r){Cx<T> sum{};terms([&](Cx<T> z,size_t q){sum=plus(sum,cm(z,x[q+r]));});y[out+r]=sum;}
            }
        }});
    }
};

struct GPUCompactSchur {
    MetalRun& m;const CompactClover& c;
    struct Half{uint32_t rows;size_t outputs;id<MTLBuffer> map;} eo,oe;
    id<MTLBuffer> links,neighbors,inv,Aodd,temp,inverted,perms,phases;
    id<MTLComputePipelineState> hopping,local,subtract,gather,scatter;
    struct Params{uint32_t rows;float beta;} params;
    GPUCompactSchur(MetalRun& m,const CompactClover& c):m(m),c(c){
        auto put=[&](const void* p,size_t n){auto b=[m.device newBufferWithLength:std::max(size_t(8),n) options:MTLResourceStorageModeShared];if(!b)throw std::runtime_error("compact Metal allocation");if(p)std::memcpy([b contents],p,n);return b;};
        links=put(c.links32.data(),c.links32.size()*8);neighbors=put(c.neighbors.data(),c.neighbors.size()*4);inv=put(c.inverseEven32.data(),c.inverseEven32.size()*8);Aodd=put(c.oddA32.data(),c.oddA32.size()*8);
        eo={c.parity.eo.p.rows,c.parity.eo.p.outputs(),put(c.parity.even.data(),c.parity.even.size()*4)};oe={eo.rows,eo.outputs,put(c.parity.odd.data(),c.parity.odd.size()*4)};temp=put(nullptr,eo.outputs*8);inverted=put(nullptr,eo.outputs*8);
        perms=put(c.perm.data(),64);auto ph=block_cast<float>(std::vector<D>(c.phase.begin(),c.phase.end()));phases=put(ph.data(),128);params={eo.rows,0};
        auto values=[MTLFunctionConstantValues new];uint32_t rhs=c.w.a.p.rhs;[values setConstantValue:&rhs type:MTLDataTypeUInt atIndex:0];
        auto pipe=[&](const char* n){NSError* e=nil;auto f=[m.library newFunctionWithName:[NSString stringWithUTF8String:n] constantValues:values error:&e];auto p=[m.device newComputePipelineStateWithFunction:f error:&e];if(!p)throw std::runtime_error([[e localizedDescription] UTF8String]);return p;};
        hopping=pipe("compact_hop_rect");local=pipe("clover_local");subtract=pipe("clover_subtract");gather=pipe("clover_gather");scatter=pipe("clover_scatter");
    }
    void hop(id<MTLComputeCommandEncoder> enc,const Half& h,id<MTLBuffer> x,id<MTLBuffer> y){struct P{uint32_t rows,time_stride,time_extent,antiperiodic;} p{h.rows,c.w.stride(3),c.w.length(3),uint32_t(c.w.antiperiodic)};[enc setComputePipelineState:hopping];size_t i=0;for(auto b:{links,neighbors,h.map,perms,phases,x,y}){[enc setBuffer:b offset:0 atIndex:i];++i;}[enc setBytes:&p length:sizeof(p) atIndex:i];[enc dispatchThreads:MTLSizeMake(h.outputs,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];[enc memoryBarrierWithScope:MTLBarrierScopeBuffers];}
    void block(id<MTLComputeCommandEncoder> enc,id<MTLBuffer> a,id<MTLBuffer> x,id<MTLBuffer> y,float beta=0){params.beta=beta;[enc setComputePipelineState:local];size_t i=0;for(auto b:{a,x,y}){[enc setBuffer:b offset:0 atIndex:i];++i;}[enc setBytes:&params length:sizeof(params) atIndex:3];[enc dispatchThreads:MTLSizeMake(eo.outputs,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];[enc memoryBarrierWithScope:MTLBarrierScopeBuffers];}
    void difference(id<MTLComputeCommandEncoder> enc,id<MTLBuffer> a,id<MTLBuffer> b,id<MTLBuffer> y){[enc setComputePipelineState:subtract];size_t i=0;for(auto p:{a,b,y}){[enc setBuffer:p offset:0 atIndex:i];++i;}[enc setBytes:&params length:sizeof(params) atIndex:3];[enc dispatchThreads:MTLSizeMake(eo.outputs,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];[enc memoryBarrierWithScope:MTLBarrierScopeBuffers];}
    void encode(id<MTLComputeCommandEncoder> enc,id<MTLBuffer> x,id<MTLBuffer> y){hop(enc,eo,x,temp);block(enc,inv,temp,inverted);hop(enc,oe,inverted,y);block(enc,Aodd,x,y,-1);}
};
