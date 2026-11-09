#pragma once
#include "compact_clover.hpp"

struct BlockMatrix {
    uint32_t bs,rows,n;int workers;double setup_ms=0;
    std::vector<uint32_t> row,col,diagonal;std::vector<D> a64,di64;std::vector<Cx<float>> a32,di32;
    BlockMatrix(uint32_t bs,uint32_t rows,int workers):bs(bs),rows(rows),n(bs*rows),workers(workers){}
    void pattern(const std::vector<std::vector<uint32_t>>& graph){row={0};diagonal.resize(rows);for(uint32_t a=0;a<rows;++a){auto v=graph[a];std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());for(auto b:v){if(a==b)diagonal[a]=uint32_t(col.size());col.push_back(b);}row.push_back(uint32_t(col.size()));}a64.resize(size_t(col.size())*bs*bs);}
    size_t position(uint32_t a,uint32_t b)const{auto it=std::lower_bound(col.begin()+row[a],col.begin()+row[a+1],b);if(it==col.begin()+row[a+1]||*it!=b)throw std::runtime_error("missing coarse edge");return it-col.begin();}
    void finish(){di64.resize(size_t(rows)*bs*bs);chunks(rows,workers,[&](int,uint32_t first,uint32_t last){for(uint32_t a=first;a<last;++a){std::vector<CD> b(bs*bs);for(size_t j=0;j<b.size();++j){auto z=a64[size_t(diagonal[a])*bs*bs+j];b[j]={z.re,z.im};}auto inv=dense_inverse(b,bs);for(size_t j=0;j<b.size();++j)di64[size_t(a)*bs*bs+j]={inv[j].real(),inv[j].imag()};}});a32=block_cast<float>(a64);di32=block_cast<float>(di64);}
    template<class T> void apply(const Cx<T>* x,Cx<T>* y,uint32_t rhs,bool diag=false)const{
        using V=SIMD<T>;const auto& a=[&]()->const std::vector<Cx<T>>&{if constexpr(std::is_same_v<T,double>)return diag?di64:a64;else return diag?di32:a32;}();
        chunks(rows,workers,[&](int,uint32_t first,uint32_t last){for(uint32_t rowi=first;rowi<last;++rowi)for(uint32_t i=0;i<bs;++i){auto each=[&](const auto& f){if(diag){for(uint32_t j=0;j<bs;++j)f(a[(size_t(rowi)*bs+i)*bs+j],rowi*bs+j);}else for(uint32_t e=row[rowi];e<row[rowi+1];++e)for(uint32_t j=0;j<bs;++j)f(a[(size_t(e)*bs+i)*bs+j],col[e]*bs+j);};uint32_t r=0;
            for(;r+V::width<=rhs;r+=V::width){auto re=V::zero(),im=V::zero();each([&](Cx<T> z,uint32_t q){auto v=V::load(x+size_t(q)*rhs+r);re=V::ma(re,v.first,z.re);re=V::ms(re,v.second,z.im);im=V::ma(im,v.second,z.re);im=V::ma(im,v.first,z.im);});V::store(y+size_t(rowi*bs+i)*rhs+r,re,im);}
            for(;r<rhs;++r){Cx<T> sum{};each([&](Cx<T> z,uint32_t q){sum=plus(sum,cm(z,x[size_t(q)*rhs+r]));});y[size_t(rowi*bs+i)*rhs+r]=sum;}
        }});
    }
    void release64(){std::vector<D>().swap(a64);std::vector<D>().swap(di64);}
};

struct Transfer {
    uint32_t rows,inbs,nvec,aggregates,coarse,grid,block;bool spin_layout;int workers;
    std::vector<uint32_t> aggregate,offsets,sites;std::vector<D> basis;std::vector<Cx<float>> basis32;
    Transfer(uint32_t rows,uint32_t inbs,uint32_t nvec,uint32_t grid,uint32_t block,bool spin,int workers,const std::vector<uint32_t>* fullsites=nullptr)
        :Transfer(rows,inbs,nvec,hypercube(grid),hypercube(block),spin,workers,fullsites){}
    Transfer(uint32_t rows,uint32_t inbs,uint32_t nvec,LatticeShape shape,LatticeShape blocks,bool spin,int workers,const std::vector<uint32_t>* fullsites=nullptr):rows(rows),inbs(inbs),nvec(nvec),grid(shape[0]),block(blocks[0]),spin_layout(spin),workers(workers){
        LatticeShape coarse_shape{};aggregates=1;
        for(int mu=0;mu<4;++mu){if(!blocks[mu]||shape[mu]%blocks[mu])throw std::runtime_error("aggregation must divide each extent");coarse_shape[mu]=shape[mu]/blocks[mu];aggregates*=coarse_shape[mu];}
        coarse=aggregates*2*nvec;aggregate.resize(rows);std::vector<std::vector<uint32_t>> lists(aggregates);
        for(uint32_t i=0;i<rows;++i){uint32_t q=fullsites?(*fullsites)[i]:i,a=0,stride=1;for(int mu=0;mu<4;++mu){a+=(q%shape[mu]/blocks[mu])*stride;q/=shape[mu];stride*=coarse_shape[mu];}aggregate[i]=a;lists[a].push_back(i);}offsets={0};for(auto& l:lists){sites.insert(sites.end(),l.begin(),l.end());offsets.push_back(uint32_t(sites.size()));}
    }
    uint32_t chirality(uint32_t q)const{return spin_layout?(q%4/2):(q/(inbs/2));}
    template<class Apply,class Diagonal> void build(Apply apply,Diagonal diagonal,int sweeps,double omega=0.65,uint32_t filter_rhs=0){
        std::mt19937_64 rng(771337+inbs);std::normal_distribution<double> rnd(0,1);std::vector<CD> constant(inbs*nvec);for(auto& z:constant)z={rnd(rng),rnd(rng)};
        basis.resize(size_t(rows)*inbs*nvec);for(size_t i=0;i<size_t(rows)*inbs;++i)for(uint32_t v=0;v<nvec;++v){auto x=constant[i%inbs*nvec+v]+.2*CD(rnd(rng),rnd(rng));basis[i*nvec+v]={x.real(),x.imag()};}
        if(filter_rhs){
            if(nvec%filter_rhs)throw std::runtime_error("basis filter tile must divide vector count");
            size_t components=size_t(rows)*inbs;
            std::vector<D> input(components*filter_rhs),a(input.size()),z(input.size());
            for(uint32_t first_vector=0;first_vector<nvec;first_vector+=filter_rhs){
                for(size_t q=0;q<components;++q)std::copy_n(basis.data()+q*nvec+first_vector,filter_rhs,input.data()+q*filter_rhs);
                for(int it=0;it<sweeps;++it){
                    apply(input.data(),a.data());diagonal(a.data(),z.data());
                    chunks(rows*inbs,workers,[&](int,uint32_t first,uint32_t last){for(size_t i=size_t(first)*filter_rhs;i<size_t(last)*filter_rhs;++i){input[i].re-=omega*z[i].re;input[i].im-=omega*z[i].im;}});
                }
                for(size_t q=0;q<components;++q)std::copy_n(input.data()+q*filter_rhs,filter_rhs,basis.data()+q*nvec+first_vector);
                std::cout<<"BASIS filtered "<<first_vector+filter_rhs<<"/"<<nvec<<std::endl;
            }
        }else{
            std::vector<D> a(basis.size()),z(basis.size());
            for(int it=0;it<sweeps;++it){apply(basis.data(),a.data());diagonal(a.data(),z.data());chunks(rows*inbs,workers,[&](int,uint32_t first,uint32_t last){for(size_t i=size_t(first)*nvec;i<size_t(last)*nvec;++i){basis[i].re-=omega*z[i].re;basis[i].im-=omega*z[i].im;}});}
        }
        chunks(aggregates,workers,[&](int,uint32_t first,uint32_t last){for(uint32_t ag=first;ag<last;++ag)for(uint32_t ch=0;ch<2;++ch)for(uint32_t v=0;v<nvec;++v){
            for(int pass=0;pass<2;++pass)for(uint32_t k=0;k<v;++k){CD dot=0;for(uint32_t p=offsets[ag];p<offsets[ag+1];++p)for(uint32_t q=0;q<inbs;++q)if(chirality(q)==ch){size_t i=(size_t(sites[p])*inbs+q)*nvec;auto a=basis[i+k],b=basis[i+v];dot+=std::conj(CD(a.re,a.im))*CD(b.re,b.im);}
                for(uint32_t p=offsets[ag];p<offsets[ag+1];++p)for(uint32_t q=0;q<inbs;++q)if(chirality(q)==ch){size_t i=(size_t(sites[p])*inbs+q)*nvec;auto b=basis[i+k];auto z=dot*CD(b.re,b.im);basis[i+v].re-=z.real();basis[i+v].im-=z.imag();}}
            double norm=0;for(uint32_t p=offsets[ag];p<offsets[ag+1];++p)for(uint32_t q=0;q<inbs;++q)if(chirality(q)==ch){auto z=basis[(size_t(sites[p])*inbs+q)*nvec+v];norm+=z.re*z.re+z.im*z.im;}if(norm<1e-24)throw std::runtime_error("rank deficient recursive basis");norm=std::sqrt(norm);
            for(uint32_t p=offsets[ag];p<offsets[ag+1];++p)for(uint32_t q=0;q<inbs;++q)if(chirality(q)==ch){auto& z=basis[(size_t(sites[p])*inbs+q)*nvec+v];z.re/=norm;z.im/=norm;}
        }});basis32=block_cast<float>(basis);
    }
    template<class T> const std::vector<Cx<T>>& weights()const{if constexpr(std::is_same_v<T,double>)return basis;else return basis32;}
    template<class T> void restrict_vector(const Cx<T>* x,Cx<T>* y,uint32_t rhs)const{const auto& b=weights<T>();using V=SIMD<T>;chunks(coarse,workers,[&](int,uint32_t first,uint32_t last){for(uint32_t out=first;out<last;++out){uint32_t ag=out/(2*nvec),ch=out%(2*nvec)/nvec,v=out%nvec,r=0;auto each=[&](const auto& f){for(uint32_t p=offsets[ag];p<offsets[ag+1];++p)for(uint32_t q=0;q<inbs;++q)if(chirality(q)==ch){size_t i=size_t(sites[p])*inbs+q;auto z=b[i*nvec+v];f(Cx<T>{z.re,-z.im},i);}};
            for(;r+V::width<=rhs;r+=V::width){auto re=V::zero(),im=V::zero();each([&](Cx<T> z,size_t i){auto a=V::load(x+i*rhs+r);re=V::ma(re,a.first,z.re);re=V::ms(re,a.second,z.im);im=V::ma(im,a.second,z.re);im=V::ma(im,a.first,z.im);});V::store(y+size_t(out)*rhs+r,re,im);}for(;r<rhs;++r){Cx<T> sum{};each([&](Cx<T> z,size_t i){sum=plus(sum,cm(z,x[i*rhs+r]));});y[size_t(out)*rhs+r]=sum;}
        }});}
    template<class T> void prolong(const Cx<T>* x,Cx<T>* y,uint32_t rhs,bool add=false)const{const auto& b=weights<T>();chunks(rows*inbs,workers,[&](int,uint32_t first,uint32_t last){for(size_t q=first;q<last;++q){uint32_t base=aggregate[q/inbs]*2*nvec+chirality(q%inbs)*nvec;for(uint32_t r=0;r<rhs;++r){Cx<T> sum{};for(uint32_t v=0;v<nvec;++v)sum=plus(sum,cm(b[q*nvec+v],x[size_t(base+v)*rhs+r]));y[q*rhs+r]=add?plus(y[q*rhs+r],sum):sum;}}});}
    void local_matrix(uint32_t site,D* out)const{std::fill(out,out+inbs*2*nvec,D{});for(uint32_t q=0;q<inbs;++q)for(uint32_t v=0;v<nvec;++v)out[q*2*nvec+chirality(q)*nvec+v]=basis[(size_t(site)*inbs+q)*nvec+v];}
    void release64(){std::vector<D>().swap(basis);}
};

inline void gemm(char ta,const D* a,const D* b,D* c,uint32_t m,uint32_t n,uint32_t k,double beta=0){D one{1,0},be{beta,0};cblas_zgemm(CblasRowMajor,ta=='C'?CblasConjTrans:CblasNoTrans,CblasNoTrans,m,n,k,reinterpret_cast<const __LAPACK_double_complex*>(&one),reinterpret_cast<const __LAPACK_double_complex*>(a),ta=='C'?m:k,reinterpret_cast<const __LAPACK_double_complex*>(b),n,reinterpret_cast<const __LAPACK_double_complex*>(&be),reinterpret_cast<__LAPACK_double_complex*>(c),n);}

inline void assemble_schur(BlockMatrix& a,const CompactClover& c,const Transfer& h){auto start=Clock::now();uint32_t bs=a.bs;size_t block=size_t(bs)*bs;std::vector<std::vector<uint32_t>> reverse(a.rows),graph(a.rows);
    for(uint32_t e=0;e<h.rows;++e){uint32_t site=c.parity.even[e];std::vector<uint32_t> ag;for(uint32_t d=0;d<8;++d)ag.push_back(h.aggregate[c.parity.index[c.neighbors[size_t(site)*8+d]]]);std::sort(ag.begin(),ag.end());ag.erase(std::unique(ag.begin(),ag.end()),ag.end());for(auto x:ag){reverse[x].push_back(e);graph[x].insert(graph[x].end(),ag.begin(),ag.end());}}
    a.pattern(graph);chunks(a.rows,a.workers,[&](int,uint32_t first,uint32_t last){std::vector<D> gram(block),pmat(12*bs),ap(12*bs),u(12*bs);
        for(uint32_t ag=first;ag<last;++ag){for(auto e:reverse[ag]){uint32_t site=c.parity.even[e];std::array<D,72> even_inverse;c.even_inverse_blocks(e,even_inverse.data());std::vector<uint32_t> list;for(uint32_t d=0;d<8;++d)list.push_back(h.aggregate[c.parity.index[c.neighbors[size_t(site)*8+d]]]);std::sort(list.begin(),list.end());list.erase(std::unique(list.begin(),list.end()),list.end());std::vector<D> t(size_t(list.size())*12*bs);
            for(uint32_t d=0;d<8;++d){uint32_t odd=c.parity.index[c.neighbors[size_t(site)*8+d]],pos=uint32_t(std::lower_bound(list.begin(),list.end(),h.aggregate[odd])-list.begin());auto* out=t.data()+size_t(pos)*12*bs;
                for(uint32_t color=0;color<3;++color)for(uint32_t s=0;s<4;++s)for(uint32_t k=0;k<3;++k){auto z=c.gauge<double>(site,d,color,k);uint32_t si=c.perm[(d/2)*4+s];auto ph=c.phase[(d/2)*4+s];double sign=d%2==0?-1:1;
                    for(auto term:{std::pair<uint32_t,D>{s,z},{si,cm(z,D{sign*ph.re,sign*ph.im})}})for(uint32_t v=0;v<h.nvec;++v){auto p=h.basis[(size_t(odd)*12+k*4+term.first)*h.nvec+v];auto& y=out[(color*4+s)*bs+(term.first/2)*h.nvec+v];y=plus(y,cm(term.second,p));}}
            }
            auto* left=t.data()+size_t(std::lower_bound(list.begin(),list.end(),ag)-list.begin())*12*bs;
            for(uint32_t p=0;p<list.size();++p){auto* right=t.data()+size_t(p)*12*bs;std::fill(u.begin(),u.end(),D{});
                for(uint32_t ch=0;ch<2;++ch)for(uint32_t q=0;q<6;++q)for(uint32_t v=0;v<bs;++v){D sum{};for(uint32_t k=0;k<6;++k)sum=plus(sum,cm(even_inverse[ch*36+q*6+k],right[((k%3)*4+ch*2+k/3)*bs+v]));double sign=ch?-1:1;u[((q%3)*4+ch*2+q/3)*bs+v]={sign*sum.re,sign*sum.im};}
                gemm('C',left,u.data(),gram.data(),bs,bs,12);auto* out=a.a64.data()+a.position(ag,list[p])*block;for(uint32_t i=0;i<bs;++i)for(uint32_t j=0;j<bs;++j){double sign=i<h.nvec?-1:1;out[i*bs+j].re+=sign*gram[i*bs+j].re;out[i*bs+j].im+=sign*gram[i*bs+j].im;}}
        }
        auto* out=a.a64.data()+size_t(a.diagonal[ag])*block;for(uint32_t p=h.offsets[ag];p<h.offsets[ag+1];++p){uint32_t site=h.sites[p];h.local_matrix(site,pmat.data());std::fill(ap.begin(),ap.end(),D{});
            for(uint32_t color=0;color<3;++color)for(uint32_t s=0;s<4;++s){uint32_t ch=s/2,q=(s%2)*3+color;for(uint32_t k=0;k<6;++k)for(uint32_t v=0;v<bs;++v)ap[(color*4+s)*bs+v]=plus(ap[(color*4+s)*bs+v],cm(c.onsite_value(true,site,ch,q,k),pmat[((k%3)*4+ch*2+k/3)*bs+v]));}gemm('C',pmat.data(),ap.data(),out,bs,bs,12,1);}
    }});a.finish();a.setup_ms=seconds(start)*1000;
}

inline void project_matrix(BlockMatrix& out,const BlockMatrix& a,const Transfer& h){auto start=Clock::now();std::vector<std::vector<uint32_t>> graph(out.rows);for(uint32_t i=0;i<a.rows;++i)for(uint32_t e=a.row[i];e<a.row[i+1];++e)graph[h.aggregate[i]].push_back(h.aggregate[a.col[e]]);out.pattern(graph);
    chunks(out.rows,out.workers,[&](int,uint32_t first,uint32_t last){std::vector<D> left(a.bs*out.bs),right(left.size()),product(left.size());for(uint32_t ag=first;ag<last;++ag)for(uint32_t p=h.offsets[ag];p<h.offsets[ag+1];++p){uint32_t i=h.sites[p];h.local_matrix(i,left.data());for(uint32_t e=a.row[i];e<a.row[i+1];++e){uint32_t j=a.col[e];h.local_matrix(j,right.data());gemm('N',a.a64.data()+size_t(e)*a.bs*a.bs,right.data(),product.data(),a.bs,out.bs,a.bs);gemm('C',left.data(),product.data(),out.a64.data()+out.position(ag,h.aggregate[j])*out.bs*out.bs,out.bs,out.bs,a.bs,1);}}});out.finish();out.setup_ms=seconds(start)*1000;
}

struct CPUBlockOperator {const BlockMatrix& a;uint32_t rhs;void apply(const Cx<float>* x,Cx<float>* y){a.apply(x,y,rhs);}};
struct GPUBlockOperator {GPUBlockCoarseSolve<BlockMatrix>& a;MetalRun& m;GPUBlockOperator(GPUBlockCoarseSolve<BlockMatrix>& a):a(a),m(a.m){}void encode(id<MTLComputeCommandEncoder> e,id<MTLBuffer> x,id<MTLBuffer> y){a.encode_product(e,x,y);}void difference(id<MTLComputeCommandEncoder> e,id<MTLBuffer> b,id<MTLBuffer> ax,id<MTLBuffer> out){a.difference(e,b,ax,out);}};

struct CPUCycle {
    const Transfer& h;uint32_t rhs;int workers;bool enabled=true;uint64_t calls=0,products=0;
    std::vector<Cx<float>> residual,work,smooth,b,x;
    std::function<void(const Cx<float>*,Cx<float>*)> product,diagonal,coarse;
    CPUCycle(const Transfer& h,uint32_t rhs,int workers):h(h),rhs(rhs),workers(workers),residual(size_t(h.rows)*h.inbs*rhs),work(residual.size()),smooth(residual.size()),b(size_t(h.coarse)*rhs),x(b.size()){}
    template<class F> void elements(F f){chunks(h.rows*h.inbs,workers,[&](int,uint32_t first,uint32_t last){for(size_t j=size_t(first)*rhs;j<size_t(last)*rhs;++j)f(j);});}
    void apply(const Cx<float>* r,Cx<float>* z){++calls;diagonal(r,z);product(z,work.data());++products;elements([&](size_t j){residual[j]=minus_z(r[j],work[j]);});
        if(enabled){h.restrict_vector(residual.data(),b.data(),rhs);coarse(b.data(),x.data());h.prolong(x.data(),z,rhs,true);product(z,work.data());++products;elements([&](size_t j){residual[j]=minus_z(r[j],work[j]);});}
        diagonal(residual.data(),smooth.data());elements([&](size_t j){z[j]=plus(z[j],smooth[j]);});
    }
};

struct GPUCycle {
    const Transfer& h;MetalRun& m;uint32_t rhs;bool enabled=true;uint64_t calls=0,products=0;
    id<MTLBuffer> basis,offsets,sites,aggregate,inv,residual,work,smooth,b,x;
    std::array<id<MTLComputePipelineState>,6> pipes;
    struct Params{uint32_t rows,inbs,nvec,coarse,mode,spin;float omega;} params;
    std::function<void(id<MTLComputeCommandEncoder>,id<MTLBuffer>,id<MTLBuffer>)> product,coarse;
    std::function<void(id<MTLBuffer>,id<MTLBuffer>)> blocking_coarse;
    std::function<void()> check;
    GPUCycle(const Transfer& h,MetalRun& m,uint32_t rhs,const std::vector<Cx<float>>& inverse,uint32_t format=0):h(h),m(m),rhs(rhs){auto put=[&](const void* p,size_t n){auto b=[m.device newBufferWithLength:std::max(size_t(8),n) options:MTLResourceStorageModeShared];if(!b)throw std::runtime_error("recursive cycle allocation");if(p)std::memcpy([b contents],p,n);return b;};
        basis=put(h.basis32.data(),h.basis32.size()*8);offsets=put(h.offsets.data(),h.offsets.size()*4);sites=put(h.sites.data(),h.sites.size()*4);aggregate=put(h.aggregate.data(),h.aggregate.size()*4);inv=put(inverse.data(),inverse.size()*8);
        size_t n=size_t(h.rows)*h.inbs*rhs*8;residual=put(nullptr,n);work=put(nullptr,n);smooth=put(nullptr,n);b=put(nullptr,size_t(h.coarse)*rhs*8);x=put(nullptr,size_t(h.coarse)*rhs*8);params={h.rows,h.inbs,h.nvec,h.coarse,0,uint32_t(h.spin_layout),1};
        const char* names[]={format==2?"mg3_diagonal_ldl_solve":(format==1?"mg3_diagonal_hermitian":"mg3_diagonal"),"mg3_restrict","mg3_prolong","mg3_difference","mg3_add","mg3_copy"};auto values=[MTLFunctionConstantValues new];[values setConstantValue:&rhs type:MTLDataTypeUInt atIndex:0];for(int i=0;i<6;++i){NSError* e=nil;auto f=[m.library newFunctionWithName:[NSString stringWithUTF8String:names[i]] constantValues:values error:&e];pipes[i]=[m.device newComputePipelineStateWithFunction:f error:&e];if(!pipes[i])throw std::runtime_error([[e localizedDescription] UTF8String]);}
    }
    void dispatch(id<MTLComputeCommandEncoder> enc,int pipe,std::initializer_list<id<MTLBuffer>> buffers,size_t n){[enc setComputePipelineState:pipes[pipe]];size_t i=0;for(auto b:buffers){[enc setBuffer:b offset:0 atIndex:i];++i;}[enc setBytes:&params length:sizeof(params) atIndex:i];[enc dispatchThreads:MTLSizeMake(n,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];[enc memoryBarrierWithScope:MTLBarrierScopeBuffers];}
    void before(id<MTLComputeCommandEncoder> enc,id<MTLBuffer> r,id<MTLBuffer> z){++calls;size_t n=size_t(h.rows)*h.inbs*rhs;dispatch(enc,0,{inv,r,z},n);product(enc,z,work);++products;dispatch(enc,3,{r,work,residual},n);if(enabled)dispatch(enc,1,{residual,basis,offsets,sites,b},size_t(h.coarse)*rhs);}
    void after(id<MTLComputeCommandEncoder> enc,id<MTLBuffer> r,id<MTLBuffer> z){size_t n=size_t(h.rows)*h.inbs*rhs;if(enabled){dispatch(enc,2,{x,basis,aggregate,z},n);product(enc,z,work);++products;dispatch(enc,3,{r,work,residual},n);}dispatch(enc,0,{inv,residual,smooth},n);dispatch(enc,4,{smooth,z},n);}
    void encode(id<MTLComputeCommandEncoder> enc,id<MTLBuffer> r,id<MTLBuffer> z){before(enc,r,z);if(enabled)coarse(enc,b,x);after(enc,r,z);}
    void finish(id<MTLCommandBuffer> cmd,id<MTLComputeCommandEncoder> enc){[enc endEncoding];[cmd commit];[cmd waitUntilCompleted];if(cmd.status==MTLCommandBufferStatusError)throw std::runtime_error([[cmd.error localizedDescription] UTF8String]);}
    void apply_device(id<MTLBuffer> r,id<MTLBuffer> z){@autoreleasepool{auto cmd=[m.queue commandBuffer];auto enc=[cmd computeCommandEncoder];before(enc,r,z);finish(cmd,enc);blocking_coarse(b,x);cmd=[m.queue commandBuffer];enc=[cmd computeCommandEncoder];after(enc,r,z);finish(cmd,enc);}}
    void copy(id<MTLBuffer> a,id<MTLBuffer> b,size_t n){@autoreleasepool{auto cmd=[m.queue commandBuffer];auto enc=[cmd computeCommandEncoder];dispatch(enc,5,{a,b},n);finish(cmd,enc);}}
    void check_coarse(){if(check)check();}
};
