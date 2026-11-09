#pragma once
#include "clover_common.hpp"

// Independent transcription of Chroma QDPCloverEnv::makeClovSiteLoop.
std::array<M6,2> chroma_clover_blocks(const CloverSystem& c,uint32_t site){
    std::array<M6,2> a{};std::array<M3,6> f{};CD I(0,1);
    for(int p=0;p<6;++p)for(int j=0;j<9;++j)f[p][j]=.5*c.coefficient*c.field[site*6+p][j];
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){
        CD tl[2]={I*(f[0][i*3+j]-f[5][i*3+j]),I*(f[5][i*3+j]+f[0][i*3+j])};
        CD lo[2]={I*(f[3][i*3+j]-f[2][i*3+j])-f[1][i*3+j]-f[4][i*3+j],
                  I*(f[2][i*3+j]+f[3][i*3+j])+f[4][i*3+j]-f[1][i*3+j]};
        for(int ch=0;ch<2;++ch){a[ch][i*6+j]=tl[ch];a[ch][(i+3)*6+j+3]=-tl[ch];a[ch][(i+3)*6+j]=lo[ch];a[ch][j*6+i+3]=std::conj(lo[ch]);}
    }
    for(auto& b:a)for(int q=0;q<6;++q)b[q*6+q]+=c.w.a.mass+4;return a;
}

M6 gaussian_inverse(M6 a){
    M6 b{};for(int i=0;i<6;++i)b[i*6+i]=1;
    for(int k=0;k<6;++k){int pivot=k;for(int i=k+1;i<6;++i)if(std::abs(a[i*6+k])>std::abs(a[pivot*6+k]))pivot=i;
        if(std::abs(a[pivot*6+k])<1e-14)throw std::runtime_error("reference clover inverse singular");
        for(int j=0;j<6;++j){std::swap(a[k*6+j],a[pivot*6+j]);std::swap(b[k*6+j],b[pivot*6+j]);}
        CD d=a[k*6+k];for(int j=0;j<6;++j){a[k*6+j]/=d;b[k*6+j]/=d;}
        for(int i=0;i<6;++i)if(i!=k){CD z=a[i*6+k];for(int j=0;j<6;++j){a[i*6+j]-=z*a[k*6+j];b[i*6+j]-=z*b[k*6+j];}}
    }return b;
}

void check_abelian_field_normalization(){
    // Quantized uniform flux on a periodic 4x4 plane gives four equal leaves:
    // F_01 = diag(i*sin(theta),-i*sin(theta),0), independently fixing sign/factor.
    auto w=make_wilson(4,1,.25,0,true,1,false,true);double theta=2*std::acos(-1.0)/16;
    for(uint32_t x=0;x<w.a.p.rows;++x){uint32_t cx=x%4,cy=(x/4)%4;
        for(int mu=0;mu<4;++mu){auto& u=w.links[x*4+mu];u=m3_identity();double phase=mu==0?-theta*cy:mu==1&&cy==3?theta*4*cx:0;
            u[0]=std::polar(1.0,phase);u[4]=std::polar(1.0,-phase);}
    }
    auto f=clover_field(w);
    for(uint32_t x=0;x<w.a.p.rows;++x)for(int p=0;p<6;++p)for(int j=0;j<9;++j){CD expected=0;
        if(p==0&&j==0)expected=CD(0,std::sin(theta));if(p==0&&j==4)expected=CD(0,-std::sin(theta));
        if(std::abs(f[x*6+p][j]-expected)>3e-13)throw std::runtime_error("analytic clover field normalization failed");}
}

std::vector<D> direct_clover(const CloverSystem& c,const std::vector<D>& x){
    auto out=direct_wilson(c.w,x);uint32_t rhs=c.w.a.p.rhs;
    for(uint32_t site=0;site<c.w.a.p.rows;++site){auto a=chroma_clover_blocks(c,site);
        for(int ch=0;ch<2;++ch)for(int q=0;q<6;++q)for(uint32_t r=0;r<rhs;++r){CD sum=0;
            for(int k=0;k<6;++k){auto v=x[((size_t(site)*3+k%3)*4+ch*2+k/3)*rhs+r];CD z=a[ch][q*6+k];if(q==k)z-=c.w.a.mass+4;sum+=z*CD(v.re,v.im);}
            auto& y=out[((size_t(site)*3+q%3)*4+ch*2+q/3)*rhs+r];y.re+=sum.real();y.im+=sum.imag();
        }
    }return out;
}

void check_clover(const CloverSystem& c,MetalRun& metal,int workers){
    static bool checked_flux=false;if(!checked_flux){check_abelian_field_normalization();checked_flux=true;}
    const auto& e=c.parity;const auto& a=c.w.a;double maxerr=0,nontrivial=0;uint32_t width=12*a.p.rhs;
    for(uint32_t site=0;site<a.p.rows;++site){auto ref=chroma_clover_blocks(c,site);
        for(int ch=0;ch<2;++ch)for(int i=0;i<6;++i)for(int j=0;j<6;++j){auto v=c.fullA[(size_t(site)*2+ch)*36+i*6+j];
            maxerr=std::max(maxerr,std::abs(CD(v.re,v.im)-ref[ch][i*6+j]));
            if(std::abs(ref[ch][i*6+j]-std::conj(ref[ch][j*6+i]))>1e-12)throw std::runtime_error("clover block Hermiticity failed");
            auto z=a.local64[(size_t(site)*2+ch)*36+i*6+j];nontrivial=std::max(nontrivial,std::hypot(z.re,z.im));
        }
    }
    if(maxerr>1e-12)throw std::runtime_error("clover sigma/F convention differs from Chroma component formulas");
    if(c.w.strength>0&&c.coefficient>0&&nontrivial<1e-5)throw std::runtime_error("nonzero clover expected");
    for(size_t i=0;i<e.even.size();++i)for(int ch=0;ch<2;++ch){auto ref=gaussian_inverse(chroma_clover_blocks(c,e.even[i])[ch]);
        for(int j=0;j<36;++j){auto v=c.inverseEven[(i*2+ch)*36+j];if(std::abs(CD(v.re,v.im)-ref[j])>1e-12)throw std::runtime_error("6x6 inverse check failed");}
    }
    auto direct=direct_clover(c,a.truth);maxerr=0;
    for(size_t j=0;j<direct.size();++j)maxerr=std::max(maxerr,std::hypot(direct[j].re-a.b[j].re,direct[j].im-a.b[j].im));
    if(maxerr>2e-13)throw std::runtime_error("full clover FP64 oracle mismatch");
    // Local gauge covariance checks all four loop orientations, site addressing, and adjoints.
    auto transformed=c.w;std::vector<M3> rotation(a.p.rows);for(uint32_t x=0;x<a.p.rows;++x)rotation[x]=c.w.links[x*4];
    for(uint32_t x=0;x<a.p.rows;++x)for(uint32_t mu=0;mu<4;++mu)transformed.links[x*4+mu]=m3_mul(m3_mul(rotation[x],c.w.links[x*4+mu]),m3_adj(rotation[c.w.neighbor(x,mu,true)]));
    auto ff=clover_field(transformed);
    for(uint32_t x=0;x<a.p.rows;++x)for(int p=0;p<6;++p){auto ref=m3_mul(m3_mul(rotation[x],c.field[x*6+p]),m3_adj(rotation[x]));
        for(int j=0;j<9;++j)if(std::abs(ff[x*6+p][j]-ref[j])>2e-13)throw std::runtime_error("clover field gauge covariance failed");}
    auto schur_reference=[&](const std::vector<D>& odd){
        std::vector<D> full(a.p.outputs()),even_only(full.size());
        for(size_t i=0;i<e.odd.size();++i)for(uint32_t k=0;k<width;++k)full[size_t(e.odd[i])*width+k]=odd[i*width+k];
        auto first=direct_clover(c,full);
        for(size_t i=0;i<e.even.size();++i){auto ref=chroma_clover_blocks(c,e.even[i]);
            for(int ch=0;ch<2;++ch){auto inv=gaussian_inverse(ref[ch]);
                for(int q=0;q<6;++q)for(uint32_t r=0;r<a.p.rhs;++r){CD sum=0;
                    for(int k=0;k<6;++k){auto v=first[((size_t(e.even[i])*3+k%3)*4+ch*2+k/3)*a.p.rhs+r];sum+=inv[q*6+k]*CD(v.re,v.im);}
                    even_only[((size_t(e.even[i])*3+q%3)*4+ch*2+q/3)*a.p.rhs+r]={sum.real(),sum.imag()};
                }
            }
        }
        auto second=direct_clover(c,even_only);std::vector<D> out(odd.size());
        for(size_t i=0;i<e.odd.size();++i)for(uint32_t k=0;k<width;++k){size_t j=size_t(e.odd[i])*width+k;out[i*width+k]=minus_z(first[j],second[j]);}
        return out;
    };
    std::vector<D> xo(e.oe.p.outputs()),ys(xo.size());for(size_t i=0;i<e.odd.size();++i)for(uint32_t k=0;k<width;++k)xo[i*width+k]=a.truth[size_t(e.odd[i])*width+k];
    auto expected=schur_reference(xo);CPUCloverSchur<double> cpu(c,workers);cpu.apply(xo.data(),ys.data());maxerr=0;
    for(size_t j=0;j<ys.size();++j)maxerr=std::max(maxerr,std::hypot(ys[j].re-expected[j].re,ys[j].im-expected[j].im));
    if(maxerr>3e-13)throw std::runtime_error("odd clover Schur composition mismatch");
    GPUCloverSchur gpu(metal,c);auto in=gpu.eo.buffers[3],out=gpu.eo.buffers[4];auto* xx=static_cast<Cx<float>*>([in contents]);auto* yy=static_cast<Cx<float>*>([out contents]);
    for(size_t j=0;j<xo.size();++j){xx[j]={float(xo[j].re),float(xo[j].im)};xo[j]={xx[j].re,xx[j].im};}
    expected=schur_reference(xo);auto cmd=[metal.queue commandBuffer];auto enc=[cmd computeCommandEncoder];gpu.encode(enc,in,out);[enc endEncoding];[cmd commit];[cmd waitUntilCompleted];
    if(cmd.status==MTLCommandBufferStatusError)throw std::runtime_error([[cmd.error localizedDescription] UTF8String]);maxerr=0;
    for(size_t j=0;j<xo.size();++j)maxerr=std::max(maxerr,std::hypot(yy[j].re-expected[j].re,yy[j].im-expected[j].im));
    if(maxerr>3e-5)throw std::runtime_error("Metal clover Schur versus direct oracle mismatch");
    std::cout<<"Clover convention/two-6x6/inverse/gauge-covariance/full-operator/odd-Schur CPU+Metal checks PASS L="<<c.w.extent<<" RHS="<<a.p.rhs<<" csw="<<c.coefficient<<std::endl;
}

void check_device_elimination(const CloverSystem& c,GPUCloverSchur& op,GPUBiCGStab& inner,int workers){
    DeviceCloverCorrection device(c,op,inner,workers);device.prepare(c.w.a.b);const auto& e=c.parity;uint32_t width=12*c.w.a.p.rhs;
    std::vector<D> be(e.eo.p.outputs()),bo(be.size()),temp(be.size()),aux(be.size());
    for(size_t i=0;i<e.even.size();++i)for(uint32_t k=0;k<width;++k){be[i*width+k]=c.w.a.b[size_t(e.even[i])*width+k];bo[i*width+k]=c.w.a.b[size_t(e.odd[i])*width+k];}
    local_blocks(c.inverseEven,be.data(),temp.data(),e.eo.p.rows,e.eo.p.rhs,workers);apply_cpu(e.oe,temp.data(),aux.data(),workers);
    auto* g=static_cast<Cx<float>*>([device.reduced contents]);double err=0;
    for(size_t j=0;j<be.size();++j)err=std::max(err,std::hypot(g[j].re-(bo[j].re-aux[j].re),g[j].im-(bo[j].im-aux[j].im)));
    if(err>3e-5)throw std::runtime_error("device clover RHS reduction mismatch");
    auto* dodd=static_cast<Cx<float>*>([inner.delta contents]);
    for(size_t i=0;i<e.odd.size();++i)for(uint32_t k=0;k<width;++k){auto z=c.w.a.truth[size_t(e.odd[i])*width+k];dodd[i*width+k]={float(z.re),float(z.im)};}
    device.reconstruct();auto* result=device.correction();err=0;
    for(size_t j=0;j<c.w.a.truth.size();++j)err=std::max(err,std::hypot(result[j].re-c.w.a.truth[j].re,result[j].im-c.w.a.truth[j].im));
    if(err>3e-5)throw std::runtime_error("device clover even reconstruction/scatter mismatch");
    std::cout<<"Metal FP32 clover RHS reduction and even reconstruction checks PASS"<<std::endl;
}
