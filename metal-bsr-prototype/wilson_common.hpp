#pragma once
#include "bicgstab.hpp"
#include "lattice_geometry.hpp"
using CD=std::complex<double>;
using M3=std::array<CD,9>;
using M4=std::array<CD,16>;

std::array<M4,4> gamma_matrices(){
    std::array<M4,4> g{};CD I(0,1);
    g[0][3]=-I;g[0][6]=-I;g[0][9]=I;g[0][12]=I;
    g[1][3]=-1;g[1][6]=1;g[1][9]=1;g[1][12]=-1;
    g[2][2]=-I;g[2][7]=I;g[2][8]=I;g[2][13]=-I;
    g[3][2]=g[3][7]=g[3][8]=g[3][13]=1;
    return g;
}

struct WilsonSystem{
    System a;uint32_t extent;double strength;bool antiperiodic;
    LatticeShape dimensions{};
    uint32_t length(uint32_t mu)const{return dimensions[mu]?dimensions[mu]:extent;}
    uint32_t stride(uint32_t mu)const{uint32_t s=1;for(uint32_t j=0;j<mu;++j)s*=length(j);return s;}
    std::vector<M3> links;const M3* mapped_links=nullptr;const M3& link_at(size_t i)const{return mapped_links?mapped_links[i]:links[i];}
    std::array<M4,4> gamma=gamma_matrices();
    uint32_t neighbor(uint32_t i,uint32_t mu,bool forward)const{
        uint32_t step=stride(mu),n=length(mu),c=(i/step)%n;
        int delta=forward?(c+1==n?-int(n-1):1):(c==0?int(n-1):-1);
        return uint32_t(int64_t(i)+int64_t(delta)*step);
    }
    int bc(uint32_t i,uint32_t mu,bool forward)const{
        if(!antiperiodic||mu!=3)return 1;uint32_t c=i/stride(3);
        return ((forward&&c==length(3)-1)||(!forward&&c==0))?-1:1;
    }
    CD link(uint32_t i,uint32_t mu,bool forward,uint32_t c,uint32_t k)const{
        return forward?link_at(i*4+mu)[c*3+k]:std::conj(link_at(neighbor(i,mu,false)*4+mu)[k*3+c]);
    }
};

WilsonSystem make_wilson(uint32_t l,uint32_t rhs,double mass,double strength,bool apbc,int workers,bool zero_rhs=false,bool qdp_basis=false){
    WilsonSystem w;w.extent=l;w.strength=strength;w.antiperiodic=apbc;
    if(qdp_basis)for(auto& g:w.gamma)for(auto& z:g)z=std::conj(z); // QDP Gamma(1,2,4,8)
    auto& a=w.a;auto& p=a.p;a.mass=mass;p.rows=l*l*l*l;p.rhs=rhs;p.colors=3;p.spins=4;
    std::mt19937_64 rng(77347);std::uniform_real_distribution<double> angle(-strength,strength),phase(-std::acos(-1.0),std::acos(-1.0));
    w.links.resize(p.rows*4);
    for(auto& u:w.links){
        u={};u[0]=u[4]=u[8]=1;
        for(auto pair: {std::pair<int,int>{0,1},{1,2},{0,2}}){
            double theta=angle(rng),phi=phase(rng);CD s=std::sin(theta)*std::polar(1.0,phi);double c=std::cos(theta);
            for(int k=0;k<3;++k){CD x=u[pair.first*3+k],y=u[pair.second*3+k];u[pair.first*3+k]=c*x+s*y;u[pair.second*3+k]=-std::conj(s)*x+c*y;}
        }
        double t0=angle(rng),t1=angle(rng);CD d[3]={std::polar(1.0,t0),std::polar(1.0,t1),std::polar(1.0,-t0-t1)};
        for(int c=0;c<3;++c)for(int k=0;k<3;++k)u[c*3+k]*=d[c];
    }
    p.row.push_back(0);
    for(uint32_t i=0;i<p.rows;++i){
        auto append=[&](uint32_t col,const M3& color,const std::array<uint32_t,4>& perm,const std::array<CD,4>& ph){
            p.col.push_back(col);for(auto z:color){a.ad.push_back({z.real(),z.imag()});a.af.push_back({float(z.real()),float(z.imag())});p.a.push_back({float(z.real()),float(z.imag())});}
            for(int s=0;s<4;++s){p.perm.push_back(perm[s]);p.phase.push_back({float(ph[s].real()),float(ph[s].imag())});}
        };
        M3 diagonal{};diagonal[0]=diagonal[4]=diagonal[8]=mass+4;
        append(i,diagonal,{0,1,2,3},{1,1,1,1});
        for(uint32_t mu=0;mu<4;++mu)for(bool forward:{true,false}){
            M3 color;for(uint32_t c=0;c<3;++c)for(uint32_t k=0;k<3;++k)color[c*3+k]=-0.5*w.link(i,mu,forward,c,k);
            int boundary=w.bc(i,mu,forward);std::array<CD,4> identity{CD(boundary),CD(boundary),CD(boundary),CD(boundary)};
            append(w.neighbor(i,mu,forward),color,{0,1,2,3},identity);
            std::array<uint32_t,4> perm{};std::array<CD,4> ph{};
            for(uint32_t s=0;s<4;++s)for(uint32_t si=0;si<4;++si)if(std::abs(w.gamma[mu][s*4+si])>0){
                perm[s]=si;ph[s]=double(boundary*(forward?-1:1))*w.gamma[mu][s*4+si];
            }
            append(w.neighbor(i,mu,forward),color,perm,ph);
        }
        p.row.push_back(uint32_t(p.col.size()));
    }
    std::uniform_real_distribution<double> values(-1,1);p.x.resize(p.outputs());a.truth.resize(p.outputs());a.b.resize(p.outputs());
    for(size_t j=0;j<p.outputs();++j){a.truth[j]={values(rng),values(rng)};p.x[j]={float(values(rng)),float(values(rng))};}
    if(zero_rhs&&rhs>1)for(size_t j=0;j<p.rows*12;++j)a.truth[j*rhs+rhs-1]={};
    apply_cpu(a,a.truth.data(),a.b.data(),workers);a.bnorm.assign(rhs,0);
    for(size_t j=0;j<p.outputs();++j)a.bnorm[j%rhs]+=a.b[j].re*a.b[j].re+a.b[j].im*a.b[j].im;
    for(auto& x:a.bnorm)x=std::sqrt(x);
    return w;
}

// Direct four-spin projector formula, independent of the 17-term encoding.
std::vector<D> direct_wilson(const WilsonSystem& w,const std::vector<D>& x){
    const auto& p=w.a.p;std::vector<D> out(x.size());
    auto read=[&](uint32_t i,uint32_t c,uint32_t s,uint32_t r){auto z=x[((size_t(i)*3+c)*4+s)*p.rhs+r];return CD(z.re,z.im);};
    for(uint32_t i=0;i<p.rows;++i)for(uint32_t c=0;c<3;++c)for(uint32_t s=0;s<4;++s)for(uint32_t r=0;r<p.rhs;++r){
        CD z=(w.a.mass+4)*read(i,c,s,r);
        for(uint32_t mu=0;mu<4;++mu)for(bool f:{true,false})for(uint32_t k=0;k<3;++k){
            uint32_t n=w.neighbor(i,mu,f);CD proj=read(n,k,s,r);
            for(uint32_t si=0;si<4;++si)proj+=double(f?-1:1)*w.gamma[mu][s*4+si]*read(n,k,si,r);
            z-=0.5*double(w.bc(i,mu,f))*w.link(i,mu,f,c,k)*proj;
        }
        out[((size_t(i)*3+c)*4+s)*p.rhs+r]={z.real(),z.imag()};
    }
    return out;
}

void validate_wilson(const WilsonSystem& w,MetalRun& metal,int workers){
    const auto& a=w.a;const auto& p=a.p;
    for(int mu=0;mu<4;++mu)for(int nu=0;nu<4;++nu)for(int s=0;s<4;++s)for(int t=0;t<4;++t){
        CD anticom=0;for(int k=0;k<4;++k)anticom+=w.gamma[mu][s*4+k]*w.gamma[nu][k*4+t]+w.gamma[nu][s*4+k]*w.gamma[mu][k*4+t];
        if(std::abs(anticom-CD(mu==nu&&s==t?2:0))>1e-13||std::abs(w.gamma[mu][s*4+t]-std::conj(w.gamma[mu][t*4+s]))>1e-13)
            throw std::runtime_error("gamma algebra failed");
    }
    for(const auto& u:w.links){
        CD det=u[0]*(u[4]*u[8]-u[5]*u[7])-u[1]*(u[3]*u[8]-u[5]*u[6])+u[2]*(u[3]*u[7]-u[4]*u[6]);
        if(std::abs(det-CD(1))>1e-12)throw std::runtime_error("SU3 determinant failed");
        for(int i=0;i<3;++i)for(int j=0;j<3;++j){CD z=0;for(int k=0;k<3;++k)z+=std::conj(u[k*3+i])*u[k*3+j];
            if(std::abs(z-CD(i==j?1:0))>1e-12)throw std::runtime_error("SU3 unitarity failed");}
    }
    auto oracle=direct_wilson(w,a.truth);double worst=0;
    for(size_t j=0;j<oracle.size();++j)worst=std::max(worst,std::hypot(oracle[j].re-a.b[j].re,oracle[j].im-a.b[j].im));
    if(worst>2e-13)throw std::runtime_error("direct Wilson versus encoded FP64 mismatch");
    // Build the small encoded matrix to test gamma5 Hermiticity entry by entry.
    size_t dim=p.rows*12;std::vector<CD> dense(dim*dim);
    for(uint32_t i=0;i<p.rows;++i)for(uint32_t e=p.row[i];e<p.row[i+1];++e)
    for(uint32_t c=0;c<3;++c)for(uint32_t s=0;s<4;++s)for(uint32_t k=0;k<3;++k){
        auto av=a.ad[e*9+c*3+k];auto ph=p.phase[e*4+s];
        dense[(i*12+c*4+s)*dim+p.col[e]*12+k*4+p.perm[e*4+s]]+=CD(av.re,av.im)*CD(ph.re,ph.im);
    }
    double nonhermitian=0;
    for(size_t i=0;i<dim;++i)for(size_t j=0;j<dim;++j){int sign=((i%4<2)==(j%4<2))?1:-1;
        if(std::abs(dense[i*dim+j]-double(sign)*std::conj(dense[j*dim+i]))>1e-12)throw std::runtime_error("gamma5 Hermiticity failed");
        nonhermitian=std::max(nonhermitian,std::abs(dense[i*dim+j]-std::conj(dense[j*dim+i])));
    }
    if(w.strength>0&&nonhermitian<1e-3)throw std::runtime_error("Wilson test unexpectedly Hermitian");
    std::vector<D> x(p.outputs());for(size_t j=0;j<x.size();++j)x[j]={p.x[j].re,p.x[j].im};
    auto expected=direct_wilson(w,x);GPUProblem gpu(metal,p);gpu.run(1);auto* y=gpu.output();worst=0;
    for(size_t j=0;j<x.size();++j)worst=std::max(worst,std::hypot(y[j].re-expected[j].re,y[j].im-expected[j].im));
    if(worst>2e-5)throw std::runtime_error("Metal versus direct Wilson mismatch");
    if(w.strength==0){
        // Free-field plane wave: D(p) = m + sum(1-cos p_mu) + i sum gamma_mu sin p_mu.
        const double pi=std::acos(-1.0);double momentum[4]={2*pi/w.extent,0,0,(w.antiperiodic?pi:0)/w.extent};
        CD spin[4]={{1,.2},{-.3,.1},{.4,-.7},{.8,.5}};std::vector<D> wave(p.outputs()),pred(p.outputs());
        for(uint32_t i=0;i<p.rows;++i){uint32_t q=i;double theta=0,diag=a.mass;
            for(int mu=0;mu<4;++mu){theta+=(q%w.extent)*momentum[mu];q/=w.extent;diag+=1-std::cos(momentum[mu]);}
            CD phase=std::polar(1.0,theta);
            for(uint32_t c=0;c<3;++c)for(uint32_t s=0;s<4;++s)for(uint32_t r=0;r<p.rhs;++r){
                CD v=phase*spin[s],z=diag*spin[s];for(int mu=0;mu<4;++mu)for(int si=0;si<4;++si)z+=CD(0,std::sin(momentum[mu]))*w.gamma[mu][s*4+si]*spin[si];
                z*=phase;size_t idx=((size_t(i)*3+c)*4+s)*p.rhs+r;wave[idx]={v.real(),v.imag()};pred[idx]={z.real(),z.imag()};
            }
        }
        std::vector<D> got(p.outputs());apply_cpu(a,wave.data(),got.data(),workers);worst=0;
        for(size_t j=0;j<got.size();++j)worst=std::max(worst,std::hypot(got[j].re-pred[j].re,got[j].im-pred[j].im));
        if(worst>1e-12)throw std::runtime_error("free-field Fourier symbol failed");
    }
    std::cout<<"Wilson structural/direct/Fourier checks PASS L="<<w.extent<<" RHS="<<p.rhs<<" gauge="<<w.strength<<" APBC="<<w.antiperiodic<<std::endl;
}

