#include "sparse_rhs.hpp"

using Vec=std::vector<D>;
struct SpectralProbe {
    CloverSystem c;int workers;double mass=0;Vec temporary;uint64_t applications=0;
    SpectralProbe(uint32_t l,int w,double strength,bool apbc):c(l,1,0,strength,1,apbc,w),workers(w),temporary(c.w.a.p.outputs()){
        std::vector<Z>().swap(c.w.a.p.x);std::vector<Z>().swap(c.parity.eo.p.x);std::vector<Z>().swap(c.parity.oe.p.x);
    }
    void gamma5(Vec& x){chunks(c.w.a.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t j=first;j<last;++j)if(j%4>=2){x[j].re=-x[j].re;x[j].im=-x[j].im;}});}
    void hermitian(const Vec& x,Vec& y){apply_cpu(c.w.a,x.data(),y.data(),workers);++applications;
        chunks(c.w.a.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t j=first;j<last;++j){y[j].re+=mass*x[j].re;y[j].im+=mass*x[j].im;if(j%4>=2){y[j].re=-y[j].re;y[j].im=-y[j].im;}}});}
    void square(const Vec& x,Vec& y){hermitian(x,temporary);hermitian(temporary,y);}
    CD dot(const Vec& x,const Vec& y){std::vector<CD> part(workers);chunks(c.w.a.p.rows*12,workers,[&](int w,uint32_t first,uint32_t last){CD sum=0;for(size_t j=first;j<last;++j)sum+=std::conj(CD(x[j].re,x[j].im))*CD(y[j].re,y[j].im);part[w]=sum;});CD sum=0;for(auto v:part)sum+=v;return sum;}
    void axpy(Vec& x,const Vec& y,CD alpha){chunks(c.w.a.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t j=first;j<last;++j){auto z=alpha*CD(y[j].re,y[j].im);x[j].re+=z.real();x[j].im+=z.imag();}});}
    void scale(Vec& x,double alpha){chunks(c.w.a.p.rows*12,workers,[&](int,uint32_t first,uint32_t last){for(size_t j=first;j<last;++j){x[j].re*=alpha;x[j].im*=alpha;}});}
};
struct Ritz {double lambda=0,bound=0;std::vector<double> weights;};
Ritz smallest(const std::vector<double>& alpha,const std::vector<double>& beta){
    __LAPACK_int n=__LAPACK_int(alpha.size()),ld=n,info=0;char vectors='V';std::vector<double> d=alpha,e(beta.begin(),beta.begin()+n-1),z(size_t(n)*n),work(std::max(1,2*int(n)-2));
    e.resize(n);dstev_(&vectors,&n,d.data(),e.data(),z.data(),&ld,work.data(),&info);if(info)throw std::runtime_error("tridiagonal eigensolver failed");
    return {d[0],std::abs(beta.back()*z[n-1]),std::vector<double>(z.begin(),z.begin()+n)};
}
struct Gap {double mass,sigma,lambda,ritz_bound,true_residual,m_stationary,ms;int steps;uint64_t applications;};
Gap estimate_gap(SpectralProbe& p,double mass,int limit,double tolerance,int seed,bool free_check,bool smooth_start){p.mass=mass;p.applications=0;auto start=Clock::now();size_t n=p.c.w.a.p.outputs();Vec v(n),w(n),best(n),check(n);std::mt19937_64 rng(seed);std::normal_distribution<double> normal(0,1);
    std::array<CD,24> constants{};if(smooth_start)for(auto& z:constants)z={normal(rng),normal(rng)};
    for(size_t j=0;j<n;++j){CD z(normal(rng),normal(rng));if(smooth_start){double theta=std::acos(-1.)*(j/12/(p.c.w.extent*p.c.w.extent*p.c.w.extent))/p.c.w.extent;z=.02*z+std::cos(theta)*constants[j%12]+std::sin(theta)*constants[12+j%12];}v[j]={z.real(),z.imag()};}p.scale(v,1/std::sqrt(p.dot(v,v).real()));
    std::vector<Vec> basis;basis.reserve(limit);std::vector<double> alpha,beta;Ritz r;
    for(int k=0;k<limit;++k){basis.push_back(v);p.square(v,w);auto a=p.dot(v,w);if(std::abs(a.imag())>1e-10*std::max(1.,std::abs(a.real())))throw std::runtime_error("H^2 not Hermitian");alpha.push_back(a.real());p.axpy(w,v,-a.real());if(k)p.axpy(w,basis[k-1],-beta[k-1]);
        // Full two-pass reorthogonalization; no unverified three-term-only gap estimate.
        for(int pass=0;pass<2;++pass)for(auto& b:basis){auto h=p.dot(b,w);p.axpy(w,b,-h);}
        double b=std::sqrt(std::max(0.,p.dot(w,w).real()));beta.push_back(b);
        if((k+1)%25==0||k+1==limit||b<1e-14){r=smallest(alpha,beta);std::cout<<"m="<<mass<<" k="<<k+1<<" lambda="<<r.lambda<<" ritz_bound="<<r.bound<<std::endl;
            if(r.bound<tolerance||b<1e-14)break;}
        if(!(b>0)||!std::isfinite(b))throw std::runtime_error("Lanczos breakdown");v=w;p.scale(v,1/b);
    }
    std::fill(best.begin(),best.end(),D{});for(size_t k=0;k<r.weights.size();++k)p.axpy(best,basis[k],r.weights[k]);p.scale(best,1/std::sqrt(p.dot(best,best).real()));p.square(best,check);double lambda=p.dot(best,check).real();p.axpy(check,best,-lambda);double residual=std::sqrt(p.dot(check,check).real());p.hermitian(best,check);double sigma=std::sqrt(p.dot(check,check).real());p.gamma5(check);double stationary=mass-p.dot(best,check).real();
    if(!std::isfinite(residual)||!(lambda>0)||std::abs(sigma*sigma-lambda)>1e-10)throw std::runtime_error("invalid singular gap estimate");
    if(free_check){double expected=1e100;uint32_t l=p.c.w.extent;for(uint32_t x=0;x<l*l*l*l;++x){uint32_t q=x;double a=mass,b=0;for(int mu=0;mu<4;++mu){double angle=2*std::acos(-1.)*(q%l+(mu==3?.5:0.))/l;q/=l;a+=1-std::cos(angle);b+=std::sin(angle)*std::sin(angle);}expected=std::min(expected,std::sqrt(a*a+b));}if(std::abs(sigma-expected)>1e-7)throw std::runtime_error("free-field spectral gap mismatch");}
    std::cout<<"GAP mass="<<mass<<" sigma="<<sigma<<" true_H2_residual="<<residual<<" stationary_m="<<stationary<<" steps="<<alpha.size()<<std::endl;
    return {mass,sigma,lambda,r.bound,residual,stationary,seconds(start)*1000,int(alpha.size()),p.applications};
}
int main(int argc,char** argv){@autoreleasepool{try{uint32_t l=16;int workers=8,limit=250,seed=12371;double tolerance=1e-7;std::vector<double> masses{.1,0,-.1,-.2};bool self=false,smooth_start=false;std::string csv;
    for(int i=1;i<argc;++i){std::string a=argv[i];auto next=[&](){if(++i>=argc)throw std::runtime_error("missing argument");return std::string(argv[i]);};if(a=="--size")l=uint32_t(std::stoul(next()));else if(a=="--steps")limit=std::stoi(next());else if(a=="--seed")seed=std::stoi(next());else if(a=="--tol")tolerance=std::stod(next());else if(a=="--mass"){masses.clear();std::stringstream in(next());std::string value;while(std::getline(in,value,','))masses.push_back(std::stod(value));}else if(a=="--csv")csv=next();else if(a=="--self-test")self=true;else if(a=="--smooth-start")smooth_start=true;else throw std::runtime_error("unknown option");}
    if(self){l=4;limit=100;masses={.1,-.1,-.25};}if(l<2||l>16||l%2||limit<25||limit>500||!(tolerance>0&&tolerance<1e-3))throw std::runtime_error("invalid options");for(auto m:masses)if(!(m>-1&&m<1))throw std::runtime_error("invalid mass");SpectralProbe probe(l,workers,self?0:.25,true);
    std::ofstream out;if(!csv.empty()){out.open(csv);out<<"extent,mass,sigma_min_estimate,lambda_min_estimate,ritz_bound,true_H2_residual,stationary_mass,steps,operator_applications,probe_ms,seed,smooth_start,converged\n";}
    for(auto m:masses){auto g=estimate_gap(probe,m,limit,tolerance,seed,self,smooth_start);if(out){out<<std::setprecision(14)<<l<<','<<m<<','<<g.sigma<<','<<g.lambda<<','<<g.ritz_bound<<','<<g.true_residual<<','<<g.m_stationary<<','<<g.steps<<','<<g.applications<<','<<g.ms<<','<<seed<<','<<smooth_start<<','<<(g.true_residual<=tolerance)<<'\n';out.flush();}}
    std::cout<<"DONE: FP64 Rayleigh diagnostics; convergence is recorded explicitly."<<std::endl;return 0;
}catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}}}
