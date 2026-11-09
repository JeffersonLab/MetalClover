struct MGParams { uint rows,coarse,nvec,mode; float omega; };
kernel void mg_scale_add(device const float2* x [[buffer(0)]],device float2* y [[buffer(1)]],constant MGParams& p [[buffer(2)]],uint gid [[thread_position_in_grid]]){
    if(gid<p.rows*12*RHS)y[gid]=p.mode?fma(float2(p.omega),x[gid],y[gid]):p.omega*x[gid];
}
kernel void mg_restrict(device const float2* x [[buffer(0)]],device const float2* basis [[buffer(1)]],device const uint* offsets [[buffer(2)]],
                        device const uint* sites [[buffer(3)]],device float2* y [[buffer(4)]],constant MGParams& p [[buffer(5)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.coarse*RHS)return;uint r=gid%RHS,row=gid/RHS,a=row/(2*p.nvec),ch=row%(2*p.nvec)/p.nvec,v=row%p.nvec;float2 sum=0;
    for(uint j=offsets[a];j<offsets[a+1];++j)for(uint c=0;c<3;++c)for(uint s=2*ch;s<2*ch+2;++s){uint q=sites[j]*12+c*4+s;
        float2 z=x[q*RHS+r],w=basis[q*p.nvec+v];sum.x=fma(w.x,z.x,sum.x);sum.x=fma(w.y,z.y,sum.x);sum.y=fma(w.x,z.y,sum.y);sum.y=fma(-w.y,z.x,sum.y);}
    y[gid]=sum;
}
kernel void mg_coarse(device const float2* inverse [[buffer(0)]],device const float2* rhs [[buffer(1)]],device float2* y [[buffer(2)]],
                      constant MGParams& p [[buffer(3)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.coarse*RHS)return;uint row=gid/RHS,r=gid%RHS;float2 sum=0;
    for(uint k=0;k<p.coarse;++k){float2 a=inverse[row*p.coarse+k],b=rhs[k*RHS+r];sum.x=fma(a.x,b.x,sum.x);sum.x=fma(-a.y,b.y,sum.x);sum.y=fma(a.x,b.y,sum.y);sum.y=fma(a.y,b.x,sum.y);}y[gid]=sum;
}
kernel void mg_prolong(device const float2* x [[buffer(0)]],device const float2* basis [[buffer(1)]],device const uint* aggregate [[buffer(2)]],
                       device float2* y [[buffer(3)]],constant MGParams& p [[buffer(4)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.rows*12*RHS)return;uint q=gid/RHS,r=gid%RHS,base=aggregate[q/12]*2*p.nvec+(q%4/2)*p.nvec;float2 sum=0;
    for(uint v=0;v<p.nvec;++v)sum+=cmul(basis[q*p.nvec+v],x[(base+v)*RHS+r]);y[gid]=p.mode?y[gid]+sum:sum;
}
struct GMParams { uint components,parts,k,stage; };
kernel void gm_init(device const float2* b [[buffer(0)]],device float2* saved [[buffer(1)]],device float2* r [[buffer(2)]],device float2* delta [[buffer(3)]],
                    constant GMParams& p [[buffer(4)]],uint gid [[thread_position_in_grid]]){if(gid<p.components*RHS){saved[gid]=r[gid]=b[gid];delta[gid]=0;}}
kernel void gm_store(device const float2* x [[buffer(0)]],device float2* stack [[buffer(1)]],constant GMParams& p [[buffer(2)]],uint gid [[thread_position_in_grid]]){
    if(gid<p.components*RHS){if(p.stage==0)stack[p.k*p.components*RHS+gid]=x[gid];}}
kernel void gm_extract(device const float2* stack [[buffer(0)]],device float2* x [[buffer(1)]],constant GMParams& p [[buffer(2)]],uint gid [[thread_position_in_grid]]){
    if(gid<p.components*RHS)x[gid]=stack[p.k*p.components*RHS+gid];}
kernel void gm_normalize(device const float2* x [[buffer(0)]],device float2* stack [[buffer(1)]],device const float* factor [[buffer(2)]],
                         constant GMParams& p [[buffer(3)]],uint gid [[thread_position_in_grid]]){
    if(gid<p.components*RHS)stack[p.k*p.components*RHS+gid]=x[gid]*factor[gid%RHS];}
kernel void gm_dots(device const float2* stack [[buffer(0)]],device const float2* w [[buffer(1)]],device float2* partial [[buffer(2)]],
                    constant GMParams& p [[buffer(3)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.k*p.parts*RHS)return;uint r=gid%RHS,part=gid/RHS%p.parts,k=gid/(RHS*p.parts);float2 sum=0;
    for(uint j=part*256;j<min((part+1)*256,p.components);++j){float2 a=stack[k*p.components*RHS+j*RHS+r],b=w[j*RHS+r];
        sum.x=fma(a.x,b.x,sum.x);sum.x=fma(a.y,b.y,sum.x);sum.y=fma(a.x,b.y,sum.y);sum.y=fma(-a.y,b.x,sum.y);}partial[gid]=sum;
}
kernel void gm_reduce(device const float2* partial [[buffer(0)]],device float2* coefficients [[buffer(1)]],device float2* h [[buffer(2)]],
                      constant GMParams& p [[buffer(3)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.k*RHS)return;uint k=gid/RHS,r=gid%RHS;float2 sum=0;for(uint j=0;j<p.parts;++j)sum+=partial[(k*p.parts+j)*RHS+r];
    coefficients[gid]=sum;h[gid]=p.stage==0?sum:h[gid]+sum;
}
kernel void gm_orthogonalize(device const float2* stack [[buffer(0)]],device const float2* h [[buffer(1)]],device float2* w [[buffer(2)]],
                             constant GMParams& p [[buffer(3)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.components*RHS)return;float2 z=w[gid];for(uint k=0;k<p.k;++k)z-=cmul(stack[k*p.components*RHS+gid],h[k*RHS+gid%RHS]);w[gid]=z;
}
kernel void gm_norm_partial(device const float2* x [[buffer(0)]],device float* partial [[buffer(1)]],constant GMParams& p [[buffer(2)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.parts*RHS)return;uint r=gid%RHS,part=gid/RHS;float sum=0;for(uint j=part*256;j<min((part+1)*256,p.components);++j){float2 v=x[j*RHS+r];sum=fma(v.x,v.x,sum);sum=fma(v.y,v.y,sum);}partial[gid]=sum;
}
kernel void gm_norm_reduce(device const float* partial [[buffer(0)]],device float* norms [[buffer(1)]],constant GMParams& p [[buffer(2)]],uint r [[thread_position_in_grid]]){
    if(r<RHS){float sum=0;for(uint j=0;j<p.parts;++j)sum+=partial[j*RHS+r];norms[r]=sum;}}
kernel void gm_combine(device const float2* zstack [[buffer(0)]],device const float2* y [[buffer(1)]],device float2* delta [[buffer(2)]],
                       constant GMParams& p [[buffer(3)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.components*RHS)return;float2 sum=delta[gid];for(uint k=0;k<p.k;++k)sum+=cmul(zstack[k*p.components*RHS+gid],y[k*RHS+gid%RHS]);delta[gid]=sum;
}
