struct SCParams {uint n,bs,parts,mode;};
kernel void sc_apply(device const uint* row [[buffer(0)]],device const uint* col [[buffer(1)]],device const float2* a [[buffer(2)]],device const float2* x [[buffer(3)]],device float2* y [[buffer(4)]],device const float2* scalars [[buffer(5)]],constant SCParams& p [[buffer(6)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.n*RHS)return;uint r=gid%RHS,q=gid/RHS,ar=q/p.bs,i=q%p.bs;float2 sum=0;
    if(p.mode && scalars[6*RHS+r].x==0){y[gid]=0;return;}
    for(uint e=row[ar];e<row[ar+1];++e)for(uint j=0;j<p.bs;++j){float2 z=a[(e*p.bs+i)*p.bs+j],v=x[(col[e]*p.bs+j)*RHS+r];sum.x=fma(z.x,v.x,sum.x);sum.x=fma(-z.y,v.y,sum.x);sum.y=fma(z.x,v.y,sum.y);sum.y=fma(z.y,v.x,sum.y);}y[gid]=sum;
}
kernel void sc_diagonal(device const float2* a [[buffer(0)]],device const float2* x [[buffer(1)]],device float2* y [[buffer(2)]],constant SCParams& p [[buffer(3)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.n*RHS)return;uint r=gid%RHS,q=gid/RHS,ar=q/p.bs,i=q%p.bs;float2 sum=0;
    for(uint j=0;j<p.bs;++j)sum+=cmul(a[(ar*p.bs+i)*p.bs+j],x[(ar*p.bs+j)*RHS+r]);y[gid]=sum;
}
kernel void sc_difference(device const float2* b [[buffer(0)]],device float2* ax [[buffer(1)]],constant SCParams& p [[buffer(2)]],uint gid [[thread_position_in_grid]]){if(gid<p.n*RHS)ax[gid]=b[gid]-ax[gid];}
kernel void sc_diagnostics(device const float2* partial [[buffer(0)]],device const float2* second [[buffer(1)]],device float2* scalars [[buffer(2)]],constant SCParams& p [[buffer(3)]],uint r [[thread_position_in_grid]]){
    if(r>=RHS)return;float a=0,b=0;for(uint k=0;k<p.parts;++k){a+=partial[k*RHS+r].x;b+=second[k*RHS+r].x;}scalars[10*RHS+r]=float2(a,b);
}
kernel void sc_copy(device const float2* x [[buffer(0)]],device float2* y [[buffer(1)]],constant SCParams& p [[buffer(2)]],uint gid [[thread_position_in_grid]]){if(gid<p.n*RHS)y[gid]=x[gid];}
