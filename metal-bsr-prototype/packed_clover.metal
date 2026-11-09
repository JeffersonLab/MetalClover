inline float2 hermitian6(device const float* a,uint i,uint j){if(i==j)return float2(a[i],0);bool lower=i>j;uint r=lower?i:j,c=lower?j:i,p=6+r*(r-1)+2*c;return float2(a[p],lower?a[p+1]:-a[p+1]);}
kernel void clover_local_hermitian(device const float* a [[buffer(0)]],device const float2* x [[buffer(1)]],device float2* y [[buffer(2)]],constant CloverParams& p [[buffer(3)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.rows*12*RHS)return;uint r=gid%RHS,s=(gid/RHS)%4,c=(gid/(RHS*4))%3,site=gid/(12*RHS),ch=s/2,q=(s%2)*3+c;float2 acc=0;
    for(uint k=0;k<6;++k){float2 z=hermitian6(a+(site*2+ch)*36,q,k),v=x[((site*3+k%3)*4+ch*2+k/3)*RHS+r];acc.x=fma(z.x,v.x,acc.x);acc.x=fma(-z.y,v.y,acc.x);acc.y=fma(z.x,v.y,acc.y);acc.y=fma(z.y,v.x,acc.y);}if(p.beta!=0)acc=fma(float2(p.beta),y[gid],acc);y[gid]=acc;
}
kernel void mg3_diagonal_hermitian(device const float* a [[buffer(0)]],device const float2* x [[buffer(1)]],device float2* y [[buffer(2)]],constant MG3Params& p [[buffer(3)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.rows*12*RHS)return;uint r=gid%RHS,q=gid/RHS,site=q/12,i=q%12,ch=i%4/2,row=(i%2)*3+i/4;float2 sum=0;
    for(uint k=0;k<6;++k)sum+=cmul(hermitian6(a+(site*2+ch)*36,row,k),x[(site*12+(k%3)*4+ch*2+k/3)*RHS+r]);y[gid]=sum;
}

inline float2 lower6(device const float* f,uint i,uint j){uint p=6+i*(i-1)+2*j;return float2(f[p],f[p+1]);}
inline float2 ldl6_value(device const float* f,device const float2* x,uint site,uint ch,uint r,uint output,bool solve){
    float2 v[6],w[6];for(uint i=0;i<6;++i)v[i]=x[(site*12+(i%3)*4+ch*2+i/3)*RHS+r];
    if(solve){
        for(uint i=0;i<6;++i){float2 z=v[i];for(uint k=0;k<i;++k)z-=cmul(lower6(f,i,k),w[k]);w[i]=z;}
        for(uint i=0;i<6;++i)w[i]/=f[i];
        for(int i=5;i>=0;--i){float2 z=w[i];for(uint k=uint(i)+1;k<6;++k){float2 l=lower6(f,k,uint(i));l.y=-l.y;z-=cmul(l,v[k]);}v[i]=z;}
    }else{
        for(uint i=0;i<6;++i){float2 z=v[i];for(uint k=i+1;k<6;++k){float2 l=lower6(f,k,i);l.y=-l.y;z+=cmul(l,v[k]);}w[i]=f[i]*z;}
        float2 z=w[output];for(uint k=0;k<output;++k)z+=cmul(lower6(f,output,k),w[k]);return z;
    }return v[output];
}
kernel void clover_ldl_multiply(device const float* a [[buffer(0)]],device const float2* x [[buffer(1)]],device float2* y [[buffer(2)]],constant CloverParams& p [[buffer(3)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.rows*12*RHS)return;uint r=gid%RHS,s=(gid/RHS)%4,c=(gid/(RHS*4))%3,site=gid/(12*RHS),ch=s/2,q=(s%2)*3+c;float2 z=ldl6_value(a+(site*2+ch)*36,x,site,ch,r,q,false);if(p.beta!=0)z=fma(float2(p.beta),y[gid],z);y[gid]=z;
}
kernel void clover_ldl_solve(device const float* a [[buffer(0)]],device const float2* x [[buffer(1)]],device float2* y [[buffer(2)]],constant CloverParams& p [[buffer(3)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.rows*12*RHS)return;uint r=gid%RHS,s=(gid/RHS)%4,c=(gid/(RHS*4))%3,site=gid/(12*RHS),ch=s/2,q=(s%2)*3+c;y[gid]=ldl6_value(a+(site*2+ch)*36,x,site,ch,r,q,true);
}
kernel void mg3_diagonal_ldl_solve(device const float* a [[buffer(0)]],device const float2* x [[buffer(1)]],device float2* y [[buffer(2)]],constant MG3Params& p [[buffer(3)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.rows*12*RHS)return;uint r=gid%RHS,q=gid/RHS,site=q/12,i=q%12,ch=i%4/2,row=(i%2)*3+i/4;y[gid]=ldl6_value(a+(site*2+ch)*36,x,site,ch,r,row,true);
}
