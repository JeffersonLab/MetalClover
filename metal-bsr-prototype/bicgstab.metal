struct BiParams { uint components,parts,limit; float tolerance2; uint stage; };
kernel void bi_initial(device const float2* input [[buffer(0)]],device float2* r [[buffer(1)]],device float2* shadow [[buffer(2)]],
                       device float2* direction [[buffer(3)]],device float2* v [[buffer(4)]],device float2* delta [[buffer(5)]],
                       constant BiParams& p [[buffer(6)]],uint gid [[thread_position_in_grid]]){
    if(gid>=p.components*RHS)return;r[gid]=shadow[gid]=input[gid];direction[gid]=v[gid]=delta[gid]=float2(0);
}
// Scalar slots: rho, previous_rho, alpha, omega, beta, initial_norm,
// active, iterations, breakdown, current_norm. Real quantities use .x.
inline float2 zdiv(float2 a,float2 b) {
    if(abs(b.x)>=abs(b.y)) {float q=b.y/b.x,d=b.x+b.y*q;return float2((a.x+a.y*q)/d,(a.y-a.x*q)/d);}
    float q=b.x/b.y,d=b.y+b.x*q;return float2((a.x*q+a.y)/d,(a.y*q-a.x)/d);
}
inline bool valid_den(float2 z){return all(isfinite(z)) && max(abs(z.x),abs(z.y))>0;}
kernel void bi_partial(device const float2* v [[buffer(0)]],device const float2* w [[buffer(1)]],
                       device float2* partial [[buffer(2)]],constant BiParams& p [[buffer(3)]],
                       uint gid [[thread_position_in_grid]]) {
    if(gid>=RHS*p.parts)return;
    uint r=gid%RHS,part=gid/RHS;float2 sum=0;
    for(uint j=part*256;j<min((part+1)*256,p.components);++j){
        float2 a=v[j*RHS+r],b=w[j*RHS+r];
        sum.x=fma(a.x,b.x,sum.x);sum.x=fma(a.y,b.y,sum.x);
        sum.y=fma(a.x,b.y,sum.y);sum.y=fma(-a.y,b.x,sum.y);
    }partial[gid]=sum;
}
kernel void bi_scalar(device const float2* partial [[buffer(0)]],device const float2* second [[buffer(1)]],
                      device float2* q [[buffer(2)]],constant BiParams& p [[buffer(3)]],uint r [[thread_position_in_grid]]) {
    if(r>=RHS)return;
    float2 z=0,w=0;for(uint k=0;k<p.parts;++k){z+=partial[k*RHS+r];if(p.stage==4)w+=second[k*RHS+r];}
    if(p.stage==0){
        for(uint k=0;k<10;++k)q[k*RHS+r]=0;
        q[RHS+r]=q[2*RHS+r]=q[3*RHS+r]=float2(1,0);
        q[5*RHS+r]=q[9*RHS+r]=float2(z.x,0);q[6*RHS+r]=float2(z.x>0?1:0,0);return;
    }
    if(q[6*RHS+r].x==0)return;
    bool bad=false;
    if(p.stage==1){
        bad=!valid_den(z)||!valid_den(q[RHS+r])||!valid_den(q[3*RHS+r]);
        if(!bad){q[r]=z;q[4*RHS+r]=cmul(zdiv(z,q[RHS+r]),zdiv(q[2*RHS+r],q[3*RHS+r]));}
    }else if(p.stage==2){bad=!valid_den(z);if(!bad)q[2*RHS+r]=zdiv(q[r],z);}
    else if(p.stage==3){
        bad=!isfinite(z.x);q[9*RHS+r]=float2(z.x,0);
        if(!bad && z.x<=p.tolerance2*q[5*RHS+r].x){q[6*RHS+r]=0;q[7*RHS+r].x+=1;}
    }else if(p.stage==4){
        bad=!(w.x>0)||!isfinite(w.x);if(!bad){q[3*RHS+r]=z/w.x;bad=!valid_den(q[3*RHS+r]);}
    }else if(p.stage==5){
        bad=!isfinite(z.x);q[9*RHS+r]=float2(z.x,0);q[7*RHS+r].x+=1;q[RHS+r]=q[r];
        if(z.x<=p.tolerance2*q[5*RHS+r].x || q[7*RHS+r].x>=p.limit)q[6*RHS+r]=0;
    }
    if(bad){q[8*RHS+r]=float2(1,0);q[6*RHS+r]=0;}
}
kernel void bi_direction(device const float2* r [[buffer(0)]],device float2* direction [[buffer(1)]],
                         device const float2* v [[buffer(2)]],device const float2* q [[buffer(3)]],
                         constant BiParams& p [[buffer(4)]],uint gid [[thread_position_in_grid]]) {
    if(gid>=p.components*RHS)return;uint k=gid%RHS;
    direction[gid]=q[6*RHS+k].x>0?r[gid]+cmul(q[4*RHS+k],direction[gid]-cmul(q[3*RHS+k],v[gid])):float2(0);
}
kernel void bi_update(device float2* delta [[buffer(0)]],device float2* r [[buffer(1)]],
                      device const float2* direction [[buffer(2)]],device const float2* v [[buffer(3)]],
                      device const float2* t [[buffer(4)]],device const float2* q [[buffer(5)]],
                      constant BiParams& p [[buffer(6)]],uint gid [[thread_position_in_grid]]) {
    if(gid>=p.components*RHS)return;uint k=gid%RHS;if(q[6*RHS+k].x==0)return;
    if(p.stage==0){delta[gid]+=cmul(q[2*RHS+k],direction[gid]);r[gid]-=cmul(q[2*RHS+k],v[gid]);}
    else {delta[gid]+=cmul(q[3*RHS+k],r[gid]);r[gid]-=cmul(q[3*RHS+k],t[gid]);}
}
