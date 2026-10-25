// ===== SIM-START =====
const SIM = (() => {
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
function hashStr(s){let h=1779033703^s.length;for(let i=0;i<s.length;i++){h=Math.imul(h^s.charCodeAt(i),3432918353);h=h<<13|h>>>19;}h=Math.imul(h^h>>>16,2246822507);h=Math.imul(h^h>>>13,3266489909);return(h^h>>>16)>>>0;}
function gauss(r){let u=0;while(u===0)u=r();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*r());}

// Tileable value noise on the torus, normalized to [0,1].
function torusNoise(W,H,rand,base,oct,persist=0.5){
  const out=new Float32Array(W*H);let amp=1;
  for(let o=0;o<oct;o++){
    const gx=base<<o, gy=Math.max(1,Math.round(gx*H/W));
    const lat=new Float32Array(gx*gy);for(let k=0;k<lat.length;k++)lat[k]=rand();
    for(let y=0;y<H;y++){
      const fy=y*gy/H,y0=Math.floor(fy),ty=fy-y0,sy=ty*ty*(3-2*ty),ya=(y0%gy)*gx,yb=((y0+1)%gy)*gx;
      for(let x=0;x<W;x++){
        const fx=x*gx/W,x0=Math.floor(fx),tx=fx-x0,sx=tx*tx*(3-2*tx),xa=x0%gx,xb=(x0+1)%gx;
        const a=lat[ya+xa],b=lat[ya+xb],c=lat[yb+xa],d=lat[yb+xb];
        const top=a+(b-a)*sx,bot=c+(d-c)*sx;
        out[y*W+x]+=amp*(top+(bot-top)*sy);
      }
    }
    amp*=persist;
  }
  let mn=Infinity,mx=-Infinity;for(const v of out){if(v<mn)mn=v;if(v>mx)mx=v;}
  for(let i=0;i<out.length;i++)out[i]=(out[i]-mn)/(mx-mn||1);
  return out;
}

// Tunable constants (km, Myr, cm/yr)
const DX=[1,-1,0,0],DY=[0,0,1,-1];
const C={
  OCEAN_T:7, CONT_BASE:33, ISO:0.16, CAP:75,
  ARC:0.6, ARCW:[0,0.5,0.3,0.2],
  COLL:0.5, COLLW:[0.35,0.3,0.2,0.15],
  KE:0.03, KD:0.03, VOLC_DECAY:0.02, SED_CONV:4, ARC_CONV:13, COLLAPSE:48,
  HOT:0.7, MERGE:180, SPLIT_TAU:40, VMAX:8, VMIN:0.5
};

class World{
  constructor(o){
    this.W=o.W;this.H=o.H;const N=this.N=o.W*o.H;
    this.rand=mulberry32(hashStr(String(o.seed)+"|"+o.W+"|"+o.nPlates+"|"+o.mode));
    this.cellKm=40000/o.W;
    this.dt=Math.min(1,0.95*this.cellKm/(C.VMAX*10));
    this.vf=10/this.cellKm; // cm/yr -> cells/Myr
    const F=Float32Array;
    this.owner=new Int16Array(N);this.thick=new F(N);this.cont=new Uint8Array(N);this.age=new F(N);this.sed=new F(N);
    this.nOwner=new Int16Array(N);this.nThick=new F(N);this.nCont=new Uint8Array(N);this.nAge=new F(N);this.nSed=new F(N);
    this.arcAdd=new F(N);this.arcPlate=new Int16Array(N);this.trench=new F(N);this.elev=new F(N);this.dT=new F(N);
    this.plates=[];this.pairs=new Map();this.events=[];this.t=0;this.steps=0;this.targetPlates=o.nPlates;
    this.init(o);
  }
  addPlate(u,v){const id=this.plates.length;const p={id,u,v,ox:0,oy:0,sx:0,sy:0,area:0,contArea:0,maxArea:0,cx:0,cy:0,alive:true,hue:(id*137.508+20)%360};this.plates.push(p);return p;}
  log(type,text,x,y){this.events.push({t:this.t,type,text,x,y});if(this.events.length>200)this.events.shift();}
  init(o){
    const {W,H,N,rand:r}=this;
    const nC=torusNoise(W,H,r,4,6), nA=torusNoise(W,H,r,3,4);
    const wx=torusNoise(W,H,r,4,4), wy=torusNoise(W,H,r,4,4);
    this.rough=torusNoise(W,H,r,3,5);const nR=torusNoise(W,H,r,5,5);
    const field=new Float32Array(N);
    if(o.mode==="super"){
      const cx=r()*W,cy=r()*H,sig=0.16*W;
      for(let y=0;y<H;y++)for(let x=0;x<W;x++){
        let dx=Math.abs(x-cx);dx=Math.min(dx,W-dx);let dy=Math.abs(y-cy);dy=Math.min(dy,H-dy);
        const i=y*W+x;field[i]=0.62*Math.exp(-(dx*dx+dy*dy)/(2*sig*sig))+0.38*nC[i];
      }
    } else field.set(nC);
    const sorted=Float32Array.from(field).sort();const thr=sorted[Math.floor(N*0.7)];
    for(let i=0;i<N;i++){
      if(field[i]>thr){this.cont[i]=1;this.thick[i]=C.CONT_BASE-1+4*Math.min(1,(field[i]-thr)/0.15)+14*Math.pow(1-Math.abs(2*nR[i]-1),3);this.age[i]=0;}
      else{this.cont[i]=0;this.thick[i]=C.OCEAN_T;this.age[i]=10+170*nA[i];}
    }
    // weighted, warped Voronoi plates
    const centers=[];
    for(let k=0;k<o.nPlates;k++){
      const ang=r()*Math.PI*2,mag=1+3.5*r();
      centers.push({x:r()*W,y:r()*H,w:0.55+0.9*r()});
      this.addPlate(Math.cos(ang)*mag,Math.sin(ang)*mag);
    }
    const warp=0.12*W;
    for(let y=0;y<H;y++)for(let x=0;x<W;x++){
      const i=y*W+x,px=x+(wx[i]-0.5)*warp,py=y+(wy[i]-0.5)*warp;
      let best=0,bd=Infinity;
      for(let k=0;k<centers.length;k++){
        const c=centers[k];let dx=Math.abs(px-c.x)%W;dx=Math.min(dx,W-dx);let dy=Math.abs(py-c.y)%H;dy=Math.min(dy,H-dy);
        const d=Math.sqrt(dx*dx+dy*dy)/c.w;if(d<bd){bd=d;best=k;}
      }
      this.owner[i]=best;
    }
    this.hotspots=[];const nh=3+Math.floor(r()*4);
    for(let k=0;k<nh;k++)this.hotspots.push({x:Math.floor(r()*W),y:Math.floor(r()*H),wasUp:false,last:-1e9});
    this.updateStats(true);this.computeElev();
  }
  wrapIdx(x,y){const W=this.W,H=this.H;x=((x%W)+W)%W;y=((y%H)+H)%H;return y*W+x;}
  collide(j,i,p){
    const q=this.nOwner[j],ic=this.cont[i],oc=this.nCont[j];
    let inWins;
    if(ic!==oc)inWins=ic===1;
    else if(ic===1)inWins=this.thick[i]>this.nThick[j]||(this.thick[i]===this.nThick[j]&&p<q);
    else inWins=this.age[i]<this.nAge[j];
    let wP,lP,lThick,lCont,lSed;
    if(inWins){
      lThick=this.nThick[j];lCont=oc;lSed=this.nSed[j];lP=q;wP=p;
      this.nOwner[j]=p;this.nThick[j]=this.thick[i];this.nCont[j]=ic;this.nAge[j]=this.age[i];this.nSed[j]=this.sed[i];
    } else {lThick=this.thick[i];lCont=ic;lSed=this.sed[i];lP=p;wP=q;}
    const PW=this.plates[wP],PL=this.plates[lP];
    let rx=PL.sx-PW.sx,ry=PL.sy-PW.sy;if(rx===0&&ry===0)return;
    const W=this.W,x=j%W,y=(j/W)|0;
    if(!lCont){
      const amt=C.ARC+0.4*lSed;
      for(let k=1;k<4;k++){const c=this.wrapIdx(x+rx*k,y+ry*k);this.arcAdd[c]+=amt*C.ARCW[k];this.arcPlate[c]=wP;}
      this.trench[this.wrapIdx(x-rx,y-ry)]=1;
    } else {
      const amt=C.COLL*lThick;
      for(let k=0;k<4;k++){const c=this.wrapIdx(x+rx*k,y+ry*k);this.arcAdd[c]+=amt*C.COLLW[k];this.arcPlate[c]=wP;}
      const key=Math.min(wP,lP)*4096+Math.max(wP,lP);this.pairs.set(key,(this.pairs.get(key)||0)+1);
    }
  }
  step(){
    const {W,H,N,dt,vf,rand:r,plates}=this;
    // 1. plate shifts
    for(const P of plates){if(!P.alive)continue;
      P.ox+=P.u*vf*dt;P.oy+=P.v*vf*dt;P.sx=0;P.sy=0;
      if(P.ox>=1){P.sx=1;P.ox-=1;}else if(P.ox<=-1){P.sx=-1;P.ox+=1;}
      if(P.oy>=1){P.sy=1;P.oy-=1;}else if(P.oy<=-1){P.sy=-1;P.oy+=1;}
    }
    // 2. advect cells; resolve conflicts
    const own=this.owner,th=this.thick,co=this.cont,ag=this.age,se=this.sed;
    const no=this.nOwner,nt=this.nThick,nc=this.nCont,na=this.nAge,ns=this.nSed;
    no.fill(-1);this.arcAdd.fill(0);
    for(let y=0;y<H;y++)for(let x=0;x<W;x++){
      const i=y*W+x,p=own[i],P=plates[p];
      let nx=x+P.sx,ny=y+P.sy;if(nx<0)nx+=W;else if(nx>=W)nx-=W;if(ny<0)ny+=H;else if(ny>=H)ny-=H;
      const j=ny*W+nx;
      if(no[j]<0){no[j]=p;nt[j]=th[i];nc[j]=co[i];na[j]=ag[i];ns[j]=se[i];}
      else this.collide(j,i,p);
    }
    // 3. fill gaps with new oceanic crust (ridges)
    let gaps=[];for(let i=0;i<N;i++)if(no[i]<0)gaps.push(i);
    for(let pass=0;pass<20&&gaps.length;pass++){
      const rest=[];
      for(const i of gaps){
        const x=i%W,y=(i/W)|0;
        const nb=[y*W+(x+1)%W,y*W+(x+W-1)%W,((y+1)%H)*W+x,((y+H-1)%H)*W+x];
        const off=Math.floor(r()*4);let got=-1,any=-1;
        // prefer the plate that just moved away from this gap (its trailing edge)
        for(let k=0;k<4;k++){const kk=(k+off)&3,c=nb[kk];if(no[c]<0)continue;if(any<0)any=no[c];
          const P=plates[no[c]];if(P.sx*DX[kk]+P.sy*DY[kk]>0){got=no[c];break;}}
        if(got<0)got=any;
        if(got<0){rest.push(i);continue;}
        no[i]=got;nt[i]=C.OCEAN_T;nc[i]=0;na[i]=0;ns[i]=0;
      }
      gaps=rest;
    }
    // 4. arc and collision thickening
    for(let i=0;i<N;i++){const a=this.arcAdd[i];if(a>0&&no[i]===this.arcPlate[i])nt[i]+=a;}
    // swap
    [this.owner,this.nOwner]=[no,own];[this.thick,this.nThick]=[nt,th];[this.cont,this.nCont]=[nc,co];[this.age,this.nAge]=[na,ag];[this.sed,this.nSed]=[ns,se];
    const T=this.thick,K=this.cont,A=this.age,S=this.sed,O=this.owner,tr=this.trench;
    for(let i=0;i<N;i++){A[i]+=dt;tr[i]*=0.85;}
    // 5. hotspots
    for(const h of this.hotspots){
      const i=h.y*W+h.x;
      if(K[i]){T[i]+=0.1*dt;}else{
        T[i]+=C.HOT*dt;
        T[this.wrapIdx(h.x+1,h.y)]+=0.3*C.HOT*dt;T[this.wrapIdx(h.x-1,h.y)]+=0.3*C.HOT*dt;
        T[this.wrapIdx(h.x,h.y+1)]+=0.3*C.HOT*dt;T[this.wrapIdx(h.x,h.y-1)]+=0.3*C.HOT*dt;
      }
    }
    // 6. erosion, diffusion, sediment
    this.computeElev();
    const E=this.elev,dT=this.dT;dT.fill(0);
    const kd=C.KD*dt;
    for(let y=0;y<H;y++){const yu=((y+H-1)%H)*W,yd=((y+1)%H)*W,yc=y*W;
      for(let x=0;x<W;x++){
        const i=yc+x,xr=(x+1)%W,xl=(x+W-1)%W;
        if(K[i]){
          // diffusion with right and down continental neighbours (each pair once)
          const jr=yc+xr,jd=yd+x;
          if(K[jr]){const f=kd*(T[jr]-T[i]);dT[i]+=f;dT[jr]-=f;}
          if(K[jd]){const f=kd*(T[jd]-T[i]);dT[i]+=f;dT[jd]-=f;}
          const e=E[i];
          if(e>0){
            const loss=C.KE*e*dt;
            let lo=-1,le=e;const nb=[jr,yc+xl,jd,yu+x];
            for(const c of nb)if(E[c]<le){le=E[c];lo=c;}
            if(lo>=0){dT[i]-=loss;if(K[lo])dT[lo]+=loss;else S[lo]+=loss;}
          }
        } else if(T[i]>C.OCEAN_T){dT[i]-=(T[i]-C.OCEAN_T)*C.VOLC_DECAY*dt;}
      }
    }
    for(let i=0;i<N;i++){
      let t=T[i]+dT[i];
      if(K[i]){if(t>C.CAP)t=C.CAP;if(t<10)t=10;}
      else if(S[i]>C.SED_CONV){K[i]=1;t=Math.max(t,24);S[i]=0;A[i]=0;}
      else if(t>C.ARC_CONV&&this.arcAdd[i]>0){K[i]=1;S[i]=0;A[i]=0;}
      T[i]=t;
    }
    // orogenic collapse: very thick crust spreads onto adjacent seafloor of the same plate
    for(let y=0;y<H;y++)for(let x=0;x<W;x++){
      const i=y*W+x;if(!K[i]||T[i]<C.COLLAPSE||r()>0.15*dt)continue;
      const nb=[y*W+(x+1)%W,y*W+(x+W-1)%W,((y+1)%H)*W+x,((y+H-1)%H)*W+x];
      const c=nb[Math.floor(r()*4)];
      if(!K[c]&&O[c]===O[i]){K[c]=1;T[c]=28;A[c]=0;S[c]=0;T[i]-=8;}
    }
    this.t+=dt;this.steps++;
    this.computeElev();
    // 7. bookkeeping and discrete events
    this.updateStats(this.steps%4===0);
    this.hotspotEvents();
    this.plateEvents();
  }
  computeElev(){
    const {N}=this,T=this.thick,K=this.cont,A=this.age,S=this.sed,E=this.elev;
    for(let i=0;i<N;i++){
      if(K[i])E[i]=C.ISO*(T[i]-C.CONT_BASE);
      else{let e=-2.6-0.33*Math.sqrt(A[i])+0.5*(T[i]-C.OCEAN_T)+0.3*S[i];E[i]=e<-7?-7:e;}
    }
  }
  updateStats(centroids){
    const {W,H,N,plates}=this,O=this.owner,K=this.cont;
    for(const P of plates){P.area=0;P.contArea=0;P._sx=0;P._cx=0;P._sy=0;P._cy=0;}
    let land=0,cont=0;const E=this.elev;
    const tx=2*Math.PI/W,ty=2*Math.PI/H;
    for(let y=0;y<H;y++){const sy=Math.sin(y*ty),cy=Math.cos(y*ty);
      for(let x=0;x<W;x++){const i=y*W+x,P=plates[O[i]];P.area++;if(K[i]){P.contArea++;cont++;}
        if(E[i]>0)land++;
        if(centroids){P._sx+=Math.sin(x*tx);P._cx+=Math.cos(x*tx);P._sy+=sy;P._cy+=cy;}
      }}
    for(const P of plates){if(!P.alive)continue;
      if(P.area>P.maxArea)P.maxArea=P.area;
      if(centroids&&P.area>0){P.cx=((Math.atan2(P._sx,P._cx)/tx)+W)%W;P.cy=((Math.atan2(P._sy,P._cy)/ty)+H)%H;}
    }
    this.landFrac=land/N;this.contFrac=cont/N;
  }
  hotspotEvents(){
    const E=this.elev,W=this.W;
    for(const h of this.hotspots){
      const up=E[h.y*W+h.x]>0;
      if(up&&!h.wasUp&&!this.cont[h.y*W+h.x]&&this.t-h.last>15){this.log("island","A volcanic island emerges over a hotspot.",h.x,h.y);h.last=this.t;}
      h.wasUp=up;
    }
  }
  plateEvents(){
    const {N,dt,rand:r,plates}=this;
    // consumed plates
    for(const P of plates){if(P.alive&&P.area===0){P.alive=false;if(P.maxArea>N*0.01)this.log("lost",`Plate ${P.id} is fully subducted.`,P.cx,P.cy);}}
    // velocity random walk
    for(const P of plates){if(!P.alive)continue;
      P.u+=gauss(r)*0.08*Math.sqrt(dt);P.v+=gauss(r)*0.08*Math.sqrt(dt);
      const s=Math.hypot(P.u,P.v);
      if(s>C.VMAX){P.u*=C.VMAX/s;P.v*=C.VMAX/s;}
      else if(s<C.VMIN){const a=r()*Math.PI*2;P.u=Math.cos(a)*C.VMIN;P.v=Math.sin(a)*C.VMIN;}
    }
    // sutures: sustained continental collision welds plates
    const decay=Math.pow(0.97,dt);let mergeKey=-1;
    for(const [k,v] of this.pairs){const nv=v*decay;if(nv<0.5)this.pairs.delete(k);else{this.pairs.set(k,nv);if(nv>C.MERGE&&mergeKey<0)mergeKey=k;}}
    if(mergeKey>=0){
      const a=plates[Math.floor(mergeKey/4096)],b=plates[mergeKey%4096];
      if(a.alive&&b.alive){
        const [Wn,L]=a.area>=b.area?[a,b]:[b,a];
        const O=this.owner;for(let i=0;i<N;i++)if(O[i]===L.id)O[i]=Wn.id;
        const tot=Wn.area+L.area;Wn.u=(Wn.u*Wn.area+L.u*L.area)/tot;Wn.v=(Wn.v*Wn.area+L.v*L.area)/tot;
        Wn.area=tot;Wn.contArea+=L.contArea;L.alive=false;L.area=0;
        this.log("suture",`Plates ${Wn.id} and ${L.id} collide and weld into one; a mountain belt marks the seam.`,Wn.cx,Wn.cy);
      }
      for(const k of [...this.pairs.keys()]){const x=Math.floor(k/4096),y=k%4096;if(x===a.id||y===a.id||x===b.id||y===b.id)this.pairs.delete(k);}
    }
    // rifts
    const alive=plates.filter(p=>p.alive);
    const maxPlates=this.targetPlates+3;
    const rate=(alive.length<this.targetPlates?2.5:1)/C.SPLIT_TAU;
    if(alive.length<maxPlates&&r()<rate*dt){
      const cands=alive.filter(p=>p.area>N*0.04);
      if(cands.length){
        let tot=0;for(const p of cands)tot+=Math.pow(p.area,1.5);
        let pick=r()*tot,P=cands[0];for(const p of cands){pick-=Math.pow(p.area,1.5);if(pick<=0){P=p;break;}}
        this.split(P);
      }
    }
  }
  split(P){
    const {W,H,N,rand:r}=this,O=this.owner,K=this.cont;
    const ang=r()*Math.PI*2,nx=Math.cos(ang),ny=Math.sin(ang);
    const Q=this.addPlate(P.u,P.v);
    const cells=[];let side=0,contA=0,contB=0;
    for(let i=0;i<N;i++){if(O[i]!==P.id)continue;
      const x=i%W,y=(i/W)|0;let dx=x-P.cx;if(dx>W/2)dx-=W;if(dx<-W/2)dx+=W;let dy=y-P.cy;if(dy>H/2)dy-=H;if(dy<-H/2)dy+=H;
      const d=dx*nx+dy*ny+(this.rough[i]-0.5)*30;
      if(d>0){cells.push(i);if(K[i])contB++;}else if(K[i])contA++;
    }
    if(cells.length<60||P.area-cells.length<60){Q.alive=false;return;}
    for(const i of cells)O[i]=Q.id;
    const s=1.5+2.5*r();
    Q.u=P.u+nx*s/2;Q.v=P.v+ny*s/2;P.u-=nx*s/2;P.v-=ny*s/2;
    Q.area=cells.length;P.area-=cells.length;
    const kind=(contA>40&&contB>40)?"Continental rift":"Rift";
    this.log("rift",`${kind}: plate ${P.id} splits, and the new plate ${Q.id} drifts away.`,P.cx,P.cy);
    this.updateStats(true);
  }
}
return {World,C};
})();
// ===== SIM-END =====
