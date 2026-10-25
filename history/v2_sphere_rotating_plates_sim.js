// ===== SIM-START =====
const SIM = (() => {
const R_KM=6371;
const ROCK={NONE:0,OCEAN:1,HOTSPOT:2,ARC:3,CRATON:4,OROGEN:5,SHELF:6,FORELAND:7};
const C={OCEAN_T:7,CONT_BASE:33,ISO:0.16,CAP:75,ARC:0.6,ARCW:[0,0.5,0.3,0.2],COLL:0.5,COLLW:[0.35,0.3,0.2,0.15],
  KE:0.03,KD:0.03,VOLC_DECAY:0.02,SED_CONV:5,ARC_CONV:17,COLLAPSE:48,HOT:0.7,CONT_TARGET:0.29,MERGE:180,COLL_LOG:30,SPLIT_TAU:40,VMAX:8,VMIN:0.5};
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
function hashStr(s){let h=1779033703^s.length;for(let i=0;i<s.length;i++){h=Math.imul(h^s.charCodeAt(i),3432918353);h=h<<13|h>>>19;}h=Math.imul(h^h>>>16,2246822507);h=Math.imul(h^h>>>13,3266489909);return(h^h>>>16)>>>0;}
function gauss(r){let u=0;while(u===0)u=r();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*r());}
function randUnit(r){const z=2*r()-1,p=2*Math.PI*r(),s=Math.sqrt(1-z*z);return[s*Math.cos(p),s*Math.sin(p),z];}
function norm3(v){const l=Math.hypot(v[0],v[1],v[2])||1;return[v[0]/l,v[1]/l,v[2]/l];}
function cross(a,b){return[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];}
function loc(x,y,z){const lat=Math.asin(Math.max(-1,Math.min(1,z)))*180/Math.PI,lon=Math.atan2(y,x)*180/Math.PI;
  return `${Math.abs(lat).toFixed(0)}°${lat>=0?"N":"S"} ${Math.abs(lon).toFixed(0)}°${lon>=0?"E":"W"}`;}

// 3D value noise (fractal), evaluated on the unit sphere
function makeNoise3(seed){
  const h=(ix,iy,iz)=>{let n=(Math.imul(ix,374761393)+Math.imul(iy,668265263)+Math.imul(iz,1440662683)+seed)|0;n=Math.imul(n^(n>>>13),1274126177);n^=n>>>16;return(n>>>0)/4294967296;};
  const sm=t=>t*t*(3-2*t),l=(a,b,t)=>a+(b-a)*t;
  const vn=(x,y,z)=>{const ix=Math.floor(x),iy=Math.floor(y),iz=Math.floor(z),fx=sm(x-ix),fy=sm(y-iy),fz=sm(z-iz);
    return l(l(l(h(ix,iy,iz),h(ix+1,iy,iz),fx),l(h(ix,iy+1,iz),h(ix+1,iy+1,iz),fx),fy),
             l(l(h(ix,iy,iz+1),h(ix+1,iy,iz+1),fx),l(h(ix,iy+1,iz+1),h(ix+1,iy+1,iz+1),fx),fy),fz);};
  return (x,y,z,freq,oct)=>{let a=1,s=0,t=0;for(let o=0;o<oct;o++){s+=a*vn(x*freq+o*31.7,y*freq+o*11.3,z*freq+o*7.9);t+=a;a*=0.5;freq*=2;}return s/t;};
}

// Equal-angle cube-sphere grid: 6 faces of n x n cells, O(1) nearest-cell lookup
function makeLookup(n){const s=n*0.5,k=4/Math.PI,nn=n*n,m=n-1;
  return function(x,y,z){const ax=x<0?-x:x,ay=y<0?-y:y,az=z<0?-z:z;let f,a,b;
    if(ax>=ay&&ax>=az){f=x>0?0:1;a=y/ax;b=z/ax;}else if(ay>=az){f=y>0?2:3;a=x/ay;b=z/ay;}else{f=z>0?4:5;a=x/az;b=y/az;}
    let i=Math.floor((Math.atan(a)*k+1)*s),j=Math.floor((Math.atan(b)*k+1)*s);
    if(i>m)i=m;else if(i<0)i=0;if(j>m)j=m;else if(j<0)j=0;return f*nn+j*n+i;};}
class Cube{
  constructor(n){
    this.n=n;const N=this.N=6*n*n;this.pos=new Float64Array(3*N);this.nb=new Int32Array(4*N);
    const q=Math.PI/4,t=[0,0,0];
    const dir=(f,a,b)=>{let x,y,z;switch(f){case 0:x=1;y=a;z=b;break;case 1:x=-1;y=a;z=b;break;case 2:x=a;y=1;z=b;break;case 3:x=a;y=-1;z=b;break;case 4:x=a;y=b;z=1;break;default:x=a;y=b;z=-1;}
      const L=Math.hypot(x,y,z);t[0]=x/L;t[1]=y/L;t[2]=z/L;};
    const lk=this.lookup=makeLookup(n);
    for(let f=0;f<6;f++)for(let j=0;j<n;j++)for(let i=0;i<n;i++){
      const id=f*n*n+j*n+i;dir(f,Math.tan(((i+.5)/n*2-1)*q),Math.tan(((j+.5)/n*2-1)*q));
      this.pos[3*id]=t[0];this.pos[3*id+1]=t[1];this.pos[3*id+2]=t[2];
      const D=[[1,0],[-1,0],[0,1],[0,-1]];
      for(let k=0;k<4;k++){dir(f,Math.tan(((i+.5+D[k][0])/n*2-1)*q),Math.tan(((j+.5+D[k][1])/n*2-1)*q));this.nb[4*id+k]=lk(t[0],t[1],t[2]);}
    }
  }
}
const cubes={};function getCube(n){return cubes[n]||(cubes[n]=new Cube(n));}
function fieldNoise(cube,nz,freq,oct,off){const N=cube.N,p=cube.pos,out=new Float32Array(N);let mn=1e9,mx=-1e9;
  for(let i=0;i<N;i++){const v=nz(p[3*i]+off,p[3*i+1]-off,p[3*i+2]+off*0.7,freq,oct);out[i]=v;if(v<mn)mn=v;if(v>mx)mx=v;}
  for(let i=0;i<N;i++)out[i]=(out[i]-mn)/(mx-mn||1);return out;}

class World{
  constructor(o){
    const cube=this.cube=getCube(o.n);const N=this.N=cube.N;this.n=o.n;
    const sn=hashStr(String(o.seed)+"|"+o.n+"|"+o.nPlates+"|"+o.mode);
    this.rand=mulberry32(sn);this.noise=makeNoise3((sn^0x9e3779b9)|0);
    this.cellRad=Math.PI/2/o.n;this.cellKm=this.cellRad*R_KM;
    this.wmax=C.VMAX*10/R_KM;this.wmin=C.VMIN*10/R_KM;
    this.dt=Math.min(1,0.7*this.cellRad/this.wmax);
    const F=Float32Array;
    this.owner=new Int16Array(N);this.nOwner=new Int16Array(N);this.bodyIdx=new Int32Array(N);
    this.wT=new F(N);this.wB=new F(N);this.wS=new F(N);this.wRA=new F(N);this.wR=new Uint8Array(N);this.wK=new Uint8Array(N);
    this.elev=new F(N);this.dT=new F(N);this.arcAdd=new F(N);this.arcPlate=new Int16Array(N);this.collAdd=new F(N);this.collPlate=new Int16Array(N);
    this.trench=new F(N);this.gaps=new Int32Array(N);this.rest=new Int32Array(N);
    this.plates=[];this.pairs=new Map();this.collLogged=new Set();this.events=[];this.nEvents=0;this.t=0;this.steps=0;this.targetPlates=o.nPlates;
    this.init(o);
  }
  addPlate(W){const N=this.N,id=this.plates.length;
    const P={id,Q:new Float64Array([1,0,0,0,1,0,0,0,1]),W:W.slice(),T:new Float32Array(N),B:new Float32Array(N),S:new Float32Array(N),RA:new Float32Array(N),
      R:new Uint8Array(N),K:new Uint8Array(N),clr:new Int32Array(N).fill(-10),area:0,contArea:0,maxArea:0,c:[1,0,0],alive:true,hue:(id*137.508+20)%360};
    this.plates.push(P);return P;}
  kill(P){P.alive=false;P.T=P.B=P.S=P.RA=P.R=P.K=P.clr=null;P.area=0;}
  log(type,text,v){this.nEvents++;this.events.push({n:this.nEvents,t:this.t,type,text,where:v?loc(v[0],v[1],v[2]):""});if(this.events.length>200)this.events.shift();}
  init(o){
    const {N,cube,rand:r,noise:nz}=this,p=cube.pos;
    const nC=fieldNoise(cube,nz,1.6,6,0),nA=fieldNoise(cube,nz,2.3,4,3.3),nR=fieldNoise(cube,nz,2.6,5,7.1),nG=fieldNoise(cube,nz,2.7,3,11.9);
    const wx=fieldNoise(cube,nz,1.8,4,17.3),wy=fieldNoise(cube,nz,1.8,4,23.9),wz=fieldNoise(cube,nz,1.8,4,29.1);
    const field=new Float32Array(N);
    if(o.mode==="super"){const c=randUnit(r),sig=0.6;
      for(let i=0;i<N;i++){const d=Math.acos(Math.max(-1,Math.min(1,p[3*i]*c[0]+p[3*i+1]*c[1]+p[3*i+2]*c[2])));field[i]=0.62*Math.exp(-d*d/(2*sig*sig))+0.38*nC[i];}
    } else field.set(nC);
    const thr=Float32Array.from(field).sort()[Math.floor(N*0.7)];
    const {wT,wB,wS,wRA,wR,wK}=this;
    for(let i=0;i<N;i++){
      wS[i]=0;
      if(field[i]>thr){
        const ridge=Math.pow(1-Math.abs(2*nR[i]-1),3);wK[i]=1;wB[i]=-2000;
        wT[i]=C.CONT_BASE-1+4*Math.min(1,(field[i]-thr)/0.15)+14*ridge;
        if(ridge>0.62){wR[i]=ROCK.OROGEN;wRA[i]=-(300+700*nG[i]);}
        else if(field[i]-thr<0.025){wR[i]=ROCK.SHELF;wRA[i]=-(30+250*nG[i]);}
        else{wR[i]=ROCK.CRATON;wRA[i]=-(1200+1800*nG[i]);}
      } else {wK[i]=0;wT[i]=C.OCEAN_T;wB[i]=-(10+170*nA[i]);wR[i]=ROCK.OCEAN;wRA[i]=wB[i];}
    }
    const centers=[];
    for(let k=0;k<o.nPlates;k++){const ax=randUnit(r),sp=(1+3.5*r())*10/R_KM;centers.push({c:randUnit(r),w:0.55+0.9*r()});this.addPlate([ax[0]*sp,ax[1]*sp,ax[2]*sp]);}
    for(let i=0;i<N;i++){
      const v=norm3([p[3*i]+0.45*(wx[i]-0.5),p[3*i+1]+0.45*(wy[i]-0.5),p[3*i+2]+0.45*(wz[i]-0.5)]);
      let best=0,bd=1e9;
      for(let k=0;k<centers.length;k++){const c=centers[k].c;const d=Math.acos(Math.max(-1,Math.min(1,v[0]*c[0]+v[1]*c[1]+v[2]*c[2])))/centers[k].w;if(d<bd){bd=d;best=k;}}
      this.owner[i]=best;this.bodyIdx[i]=i;
    }
    this.scatter();
    this.hotspots=[];const nh=3+Math.floor(r()*4);
    for(let k=0;k<nh;k++){const v=randUnit(r);this.hotspots.push({v,i:cube.lookup(v[0],v[1],v[2]),wasUp:false,last:-1e9});}
    this.computeElev();this.updateStats(true);
  }
  scatter(){const {N,plates,owner,bodyIdx,wT,wB,wS,wRA,wR,wK}=this;
    for(let i=0;i<N;i++){const P=plates[owner[i]],k=bodyIdx[i];P.T[k]=wT[i];P.B[k]=wB[i];P.S[k]=wS[i];P.RA[k]=wRA[i];P.R[k]=wR[i];P.K[k]=wK[i];}}
  gather(){const {N,plates,owner,bodyIdx,wT,wB,wS,wRA,wR,wK}=this;
    for(let i=0;i<N;i++){const P=plates[owner[i]],k=bodyIdx[i];wT[i]=P.T[k];wB[i]=P.B[k];wS[i]=P.S[k];wRA[i]=P.RA[k];wR[i]=P.R[k];wK[i]=P.K[k];}}
  copyCell(Ps,ks,Pd,kd){Pd.T[kd]=Ps.T[ks];Pd.B[kd]=Ps.B[ks];Pd.S[kd]=Ps.S[ks];Pd.RA[kd]=Ps.RA[ks];Pd.R[kd]=Ps.R[ks];Pd.K[kd]=Ps.K[ks];}
  deposit(arr,parr,plate,x,y,z,ux,uy,uz,amt,W,k0){
    const d=this.cellRad,lk=this.cube.lookup;
    for(let k=k0;k<W.length;k++){const c=lk(x+ux*d*k,y+uy*d*k,z+uz*d*k);arr[c]+=amt*W[k];parr[c]=plate;}
  }
  resolve(x,y,z,p1,k1,p2,k2){
    const P1=this.plates[p1],P2=this.plates[p2],c1=P1.K[k1],c2=P2.K[k2];
    let w;if(c1!==c2)w=c1?0:1;else if(c1)w=P1.T[k1]>=P2.T[k2]?0:1;else w=P1.B[k1]>=P2.B[k2]?0:1;
    const Pw=w===0?P1:P2,Pl=w===0?P2:P1,kl=w===0?k2:k1;
    const dx=Pl.W[0]-Pw.W[0],dy=Pl.W[1]-Pw.W[1],dz=Pl.W[2]-Pw.W[2];
    let rx=dy*z-dz*y,ry=dz*x-dx*z,rz=dx*y-dy*x;const len=Math.hypot(rx,ry,rz);
    if(len>1e-9){rx/=len;ry/=len;rz/=len;
      if(!Pl.K[kl]){
        this.deposit(this.arcAdd,this.arcPlate,Pw.id,x,y,z,rx,ry,rz,C.ARC+0.4*Pl.S[kl],C.ARCW,1);
        const d=this.cellRad;this.trench[this.cube.lookup(x-rx*d,y-ry*d,z-rz*d)]=1;
      } else {
        this.deposit(this.collAdd,this.collPlate,Pw.id,x,y,z,rx,ry,rz,C.COLL*Pl.T[kl],C.COLLW,0);
        const key=Math.min(Pw.id,Pl.id)*4096+Math.max(Pw.id,Pl.id);this.pairs.set(key,(this.pairs.get(key)||0)+1);
        if(!this.pairLoc)this.pairLoc=new Map();this.pairLoc.set(key,[x,y,z]);
      }
    }
    Pl.T[kl]=-Pl.T[kl];Pl.clr[kl]=this.steps; // negative = subducted this step, restorable
    return w;
  }
  step(){
    const {N,cube,dt,plates,owner,nOwner,bodyIdx}=this,pos=cube.pos,nb=cube.nb,lk=cube.lookup;
    // 1. rotate every plate about its Euler pole
    for(const P of plates){if(!P.alive)continue;
      const w=Math.hypot(P.W[0],P.W[1],P.W[2]);if(w<=0)continue;
      const ux=P.W[0]/w,uy=P.W[1]/w,uz=P.W[2]/w,a=w*dt,c=Math.cos(a),s=Math.sin(a),t=1-c;
      const R=[c+ux*ux*t,ux*uy*t-uz*s,ux*uz*t+uy*s, uy*ux*t+uz*s,c+uy*uy*t,uy*uz*t-ux*s, uz*ux*t-uy*s,uz*uy*t+ux*s,c+uz*uz*t];
      const Q=P.Q,M=new Float64Array(9);
      for(let r=0;r<3;r++)for(let k=0;k<3;k++)M[r*3+k]=R[r*3]*Q[k]+R[r*3+1]*Q[3+k]+R[r*3+2]*Q[6+k];
      if(this.steps%64===0){for(let r=0;r<3;r++){for(let q=0;q<r;q++){const d=M[r*3]*M[q*3]+M[r*3+1]*M[q*3+1]+M[r*3+2]*M[q*3+2];for(let k=0;k<3;k++)M[r*3+k]-=d*M[q*3+k];}
        const l=Math.hypot(M[r*3],M[r*3+1],M[r*3+2]);for(let k=0;k<3;k++)M[r*3+k]/=l;}}
      P.Q=M;
    }
    // 2. which plates cover each fixed cell now? resolve overlaps
    nOwner.fill(-1);this.arcAdd.fill(0);this.collAdd.fill(0);
    const cand=new Int16Array(5);let ng=0;
    for(let i=0;i<N;i++){
      const x=pos[3*i],y=pos[3*i+1],z=pos[3*i+2];
      let nc=0;cand[nc++]=owner[i];
      for(let k=0;k<4;k++){const q=owner[nb[4*i+k]];let dup=false;for(let m=0;m<nc;m++)if(cand[m]===q){dup=true;break;}if(!dup)cand[nc++]=q;}
      let cov=0,wP=-1,wK=-1;
      for(let m=0;m<nc;m++){const p=cand[m],P=plates[p];if(!P.alive)continue;const Q=P.Q;
        const k=lk(Q[0]*x+Q[3]*y+Q[6]*z,Q[1]*x+Q[4]*y+Q[7]*z,Q[2]*x+Q[5]*y+Q[8]*z);
        if(P.T[k]>0){if(cov===0){wP=p;wK=k;}else if(this.resolve(x,y,z,wP,wK,p,k)===1){wP=p;wK=k;}cov++;}
      }
      if(cov===0)this.gaps[ng++]=i;else{nOwner[i]=wP;bodyIdx[i]=wK;}
    }
    // 3. fill gaps: interior holes extend the surrounding plate; true gaps become new ridge crust
    let gaps=this.gaps,rest=this.rest;
    for(let pass=0;pass<12&&ng>0;pass++){
      let nr=0;
      for(let g=0;g<ng;g++){const i=gaps[g],x=pos[3*i],y=pos[3*i+1],z=pos[3*i+2];
        const p=owner[i],Po=plates[p];let kOld=-1;
        if(Po&&Po.alive){const Q=Po.Q;kOld=lk(Q[0]*x+Q[3]*y+Q[6]*z,Q[1]*x+Q[4]*y+Q[7]*z,Q[2]*x+Q[5]*y+Q[8]*z);}
        const cleared=kOld>=0&&Po.clr[kOld]>=this.steps-2;
        let best=-1,bc=0,bj=-1;
        for(let k=0;k<4;k++){const j=nb[4*i+k],q=nOwner[j];if(q<0||(cleared&&q===p))continue;let c=0;for(let m=0;m<4;m++)if(nOwner[nb[4*i+m]]===q)c++;if(c>bc){bc=c;best=q;bj=j;}}
        if(cleared&&Po.T[kOld]<0){ // duplicate sample of crust just subducted elsewhere: keep it, it is still at the surface here
          Po.T[kOld]=-Po.T[kOld];nOwner[i]=p;bodyIdx[i]=kOld;
        } else if(kOld>=0&&!cleared&&bc<3){
          Po.T[kOld]=C.OCEAN_T;Po.B[kOld]=this.t;Po.S[kOld]=0;Po.RA[kOld]=this.t;Po.R[kOld]=ROCK.OCEAN;Po.K[kOld]=0;
          nOwner[i]=p;bodyIdx[i]=kOld;
        } else if(best>=0){
          const Pq=plates[best],Q=Pq.Q;const kq=lk(Q[0]*x+Q[3]*y+Q[6]*z,Q[1]*x+Q[4]*y+Q[7]*z,Q[2]*x+Q[5]*y+Q[8]*z);
          if(Pq.T[kq]<=0){ // prefer a neighbour whose crust type matches what was here last step
            let src=bodyIdx[bj];
            for(let m=0;m<4;m++){const j=nb[4*i+m];if(nOwner[j]===best&&Pq.K[bodyIdx[j]]===this.wK[i]){src=bodyIdx[j];break;}}
            this.copyCell(Pq,src,Pq,kq);}
          nOwner[i]=best;bodyIdx[i]=kq;
        } else rest[nr++]=i;
      }
      const tmp=gaps;gaps=rest;rest=tmp;ng=nr;
    }
    for(let g=0;g<ng;g++){const i=gaps[g];nOwner[i]=owner[i]>=0&&plates[owner[i]].alive?owner[i]:0;} // pathological leftovers
    this.owner=nOwner;this.nOwner=owner;
    const O=this.owner;
    // guard any leftover that maps to empty crust
    for(let g=0;g<ng;g++){const i=gaps[g],P=plates[O[i]];const Q=P.Q,x=pos[3*i],y=pos[3*i+1],z=pos[3*i+2];
      const k=lk(Q[0]*x+Q[3]*y+Q[6]*z,Q[1]*x+Q[4]*y+Q[7]*z,Q[2]*x+Q[5]*y+Q[8]*z);
      if(P.T[k]<=0){P.T[k]=C.OCEAN_T;P.B[k]=this.t;P.S[k]=0;P.RA[k]=this.t;P.R[k]=ROCK.OCEAN;P.K[k]=0;}bodyIdx[i]=k;}
    // a cell sampled before its crust was subducted via a duplicate sample keeps that crust for now
    for(let i=0;i<N;i++){const P=plates[O[i]],k=bodyIdx[i],v=P.T[k];
      if(v<0)P.T[k]=-v;else if(v===0){P.T[k]=C.OCEAN_T;P.B[k]=this.t;P.S[k]=0;P.RA[k]=this.t;P.R[k]=ROCK.OCEAN;P.K[k]=0;}}
    // 4. surface processes on the fixed grid
    this.gather();
    const {wT,wB,wS,wRA,wR,wK,dT,arcAdd,arcPlate,collAdd,collPlate,trench}=this,t=this.t,r=this.rand;
    for(let i=0;i<N;i++){
      const a=arcAdd[i];if(a>0&&arcPlate[i]===O[i]){wT[i]+=a;if(wK[i]?a>0.25:wT[i]>C.OCEAN_T+4){wR[i]=ROCK.ARC;wRA[i]=t;}}
      const c=collAdd[i];if(c>0&&collPlate[i]===O[i]){wT[i]+=c;if(c>0.4){wR[i]=ROCK.OROGEN;wRA[i]=t;}}
      trench[i]*=0.85;
    }
    for(const h of this.hotspots){const i=h.i;
      if(wK[i]){wT[i]+=0.1*dt;wR[i]=ROCK.HOTSPOT;wRA[i]=t;}
      else{wT[i]+=C.HOT*dt;wR[i]=ROCK.HOTSPOT;wRA[i]=t;
        for(let k=0;k<4;k++){const j=nb[4*i+k];if(!wK[j]){wT[j]+=0.3*C.HOT*dt;wR[j]=ROCK.HOTSPOT;wRA[j]=t;}}}
    }
    this.computeElev();
    // crust-budget thermostat: margin growth slows when continents exceed the target fraction, speeds up below it
    const fb=Math.pow((this.contFrac||C.CONT_TARGET)/C.CONT_TARGET,6),sedConv=C.SED_CONV*fb,collapseP=0.05*dt/fb;
    const E=this.elev,kd=C.KD*dt;dT.fill(0);
    for(let i=0;i<N;i++){
      if(wK[i]){
        let f=0,lo=-1,le=E[i];
        for(let k=0;k<4;k++){const j=nb[4*i+k];if(wK[j])f+=wT[j]-wT[i];if(E[j]<le){le=E[j];lo=j;}}
        dT[i]+=kd*f;
        if(E[i]>0&&lo>=0){const loss=C.KE*E[i]*dt;dT[i]-=loss;if(wK[lo])dT[lo]+=loss;else wS[lo]+=0.5*loss;}
      } else if(wT[i]>C.OCEAN_T)dT[i]-=(wT[i]-C.OCEAN_T)*C.VOLC_DECAY*dt;
    }
    for(let i=0;i<N;i++){let v=wT[i]+dT[i];
      if(wK[i]){if(v>C.CAP)v=C.CAP;if(v<10)v=10;}
      else{
        if(wS[i]>sedConv){wK[i]=1;v=Math.max(v,24);wS[i]=0;wR[i]=ROCK.SHELF;wRA[i]=t;}
        else if(v>C.ARC_CONV*fb&&arcAdd[i]>0&&arcPlate[i]===O[i])wK[i]=1;
        if(v<5)v=5;
      }
      wT[i]=v;}
    for(let i=0;i<N;i++){
      if(!wK[i]||wT[i]<C.COLLAPSE||r()>collapseP)continue;
      const c=nb[4*i+Math.floor(r()*4)];
      if(!wK[c]&&O[c]===O[i]){wK[c]=1;wT[c]=28;wR[c]=ROCK.FORELAND;wRA[c]=t;wS[c]=0;wT[i]-=8;}
    }
    if(fb>1){ // thermostat, loss side: drowned continental margins founder back to oceanic crust
      const pz=0.1*dt*(fb-1);
      for(let i=0;i<N;i++){if(!wK[i]||E[i]>0.3||r()>pz)continue;let coast=false;for(let k=0;k<4;k++)if(!wK[nb[4*i+k]])coast=true;
        if(coast){wK[i]=0;wT[i]=C.OCEAN_T+1;wB[i]=t-60;wR[i]=ROCK.OCEAN;wRA[i]=t;wS[i]=2;}}
    }
    this.t+=dt;this.steps++;
    this.computeElev();this.scatter();
    // 5. bookkeeping and discrete events
    this.updateStats(this.steps%4===0);this.hotspotEvents();this.plateEvents();
  }
  computeElev(){const {N,wT,wK,wB,wS,elev:E,t}=this;
    for(let i=0;i<N;i++){if(wK[i])E[i]=C.ISO*(wT[i]-C.CONT_BASE);
      else{let e=-2.6-0.33*Math.sqrt(Math.max(0,t-wB[i]))+0.5*(wT[i]-C.OCEAN_T)+0.3*wS[i];E[i]=e<-7?-7:e;}}}
  updateStats(centroids){
    const {N,plates,owner,wK,wB,elev:E,t}=this,p=this.cube.pos;
    for(const P of plates){P.area=0;P.contArea=0;P._c=[0,0,0];}
    let land=0,cont=0,oa=0,on=0;
    for(let i=0;i<N;i++){const P=plates[owner[i]];P.area++;if(wK[i]){P.contArea++;cont++;}else{oa+=t-wB[i];on++;}if(E[i]>0)land++;
      if(centroids){P._c[0]+=p[3*i];P._c[1]+=p[3*i+1];P._c[2]+=p[3*i+2];}}
    for(const P of plates){if(!P.alive)continue;if(P.area>P.maxArea)P.maxArea=P.area;if(centroids&&P.area>0)P.c=norm3(P._c);}
    this.landFrac=land/N;this.contFrac=cont/N;this.oceanAge=on?oa/on:0;
  }
  hotspotEvents(){for(const h of this.hotspots){const up=this.elev[h.i]>0;
    if(up&&!h.wasUp&&!this.wK[h.i]&&this.t-h.last>15){this.log("island","A volcanic island emerges over a hotspot",h.v);h.last=this.t;}h.wasUp=up;}}
  plateEvents(){
    const {N,dt,rand:r,plates}=this;
    for(const P of plates){if(P.alive&&P.area===0){if(P.maxArea>N*0.01)this.log("lost",`Plate ${P.id} is fully subducted`,P.c);this.kill(P);}}
    const sig=0.08*10/R_KM*Math.sqrt(dt);
    for(const P of plates){if(!P.alive)continue;for(let k=0;k<3;k++)P.W[k]+=gauss(r)*sig;
      const w=Math.hypot(P.W[0],P.W[1],P.W[2]);
      if(w>this.wmax)for(let k=0;k<3;k++)P.W[k]*=this.wmax/w;
      else if(w<this.wmin){const a=randUnit(r);for(let k=0;k<3;k++)P.W[k]=a[k]*this.wmin;}}
    const decay=Math.pow(0.97,dt);let mergeKey=-1;
    for(const [k,v] of this.pairs){const nv=v*decay;if(nv<0.5){this.pairs.delete(k);continue;}this.pairs.set(k,nv);
      if(nv>C.COLL_LOG&&!this.collLogged.has(k)){this.collLogged.add(k);this.log("collide",`Continents on plates ${Math.floor(k/4096)} and ${k%4096} begin to collide; mountains rise`,this.pairLoc&&this.pairLoc.get(k));}
      if(nv>C.MERGE&&mergeKey<0)mergeKey=k;}
    if(mergeKey>=0)this.merge(plates[Math.floor(mergeKey/4096)],plates[mergeKey%4096]);
    const alive=plates.filter(p=>p.alive);
    const rate=(alive.length<this.targetPlates?2.5:1)/C.SPLIT_TAU;
    if(alive.length<this.targetPlates+3&&r()<rate*dt){
      const cands=alive.filter(p=>p.area>N*0.04);
      if(cands.length){let tot=0;for(const p of cands)tot+=Math.pow(p.area,1.5);
        let pick=r()*tot,P=cands[0];for(const p of cands){pick-=Math.pow(p.area,1.5);if(pick<=0){P=p;break;}}this.split(P);}
    }
  }
  merge(a,b){
    const purge=()=>{for(const k of [...this.pairs.keys()]){const x=Math.floor(k/4096),y=k%4096;if(x===a.id||y===a.id||x===b.id||y===b.id)this.pairs.delete(k);}};
    if(!a.alive||!b.alive){purge();return;}
    const [Wn,L]=a.area>=b.area?[a,b]:[b,a];const {N,owner,bodyIdx}=this,p=this.cube.pos,lk=this.cube.lookup,QL=L.Q,QW=Wn.Q;
    for(let k=0;k<N;k++){if(L.T[k]<=0)continue;const bx=p[3*k],by=p[3*k+1],bz=p[3*k+2];
      const x=QL[0]*bx+QL[1]*by+QL[2]*bz,y=QL[3]*bx+QL[4]*by+QL[5]*bz,z=QL[6]*bx+QL[7]*by+QL[8]*bz;
      const kk=lk(QW[0]*x+QW[3]*y+QW[6]*z,QW[1]*x+QW[4]*y+QW[7]*z,QW[2]*x+QW[5]*y+QW[8]*z);if(Wn.T[kk]<=0)this.copyCell(L,k,Wn,kk);}
    for(let i=0;i<N;i++){if(owner[i]!==L.id)continue;const x=p[3*i],y=p[3*i+1],z=p[3*i+2];
      const kk=lk(QW[0]*x+QW[3]*y+QW[6]*z,QW[1]*x+QW[4]*y+QW[7]*z,QW[2]*x+QW[5]*y+QW[8]*z);
      if(Wn.T[kk]<=0)this.copyCell(L,bodyIdx[i],Wn,kk);owner[i]=Wn.id;bodyIdx[i]=kk;}
    const tot=Wn.area+L.area;for(let k=0;k<3;k++)Wn.W[k]=(Wn.W[k]*Wn.area+L.W[k]*L.area)/tot;
    Wn.area=tot;Wn.contArea+=L.contArea;
    this.log("suture",`Plates ${Wn.id} and ${L.id} weld into one along a suture`,Wn.c);
    this.kill(L);purge();
  }
  split(P){
    const {N,owner,bodyIdx,rand:r,noise:nz}=this,p=this.cube.pos;
    const c=P.c;let n=randUnit(r);const d=n[0]*c[0]+n[1]*c[1]+n[2]*c[2];n=norm3([n[0]-d*c[0],n[1]-d*c[1],n[2]-d*c[2]]);
    const a=norm3(cross(c,n));const Qp=this.addPlate(P.W);Qp.Q.set(P.Q);const Q=P.Q;
    const moved=[];let contA=0,contB=0;
    for(let k=0;k<N;k++){if(P.T[k]<=0)continue;const bx=p[3*k],by=p[3*k+1],bz=p[3*k+2];
      const x=Q[0]*bx+Q[1]*by+Q[2]*bz,y=Q[3]*bx+Q[4]*by+Q[5]*bz,z=Q[6]*bx+Q[7]*by+Q[8]*bz;
      const s=x*n[0]+y*n[1]+z*n[2]+(nz(x+40,y,z,2.2,4)-0.5)*0.5;
      if(s>0){this.copyCell(P,k,Qp,k);P.T[k]=0;moved.push(k);if(P.K[k])contB++;}else if(P.K[k])contA++;}
    let nNew=0;for(let i=0;i<N;i++)if(owner[i]===P.id&&Qp.T[bodyIdx[i]]>0)nNew++;
    if(nNew<60||P.area-nNew<60){for(const k of moved)this.copyCell(Qp,k,P,k);this.kill(Qp);return;}
    for(let i=0;i<N;i++)if(owner[i]===P.id&&Qp.T[bodyIdx[i]]>0)owner[i]=Qp.id;
    const s=(1.5+2.5*r())*10/R_KM;for(let k=0;k<3;k++){Qp.W[k]=P.W[k]+a[k]*s/2;P.W[k]-=a[k]*s/2;}
    this.log("rift",`${contA>40&&contB>40?"Continental rift":"Rift"}: plate ${P.id} splits and plate ${Qp.id} drifts away`,c);
    this.updateStats(true);
  }
}
return {World,C,ROCK,R_KM,getCube};
})();
// ===== SIM-END =====
