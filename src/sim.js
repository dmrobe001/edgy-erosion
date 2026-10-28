// ===== SIM-START =====
// A function, so its source can be handed to a Web Worker (the marker core runs in one; see mesh.js)
function SIMmodule(){
const R_KM=6371, U=8, FS=15;
const ROCK={NONE:0,OCEAN:1,HOTSPOT:2,ARC:3,CRATON:4,OROGEN:5,SHELF:6,CLASTIC:7};
const ROCK_NAMES=["none","ocean-ridge basalt","hotspot or flood basalt","arc volcanics","cratonic basement","collisional metamorphics","marine shelf sediments","clastic sediments"];
const FK={NONE:0,THRUST:1,NORMAL:2,STRIKE:3}, FK_NAMES=["none","thrust","normal","strike-slip"];
const C={OCEAN_T:7,CONT_BASE:33,ISO:0.16,CAP:75,ARC:0.6,ARCW:[0,0.5,0.3,0.2],COLL:0.5,COLLW:[0.35,0.3,0.2,0.15],
  KE:0.03,KD:0.03,VOLC_DECAY:0.02,SED_CONV:5,ARC_CONV:17,COLLAPSE:48,HOT:0.7,CONT_TARGET:0.40,MERGE:180,COLL_LOG:30,SPLIT_TAU:40,VMAX:8,VMIN:0.5,
  UNIT_MERGE:150,SLAB_BINS:25,SLAB_DT:10,SLAB_INBOARD:3,COVER_DEPTH:1,SPEED_TAU:40,FAULT_GAP:10,FAULT_SWITCH:10};
function mulberry32(a){const f=function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};
  f.getState=()=>a;f.setState=v=>{a=v;};return f;} // state exposed so a world can be saved and restored exactly
function hashStr(s){let h=1779033703^s.length;for(let i=0;i<s.length;i++){h=Math.imul(h^s.charCodeAt(i),3432918353);h=h<<13|h>>>19;}h=Math.imul(h^h>>>16,2246822507);h=Math.imul(h^h>>>13,3266489909);return(h^h>>>16)>>>0;}
function gauss(r){let u=0;while(u===0)u=r();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*r());}
function randUnit(r){const z=2*r()-1,p=2*Math.PI*r(),s=Math.sqrt(1-z*z);return[s*Math.cos(p),s*Math.sin(p),z];}
function norm3(v){const l=Math.hypot(v[0],v[1],v[2])||1;return[v[0]/l,v[1]/l,v[2]/l];}
function cross(a,b){return[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];}
function latlon(x,y,z){return[Math.asin(Math.max(-1,Math.min(1,z)))*180/Math.PI,Math.atan2(y,x)*180/Math.PI];}
function loc(x,y,z){const [lat,lon]=latlon(x,y,z);return `${Math.abs(lat).toFixed(0)}°${lat>=0?"N":"S"} ${Math.abs(lon).toFixed(0)}°${lon>=0?"E":"W"}`;}
function elevOf(K,T,B,S,t){if(K)return C.ISO*(T-C.CONT_BASE);const e=-2.6-0.33*Math.sqrt(Math.max(0,t-B))+0.5*(T-C.OCEAN_T)+0.3*S;return e<-7?-7:e;}
function kindOf(X,S){if(Math.abs(S)>2*Math.abs(X))return FK.STRIKE;return X>=0?FK.THRUST:FK.NORMAL;}

function makeNoise3(seed){
  const h=(ix,iy,iz)=>{let n=(Math.imul(ix,374761393)+Math.imul(iy,668265263)+Math.imul(iz,1440662683)+seed)|0;n=Math.imul(n^(n>>>13),1274126177);n^=n>>>16;return(n>>>0)/4294967296;};
  const sm=t=>t*t*(3-2*t),l=(a,b,t)=>a+(b-a)*t;
  const vn=(x,y,z)=>{const ix=Math.floor(x),iy=Math.floor(y),iz=Math.floor(z),fx=sm(x-ix),fy=sm(y-iy),fz=sm(z-iz);
    return l(l(l(h(ix,iy,iz),h(ix+1,iy,iz),fx),l(h(ix,iy+1,iz),h(ix+1,iy+1,iz),fx),fy),
             l(l(h(ix,iy,iz+1),h(ix+1,iy,iz+1),fx),l(h(ix,iy+1,iz+1),h(ix+1,iy+1,iz+1),fx),fy),fz);};
  return (x,y,z,freq,oct)=>{let a=1,s=0,t=0;for(let o=0;o<oct;o++){s+=a*vn(x*freq+o*31.7,y*freq+o*11.3,z*freq+o*7.9);t+=a;a*=0.5;freq*=2;}return s/t;};
}
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
    const lk=this.lookup=makeLookup(n),D=[[1,0],[-1,0],[0,1],[0,-1]];
    {const s=n*0.5,kk=4/Math.PI,nn=n*n,m=n-1;
     this.lookupFrac=(x,y,z,out)=>{const ax=x<0?-x:x,ay=y<0?-y:y,az=z<0?-z:z;let f,a,b;
       if(ax>=ay&&ax>=az){f=x>0?0:1;a=y/ax;b=z/ax;}else if(ay>=az){f=y>0?2:3;a=x/ay;b=z/ay;}else{f=z>0?4:5;a=x/az;b=y/az;}
       const u=(Math.atan(a)*kk+1)*s,v=(Math.atan(b)*kk+1)*s;let i=Math.floor(u),j=Math.floor(v);if(i>m)i=m;else if(i<0)i=0;if(j>m)j=m;else if(j<0)j=0;
       out[0]=f*nn+j*n+i;out[1]=u-(i+0.5);out[2]=v-(j+0.5);};}
    for(let f=0;f<6;f++)for(let j=0;j<n;j++)for(let i=0;i<n;i++){
      const id=f*n*n+j*n+i;dir(f,Math.tan(((i+.5)/n*2-1)*q),Math.tan(((j+.5)/n*2-1)*q));
      this.pos[3*id]=t[0];this.pos[3*id+1]=t[1];this.pos[3*id+2]=t[2];
      for(let k=0;k<4;k++){dir(f,Math.tan(((i+.5+D[k][0])/n*2-1)*q),Math.tan(((j+.5+D[k][1])/n*2-1)*q));this.nb[4*id+k]=lk(t[0],t[1],t[2]);}
    }
    // cell areas on the unit sphere (equal-angle cube cells vary by about ±25%)
    this.area=new Float64Array(N);let tot=0;
    for(let f=0;f<6;f++)for(let j=0;j<n;j++)for(let i=0;i<n;i++){const a=Math.tan(((i+.5)/n*2-1)*q),b=Math.tan(((j+.5)/n*2-1)*q);
      const A=(1+a*a)*(1+b*b)/Math.pow(1+a*a+b*b,1.5);this.area[f*n*n+j*n+i]=A;tot+=A;}
    const Mt=new Float64Array(9);
    for(let id=0;id<N;id++){const A=this.area[id]*=4*Math.PI/tot,x=this.pos[3*id],y=this.pos[3*id+1],z=this.pos[3*id+2],v=[x,y,z];
      for(let r=0;r<3;r++)for(let c=0;c<3;c++)Mt[r*3+c]+=A*((r===c?1:0)-v[r]*v[c]);}
    const [a,b,c,d,e,f2,g,h,k]=Mt,det=a*(e*k-f2*h)-b*(d*k-f2*g)+c*(d*h-e*g);
    this.Minv=[(e*k-f2*h)/det,(c*h-b*k)/det,(b*f2-c*e)/det,(f2*g-d*k)/det,(a*k-c*g)/det,(c*d-a*f2)/det,(d*h-e*g)/det,(b*g-a*h)/det,(a*e-b*d)/det];
  }
}
const cubes={};function getCube(n){return cubes[n]||(cubes[n]=new Cube(n));}
function fieldNoise(cube,nz,freq,oct,off){const N=cube.N,p=cube.pos,out=new Float32Array(N);let mn=1e9,mx=-1e9;
  for(let i=0;i<N;i++){const v=nz(p[3*i]+off,p[3*i+1]-off,p[3*i+2]+off*0.7,freq,oct);out[i]=v;if(v<mn)mn=v;if(v>mx)mx=v;}
  for(let i=0;i<N;i++)out[i]=(out[i]-mn)/(mx-mn||1);return out;}

// Pool of rock columns: each column is a stack of up to U units (bottom first), plus an optional fault record.
class ColumnPool{
  constructor(cap){this.cap=0;this.nFree=0;this.free=new Int32Array(0);this.grow(cap);}
  grow(nc){const R=new Uint8Array(nc*U),B=new Float32Array(nc*U),T=new Float32Array(nc*U),P=new Int16Array(nc*U),cnt=new Uint8Array(nc),used=new Uint8Array(nc),mark=new Uint8Array(nc),flt=new Int32Array(nc).fill(-1);
    if(this.cap){R.set(this.R);B.set(this.B);T.set(this.T);P.set(this.P);cnt.set(this.cnt);used.set(this.used);flt.set(this.flt);}
    const free=new Int32Array(nc);free.set(this.free.subarray(0,this.nFree));let nf=this.nFree;for(let c=nc-1;c>=this.cap;c--)free[nf++]=c;
    Object.assign(this,{R,B,T,P,cnt,used,mark,flt,free,nFree:nf,cap:nc});}
  alloc(){if(!this.nFree)this.grow(Math.ceil(this.cap*1.5));const c=this.free[--this.nFree];this.used[c]=1;this.cnt[c]=0;this.flt[c]=-1;return c;}
  release(c){this.used[c]=0;this.free[this.nFree++]=c;}
  make(rock,t,th,pl){const c=this.alloc();this.R[c*U]=rock;this.B[c*U]=t;this.T[c*U]=th;this.P[c*U]=pl;this.cnt[c]=1;return c;}
  push(c,rock,t,th,pl){if(th<=0)return;const o=c*U;
    while(this.cnt[c]>1&&this.T[o+this.cnt[c]-1]<0.05&&this.R[o+this.cnt[c]-1]!==rock){const k=o+this.cnt[c]-1;this.T[k-1]+=this.T[k];this.cnt[c]--;} // absorb veneers
    const n=this.cnt[c];
    if(n>0&&this.R[o+n-1]===rock&&this.P[o+n-1]===pl&&t-this.B[o+n-1]<C.UNIT_MERGE){this.T[o+n-1]+=th;return;}
    if(n===U){this.T[o]+=this.T[o+1];for(let k=1;k<U-1;k++){this.R[o+k]=this.R[o+k+1];this.B[o+k]=this.B[o+k+1];this.T[o+k]=this.T[o+k+1];this.P[o+k]=this.P[o+k+1];}this.cnt[c]=U-1;}
    const m=this.cnt[c];this.R[o+m]=rock;this.B[o+m]=t;this.T[o+m]=th;this.P[o+m]=pl;this.cnt[c]=m+1;}
  erode(c,dz){const o=c*U;let n=this.cnt[c];
    while(dz>0&&n>0){const k=o+n-1;if(n===1){this.T[k]=Math.max(0.5,this.T[k]-dz);break;}if(this.T[k]>dz){this.T[k]-=dz;dz=0;}else{dz-=this.T[k];n--;}}
    this.cnt[c]=n;}
  sum(c){const o=c*U;let s=0;for(let k=0;k<this.cnt[c];k++)s+=this.T[o+k];return s;}
  reconcile(c,target){const d=target-this.sum(c),o=c*U;
    if(d>0)this.T[o]+=d;else if(d<0){const take=Math.min(-d,Math.max(0,this.T[o]-0.5));this.T[o]-=take;const rem=-d-take;if(rem>1e-6)this.erode(c,rem);}}
  clone(c){const d=this.alloc(),o=c*U,p=d*U,n=this.cnt[c];for(let k=0;k<n;k++){this.R[p+k]=this.R[o+k];this.B[p+k]=this.B[o+k];this.T[p+k]=this.T[o+k];this.P[p+k]=this.P[o+k];}this.cnt[d]=n;return d;}
  top(c){const o=c*U;for(let k=this.cnt[c]-1;k>0;k--)if(this.T[o+k]>0.05)return o+k;return o;}
  units(c){const o=c*U,a=[];for(let k=0;k<this.cnt[c];k++)a.push({rock:this.R[o+k],born:this.B[o+k],th:this.T[o+k],plate:this.P[o+k]});return a;}
}
// Pool of fault records: [switch-counter, previous episode(7), current episode(7)]; episode = t0,t1,across,along,strike xyz (plate frame)
class FaultPool{
  constructor(cap){this.cap=0;this.nFree=0;this.free=new Int32Array(0);this.grow(cap);}
  grow(nc){const D=new Float32Array(nc*FS).fill(NaN),used=new Uint8Array(nc),mark=new Uint8Array(nc);if(this.cap){D.set(this.D);used.set(this.used);}
    const free=new Int32Array(nc);free.set(this.free.subarray(0,this.nFree));let nf=this.nFree;for(let c=nc-1;c>=this.cap;c--)free[nf++]=c;
    Object.assign(this,{D,used,mark,free,nFree:nf,cap:nc});}
  alloc(){if(!this.nFree)this.grow(Math.ceil(this.cap*1.5));const f=this.free[--this.nFree];this.used[f]=1;const o=f*FS;this.D[o]=0;for(let k=1;k<FS;k++)this.D[o+k]=NaN;return f;}
  release(f){this.used[f]=0;this.free[this.nFree++]=f;}
  clone(f){const g=this.alloc();this.D.copyWithin(g*FS,f*FS,f*FS+FS);return g;}
}

class World{
  constructor(o){
    const cube=this.cube=getCube(o.n);const N=this.N=cube.N;this.n=o.n;this.opts=o;
    const sn=hashStr(String(o.seed)+"|"+o.n+"|"+o.nPlates+"|"+o.mode);
    this.rand=mulberry32(sn);this.noise=makeNoise3((sn^0x9e3779b9)|0);
    this.cellRad=Math.PI/2/o.n;this.cellKm=this.cellRad*R_KM;
    this.wmax=C.VMAX*10/R_KM;this.wmin=C.VMIN*10/R_KM;
    this.dt=Math.min(1,0.7*this.cellRad/this.wmax);this.contTarget=o.contTarget||C.CONT_TARGET;
    const F=Float32Array;
    this.owner=new Int16Array(N);this.nOwner=new Int16Array(N);this.bodyIdx=new Int32Array(N);
    this.wT=new F(N);this.wB=new F(N);this.wS=new F(N);this.wK=new Uint8Array(N);this.wCol=new Int32Array(N);
    this.wR=new Uint8Array(N);this.wRA=new F(N);this.bn=new F(3*N);
    this.elev=new F(N);this.dT=new F(N);this.arcAdd=new F(N);this.arcPlate=new Int16Array(N);this.collAdd=new F(N);this.collPlate=new Int16Array(N);
    this._cw=[0,0,0,0];this.trench=new F(N);
    // slabs (recording only): oceanic area consumed at trenches, weighted by the slab's age, in 10 Myr bins on the fixed
    // (mantle) grid, about 300 km beneath the overriding plate; the last SLAB_BINS bins are kept. L1 turns it into
    // dynamic topography.
    this.slab=new F(C.SLAB_BINS*N);this.slabT=new Float64Array(C.SLAB_BINS).fill(-1e9);this.gaps=new Int32Array(N);this.rest=new Int32Array(N);
    this.pool=new ColumnPool(Math.ceil(N*1.6));this.fpool=new FaultPool(Math.ceil(N*0.4));
    this.plates=[];this.pairs=new Map();this.pairLoc=new Map();this.collLogged=new Set();this.events=[];this.nEvents=0;this.t=0;this.steps=0;this.targetPlates=o.nPlates;
    this.init(o);
  }
  addPlate(W){const N=this.N,id=this.plates.length;
    const P={id,Q:new Float64Array([1,0,0,0,1,0,0,0,1]),W:W.slice(),T:new Float32Array(N),B:new Float32Array(N),S:new Float32Array(N),K:new Uint8Array(N),
      col:new Int32Array(N).fill(-1),clr:new Int32Array(N).fill(-10),seen:new Int32Array(N).fill(-1),w0:Math.max(this.wmin*3,Math.min(this.wmax*0.75,Math.hypot(W[0],W[1],W[2]))),area:0,contArea:0,maxArea:0,c:[1,0,0],alive:true,hue:(id*137.508+20)%360};
    if(this.tracking)P.pid=new Int32Array(N).fill(-1);this.plates.push(P);return P;}
  kill(P){P.alive=false;P.T=P.B=P.S=P.K=P.col=P.clr=P.seen=P.pid=null;P.area=0;}
  // Optional parcel identities for L1 (recording only): every crust cell present now gets id = plate*N + body cell,
  // carried wherever L0 moves crust (splits, welds, hole fills) and cleared where new crust forms.
  untrackParcels(){this.tracking=false;for(const P of this.plates)delete P.pid;}
  trackParcels(){const N=this.N;this.tracking=true;for(const P of this.plates){if(!P.alive)continue;P.pid=new Int32Array(N).fill(-1);for(let k=0;k<N;k++)if(P.T[k]>0)P.pid[k]=P.id*N+k;}}
  log(type,text,v){this.nEvents++;const ll=v?latlon(v[0],v[1],v[2]):null;
    this.events.push({n:this.nEvents,t:this.t,type,text,where:v?loc(v[0],v[1],v[2]):"",lat:ll?+ll[0].toFixed(2):null,lon:ll?+ll[1].toFixed(2):null});
    if(this.events.length>400)this.events.shift();}
  init(o){
    const {N,cube,rand:r,noise:nz,pool}=this,p=cube.pos;
    const nC=fieldNoise(cube,nz,1.6,6,0),nA=fieldNoise(cube,nz,2.3,4,3.3),nR=fieldNoise(cube,nz,2.6,5,7.1),nG=fieldNoise(cube,nz,2.7,3,11.9),nS=fieldNoise(cube,nz,2.0,4,37.7);
    const wx=fieldNoise(cube,nz,1.8,4,17.3),wy=fieldNoise(cube,nz,1.8,4,23.9),wz=fieldNoise(cube,nz,1.8,4,29.1);
    const field=new Float32Array(N);
    if(o.mode==="super"){const c=randUnit(r),sig=0.6;
      for(let i=0;i<N;i++){const d=Math.acos(Math.max(-1,Math.min(1,p[3*i]*c[0]+p[3*i+1]*c[1]+p[3*i+2]*c[2])));field[i]=0.62*Math.exp(-d*d/(2*sig*sig))+0.38*nC[i];}
    } else field.set(nC);
    const thr=Float32Array.from(field).sort()[Math.floor(N*(1-this.contTarget))];
    const {wT,wB,wS,wK,wCol}=this;
    for(let i=0;i<N;i++){
      wS[i]=0;
      if(field[i]>thr){
        const ridge=Math.pow(1-Math.abs(2*nR[i]-1),3);wK[i]=1;wB[i]=-2000;
        const T=C.CONT_BASE-1+4*Math.min(1,(field[i]-thr)/0.15)+14*ridge;wT[i]=T;
        const cAge=-(1200+1800*nG[i]);
        if(ridge>0.62){wCol[i]=pool.make(ROCK.CRATON,cAge,T-8,0);pool.push(wCol[i],ROCK.OROGEN,-(300+700*nG[i]),8,0);}
        else if(field[i]-thr<0.025){wCol[i]=pool.make(ROCK.CRATON,cAge,T-3,0);pool.push(wCol[i],ROCK.SHELF,-(30+250*nG[i]),3,0);}
        else{wCol[i]=pool.make(ROCK.CRATON,cAge,nS[i]>0.68?T-2:T,0);if(nS[i]>0.68)pool.push(wCol[i],ROCK.SHELF,-(400+500*nG[i]),2,0);}
      } else {wK[i]=0;wT[i]=C.OCEAN_T;wB[i]=-(10+170*nA[i]);wCol[i]=pool.make(ROCK.OCEAN,wB[i],C.OCEAN_T,0);}
    }
    const centers=[];
    for(let k=0;k<o.nPlates;k++){const ax=randUnit(r),sp=(1+3.5*r())*10/R_KM;centers.push({c:randUnit(r),w:0.55+0.9*r()});this.addPlate([ax[0]*sp,ax[1]*sp,ax[2]*sp]);}
    for(let i=0;i<N;i++){
      const v=norm3([p[3*i]+0.45*(wx[i]-0.5),p[3*i+1]+0.45*(wy[i]-0.5),p[3*i+2]+0.45*(wz[i]-0.5)]);
      let best=0,bd=1e9;
      for(let k=0;k<centers.length;k++){const c=centers[k].c;const d=Math.acos(Math.max(-1,Math.min(1,v[0]*c[0]+v[1]*c[1]+v[2]*c[2])))/centers[k].w;if(d<bd){bd=d;best=k;}}
      this.owner[i]=best;this.bodyIdx[i]=i;
      const c=wCol[i];for(let k=0;k<pool.cnt[c];k++)pool.P[c*U+k]=best;
    }
    this.scatter();
    if(o.nnr!==false)this.removeNetRotation();
    this.hotspots=[];const nh=3+Math.floor(r()*4);
    for(let k=0;k<nh;k++){const v=randUnit(r);this.hotspots.push({v,i:cube.lookup(v[0],v[1],v[2]),wasUp:false,last:-1e9});}
    this.computeElev();this.surface();this.updateStats(true);
  }
  scatter(){const {N,plates,owner,bodyIdx,wT,wB,wS,wK,wCol}=this;
    for(let i=0;i<N;i++){const P=plates[owner[i]],k=bodyIdx[i];P.T[k]=wT[i];P.B[k]=wB[i];P.S[k]=wS[i];P.K[k]=wK[i];P.col[k]=wCol[i];}}
  gather(){const {N,plates,owner,bodyIdx,wT,wB,wS,wK,wCol,pool}=this;
    for(let i=0;i<N;i++){const P=plates[owner[i]],k=bodyIdx[i];wT[i]=P.T[k];wB[i]=P.B[k];wS[i]=P.S[k];wK[i]=P.K[k];
      let c=P.col[k];if(c<0){c=P.col[k]=pool.make(wK[i]?ROCK.CRATON:ROCK.OCEAN,wB[i],wT[i],P.id);}wCol[i]=c;}}
  cloneCol(c,M){if(c<0)return -1;const pool=this.pool,d=pool.clone(c),f=pool.flt[c];
    if(f>=0){const g=this.fpool.clone(f);pool.flt[d]=g;
      if(M){const D=this.fpool.D;for(const o of [g*FS+5,g*FS+12]){const x=D[o],y=D[o+1],z=D[o+2];if(isNaN(x))continue;
        D[o]=M[0]*x+M[1]*y+M[2]*z;D[o+1]=M[3]*x+M[4]*y+M[5]*z;D[o+2]=M[6]*x+M[7]*y+M[8]*z;}}}
    return d;}
  copyCell(Ps,ks,Pd,kd,M){Pd.T[kd]=Ps.T[ks];Pd.B[kd]=Ps.B[ks];Pd.S[kd]=Ps.S[ks];Pd.K[kd]=Ps.K[ks];Pd.col[kd]=this.cloneCol(Ps.col[ks],M);if(Pd.pid)Pd.pid[kd]=Ps.pid?Ps.pid[ks]:-1;}
  newOcean(P,k){if(P.pid)P.pid[k]=-1;P.T[k]=C.OCEAN_T;P.B[k]=this.t;P.S[k]=0;P.K[k]=0;P.col[k]=this.pool.make(ROCK.OCEAN,this.t,C.OCEAN_T,P.id);}
  deposit(arr,parr,plate,x,y,z,ux,uy,uz,amt,W,k0){
    const d=this.cellRad,lk=this.cube.lookup;
    for(let k=k0;k<W.length;k++){const c=lk(x+ux*d*k,y+uy*d*k,z+uz*d*k);arr[c]+=amt*W[k];parr[c]=plate;}
  }
  resolve(i,x,y,z,p1,k1,p2,k2){
    const P1=this.plates[p1],P2=this.plates[p2],c1=P1.K[k1],c2=P2.K[k2];
    let w;if(c1!==c2)w=c1?0:1;else if(c1)w=P1.T[k1]>=P2.T[k2]?0:1;else w=P1.B[k1]>=P2.B[k2]?0:1;
    const Pw=w===0?P1:P2,Pl=w===0?P2:P1,kl=w===0?k2:k1;
    const dx=Pl.W[0]-Pw.W[0],dy=Pl.W[1]-Pw.W[1],dz=Pl.W[2]-Pw.W[2];
    let rx=dy*z-dz*y,ry=dz*x-dx*z,rz=dx*y-dy*x;const len=Math.hypot(rx,ry,rz);
    if(len>1e-9){rx/=len;ry/=len;rz/=len;
      if(!Pl.K[kl]){
        this.deposit(this.arcAdd,this.arcPlate,Pw.id,x,y,z,rx,ry,rz,C.ARC+0.4*Pl.S[kl],C.ARCW,1);
        const d=this.cellRad;this.trench[this.cube.lookup(x-rx*d,y-ry*d,z-rz*d)]=1;
        {const age=this.t-Pl.B[kl],c=this.cube.lookup(x+rx*d*C.SLAB_INBOARD,y+ry*d*C.SLAB_INBOARD,z+rz*d*C.SLAB_INBOARD);this.slab[this.slabQ*this.N+c]+=this.cube.area[i]*R_KM*R_KM*Math.min(1,Math.sqrt(Math.max(0,age)/80));}
      } else {
        // thickening scales with the share of motion across the boundary: head-on stacks crust, sideways sliding barely does
        const bn=this.bn,nl=Math.hypot(bn[3*i],bn[3*i+1],bn[3*i+2]);
        const f=nl>0.5?Math.abs(rx*bn[3*i]+ry*bn[3*i+1]+rz*bn[3*i+2])/nl:1;
        // all of the consumed continental crust is kept: head-on, it underthrusts and thickens the winner over a few cells inland;
        // sideways, it piles up at the boundary cell itself
        const amt=Pl.T[kl],W0=C.COLLW;for(let k=0;k<4;k++)this._cw[k]=f*W0[k]+(k===0?1-f:0);
        this.deposit(this.collAdd,this.collPlate,Pw.id,x,y,z,rx,ry,rz,amt,this._cw,0);
        const key=Math.min(Pw.id,Pl.id)*4096+Math.max(Pw.id,Pl.id);this.pairs.set(key,(this.pairs.get(key)||0)+1);this.pairLoc.set(key,[x,y,z]);
      }
    }
    Pl.T[kl]=-Pl.T[kl];Pl.clr[kl]=this.steps;
    return w;
  }
  step(){
    const {N,cube,dt,plates,owner,nOwner,bodyIdx}=this,pos=cube.pos,nb=cube.nb,lk=cube.lookup;
    // 1. rotate plates about their Euler poles
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
    // 2. coverage and overlaps
    nOwner.fill(-1);this.arcAdd.fill(0);this.collAdd.fill(0);
    {const b=Math.floor(this.t/C.SLAB_DT),q=((b%C.SLAB_BINS)+C.SLAB_BINS)%C.SLAB_BINS;if(this.slabT[q]!==b*C.SLAB_DT){this.slabT[q]=b*C.SLAB_DT;this.slab.fill(0,q*N,(q+1)*N);}this.slabQ=q;}
    const cand=new Int16Array(5);let ng=0;
    for(let i=0;i<N;i++){
      const x=pos[3*i],y=pos[3*i+1],z=pos[3*i+2];
      let nc=0;cand[nc++]=owner[i];
      for(let k=0;k<4;k++){const q=owner[nb[4*i+k]];let dup=false;for(let m=0;m<nc;m++)if(cand[m]===q){dup=true;break;}if(!dup)cand[nc++]=q;}
      let cov=0,wP=-1,wK=-1;
      for(let m=0;m<nc;m++){const p=cand[m],P=plates[p];if(!P.alive)continue;const Q=P.Q;
        const k=lk(Q[0]*x+Q[3]*y+Q[6]*z,Q[1]*x+Q[4]*y+Q[7]*z,Q[2]*x+Q[5]*y+Q[8]*z);
        if(P.T[k]>0){if(cov===0){wP=p;wK=k;}else if(this.resolve(i,x,y,z,wP,wK,p,k)===1){wP=p;wK=k;}cov++;}
      }
      if(cov===0)this.gaps[ng++]=i;else{nOwner[i]=wP;bodyIdx[i]=wK;}
    }
    // 3. gaps: duplicate samples keep their crust; interior holes extend the plate; true gaps become ridge crust
    let gaps=this.gaps,rest=this.rest;
    for(let pass=0;pass<12&&ng>0;pass++){
      let nr=0;
      for(let g=0;g<ng;g++){const i=gaps[g],x=pos[3*i],y=pos[3*i+1],z=pos[3*i+2];
        const p=owner[i],Po=plates[p];let kOld=-1;
        if(Po&&Po.alive){const Q=Po.Q;kOld=lk(Q[0]*x+Q[3]*y+Q[6]*z,Q[1]*x+Q[4]*y+Q[7]*z,Q[2]*x+Q[5]*y+Q[8]*z);}
        const cleared=kOld>=0&&Po.clr[kOld]>=this.steps-2;
        let best=-1,bc=0,bj=-1;
        for(let k=0;k<4;k++){const j=nb[4*i+k],q=nOwner[j];if(q<0||(cleared&&q===p))continue;let c=0;for(let m=0;m<4;m++)if(nOwner[nb[4*i+m]]===q)c++;if(c>bc){bc=c;best=q;bj=j;}}
        let holeSrc=-1;
        if(kOld>=0&&!cleared){let nf=0,src=-1;for(let m=0;m<4;m++){const kb=nb[4*kOld+m];if(Po.T[kb]>0){nf++;if(src<0||Po.K[kb]===this.wK[i])src=kb;}}if(nf>=3)holeSrc=src;}
        if(cleared&&Po.T[kOld]<0){Po.T[kOld]=-Po.T[kOld];nOwner[i]=p;bodyIdx[i]=kOld;}
        else if(holeSrc>=0){this.copyCell(Po,holeSrc,Po,kOld);nOwner[i]=p;bodyIdx[i]=kOld;}
        else if(kOld>=0&&!cleared&&bc<3){this.newOcean(Po,kOld);nOwner[i]=p;bodyIdx[i]=kOld;}
        else if(best>=0){
          const Pq=plates[best],Q=Pq.Q;const kq=lk(Q[0]*x+Q[3]*y+Q[6]*z,Q[1]*x+Q[4]*y+Q[7]*z,Q[2]*x+Q[5]*y+Q[8]*z);
          if(Pq.T[kq]<=0){let src=bodyIdx[bj];
            for(let m=0;m<4;m++){const j=nb[4*i+m];if(nOwner[j]===best&&Pq.K[bodyIdx[j]]===this.wK[i]){src=bodyIdx[j];break;}}
            this.copyCell(Pq,src,Pq,kq);}
          nOwner[i]=best;bodyIdx[i]=kq;
        } else rest[nr++]=i;
      }
      const tmp=gaps;gaps=rest;rest=tmp;ng=nr;
    }
    for(let g=0;g<ng;g++){const i=gaps[g];nOwner[i]=owner[i]>=0&&plates[owner[i]].alive?owner[i]:0;}
    this.owner=nOwner;this.nOwner=owner;
    const O=this.owner;
    for(let g=0;g<ng;g++){const i=gaps[g],P=plates[O[i]];const Q=P.Q,x=pos[3*i],y=pos[3*i+1],z=pos[3*i+2];
      const k=lk(Q[0]*x+Q[3]*y+Q[6]*z,Q[1]*x+Q[4]*y+Q[7]*z,Q[2]*x+Q[5]*y+Q[8]*z);if(P.T[k]<=0)this.newOcean(P,k);bodyIdx[i]=k;}
    for(let i=0;i<N;i++){const P=plates[O[i]],k=bodyIdx[i],v=P.T[k];if(v<0)P.T[k]=-v;else if(v===0)this.newOcean(P,k);}
    // 4. surface processes on the fixed grid; rock units are pushed or eroded as they happen
    this.gather();
    const {wT,wB,wS,wK,wCol,dT,arcAdd,arcPlate,collAdd,collPlate,trench,pool}=this,t=this.t,r=this.rand;
    // optional record of tectonic thickening per fixed cell (arcs, collisions, hotspots, collapse), read by L1
    const rec=this.tecRec||null; // new crust from arcs and hotspots (what L1's sheet cannot make itself)
    for(let i=0;i<N;i++){
      const a=arcAdd[i];if(a>0&&arcPlate[i]===O[i]){wT[i]+=a;if(rec)rec[i]+=a;if(wK[i]?a>0.25:wT[i]>C.OCEAN_T+4)pool.push(wCol[i],ROCK.ARC,t,a,O[i]);}
      const c=collAdd[i];if(c>0&&collPlate[i]===O[i]){wT[i]+=c;if(!wK[i]&&wT[i]>C.ARC_CONV)wK[i]=1; // continental slices pushed onto seafloor accrete as continent
        if(c>0.4)pool.push(wCol[i],ROCK.OROGEN,t,c*0.4,O[i]);}
      trench[i]*=0.85;
    }
    for(const h of this.hotspots){const i=h.i;
      if(wK[i]){wT[i]+=0.1*dt;if(rec)rec[i]+=0.1*dt;pool.push(wCol[i],ROCK.HOTSPOT,t,0.1*dt,O[i]);}
      else{wT[i]+=C.HOT*dt;if(rec)rec[i]+=C.HOT*dt;pool.push(wCol[i],ROCK.HOTSPOT,t,C.HOT*dt,O[i]);
        for(let k=0;k<4;k++){const j=nb[4*i+k];if(!wK[j]){wT[j]+=0.3*C.HOT*dt;if(rec)rec[j]+=0.3*C.HOT*dt;pool.push(wCol[j],ROCK.HOTSPOT,t,0.3*C.HOT*dt,O[j]);}}}
    }
    this.computeElev();
    const fb=Math.pow((this.contFrac||this.contTarget)/this.contTarget,6),sedConv=C.SED_CONV*fb,collapseP=0.05*dt/fb;
    const E=this.elev,kd=C.KD*dt;dT.fill(0);
    for(let i=0;i<N;i++){
      if(wK[i]){
        let f=0,drop=0;
        for(let k=0;k<4;k++){const j=nb[4*i+k];if(wK[j])f+=wT[j]-wT[i];if(E[j]<E[i])drop+=E[i]-E[j];}
        dT[i]+=kd*f;
        if(E[i]>0&&drop>0){const loss=C.KE*E[i]*dt;dT[i]-=loss;pool.erode(wCol[i],loss);
          // spread sediment over every lower neighbour in proportion to its drop (multiple-flow routing)
          for(let k=0;k<4;k++){const j=nb[4*i+k];if(E[j]>=E[i])continue;const sh=loss*(E[i]-E[j])/drop;
            if(wK[j]){dT[j]+=sh;pool.push(wCol[j],ROCK.CLASTIC,t,sh,O[j]);}else wS[j]+=0.5*sh;}}
      } else if(wT[i]>C.OCEAN_T){const d=(wT[i]-C.OCEAN_T)*C.VOLC_DECAY*dt;dT[i]-=d;pool.erode(wCol[i],d);}
    }
    for(let i=0;i<N;i++){let v=wT[i]+dT[i];
      if(wK[i]){if(v>C.CAP)v=C.CAP;if(v<10)v=10;}
      else{
        if(wS[i]>sedConv){const add=Math.max(0.5,24-v);wK[i]=1;v+=add;pool.push(wCol[i],ROCK.SHELF,t,add,O[i]);wS[i]=0;}
        else if(v>C.ARC_CONV*fb&&arcAdd[i]>0&&arcPlate[i]===O[i])wK[i]=1;
        if(v<5)v=5;if(v>30)v=30;
      }
      wT[i]=v;}
    for(let i=0;i<N;i++){
      if(!wK[i]||wT[i]<C.COLLAPSE||r()>collapseP)continue;
      const c=nb[4*i+Math.floor(r()*4)];
      if(!wK[c]&&O[c]===O[i]){const base=pool.R[wCol[i]*U];wK[c]=1;wT[c]=28;wS[c]=0;wCol[c]=pool.make(base,t,26,O[i]);pool.push(wCol[c],ROCK.CLASTIC,t,2,O[i]);wT[i]-=8;}
    }
    if(fb>1){const pz=0.1*dt*(fb-1);
      for(let i=0;i<N;i++){if(!wK[i]||E[i]>0.3||r()>pz)continue;let coast=false;for(let k=0;k<4;k++)if(!wK[nb[4*i+k]])coast=true;
        if(coast){wK[i]=0;wT[i]=C.OCEAN_T+1;wB[i]=t-60;wS[i]=2;wCol[i]=pool.make(ROCK.OCEAN,t-60,wT[i],O[i]);}}
    }
    for(let i=0;i<N;i++)pool.reconcile(wCol[i],wT[i]);
    this.recordFaults();
    this.t+=dt;this.steps++;
    this.computeElev();this.surface();this.scatter();this.syncUnsampled();
    // 5. bookkeeping and discrete events
    this.updateStats(this.steps%4===0);this.hotspotEvents();this.plateEvents();
    if(this.steps%32===0)this.gc();
  }
  // Plate-grid cells that no fixed cell sampled this step miss every surface process, so they go stale:
  // an old ocean cell can sit inside a continent that grew around it, and pop up as a pinhole when
  // rotation brings it under a fixed cell. Each step, an unsampled cell surrounded by sampled
  // neighbours takes their crust type (copying a neighbour's column if it differs) and their mean thickness.
  syncUnsampled(){
    const {N,plates,owner,bodyIdx,steps,pool}=this,nb=this.cube.nb;
    for(let i=0;i<N;i++)plates[owner[i]].seen[bodyIdx[i]]=steps;
    for(const P of plates){if(!P.alive)continue;const T=P.T,K=P.K,seen=P.seen;
      for(let k=0;k<N;k++){if(T[k]<=0||seen[k]===steps)continue;
        let ns=0,nk=0,sumT=0,src1=-1,src0=-1;
        for(let m=0;m<4;m++){const kb=nb[4*k+m];if(seen[kb]!==steps||T[kb]<=0)continue;ns++;sumT+=T[kb];if(K[kb]){nk++;src1=kb;}else src0=kb;}
        if(ns<3)continue;
        const major=nk*2>ns?1:0;
        if(K[k]!==major&&(major?nk:ns-nk)>=3){const src=major?src1:src0;const B=P.B[k];this.copyCell(P,src,P,k);if(!major)P.B[k]=P.B[src];}
        else{T[k]=sumT/ns;if(P.col[k]>=0)pool.reconcile(P.col[k],T[k]);}
      }}
  }
  // Boundary kinematics: across (convergence +) and along slip, accumulated into the fault record of the crust parcel
  recordFaults(){
    const {N,cube,plates,owner:O,wCol,pool,fpool,bn,dt,t}=this,pos=cube.pos,nb=cube.nb;
    for(let i=0;i<N;i++){
      const p=O[i],P=plates[p],x=pos[3*i],y=pos[3*i+1],z=pos[3*i+2];
      let cnt=0,ac=0,al=0,nx=0,ny=0,nz=0;
      for(let k=0;k<4;k++){const j=nb[4*i+k],q=O[j];if(q===p)continue;const Qp=plates[q];
        let ex=pos[3*j]-x,ey=pos[3*j+1]-y,ez=pos[3*j+2]-z;const d=ex*x+ey*y+ez*z;ex-=d*x;ey-=d*y;ez-=d*z;const el=Math.hypot(ex,ey,ez)||1;ex/=el;ey/=el;ez/=el;
        const dx=Qp.W[0]-P.W[0],dy=Qp.W[1]-P.W[1],dz=Qp.W[2]-P.W[2];const vx=dy*z-dz*y,vy=dz*x-dx*z,vz=dx*y-dy*x;
        const tx=y*ez-z*ey,ty=z*ex-x*ez,tz=x*ey-y*ex;
        ac+=-(vx*ex+vy*ey+vz*ez);al+=vx*tx+vy*ty+vz*tz;nx+=ex;ny+=ey;nz+=ez;cnt++;}
      if(!cnt){bn[3*i]=bn[3*i+1]=bn[3*i+2]=0;continue;}
      const nl=Math.hypot(nx,ny,nz);if(nl<1e-6){bn[3*i]=bn[3*i+1]=bn[3*i+2]=0;continue;}
      nx/=nl;ny/=nl;nz/=nl;bn[3*i]=nx;bn[3*i+1]=ny;bn[3*i+2]=nz;
      ac=ac/cnt*R_KM;al=al/cnt*R_KM; // km/Myr
      const c=wCol[i];let f=pool.flt[c];if(f<0){f=fpool.alloc();pool.flt[c]=f;}
      const D=fpool.D,o=f*FS,e=o+8,ki=kindOf(ac,al);
      if(isNaN(D[e])){D[e]=t;D[e+1]=t;D[e+2]=0;D[e+3]=0;D[e+4]=0;D[e+5]=0;D[e+6]=0;}
      else{let shift=t-D[e+1]>C.FAULT_GAP;
        if(!shift){if(kindOf(D[e+2],D[e+3])!==ki&&Math.hypot(D[e+2],D[e+3])>5){D[o]+=dt;if(D[o]>C.FAULT_SWITCH)shift=true;}else D[o]=0;}
        if(shift){for(let k=0;k<7;k++)D[o+1+k]=D[e+k];D[e]=t;D[e+1]=t;D[e+2]=0;D[e+3]=0;D[e+4]=0;D[e+5]=0;D[e+6]=0;D[o]=0;}}
      // strike = boundary tangent, stored in the plate frame so it travels with the crust
      const sx=y*nz-z*ny,sy=z*nx-x*nz,sz=x*ny-y*nx,Q=P.Q;
      let bx=Q[0]*sx+Q[3]*sy+Q[6]*sz,by=Q[1]*sx+Q[4]*sy+Q[7]*sz,bz=Q[2]*sx+Q[5]*sy+Q[8]*sz;
      if(bx*D[e+4]+by*D[e+5]+bz*D[e+6]<0){bx=-bx;by=-by;bz=-bz;}
      const m=Math.hypot(ac,al)*dt;
      D[e+1]=t;D[e+2]+=ac*dt;D[e+3]+=al*dt;D[e+4]+=bx*m;D[e+5]+=by*m;D[e+6]+=bz*m;
    }
  }
  // net rotation of the whole lithosphere relative to the fixed (mantle) frame:
  // omega_net = Minv * sum_i (I - x x^T) Omega_owner(i) dA_i ; zero means no net rotation
  netRotation(){const {N,cube,plates,owner}=this,p=cube.pos,A=cube.area,M=cube.Minv;let sx=0,sy=0,sz=0;
    for(let i=0;i<N;i++){const W=plates[owner[i]].W,x=p[3*i],y=p[3*i+1],z=p[3*i+2],d=W[0]*x+W[1]*y+W[2]*z,a=A[i];sx+=a*(W[0]-d*x);sy+=a*(W[1]-d*y);sz+=a*(W[2]-d*z);}
    return[M[0]*sx+M[1]*sy+M[2]*sz,M[3]*sx+M[4]*sy+M[5]*sz,M[6]*sx+M[7]*sy+M[8]*sz];}
  removeNetRotation(){const d=this.netRotation();for(const P of this.plates)if(P.alive)for(let k=0;k<3;k++)P.W[k]-=d[k];}
  // displayed surface rock: the rock type with the most thickness in the top kilometre (map generalisation, not the topmost veneer)
  surface(){const {N,wCol,wR,wRA,pool}=this,acc=new Float32Array(8),age=new Float32Array(8);
    for(let i=0;i<N;i++){const c=wCol[i],o=c*U;acc.fill(0);let rem=C.COVER_DEPTH;
      for(let k=pool.cnt[c]-1;k>=0&&rem>0;k--){const th=Math.min(pool.T[o+k],rem),r=pool.R[o+k];if(acc[r]===0)age[r]=pool.B[o+k];acc[r]+=th;rem-=th;}
      let best=pool.R[o],bv=-1;for(let r=1;r<8;r++)if(acc[r]>bv){bv=acc[r];best=r;}wR[i]=best;wRA[i]=age[best];}}
  gc(){
    const {pool,fpool,plates,steps}=this;pool.mark.fill(0);fpool.mark.fill(0);
    for(const P of plates){if(!P.alive)continue;const T=P.T,col=P.col,clr=P.clr;
      for(let k=0;k<this.N;k++){const c=col[k];if(c<0)continue;
        if(T[k]>0||(T[k]<0&&clr[k]>=steps-2))pool.mark[c]=1;else{col[k]=-1;if(T[k]<0)T[k]=0;}}}
    for(let c=0;c<pool.cap;c++){if(!pool.used[c])continue;if(!pool.mark[c]){const f=pool.flt[c];pool.release(c);if(f>=0){fpool.release(f);pool.flt[c]=-1;}}else if(pool.flt[c]>=0)fpool.mark[pool.flt[c]]=1;}
    for(let f=0;f<fpool.cap;f++)if(fpool.used[f]&&!fpool.mark[f])fpool.release(f);
  }
  computeElev(){const {N,wT,wK,wB,wS,elev:E,t}=this;for(let i=0;i<N;i++)E[i]=elevOf(wK[i],wT[i],wB[i],wS[i],t);}
  updateStats(centroids){
    const {N,plates,owner,wK,wB,wCol,pool,elev:E,t}=this,p=this.cube.pos;
    for(const P of plates){P.area=0;P.contArea=0;P._c=[0,0,0];}
    let land=0,cont=0,oa=0,on=0,units=0;
    for(let i=0;i<N;i++){const P=plates[owner[i]];P.area++;if(wK[i]){P.contArea++;cont++;units+=pool.cnt[wCol[i]];}else{oa+=t-wB[i];on++;}if(E[i]>0)land++;
      if(centroids){P._c[0]+=p[3*i];P._c[1]+=p[3*i+1];P._c[2]+=p[3*i+2];}}
    for(const P of plates){if(!P.alive)continue;if(P.area>P.maxArea)P.maxArea=P.area;if(centroids&&P.area>0)P.c=norm3(P._c);}
    this.landFrac=land/N;this.contFrac=cont/N;this.oceanAge=on?oa/on:0;this.meanUnits=cont?units/cont:0;
    // speed: the root mean square over the crust's area of |omega x r| (cm/yr)
    if(centroids){this.netRot=Math.hypot(...this.netRotation())*180/Math.PI;const A=this.cube.area;let ss=0,sw=0;
      for(let i=0;i<N;i++){const W=plates[owner[i]].W,d=W[0]*p[3*i]+W[1]*p[3*i+1]+W[2]*p[3*i+2];ss+=A[i]*(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]-d*d);sw+=A[i];}
      this.rmsSpeed=Math.sqrt(ss/(sw||1))*R_KM/10;}
  }
  hotspotEvents(){for(const h of this.hotspots){const up=this.elev[h.i]>0;
    if(up&&!h.wasUp&&!this.wK[h.i]&&this.t-h.last>15){this.log("island","A volcanic island emerges over a hotspot",h.v);h.last=this.t;}h.wasUp=up;}}
  plateEvents(){
    const {N,dt,rand:r,plates}=this;
    for(const P of plates){if(P.alive&&P.area===0){if(P.maxArea>N*0.01)this.log("lost",`Plate ${P.id} is fully subducted`,P.c);this.kill(P);}}
    this.walkVelocities();
    const decay=Math.pow(0.97,dt);let mergeKey=-1;
    for(const [k,v] of this.pairs){const nv=v*decay;if(nv<0.5){this.pairs.delete(k);continue;}this.pairs.set(k,nv);
      if(nv>C.COLL_LOG&&!this.collLogged.has(k)){this.collLogged.add(k);this.log("collide",`Continents on plates ${Math.floor(k/4096)} and ${k%4096} begin to collide; mountains rise`,this.pairLoc.get(k));}
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
  // angular velocities: random walk in direction, speed relaxes toward each plate's preferred speed
  // (an unconstrained 3D random walk drifts outward in magnitude until it hits the cap)
  walkVelocities(){const {dt,rand:r,plates}=this;
    const sig=0.08*10/R_KM*Math.sqrt(dt),relax=Math.min(1,dt/C.SPEED_TAU);
    for(const P of plates){if(!P.alive)continue;for(let k=0;k<3;k++)P.W[k]+=gauss(r)*sig;
      let w=Math.hypot(P.W[0],P.W[1],P.W[2]);
      if(w>1e-12){const f=1+(P.w0/w-1)*relax;for(let k=0;k<3;k++)P.W[k]*=f;w*=f;}
      if(w>this.wmax)for(let k=0;k<3;k++)P.W[k]*=this.wmax/w;
      else if(w<this.wmin){const a=randUnit(r);for(let k=0;k<3;k++)P.W[k]=a[k]*this.wmin;}}
    if(this.opts.nnr!==false)this.removeNetRotation();}
  merge(a,b){
    const purge=()=>{for(const k of [...this.pairs.keys()]){const x=Math.floor(k/4096),y=k%4096;if(x===a.id||y===a.id||x===b.id||y===b.id)this.pairs.delete(k);}};
    if(!a.alive||!b.alive){purge();return;}
    const [Wn,L]=a.area>=b.area?[a,b]:[b,a];const {N,owner,bodyIdx}=this,p=this.cube.pos,lk=this.cube.lookup,QL=L.Q,QW=Wn.Q;
    const M=new Float64Array(9);for(let r=0;r<3;r++)for(let c=0;c<3;c++)M[r*3+c]=QW[r]*QL[c]+QW[3+r]*QL[3+c]+QW[6+r]*QL[6+c]; // QW^T QL
    // optional provenance for L1's lookahead: which absorbed-plate body cell filled each surviving-plate cell
    const from=this.mergeLog?new Int32Array(N).fill(-1):null;if(from)this.mergeLog.push({from:L.id,into:Wn.id,t:this.t,map:from});
    // pull: every empty cell of the surviving plate's grid asks which absorbed-plate cell lies under it.
    // (Pushing each absorbed cell to its nearest target is many-to-one and leaves a scatter of empty cells.)
    // Cells of the surviving plate's grid that now sit on the absorbed plate's territory are overwritten too:
    // they are stale leftovers (never sampled while the other plate covered them), not real crust.
    for(let kk=0;kk<N;kk++){const bx=p[3*kk],by=p[3*kk+1],bz=p[3*kk+2];
      const x=QW[0]*bx+QW[1]*by+QW[2]*bz,y=QW[3]*bx+QW[4]*by+QW[5]*bz,z=QW[6]*bx+QW[7]*by+QW[8]*bz;
      if(Wn.T[kk]>0&&owner[lk(x,y,z)]!==L.id)continue;
      const k=lk(QL[0]*x+QL[3]*y+QL[6]*z,QL[1]*x+QL[4]*y+QL[7]*z,QL[2]*x+QL[5]*y+QL[8]*z);if(L.T[k]>0){this.copyCell(L,k,Wn,kk,M);if(from)from[kk]=k;}}
    for(let i=0;i<N;i++){if(owner[i]!==L.id)continue;const x=p[3*i],y=p[3*i+1],z=p[3*i+2];
      const kk=lk(QW[0]*x+QW[3]*y+QW[6]*z,QW[1]*x+QW[4]*y+QW[7]*z,QW[2]*x+QW[5]*y+QW[8]*z);
      if(Wn.T[kk]<=0){this.copyCell(L,bodyIdx[i],Wn,kk,M);if(from)from[kk]=bodyIdx[i];}owner[i]=Wn.id;bodyIdx[i]=kk;}
    const tot=Wn.area+L.area;for(let k=0;k<3;k++)Wn.W[k]=(Wn.W[k]*Wn.area+L.W[k]*L.area)/tot;
    Wn.area=tot;Wn.contArea+=L.contArea;
    this.log("suture",`Plates ${Wn.id} and ${L.id} weld into one along a suture`,Wn.c);
    L.mergedInto=Wn.id;this.kill(L);purge();this.gather();this.gc();
  }
  split(P){
    const {N,owner,bodyIdx,rand:r,noise:nz}=this,p=this.cube.pos;
    const c=P.c;let n=randUnit(r);const d=n[0]*c[0]+n[1]*c[1]+n[2]*c[2];n=norm3([n[0]-d*c[0],n[1]-d*c[1],n[2]-d*c[2]]);
    const a=norm3(cross(c,n));const Qp=this.addPlate(P.W);Qp.Q.set(P.Q);Qp.parent=P.id;const Q=P.Q;
    const moved=[];let contA=0,contB=0;
    for(let k=0;k<N;k++){if(P.T[k]<=0)continue;const bx=p[3*k],by=p[3*k+1],bz=p[3*k+2];
      const x=Q[0]*bx+Q[1]*by+Q[2]*bz,y=Q[3]*bx+Q[4]*by+Q[5]*bz,z=Q[6]*bx+Q[7]*by+Q[8]*bz;
      const s=x*n[0]+y*n[1]+z*n[2]+(nz(x+40,y,z,2.2,4)-0.5)*0.5;
      if(s>0){Qp.T[k]=P.T[k];Qp.B[k]=P.B[k];Qp.S[k]=P.S[k];Qp.K[k]=P.K[k];Qp.col[k]=P.col[k];if(P.pid){Qp.pid[k]=P.pid[k];P.pid[k]=-1;}P.T[k]=0;P.col[k]=-1;moved.push(k);if(Qp.K[k])contB++;}
      else if(P.K[k])contA++;}
    let nNew=0;for(let i=0;i<N;i++)if(owner[i]===P.id&&Qp.T[bodyIdx[i]]>0)nNew++;
    if(nNew<60||P.area-nNew<60){for(const k of moved){P.T[k]=Qp.T[k];P.B[k]=Qp.B[k];P.S[k]=Qp.S[k];P.K[k]=Qp.K[k];P.col[k]=Qp.col[k];if(P.pid)P.pid[k]=Qp.pid[k];}this.kill(Qp);return;}
    for(let i=0;i<N;i++)if(owner[i]===P.id&&Qp.T[bodyIdx[i]]>0)owner[i]=Qp.id;
    const s=(1.5+2.5*r())*10/R_KM;for(let k=0;k<3;k++){Qp.W[k]=P.W[k]+a[k]*s/2;P.W[k]-=a[k]*s/2;}
    this.log("rift",`${contA>40&&contB>40?"Continental rift":"Rift"}: plate ${P.id} splits and plate ${Qp.id} drifts away`,c);
    this.updateStats(true);
  }
  // Save and restore the whole state, so the future can be previewed (L1 looks ahead for sutures) and then
  // replayed identically. The cube is shared and the noise function is stateless; the generator's state is saved.
  snapshot(){const d={};for(const k of Object.keys(this))if(k!=="cube"&&k!=="rand"&&k!=="noise")d[k]=this[k];
    return{data:structuredClone(d),rand:this.rand.getState()};}
  restore(snap){const d=structuredClone(snap.data);Object.setPrototypeOf(d.pool,ColumnPool.prototype);Object.setPrototypeOf(d.fpool,FaultPool.prototype);
    for(const k of Object.keys(this))if(k!=="cube"&&k!=="rand"&&k!=="noise"&&!(k in d))delete this[k];
    Object.assign(this,d);this.rand.setState(snap.rand);}
  // ---------- inspection and L1 handoff ----------
  velocityAt(i){const P=this.plates[this.owner[i]],p=this.cube.pos,x=p[3*i],y=p[3*i+1],z=p[3*i+2];
    const v=[P.W[1]*z-P.W[2]*y,P.W[2]*x-P.W[0]*z,P.W[0]*y-P.W[1]*x];
    const east=norm3([-y,x,0]),north=cross([x,y,z],east);const ve=(v[0]*east[0]+v[1]*east[1]+v[2]*east[2])*R_KM/10,vn=(v[0]*north[0]+v[1]*north[1]+v[2]*north[2])*R_KM/10;
    return{speed:Math.hypot(ve,vn),bearing:((Math.atan2(ve,vn)*180/Math.PI)+360)%360};}
  faultEpisodes(i){const f=this.pool.flt[this.wCol[i]];if(f<0)return[];const D=this.fpool.D,out=[];
    for(const o of [f*FS+1,f*FS+8]){if(isNaN(D[o]))continue;const X=D[o+2],S=D[o+3];
      out.push({kind:FK_NAMES[kindOf(X,S)],t0:D[o],t1:D[o+1],across:X,along:S,sense:S<0?"right-lateral":"left-lateral",active:this.t-D[o+1]<3});}
    return out;}
  // what the inspector shows for a cell: its rock column (bottom first) and fault episodes
  columnInfo(i){return{units:this.pool.units(this.wCol[i]),eps:this.faultEpisodes(i)};}
  faultCellInfo(i){const f=this.pool.flt[this.wCol[i]];if(f<0)return null;const D=this.fpool.D,o=f*FS+8,op=f*FS+1;if(isNaN(D[o]))return null;
    let tot=Math.hypot(D[o+2],D[o+3]);if(!isNaN(D[op]))tot+=Math.hypot(D[op+2],D[op+3]);
    return{kind:kindOf(D[o+2],D[o+3]),slip:tot,active:this.t-D[o+1]<3};}
  extractFaults(minSlip){
    minSlip=minSlip||60;
    const {N,cube,pool,fpool,owner,wCol,wK,plates,t}=this,pos=cube.pos,nb=cube.nb,D=fpool.D;
    const kind=new Uint8Array(N),str=new Float32Array(3*N),X=new Float32Array(N),S=new Float32Array(N),act=new Uint8Array(N),t0=new Float32Array(N),t1=new Float32Array(N);
    for(let i=0;i<N;i++){if(!wK[i])continue;const f=pool.flt[wCol[i]];if(f<0)continue;const o=f*FS+8,op=f*FS+1;if(isNaN(D[o]))continue;
      let tot=Math.hypot(D[o+2],D[o+3]);const hasP=!isNaN(D[op]);if(hasP)tot+=Math.hypot(D[op+2],D[op+3]);if(tot<minSlip)continue;
      const Q=plates[owner[i]].Q,bx=D[o+4],by=D[o+5],bz=D[o+6];
      let sx=Q[0]*bx+Q[1]*by+Q[2]*bz,sy=Q[3]*bx+Q[4]*by+Q[5]*bz,sz=Q[6]*bx+Q[7]*by+Q[8]*bz;const sl=Math.hypot(sx,sy,sz);if(sl<1e-9)continue;
      kind[i]=kindOf(D[o+2],D[o+3]);str[3*i]=sx/sl;str[3*i+1]=sy/sl;str[3*i+2]=sz/sl;X[i]=D[o+2];S[i]=D[o+3];act[i]=t-D[o+1]<3?1:0;t0[i]=hasP?D[op]:D[o];t1[i]=D[o+1];}
    const sdot=(a,b)=>Math.abs(str[3*a]*str[3*b]+str[3*a+1]*str[3*b+1]+str[3*a+2]*str[3*b+2]);
    const comp=new Int32Array(N).fill(-1),comps=[];
    for(let i=0;i<N;i++){if(!kind[i]||comp[i]>=0)continue;const id=comps.length,cells=[i];comp[i]=id;
      for(let h=0;h<cells.length;h++){const a=cells[h];for(let k=0;k<4;k++){const j=nb[4*a+k];if(kind[j]&&comp[j]<0&&sdot(a,j)>0.7){comp[j]=id;cells.push(j);}}}
      comps.push(cells);}
    const owned=new Int32Array(N).fill(-1),faults=[],cellFault=new Int32Array(N).fill(-1);
    const P3=i=>[pos[3*i],pos[3*i+1],pos[3*i+2]];
    for(let ci=0;ci<comps.length;ci++){const cells=comps[ci];if(cells.length<4)continue;
      let guard=0;
      while(guard++<50){
        const left=cells.filter(c=>owned[c]<0);if(left.length<3)break;
        let cx=0,cy=0,cz=0;for(const c of left){cx+=pos[3*c];cy+=pos[3*c+1];cz+=pos[3*c+2];}
        let start=left[0],bd=-1;for(const c of left){const d=(pos[3*c]-cx/left.length)**2+(pos[3*c+1]-cy/left.length)**2+(pos[3*c+2]-cz/left.length)**2;if(d>bd){bd=d;start=c;}}
        const fid=faults.length,path=[start];
        const claim=c=>{if(owned[c]<0)owned[c]=fid;for(let k=0;k<4;k++){const j=nb[4*c+k];if(comp[j]===ci&&owned[j]<0)owned[j]=fid;}};
        let dir=[str[3*start],str[3*start+1],str[3*start+2]];
        const toC=[cx/left.length-pos[3*start],cy/left.length-pos[3*start+1],cz/left.length-pos[3*start+2]];
        if(dir[0]*toC[0]+dir[1]*toC[1]+dir[2]*toC[2]<0)dir=[-dir[0],-dir[1],-dir[2]];
        claim(start);let cur=start;
        for(let it=0;it<2000;it++){
          const seen=new Map([[cur,0]]),q=[cur];
          for(let h=0;h<q.length;h++){const a=q[h],dp=seen.get(a);if(dp>=3)continue;for(let k=0;k<4;k++){const j=nb[4*a+k];if(comp[j]===ci&&!seen.has(j)){seen.set(j,dp+1);q.push(j);}}}
          let best=-1,bs=0;
          for(const [c] of seen){if(c===cur||(owned[c]>=0&&owned[c]!==fid)||path.includes(c))continue;
            const dx=pos[3*c]-pos[3*cur],dy=pos[3*c+1]-pos[3*cur+1],dz=pos[3*c+2]-pos[3*cur+2],dist=Math.hypot(dx,dy,dz);
            const al=(dx*dir[0]+dy*dir[1]+dz*dir[2])/dist;if(al<0.45)continue;const sc=dist*al;if(sc>bs){bs=sc;best=c;}}
          if(best<0)break;
          const nd=norm3([pos[3*best]-pos[3*cur],pos[3*best+1]-pos[3*cur+1],pos[3*best+2]-pos[3*cur+2]]);
          let sb=[str[3*best],str[3*best+1],str[3*best+2]];if(sb[0]*nd[0]+sb[1]*nd[1]+sb[2]*nd[2]<0)sb=[-sb[0],-sb[1],-sb[2]];
          dir=norm3([nd[0]*0.6+sb[0]*0.4,nd[1]*0.6+sb[1]*0.4,nd[2]*0.6+sb[2]*0.4]);
          path.push(best);claim(best);
          for(const [c,dp] of seen)if(dp<=1&&owned[c]<0)owned[c]=fid;
          cur=best;
        }
        if(path.length<2){owned[start]=-2;continue;}
        // smooth the trace a little
        const pts=path.map(P3),sm=pts.map((p,k)=>{if(k===0||k===pts.length-1)return p;const a=pts[k-1],b=pts[k+1];return norm3([a[0]+2*p[0]+b[0],a[1]+2*p[1]+b[1],a[2]+2*p[2]+b[2]]);});
        let len=0;for(let k=1;k<sm.length;k++)len+=Math.hypot(sm[k][0]-sm[k-1][0],sm[k][1]-sm[k-1][1],sm[k][2]-sm[k-1][2])*R_KM;
        faults.push({id:fid,kind:kind[start],trace:sm,len,cells:[]});
      }
      for(const c of cells){const f=owned[c];if(f>=0&&faults[f]){faults[f].cells.push(c);cellFault[c]=f;}}
    }
    const out=[];
    for(const F of faults){if(!F.cells.length)continue;let sX=0,sS=0,na=0,first=Infinity,last=-Infinity;
      for(const c of F.cells){sX+=X[c];sS+=S[c];na+=act[c];if(t0[c]<first)first=t0[c];if(t1[c]>last)last=t1[c];}
      const n=F.cells.length;if(n<4||F.len<this.cellKm*3)continue;const fk=kindOf(sX,sS);
      out.push({id:F.id,kind:fk,kindName:FK_NAMES[fk],sense:fk===FK.STRIKE?(sS<0?"right-lateral":"left-lateral"):null,active:na>n/2,
        startMa:t-first,endMa:Math.max(0,t-last),slipAcrossKm:sX/n,slipAlongKm:sS/n,lengthKm:F.len,cells:n,trace:F.trace});}
    const remap=new Map(out.map((f,k)=>[f.id,k]));out.forEach((f,k)=>f.id=k);
    for(let i=0;i<N;i++)if(cellFault[i]>=0)cellFault[i]=remap.has(cellFault[i])?remap.get(cellFault[i]):-1;
    return{t,faults:out,cellFault};
  }
  // Correlate rock units between neighbouring columns into formations, group columns into domains,
  // and build band-limited fields (smoothed elevation, smoothed formation thickness) for the next level.
  buildL1(){
    const {N,cube,pool,wCol,wK,owner,elev,t}=this,nb=cube.nb,R=pool.R,B=pool.B,T=pool.T,PP=pool.P,cnt=pool.cnt;
    const parent=new Int32Array(N*U).fill(-1),dpar=new Int32Array(N);
    for(let i=0;i<N;i++){dpar[i]=i;if(!wK[i])continue;const n=cnt[wCol[i]];for(let k=0;k<n;k++)parent[i*U+k]=i*U+k;}
    const find=(a,P)=>{while(P[a]!==a){P[a]=P[P[a]];a=P[a];}return a;};
    const uni=(a,b,P)=>{a=find(a,P);b=find(b,P);if(a!==b)P[a]=b;};
    const tot=new Float32Array(N);for(let i=0;i<N;i++)if(wK[i])tot[i]=pool.sum(wCol[i]);
    for(let i=0;i<N;i++){if(!wK[i])continue;const ci=wCol[i],oi=ci*U,ni=cnt[ci];
      for(let q=0;q<4;q++){const j=nb[4*i+q];if(!wK[j]||j===i)continue;const cj=wCol[j],oj=cj*U,nj=cnt[cj];let p=0,m=0;
        for(let a=0;a<ni;a++)for(let b=p;b<nj;b++){
          const ra=oi+a,rb=oj+b;if(R[ra]!==R[rb]||PP[ra]!==PP[rb])continue;
          const tol=Math.max(40,0.1*(t-Math.min(B[ra],B[rb])));if(Math.abs(B[ra]-B[rb])>tol)continue;
          uni(i*U+a,j*U+b,parent);m+=Math.min(T[ra],T[rb]);p=b+1;break;}
        if(owner[i]===owner[j]&&m>=0.4*Math.max(tot[i],tot[j]))uni(i,j,dpar);}}
    // formations
    const fidOf=new Map(),forms=[],fId=new Int32Array(N*U).fill(-1),fTh=new Float32Array(N*U);
    for(let i=0;i<N;i++){if(!wK[i])continue;const c=wCol[i],o=c*U;
      for(let k=0;k<cnt[c];k++){const r=find(i*U+k,parent);let f=fidOf.get(r);
        if(f===undefined){f=forms.length;fidOf.set(r,f);forms.push({id:f,rock:R[o+k],plate:PP[o+k],bornMin:Infinity,bornMax:-Infinity,bSum:0,thSum:0,cells:0,basal:0});}
        const F=forms[f],th=T[o+k];F.bornMin=Math.min(F.bornMin,B[o+k]);F.bornMax=Math.max(F.bornMax,B[o+k]);F.bSum+=B[o+k]*th;F.thSum+=th;F.cells++;if(k===0)F.basal++;
        fId[i*U+k]=f;fTh[i*U+k]=th;}}
    const keys=new Float64Array(forms.length),rockOf=new Uint8Array(forms.length);
    for(const F of forms){F.born=F.thSum>0?F.bSum/F.thSum:F.bornMin;keys[F.id]=(F.basal>F.cells/2?-1e7:0)+F.born;rockOf[F.id]=F.rock;}
    // domains
    const domain=new Int32Array(N).fill(-1),dmap=new Map();let nd=0;
    for(let i=0;i<N;i++){if(!wK[i])continue;const r=find(i,dpar);let d=dmap.get(r);if(d===undefined){d=nd++;dmap.set(r,d);}domain[i]=d;}
    // absorb slivers: domains under 8 columns join the neighbouring domain they share the most edges with
    for(let pass=0;pass<3;pass++){const size=new Int32Array(nd);for(let i=0;i<N;i++)if(domain[i]>=0)size[domain[i]]++;
      for(let i=0;i<N;i++){const d=domain[i];if(d<0||size[d]>=8)continue;let best=-1,bc=0;const votes=new Map();
        for(let k=0;k<4;k++){const e=domain[nb[4*i+k]];if(e<0||e===d||size[e]<size[d])continue;const v=(votes.get(e)||0)+1;votes.set(e,v);if(v>bc){bc=v;best=e;}}
        if(best>=0){domain[i]=best;size[d]--;size[best]++;}}}
    {const remap=new Map();let k=0;for(let i=0;i<N;i++){const d=domain[i];if(d<0)continue;if(!remap.has(d))remap.set(d,k++);domain[i]=remap.get(d);}nd=k;}
    // smoothed elevation: three passes of neighbour averaging, removing structure below ~3 cells
    let Es=Float32Array.from(elev),tmp=new Float32Array(N);
    for(let pass=0;pass<3;pass++){for(let i=0;i<N;i++){let s=Es[i];for(let k=0;k<4;k++)s+=Es[nb[4*i+k]];tmp[i]=s/5;}[Es,tmp]=[tmp,Es];}
    // smoothed formation stacks: thickness averaged over neighbours in the same domain (absent = zero, so specks fade and edges pinch out)
    let cur=new Array(N);
    for(let i=0;i<N;i++){if(!wK[i]){cur[i]=null;continue;}const n=cnt[wCol[i]],m=new Map();for(let k=0;k<n;k++){const f=fId[i*U+k];m.set(f,(m.get(f)||0)+fTh[i*U+k]);}cur[i]=m;}
    for(let pass=0;pass<2;pass++){const nxt=new Array(N);
      for(let i=0;i<N;i++){if(!cur[i]){nxt[i]=null;continue;}const m=new Map();let w=1;
        for(const [f,th] of cur[i])m.set(f,th);
        for(let k=0;k<4;k++){const j=nb[4*i+k];if(!cur[j]||domain[j]!==domain[i])continue;w++;for(const [f,th] of cur[j])m.set(f,(m.get(f)||0)+th);}
        const o=new Map();for(const [f,th] of m){const v=th/w;if(v>0.02)o.set(f,v);}nxt[i]=o;}
      cur=nxt;}
    const sOff=new Int32Array(N+1);let tot2=0;for(let i=0;i<N;i++){sOff[i]=tot2;if(cur[i])tot2+=Math.min(12,cur[i].size);}sOff[N]=tot2;
    const sId=new Int32Array(tot2),sTh=new Float32Array(tot2);
    for(let i=0;i<N;i++){if(!cur[i])continue;const a=[...cur[i]].sort((x,y)=>keys[x[0]]-keys[y[0]]);const take=a.slice(Math.max(0,a.length-12));
      take.forEach(([f,th],k)=>{sId[sOff[i]+k]=f;sTh[sOff[i]+k]=th;});}
    return{steps:this.steps,t,forms,keys,rockOf,fId,fTh,domain,nDomains:nd,Es,sOff,sId,sTh};
  }
  exportHandoff(h){
    const {N,cube,pool,owner,wK,elev,plates,t}=this,p=cube.pos;
    const r1=v=>Math.round(v*10)/10,r2=v=>Math.round(v*100)/100;
    const L=h&&h.l1;
    const lat=new Array(N),lon=new Array(N),plate=new Array(N),el=new Array(N),els=new Array(N),cont=new Array(N),fault=new Array(N),dom=new Array(N),off=new Array(N+1),units=[];
    for(let i=0;i<N;i++){const ll=latlon(p[3*i],p[3*i+1],p[3*i+2]);lat[i]=r2(ll[0]);lon[i]=r2(ll[1]);plate[i]=owner[i];el[i]=r2(elev[i]);els[i]=L?r2(L.Es[i]):el[i];cont[i]=wK[i];fault[i]=h?h.cellFault[i]:-1;dom[i]=L?L.domain[i]:-1;
      off[i]=units.length/4;const c=this.wCol[i],o=c*U;for(let k=0;k<pool.cnt[c];k++)units.push(pool.R[o+k],Math.round(t-pool.B[o+k]),r1(pool.T[o+k]),L&&wK[i]?L.fId[i*U+k]:-1);}
    off[N]=units.length/4;
    return{
      meta:{generator:"L0 plate tectonics sandbox",seed:this.opts.seed,settings:{resolution:this.n,initialPlates:this.opts.nPlates,start:this.opts.mode},
        snapshotMyr:r1(t),units:{distance:"km",time:"Myr before snapshot (ma)",elevation:"km",velocity:"cm/yr"},
        grid:{type:"equal-angle cube-sphere",facesCells:this.n,cells:N,radiusKm:R_KM,cellIndex:"face*n*n + j*n + i"},
        rockTypes:ROCK_NAMES,faultKinds:FK_NAMES,
        notes:"Columns list units bottom to top as [rock, formed_ma, thickness_km, formation_id] (formation -1 on oceanic crust). Formations correlate units across neighbouring continental columns; interpolate each formation's thickness only within a domain. smoothedElevationKm is band-limited to about three cells: interpolate it, threshold at your own resolution, and add your own detail. Faults are extracted from continental crust only. Plate rotations are Euler poles; velocity at a point is omega x position."},
      plates:plates.filter(P=>P.alive).map(P=>{const w=Math.hypot(...P.W),ll=latlon(P.W[0]/w,P.W[1]/w,P.W[2]/w);
        return{id:P.id,eulerPoleLat:r2(ll[0]),eulerPoleLon:r2(ll[1]),rateDegPerMyr:+(w*180/Math.PI).toFixed(4),maxSpeedCmPerYr:r2(w*R_KM/10),areaFraction:+(P.area/N).toFixed(4),continentalFraction:+(P.contArea/N).toFixed(4)};}),
      events:this.events.map(e=>({ma:r1(t-e.t),type:e.type,text:e.text,lat:e.lat,lon:e.lon})),
      faults:h?h.faults.map(f=>({id:f.id,kind:f.kindName,sense:f.sense,active:f.active,startMa:r1(f.startMa),endMa:r1(f.endMa),slipAcrossKm:r1(f.slipAcrossKm),slipAlongKm:r1(f.slipAlongKm),lengthKm:Math.round(f.lengthKm),
        trace:f.trace.map(v=>{const ll=latlon(v[0],v[1],v[2]);return[r2(ll[0]),r2(ll[1])];})})):[],
      formations:L?L.forms.map(F=>({id:F.id,rock:ROCK_NAMES[F.rock],formedFromMa:r1(t-F.bornMin),formedToMa:r1(t-F.bornMax),originPlate:F.plate,columns:F.cells,meanThicknessKm:r2(F.thSum/F.cells)})):[],
      domains:L?L.nDomains:0,
      cells:{lat,lon,plate,elevationKm:el,smoothedElevationKm:els,continental:cont,domain:dom,fault,columnOffsets:off,columnUnits:units}
    };
  }
}
return {World,C,ROCK,ROCK_NAMES,FK,FK_NAMES,R_KM,getCube,U,FS,makeNoise3,makeLookup,elevOf,hashStr,
  ColumnPool,FaultPool,kindOf,gauss,randUnit,norm3,cross,latlon,loc,mulberry32}; // the last row is shared with the marker core (mesh.js)
}
const SIM=SIMmodule();
// ===== SIM-END =====
