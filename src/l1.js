// ===== L1-START =====
// L1: one landmass at the present, refined, with the last few tens of Myr of its tectonic history replayed as a
// thin viscous sheet: where plates converge, its crust shortens and thickens into mountain belts and plateaus at L1
// resolution, carried by markers with L0's rock columns. A short frozen erosion pass then settles its drainage
// basins, the tiles for L2.
const L1 = (() => {
const {R_KM,C,makeLookup,elevOf,makeNoise3,hashStr}=SIM;
const ROOT=1-C.ISO, EPS=1e-6; // share of a thickness change taken by the root under full isostasy; priority-flood step (km)

// ---------- generic graph tools (also used on a planar grid by the tests) ----------
function Heap(cap){this.k=new Float64Array(cap);this.v=new Int32Array(cap);this.n=0;}
Heap.prototype.push=function(v,key){const K=this.k,V=this.v;let i=this.n++;while(i>0){const p=(i-1)>>1;if(K[p]<=key)break;K[i]=K[p];V[i]=V[p];i=p;}K[i]=key;V[i]=v;};
Heap.prototype.pop=function(){const K=this.k,V=this.v,top=V[0],n=--this.n;
  if(n>0){const key=K[n],v=V[n];let i=0;for(;;){let c=2*i+1;if(c>=n)break;if(c+1<n&&K[c+1]<K[c])c++;if(K[c]>=key)break;K[i]=K[c];V[i]=V[c];i=c;}K[i]=key;V[i]=v;}
  return top;};
function scratch(nC){return{f:new Float64Array(nC),closed:new Uint8Array(nC),rec:new Int32Array(nC),rdist:new Float32Array(nC),stack:new Int32Array(nC),
  ndon:new Int32Array(nC+1),don:new Int32Array(nC),A:new Float64Array(nC),hN:new Float64Array(nC),heap:new Heap(nC),pit:new Int32Array(nC),ns:0};}
// Route flow over a graph with 8 neighbour slots per cell (-1 = none), plus optional extra links X in CSR form
// ({off, adj, dist}: contacts between pieces of crust). Outlets are fixed base level. Cells flagged in S.seaEdge
// border open sea that L1 does not represent; they drain to a virtual sea-level neighbour one cell away.
// Depressions are filled with an epsilon gradient (Barnes et al. 2014, Priority-Flood+epsilon) so every
// cell drains to an outlet; receivers are steepest descent on the filled surface; area accumulates down the tree.
function route(G,e,alive,outlet,S,precip,X){
  const {nC,nbr,dist,area}=G,{f,closed,rec,rdist,stack,ndon,don,A,heap,pit}=S,sea=S.seaEdge;
  const xo=X?X.off:null,xa=X?X.adj:null,xd=X?X.dist:null;
  heap.n=0;let ph=0,pt=0;
  for(let l=0;l<nC;l++){closed[l]=0;if(!alive[l]){closed[l]=1;rec[l]=-1;}else if(outlet[l]||(sea&&sea[l])){closed[l]=1;f[l]=e[l];}}
  for(let l=0;l<nC;l++){if(!alive[l]||!closed[l])continue;let push=false;
    for(let k=0;k<8&&!push;k++){const j=nbr[8*l+k];if(j>=0&&alive[j]&&!closed[j])push=true;}
    if(!push&&xo)for(let q=xo[l];q<xo[l+1]&&!push;q++){const j=xa[q];if(alive[j]&&!closed[j])push=true;}
    if(push)heap.push(l,f[l]);}
  // Zhou et al. (2016) variant: a cell strictly higher than a finalized neighbour is final at its own height,
  // so whole hillslopes are swept uphill with a plain queue; only cells that still have lower, unresolved
  // neighbours wait in the heap. Same filled surface as the plain algorithm, far fewer heap operations.
  const up=S.up||(S.up=new Int32Array(nC));
  const visit=(c,fc,j)=>{if(j<0||closed[j])return;closed[j]=1;
      if(e[j]<=fc){f[j]=fc+EPS;pit[pt++]=j;return;}
      f[j]=e[j];let uh=0,ut=0;up[ut++]=j;
      while(uh<ut){const a=up[uh++],fa=f[a];let defer=false;
        for(let q=0;q<8;q++){const b=nbr[8*a+q];if(b<0||closed[b])continue;
          if(e[b]>fa){closed[b]=1;f[b]=e[b];up[ut++]=b;}else defer=true;}
        if(xo)for(let q=xo[a];q<xo[a+1];q++){const b=xa[q];if(closed[b])continue;
          if(e[b]>fa){closed[b]=1;f[b]=e[b];up[ut++]=b;}else defer=true;}
        if(defer)heap.push(a,fa);}};
  const flood=()=>{while(ph<pt||heap.n>0){const c=ph<pt?pit[ph++]:heap.pop(),fc=f[c];
    for(let k=0;k<8;k++)visit(c,fc,nbr[8*c+k]);
    if(xo)for(let q=xo[c];q<xo[c+1];q++)visit(c,fc,xa[q]);}};
  flood();
  // land the flood never reached is cut off from the sea (e.g. enclosed by consumed parcels): each such piece
  // becomes an inland basin draining to its lowest cell
  const sink=S.sink||(S.sink=new Uint8Array(nC));sink.fill(0);S.sinks=0;
  for(;;){let lo=-1;for(let l=0;l<nC;l++)if(!closed[l]&&(lo<0||e[l]<e[lo]))lo=l;if(lo<0)break;
    closed[lo]=1;sink[lo]=1;f[lo]=e[lo];heap.push(lo,f[lo]);S.sinks++;flood();}
  for(let l=0;l<nC;l++){if(!alive[l])continue;
    if(outlet[l]||sink[l]||(sea&&sea[l])){rec[l]=l;rdist[l]=sea&&sea[l]?Math.sqrt(area[l]):0;continue;}
    let best=l,bs=0,bd=0;
    for(let k=0;k<8;k++){const j=nbr[8*l+k];if(j<0||!alive[j])continue;const d=dist[8*l+k],s=(f[l]-f[j])/d;if(s>bs){bs=s;best=j;bd=d;}}
    if(xo)for(let q=xo[l];q<xo[l+1];q++){const j=xa[q];if(!alive[j])continue;const d=xd[q],s=(f[l]-f[j])/d;if(s>bs){bs=s;best=j;bd=d;}}
    rec[l]=best;rdist[l]=bd;}
  // stack: receivers before donors (breadth-first from each outlet)
  ndon.fill(0);for(let l=0;l<nC;l++){const r=rec[l];if(r>=0&&r!==l)ndon[r+1]++;}
  for(let l=0;l<nC;l++)ndon[l+1]+=ndon[l];
  const fill=S.fillp||(S.fillp=new Int32Array(nC));fill.set(ndon.subarray(0,nC));
  for(let l=0;l<nC;l++){const r=rec[l];if(r>=0&&r!==l)don[fill[r]++]=l;}
  let ns=0;for(let l=0;l<nC;l++)if(rec[l]===l)stack[ns++]=l;
  for(let h=0;h<ns;h++){const c=stack[h];for(let q=ndon[c];q<ndon[c+1];q++)stack[ns++]=don[q];}
  S.ns=ns;
  for(let h=0;h<ns;h++){const l=stack[h];A[l]=area[l]*(precip?precip[l]:1);}
  for(let h=ns-1;h>=0;h--){const l=stack[h],r=rec[l];if(r!==l)A[r]+=A[l];}
  return S;
}
// Implicit stream-power incision, n = 1 (Braun & Willett 2013): dh/dt = -K A^m S. Units: km, Myr, A in km^2.
// Outlets hold base level max(e, 0): sea level, or their own height for a land cell no outlet could reach.
// Sea-edge cells erode toward a virtual sea-level receiver one cell away.
// A cell lower than its receiver's new height (a filled depression) is left alone. Returns hN, the new heights.
function streamPower(e,S,K,m,dt){
  const {rec,rdist,stack,A,hN,ns}=S,sea=S.seaEdge;
  for(let h=0;h<ns;h++){const l=stack[h],r=rec[l];
    if(r===l){if(sea&&sea[l]&&e[l]>0){const F=(typeof K==="number"?K:K[l])*dt*Math.pow(A[l],m)/rdist[l];hN[l]=e[l]/(1+F);}else hN[l]=e[l]>0?e[l]:0;continue;}
    const hb=hN[r];if(e[l]<=hb){hN[l]=e[l];continue;}
    const F=(typeof K==="number"?K:K[l])*dt*Math.pow(A[l],m)/rdist[l];hN[l]=(e[l]+F*hb)/(1+F);}
  return hN;
}
function planarGrid(nx,ny,dx){const nC=nx*ny,nbr=new Int32Array(8*nC),dist=new Float32Array(8*nC),area=new Float32Array(nC).fill(dx*dx);
  const D=[[1,0],[-1,0],[0,1],[0,-1],[1,1],[-1,1],[1,-1],[-1,-1]];
  for(let y=0;y<ny;y++)for(let x=0;x<nx;x++){const l=y*nx+x;for(let k=0;k<8;k++){const X=x+D[k][0],Y=y+D[k][1];
    nbr[8*l+k]=X<0||Y<0||X>=nx||Y>=ny?-1:Y*nx+X;dist[8*l+k]=dx*Math.hypot(D[k][0],D[k][1]);}}
  return{nC,nbr,dist,area};}


// ---------- small vector helpers (3x3 rotations row-major; body -> world is Q b) ----------
function mat(A,B){const o=new Array(9);for(let r=0;r<3;r++)for(let c=0;c<3;c++)o[r*3+c]=A[r*3]*B[c]+A[r*3+1]*B[3+c]+A[r*3+2]*B[6+c];return o;}
function matT(A){return[A[0],A[3],A[6],A[1],A[4],A[7],A[2],A[5],A[8]];}
function mv(A,x,y,z,o){o[0]=A[0]*x+A[1]*y+A[2]*z;o[1]=A[3]*x+A[4]*y+A[5]*z;o[2]=A[6]*x+A[7]*y+A[8]*z;return o;}
function mtv(A,x,y,z,o){o[0]=A[0]*x+A[3]*y+A[6]*z;o[1]=A[1]*x+A[4]*y+A[7]*z;o[2]=A[2]*x+A[5]*y+A[8]*z;return o;}
function faceDir(f,a,b,t){let x,y,z;switch(f){case 0:x=1;y=a;z=b;break;case 1:x=-1;y=a;z=b;break;case 2:x=a;y=1;z=b;break;case 3:x=a;y=-1;z=b;break;case 4:x=a;y=b;z=1;break;default:x=a;y=b;z=-1;}
  const L=Math.hypot(x,y,z);t[0]=x/L;t[1]=y/L;t[2]=z/L;return t;}

// ---------- the landmass: picked at the end, its history found by replay ----------
// The object is one landmass at the present: the fixed-grid land cells connected to the picked one, plus a ring of
// submerged continental crust around them (so the coast can move). Its crust is every parcel that ends up there.
function landmassAt(world,cell){
  const N=world.N,nb=world.cube.nb,E=world.elev,K=world.wK;if(E[cell]<=0)return null;
  const inSet=new Uint8Array(N),q=[cell];inSet[cell]=1;
  for(let h=0;h<q.length;h++){const a=q[h];for(let k=0;k<4;k++){const b=nb[4*a+k];if(!inSet[b]&&E[b]>0){inSet[b]=1;q.push(b);}}}
  const land=q.length;for(let h=0;h<land;h++){const a=q[h];for(let k=0;k<4;k++){const b=nb[4*a+k];if(!inSet[b]&&K[b]){inSet[b]=2;q.push(b);}}}
  const votes=new Map();for(let h=0;h<land;h++){const p=world.owner[q[h]];votes.set(p,(votes.get(p)||0)+1);}
  let ref=-1,bv=-1;for(const [p,v] of votes)if(v>bv){bv=v;ref=p;}
  return{cells:Int32Array.from(q),land,ref};
}
// Pass 1: restore the snapshot, tag every parcel, replay to the present step and record where parcels are from
// t0 = present - window on. At the end the object's parcels are read off the picked cells; the reference frame is
// the plate carrying most of the landmass (followed back through its split ancestry); the footprint is where the
// object's crust has been during the window, in that frame. L0 ends exactly where it started (deterministic replay).
class Replay{
  constructor(world,entry,window,pick){
    this.world=world;this.entry=entry;this.pick=pick;this.endSteps=world.steps;this.endT=world.t;this.window=window;
    world.restore(entry.snap);world.trackParcels();this.idSpace=world.plates.length*world.N;this.frames=[];this.done=false;
    if(world.t>=this.endT-window-1e-9)throw new Error("the snapshot is not old enough for this window");}
  record(){const w=this.world,N=w.N,id=new Int32Array(N);
    // parcel identity on the fixed grid, negative-encoded (-2-id) where the parcel is oceanic
    for(let i=0;i<N;i++){const P=w.plates[w.owner[i]],v=P.pid?P.pid[w.bodyIdx[i]]:-1;id[i]=v>=0&&!w.wK[i]?-2-v:v;}
    const Q=new Map();for(const P of w.plates)if(P.alive)Q.set(P.id,Float64Array.from(P.Q));
    this.frames.push({steps:w.steps,t:w.t,id,Q,parent:w.plates.map(P=>P.parent)});}
  advance(maxMs){const w=this.world,t1=Date.now();
    while(!this.done&&Date.now()-t1<maxMs){
      if(!this.frames.length&&w.t>=this.endT-this.window-1e-9)this.record();
      if(w.steps>=this.endSteps){this.done=true;break;}
      w.step();if(this.frames.length)this.record();}
    return this.done;}
  // the plate whose frame is the reference at a frame: the reference plate, or its ancestor before it split off
  static refAt(fr,chain){for(const p of chain)if(fr.Q.has(p))return p;return chain[chain.length-1];}
  finish(){const w=this.world,N=w.N,last=this.frames[this.frames.length-1],obj=new Uint8Array(this.idSpace);let nObj=0;
    for(const c of this.pick.cells){const id=last.id[c];if(id>=0&&!obj[id]){obj[id]=1;nObj++;}}
    // plus the crust it consumes: a continental parcel that touches the object's crust at some frame and is gone by
    // the end (overrun in a collision, or relabelled by L0) was squeezed into it, so it belongs to the sheet
    const atEnd=new Uint8Array(this.idSpace),nb=w.cube.nb;for(let i=0;i<N;i++){const v=last.id[i];if(v>=0)atEnd[v]=1;else if(v<-1)atEnd[-2-v]=1;}
    for(const fr of this.frames){const id=fr.id,add=[];
      for(let i=0;i<N;i++){const v=id[i];if(v<0||obj[v]||atEnd[v])continue;for(let k=0;k<4;k++){const u=id[nb[4*i+k]];if(u>=0&&obj[u]){add.push(v);break;}}}
      for(const v of add)if(!obj[v]){obj[v]=1;nObj++;}}
    const chain=[this.pick.ref];for(let g=0;g<64;g++){const p=last.parent[chain[chain.length-1]];if(p===undefined)break;chain.push(p);}
    // centre: the landmass's centroid at the end, in the reference plate's body frame
    const pos=w.cube.pos,t=[0,0,0];let cx=0,cy=0,cz=0;
    for(let h=0;h<this.pick.land;h++){const c=this.pick.cells[h];cx+=pos[3*c];cy+=pos[3*c+1];cz+=pos[3*c+2];}
    const Qe=last.Q.get(Replay.refAt(last,chain));mtv(Qe,cx,cy,cz,t);const cl=Math.hypot(t[0],t[1],t[2]);const center=[t[0]/cl,t[1]/cl,t[2]/cl];
    // map axes: east and north at the centre, north taken from the world's pole at the end
    const np=mtv(Qe,0,0,1,[0,0,0]);let d=np[0]*center[0]+np[1]*center[1]+np[2]*center[2],nr=[np[0]-d*center[0],np[1]-d*center[1],np[2]-d*center[2]];
    let L=Math.hypot(nr[0],nr[1],nr[2])||1;nr=nr.map(v=>v/L);const ea=[nr[1]*center[2]-nr[2]*center[1],nr[2]*center[0]-nr[0]*center[2],nr[0]*center[1]-nr[1]*center[0]];
    const proj=new Projection(center,ea,nr);
    // footprint of the object's crust over the window
    let x0=Infinity,x1=-Infinity,y0=Infinity,y1=-Infinity;const xy=[0,0];
    for(const fr of this.frames){const Q=fr.Q.get(Replay.refAt(fr,chain));
      for(let i=0;i<N;i++){const id=fr.id[i];if(id<0||!obj[id])continue;mtv(Q,pos[3*i],pos[3*i+1],pos[3*i+2],t);proj.fwd(t,xy);
        if(xy[0]<x0)x0=xy[0];if(xy[0]>x1)x1=xy[0];if(xy[1]<y0)y0=xy[1];if(xy[1]>y1)y1=xy[1];}}
    const t0Steps=this.frames[0].steps;
    // L0 is now back at the present (the replay reached it); the sheet's pass restores and replays again
    return{entry:this.entry,endSteps:this.endSteps,t0Steps,window:this.endT-this.frames[0].t,obj,nObj,chain,proj,bbox:[x0,x1,y0,y1],pick:this.pick,idSpace:this.idSpace};}
}
// azimuthal equidistant projection about `center` (unit vector) with axes ea (east), nr (north); x, y in km
class Projection{
  constructor(c,ea,nr){this.c=c;this.ea=ea;this.nr=nr;}
  fwd(b,out){const c=this.c,ct=Math.max(-1,Math.min(1,b[0]*c[0]+b[1]*c[1]+b[2]*c[2])),th=Math.acos(ct),st=Math.sin(th);
    if(st<1e-12){out[0]=0;out[1]=0;return out;}const k=R_KM*th/st,dx=b[0]-ct*c[0],dy=b[1]-ct*c[1],dz=b[2]-ct*c[2];
    out[0]=k*(dx*this.ea[0]+dy*this.ea[1]+dz*this.ea[2]);out[1]=k*(dx*this.nr[0]+dy*this.nr[1]+dz*this.nr[2]);return out;}
  inv(x,y,out){const r=Math.hypot(x,y),th=r/R_KM,c=this.c;if(r<1e-9){out[0]=c[0];out[1]=c[1];out[2]=c[2];return out;}
    const ct=Math.cos(th),st=Math.sin(th),ux=x/r,uy=y/r;
    for(let k=0;k<3;k++)out[k]=ct*c[k]+st*(ux*this.ea[k]+uy*this.nr[k]);return out;}
}

// ---------- the thin viscous sheet ----------
// L1 models the landmass's lithosphere as a thin viscous sheet (England & McKenzie 1982; Houseman & England 1986):
// a vertically averaged layer of power-law fluid (stress exponent n), on a planar grid co-moving with the reference
// plate. Depth-integrated stress sigma_ij = 2 eta (e_ij + delta_ij e_kk) (the vertical strain rate takes up the
// horizontal divergence), in balance with buoyancy from crustal thickness and with basal drag toward the velocity of
// the plate each piece of crust rides:
//     d_j sigma_ij - kappa (u_i - U_i) = d_i Phi,   Phi = (S^2 - S_ref^2)/2
// Only thickening since t0 drives spreading (S_ref is the inherited thickness; inherited contrasts are equilibrium,
// as everywhere in L1). The edges of the crust are traction-free. Units: km, Myr; Phi in km^2, so eta is in km^2 Myr
// and kappa in Myr. Two knobs: the Argand number Ar (buoyancy against viscous stress at strain rate U0/ell, i.e.
// plateau vs narrow belt) and the coupling length ell = sqrt(4 eta / kappa) over which a plate's motion is
// transmitted into its neighbour (the width of the crumple zone).
// Crust is carried by markers (Lagrangian points, each a piece of a parcel with L0's rock column), each with its own
// thickness and area: strain multiplies thickness by exp(-div dt) and area by exp(div dt), so volume is exact. The grid
// averages marker thickness (an edge cell reads the crust's true thickness, not a coverage fraction).
// Discretization: staggered (MAC) grid, u on vertical faces, v on horizontal faces, thickness and viscosity at cell
// centres, shear at corners; the operator is the gradient of a discrete dissipation, so it is symmetric positive
// definite, and it is solved by conjugate gradients with a diagonal preconditioner and Picard iteration on viscosity.
class Sheet{
  constructor(world,plan,o){
    o=Object.assign({h:50,Ar:3,n:3,ellKm:500,U0:50,S0:35,ero5:0.5,eroExp:1.4,flexKm:40,strengthKm:3,cratonStrength:10,markers:6,picard:2,marginKm:400,cgTol:1e-4,cgMax:200},o);this.o=o;
    this.world=world;this.plan=plan;this.proj=plan.proj;this.t=0;this.steps=0;this.t0=world.t;
    const h=o.h,[bx0,bx1,by0,by1]=plan.bbox,mg=o.marginKm;
    this.x0=Math.floor((bx0-mg)/h)*h;this.y0=Math.floor((by0-mg)/h)*h;this.nx=Math.ceil((bx1+mg-this.x0)/h);this.ny=Math.ceil((by1+mg-this.y0)/h);
    const nx=this.nx,ny=this.ny,nc=nx*ny;
    // viscosity scale from the Argand number at the reference strain rate U0/ell; drag from the coupling length
    const e0=o.U0/o.ellKm;this.e0=e0;this.eta0=0.5*o.S0*o.S0/(2*o.Ar*e0);this.kappa=4*this.eta0/(o.ellKm*o.ellKm);
    // fields
    this.S=new Float64Array(nc);this.Sref=new Float64Array(nc);this.eref=new Float64Array(nc);this.str=new Float64Array(nc);this.Ux=new Float64Array(nc);this.Uy=new Float64Array(nc);
    this.wsum=new Float64Array(nc);this.act=new Uint8Array(nc);this.eta=new Float64Array(nc).fill(this.eta0);this.e=new Float64Array(nc).fill(NaN);this.eps=new Float64Array(nc);
    this.u=new Float64Array((nx+1)*ny);this.v=new Float64Array(nx*(ny+1));this.dS=new Float64Array(nc);this.eroded=new Float64Array(nc);this.rockTop=new Uint8Array(nc);
    this.budget={initial:0,arcsHotspots:0,eroded:0,left:0};this.history=[];
    this.seed();this.budget.initial=this.volume();this.deposit();this.solve();this.elevation();this.stats(0);
  }
  // markers: m x m per parcel of the object present at t0, at sub-positions of its body cell
  seed(){const w=this.world,N=w.N,n=w.n,nn=n*n,m=this.o.markers,q=Math.PI/4,t=[0,0,0],b=[0,0,0],xy=[0,0],obj=this.plan.obj,cube=w.cube,pool=w.pool,U=SIM.U;
    this.Qref=this.refQ();const Qr=this.Qref;
    // L0 elevation smoothed (three passes) on the fixed grid, sampled bilinearly at each marker: the inherited surface
    let Es=Float64Array.from(w.elev),tmp=new Float64Array(N);const nb=cube.nb;
    for(let p=0;p<3;p++){for(let i=0;i<N;i++){let s=Es[i];for(let k=0;k<4;k++)s+=Es[nb[4*i+k]];tmp[i]=s/5;}[Es,tmp]=[tmp,Es];}
    const out=[0,0,0],X=[],Y=[],V=[],S0=[],E0=[],PID=[],PL=[],COL=[],STR=[];
    // L0 relabels some parcels during the window (resampling copies a neighbour's identity), leaving pinholes in the
    // traced set: a continental parcel with at least three traced neighbours on its plate counts as traced too
    const take=new Map();
    for(const P of w.plates){if(!P.alive||!P.pid)continue;const pid=P.pid,tk=new Uint8Array(N);for(let k=0;k<N;k++)if(pid[k]>=0&&obj[pid[k]]&&P.T[k]>0&&P.K[k])tk[k]=1;
      for(let pass=0;pass<2;pass++){const add=[];for(let k=0;k<N;k++){if(tk[k]||!(P.T[k]>0)||!P.K[k])continue;let s=0;for(let d=0;d<4;d++)s+=tk[nb[4*k+d]];if(s>=3)add.push(k);}for(const k of add)tk[k]=1;}
      take.set(P.id,tk);}
    for(const P of w.plates){if(!P.alive||!P.pid)continue;const tk=take.get(P.id);
      for(let k=0;k<N;k++){const id=P.pid[k];if(!tk[k])continue;
        const f=(k/nn)|0,j=((k-f*nn)/n)|0,i=k-f*nn-j*n,area=cube.area[k]*R_KM*R_KM/(m*m);
        // the parcel's rock column, bottom to top, as fractions of its thickness
        const c=P.col[k],cnt=c>=0?pool.cnt[c]:0,tot=c>=0?pool.sum(c):0,col=[];
        for(let u=0;u<cnt;u++)col.push([pool.R[c*U+u],pool.B[c*U+u],tot>0?pool.T[c*U+u]/tot:1]);
        let craton=0;for(const [r,born,fr] of col)if(r===SIM.ROCK.CRATON&&w.t-born>500)craton+=fr;
        for(let a=0;a<m;a++)for(let bb=0;bb<m;bb++){
          faceDir(f,Math.tan(((i+(a+.5)/m)/n*2-1)*q),Math.tan(((j+(bb+.5)/m)/n*2-1)*q),b);
          mv(P.Q,b[0],b[1],b[2],t);cube.lookupFrac(t[0],t[1],t[2],out);
          const cc=out[0],tx=out[1],ty=out[2],cB=nb[4*cc+(tx>=0?0:1)],cC=nb[4*cc+(ty>=0?2:3)],cD=nb[4*cB+(ty>=0?2:3)],ax=Math.abs(tx),ay=Math.abs(ty);
          const e0=(1-ax)*(1-ay)*Es[cc]+ax*(1-ay)*Es[cB]+(1-ax)*ay*Es[cC]+ax*ay*Es[cD],s0=Math.max(10,C.CONT_BASE+e0/C.ISO);
          mtv(Qr,t[0],t[1],t[2],b);this.proj.fwd(b,xy);
          X.push(xy[0]);Y.push(xy[1]);V.push(area);S0.push(s0);E0.push(e0);PID.push(id);PL.push(P.id);COL.push(col);STR.push(1+(this.o.cratonStrength-1)*craton);}}}
    const M=this.nm=X.length;this.mx=Float64Array.from(X);this.my=Float64Array.from(Y);this.mA=Float64Array.from(V);this.mS=Float64Array.from(S0);this.mS0=Float64Array.from(S0);this.mE0=Float64Array.from(E0);
    this.mPid=Int32Array.from(PID);this.mPlate=Int32Array.from(PL);this.mCol=COL;this.mStr=Float64Array.from(STR);this.mIn=new Uint8Array(M).fill(1);this.mUx=new Float64Array(M);this.mUy=new Float64Array(M);
  }
  // the reference plate's orientation now (the reference plate, or its ancestor before it split off)
  refQ(){const w=this.world;for(const p of this.plan.chain){const P=w.plates[p];if(P&&P.alive)return Array.from(P.Q);}return Array.from(w.plates[this.plan.chain[0]].Q);}
  refW(){const w=this.world;for(const p of this.plan.chain){const P=w.plates[p];if(P&&P.alive)return P.W;}return w.plates[this.plan.chain[0]].W;}
  volume(){let v=0;for(let k=0;k<this.nm;k++)if(this.mIn[k])v+=this.mS[k]*this.mA[k];return v;}
  // which plate carries each marker's parcel now (found by identity; a parcel L0 no longer has keeps its last plate)
  locate(){const w=this.world,N=w.N,plates=w.plates,M=this.nm,where=this._where||(this._where=new Int32Array(this.plan.idSpace));where.fill(-1);
    for(const P of plates){if(!P.alive||!P.pid)continue;const pid=P.pid,T=P.T;for(let k=0;k<N;k++){const id=pid[k];if(id>=0&&id<where.length&&T[k]>0&&where[id]<0)where[id]=P.id;}}
    for(let k=0;k<M;k++){const p=where[this.mPid[k]];if(p>=0)this.mPlate[k]=p;else{let P=plates[this.mPlate[k]];for(let g=0;g<16&&!P.alive&&P.mergedInto!==undefined;g++)P=plates[P.mergedInto];this.mPlate[k]=P.id;}}}
  // each marker's plate velocity, relative to the reference plate, in the map plane (km/Myr)
  plateVelocities(){const w=this.world,Wr=this.refW(),Qr=this.Qref,proj=this.proj,b=[0,0,0],b2=[0,0,0],p0=[0,0],p1=[0,0],om=new Map(),dl=1e-3;
    for(let k=0;k<this.nm;k++){if(!this.mIn[k])continue;const p=this.mPlate[k];let o=om.get(p);
      if(!o){const P=w.plates[p],W=P.W;o=mtv(Qr,W[0]-Wr[0],W[1]-Wr[1],W[2]-Wr[2],[0,0,0]);om.set(p,o);}
      proj.inv(this.mx[k],this.my[k],b);const vx=o[1]*b[2]-o[2]*b[1],vy=o[2]*b[0]-o[0]*b[2],vz=o[0]*b[1]-o[1]*b[0];
      b2[0]=b[0]+vx*dl;b2[1]=b[1]+vy*dl;b2[2]=b[2]+vz*dl;const L=Math.hypot(b2[0],b2[1],b2[2]);b2[0]/=L;b2[1]/=L;b2[2]/=L;
      proj.fwd(b,p0);proj.fwd(b2,p1);this.mUx[k]=(p1[0]-p0[0])/dl;this.mUy[k]=(p1[1]-p0[1])/dl;}}
  // markers -> grid (bilinear weights at cell centres): thickness, reference thickness, inherited elevation, strength
  // and plate velocity as weighted averages; a cell is crust where the markers' weight covers a quarter of it
  deposit(){const {nx,ny,o,S,Sref,eref,str,Ux,Uy,wsum,act}=this,h=o.h,nc=nx*ny;S.fill(0);Sref.fill(0);eref.fill(0);str.fill(0);Ux.fill(0);Uy.fill(0);wsum.fill(0);
    const add=(k,c,w)=>{S[c]+=w*this.mS[k];Sref[c]+=w*this.mS0[k];eref[c]+=w*this.mE0[k];str[c]+=w*this.mStr[k];Ux[c]+=w*this.mUx[k];Uy[c]+=w*this.mUy[k];wsum[c]+=w;};
    for(let k=0;k<this.nm;k++){if(!this.mIn[k])continue;const gx=(this.mx[k]-this.x0)/h-0.5,gy=(this.my[k]-this.y0)/h-0.5,i=Math.floor(gx),j=Math.floor(gy),fx=gx-i,fy=gy-j;
      if(i>=0&&j>=0)add(k,j*nx+i,(1-fx)*(1-fy));if(i+1<nx&&j>=0)add(k,j*nx+i+1,fx*(1-fy));if(i>=0&&j+1<ny)add(k,(j+1)*nx+i,(1-fx)*fy);if(i+1<nx&&j+1<ny)add(k,(j+1)*nx+i+1,fx*fy);}
    for(let c=0;c<nc;c++){const ws=wsum[c];if(ws>0){S[c]/=ws;Sref[c]/=ws;eref[c]/=ws;str[c]/=ws;Ux[c]/=ws;Uy[c]/=ws;}act[c]=S[c]>2&&ws>0.25?1:0;}}
  // A u for the sheet operator (gradient of the discrete dissipation), on the staggered grid; inactive faces are 0
  applyK(u,v,ou,ov){const {nx,ny,o,act,eta}=this,h=o.h,kap=this.kappa;ou.fill(0);ov.fill(0);
    for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const c=j*nx+i;if(!act[c])continue;const et=eta[c];
      const exx=(u[j*(nx+1)+i+1]-u[j*(nx+1)+i])/h,eyy=(v[(j+1)*nx+i]-v[j*nx+i])/h,sxx=2*et*(2*exx+eyy),syy=2*et*(exx+2*eyy);
      ou[j*(nx+1)+i+1]+=h*sxx;ou[j*(nx+1)+i]-=h*sxx;ov[(j+1)*nx+i]+=h*syy;ov[j*nx+i]-=h*syy;}
    for(let j=1;j<ny;j++)for(let i=1;i<nx;i++){const a=(j-1)*nx+i-1,b=a+1,c=j*nx+i-1,d=c+1;if(!(act[a]&&act[b]&&act[c]&&act[d]))continue;
      const et=0.25*(eta[a]+eta[b]+eta[c]+eta[d]),exy=0.5*((u[j*(nx+1)+i]-u[(j-1)*(nx+1)+i])+(v[j*nx+i]-v[j*nx+i-1]))/h,g=4*h*h*et*exy/(2*h);
      ou[j*(nx+1)+i]+=g;ou[(j-1)*(nx+1)+i]-=g;ov[j*nx+i]+=g;ov[j*nx+i-1]-=g;}
    const fu=this.fu,fv=this.fv,hh=h*h*kap;
    for(let q=0;q<u.length;q++)if(fu[q])ou[q]+=hh*u[q];else ou[q]=0;
    for(let q=0;q<v.length;q++)if(fv[q])ov[q]+=hh*v[q];else ov[q]=0;}
  // velocities: Picard iterations on the power-law viscosity, each a preconditioned conjugate-gradient solve
  solve(){const {nx,ny,o,act,S,Sref,Ux,Uy,str,eta}=this,h=o.h,nc=nx*ny,nu=(nx+1)*ny,nv=nx*(ny+1),kap=this.kappa;
    const fu=this.fu=this.fu||new Uint8Array(nu),fv=this.fv=this.fv||new Uint8Array(nv);fu.fill(0);fv.fill(0);
    for(let j=0;j<ny;j++)for(let i=0;i<nx;i++)if(act[j*nx+i]){fu[j*(nx+1)+i]=1;fu[j*(nx+1)+i+1]=1;fv[j*nx+i]=1;fv[(j+1)*nx+i]=1;}
    // right-hand side: buoyancy (excess potential energy of thickened crust) and drag toward plate velocities
    const bu=new Float64Array(nu),bv=new Float64Array(nv);
    for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const c=j*nx+i;if(!act[c])continue;const ph=0.5*(S[c]*S[c]-Sref[c]*Sref[c]);
      bu[j*(nx+1)+i+1]+=h*ph;bu[j*(nx+1)+i]-=h*ph;bv[(j+1)*nx+i]+=h*ph;bv[j*nx+i]-=h*ph;}
    const avg=(a,b,fa,fb)=>fa&&fb?(a+b)/2:fa?a:b;
    for(let j=0;j<ny;j++)for(let i=0;i<=nx;i++){const q=j*(nx+1)+i;if(!fu[q])continue;const cl=i>0?j*nx+i-1:-1,cr=i<nx?j*nx+i:-1;bu[q]+=h*h*kap*avg(cl>=0?Ux[cl]:0,cr>=0?Ux[cr]:0,cl>=0&&act[cl],cr>=0&&act[cr]);}
    for(let j=0;j<=ny;j++)for(let i=0;i<nx;i++){const q=j*nx+i;if(!fv[q])continue;const cd=j>0?(j-1)*nx+i:-1,cu=j<ny?j*nx+i:-1;bv[q]+=h*h*kap*avg(cd>=0?Uy[cd]:0,cu>=0?Uy[cu]:0,cd>=0&&act[cd],cu>=0&&act[cu]);}
    let iters=0;
    for(let pic=0;pic<o.picard;pic++){
      // viscosity from the current strain rate (second invariant including the vertical), power law, strength
      const {u,v}=this,cxy=this._cxy||(this._cxy=new Float64Array((nx+1)*(ny+1))),cn=this._cn||(this._cn=new Uint8Array((nx+1)*(ny+1)));cxy.fill(0);cn.fill(0);
      for(let j=1;j<ny;j++)for(let i=1;i<nx;i++){const a=(j-1)*nx+i-1,b=a+1,c=j*nx+i-1,d=c+1;if(!(act[a]&&act[b]&&act[c]&&act[d]))continue;
        cxy[j*(nx+1)+i]=0.5*((u[j*(nx+1)+i]-u[(j-1)*(nx+1)+i])+(v[j*nx+i]-v[j*nx+i-1]))/h;cn[j*(nx+1)+i]=1;}
      for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const c=j*nx+i;if(!act[c])continue;
        const exx=(u[j*(nx+1)+i+1]-u[j*(nx+1)+i])/h,eyy=(v[(j+1)*nx+i]-v[j*nx+i])/h;let sxy=0,k=0;
        for(const q of [j*(nx+1)+i,j*(nx+1)+i+1,(j+1)*(nx+1)+i,(j+1)*(nx+1)+i+1])if(cn[q]){sxy+=cxy[q];k++;}
        const exy=k?sxy/k:0,eE=Math.sqrt(exx*exx+eyy*eyy+exx*eyy+exy*exy);this.eps[c]=eE;
        const r=Math.max(eE,0.02*this.e0)/this.e0;eta[c]=this.eta0*str[c]*Math.pow(r,1/o.n-1);}
      // preconditioned conjugate gradients
      const diagU=new Float64Array(nu),diagV=new Float64Array(nv);
      for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const c=j*nx+i;if(!act[c])continue;const d=4*eta[c];diagU[j*(nx+1)+i]+=d;diagU[j*(nx+1)+i+1]+=d;diagV[j*nx+i]+=d;diagV[(j+1)*nx+i]+=d;}
      for(let j=1;j<ny;j++)for(let i=1;i<nx;i++){const a=(j-1)*nx+i-1,b=a+1,c=j*nx+i-1,d=c+1;if(!(act[a]&&act[b]&&act[c]&&act[d]))continue;const et=0.25*(eta[a]+eta[b]+eta[c]+eta[d]);
        diagU[j*(nx+1)+i]+=et;diagU[(j-1)*(nx+1)+i]+=et;diagV[j*nx+i]+=et;diagV[j*nx+i-1]+=et;}
      for(let q=0;q<nu;q++)diagU[q]=fu[q]?1/(diagU[q]+h*h*kap):0;for(let q=0;q<nv;q++)diagV[q]=fv[q]?1/(diagV[q]+h*h*kap):0;
      const ru=new Float64Array(nu),rv=new Float64Array(nv),zu=new Float64Array(nu),zv=new Float64Array(nv),pu=new Float64Array(nu),pv=new Float64Array(nv),qu=new Float64Array(nu),qv=new Float64Array(nv);
      for(let q=0;q<nu;q++)if(!fu[q])u[q]=0;for(let q=0;q<nv;q++)if(!fv[q])v[q]=0;
      this.applyK(u,v,qu,qv);let bn=0;
      for(let q=0;q<nu;q++){ru[q]=fu[q]?bu[q]-qu[q]:0;zu[q]=ru[q]*diagU[q];pu[q]=zu[q];bn+=bu[q]*bu[q];}
      for(let q=0;q<nv;q++){rv[q]=fv[q]?bv[q]-qv[q]:0;zv[q]=rv[q]*diagV[q];pv[q]=zv[q];bn+=bv[q]*bv[q];}
      let rz=0;for(let q=0;q<nu;q++)rz+=ru[q]*zu[q];for(let q=0;q<nv;q++)rz+=rv[q]*zv[q];bn=Math.sqrt(bn)||1;
      for(let it=0;it<this.o.cgMax;it++){this.applyK(pu,pv,qu,qv);let pq=0;for(let q=0;q<nu;q++)pq+=pu[q]*qu[q];for(let q=0;q<nv;q++)pq+=pv[q]*qv[q];if(!(pq>0))break;
        const al=rz/pq;let rn=0;for(let q=0;q<nu;q++){u[q]+=al*pu[q];ru[q]-=al*qu[q];rn+=ru[q]*ru[q];}for(let q=0;q<nv;q++){v[q]+=al*pv[q];rv[q]-=al*qv[q];rn+=rv[q]*rv[q];}
        iters++;this.cgRes=Math.sqrt(rn)/bn;if(this.cgRes<this.o.cgTol)break;
        let rz2=0;for(let q=0;q<nu;q++){zu[q]=ru[q]*diagU[q];rz2+=ru[q]*zu[q];}for(let q=0;q<nv;q++){zv[q]=rv[q]*diagV[q];rz2+=rv[q]*zv[q];}
        const be=rz2/rz;rz=rz2;for(let q=0;q<nu;q++)pu[q]=zu[q]+be*pu[q];for(let q=0;q<nv;q++)pv[q]=zv[q]+be*pv[q];}}
    this.cgIters=iters;}
  // velocity at a point (bilinear on the staggered components)
  velAt(x,y,out){const {nx,ny,o,u,v}=this,h=o.h;
    let gx=(x-this.x0)/h,gy=(y-this.y0)/h-0.5,i=Math.max(0,Math.min(nx-1,Math.floor(gx))),j=Math.max(0,Math.min(ny-2,Math.floor(gy))),fx=Math.min(1,Math.max(0,gx-i)),fy=Math.min(1,Math.max(0,gy-j));
    out[0]=(1-fx)*(1-fy)*u[j*(nx+1)+i]+fx*(1-fy)*u[j*(nx+1)+i+1]+(1-fx)*fy*u[(j+1)*(nx+1)+i]+fx*fy*u[(j+1)*(nx+1)+i+1];
    gx=(x-this.x0)/h-0.5;gy=(y-this.y0)/h;i=Math.max(0,Math.min(nx-2,Math.floor(gx)));j=Math.max(0,Math.min(ny-1,Math.floor(gy)));fx=Math.min(1,Math.max(0,gx-i));fy=Math.min(1,Math.max(0,gy-j));
    out[1]=(1-fx)*(1-fy)*v[j*nx+i]+fx*(1-fy)*v[j*nx+i+1]+(1-fx)*fy*v[(j+1)*nx+i]+fx*fy*v[(j+1)*nx+i+1];return out;}
  // move markers with the sheet (midpoint rule) and strain them: the horizontal divergence there thins (or, negative,
  // thickens) the crust and grows (shrinks) the marker's area; markers leaving the grid leave the budget
  advect(dt){const vel=[0,0],v2=[0,0],{nx,ny,o,u,v}=this,h=o.h,X1=this.x0+nx*h,Y1=this.y0+ny*h,div=this._div||(this._div=new Float64Array(nx*ny));
    for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const c=j*nx+i;div[c]=this.act[c]?(u[j*(nx+1)+i+1]-u[j*(nx+1)+i]+v[(j+1)*nx+i]-v[j*nx+i])/h:0;}
    for(let k=0;k<this.nm;k++){if(!this.mIn[k])continue;const x=this.mx[k],y=this.my[k];this.velAt(x,y,vel);this.velAt(x+0.5*dt*vel[0],y+0.5*dt*vel[1],v2);
      const nxk=x+dt*v2[0],nyk=y+dt*v2[1];if(nxk<this.x0||nyk<this.y0||nxk>=X1||nyk>=Y1){this.mIn[k]=0;this.budget.left+=this.mS[k]*this.mA[k];continue;}
      const d=this.sample(div,x+0.5*dt*vel[0],y+0.5*dt*vel[1]),g=Math.exp(d*dt);this.mS[k]/=g;this.mA[k]*=g;this.mx[k]=nxk;this.my[k]=nyk;}}
  // bilinear sample of a cell-centred field (inactive cells count as zero)
  sample(f,x,y){const {nx,ny,o}=this,gx=(x-this.x0)/o.h-0.5,gy=(y-this.y0)/o.h-0.5,i=Math.max(0,Math.min(nx-2,Math.floor(gx))),j=Math.max(0,Math.min(ny-2,Math.floor(gy))),fx=Math.min(1,Math.max(0,gx-i)),fy=Math.min(1,Math.max(0,gy-j));
    return(1-fx)*(1-fy)*f[j*nx+i]+fx*(1-fy)*f[j*nx+i+1]+(1-fx)*fy*f[(j+1)*nx+i]+fx*fy*f[(j+1)*nx+i+1];}
  // Elevation: the inherited elevation plus the thickness change, minus a flexural root (Gaussian of width
  // 1.18 x flexural length over the thickness change, the plate yielding beyond strengthKm from local balance)
  elevation(){const {nx,ny,o,S,Sref,eref,act,e,dS}=this,nc=nx*ny,sig=1.18*o.flexKm/o.h;
    for(let c=0;c<nc;c++)dS[c]=act[c]?S[c]-Sref[c]:0;
    const F=this.blur(dS,act,sig),lim=o.strengthKm;
    for(let c=0;c<nc;c++){if(!act[c]){e[c]=NaN;continue;}let rt=ROOT*F[c];const air=ROOT*dS[c];if(rt<air-lim)rt=air-lim;else if(rt>air+lim)rt=air+lim;e[c]=eref[c]+dS[c]-rt;}}
  // separable Gaussian blur of f over cells with mask set (normalized by the mask's own blur)
  blur(f,mask,sig){const {nx,ny}=this,nc=nx*ny,R=Math.max(1,Math.ceil(3*sig)),w=[];for(let k=-R;k<=R;k++)w.push(Math.exp(-k*k/(2*sig*sig)));
    const a=new Float64Array(nc),b=new Float64Array(nc),ma=new Float64Array(nc),mb=new Float64Array(nc);
    for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){let s=0,sm=0;for(let k=-R;k<=R;k++){const ii=i+k;if(ii<0||ii>=nx)continue;const c=j*nx+ii;if(!mask[c])continue;s+=w[k+R]*f[c];sm+=w[k+R];}a[j*nx+i]=s;ma[j*nx+i]=sm;}
    for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){let s=0,sm=0;for(let k=-R;k<=R;k++){const jj=j+k;if(jj<0||jj>=ny)continue;const c=jj*nx+i;s+=w[k+R]*a[c];sm+=w[k+R]*ma[c];}b[j*nx+i]=sm>0?s/sm:0;}
    return b;}
  // L0's arcs and hotspots as new crust on the markers beneath them: each marker thickens by what its L0 cell gained
  addVolcanism(rec){const w=this.world,Qr=this.Qref,b=[0,0,0],t=[0,0,0];let done=0;
    for(let k=0;k<this.nm;k++){if(!this.mIn[k])continue;this.proj.inv(this.mx[k],this.my[k],b);mv(Qr,b[0],b[1],b[2],t);const a=rec[w.cube.lookup(t[0],t[1],t[2])];if(!(a>0))continue;
      const col=this.mCol[k],top=col[col.length-1],fr=a/(this.mS[k]+a);for(const u of col)u[2]*=1-fr;
      // a new top unit of arc or hotspot rock (arcs are by far the commoner source on continents)
      if(top&&top[0]===SIM.ROCK.ARC)top[2]+=fr;else col.push([SIM.ROCK.ARC,w.t,fr]);this.mS[k]+=a;done+=a*this.mA[k];}
    this.budget.arcsHotspots+=done;}
  // share a per-cell volume change among the markers in each cell, in proportion to their volume
  toMarkers(dvCell,apply,line){const {nx,ny,o}=this,h=o.h,tot=new Float64Array(nx*ny),cellOf=this._cellOf||(this._cellOf=new Int32Array(this.nm));
    for(let k=0;k<this.nm;k++){if(!this.mIn[k]){cellOf[k]=-1;continue;}const i=Math.floor((this.mx[k]-this.x0)/h),j=Math.floor((this.my[k]-this.y0)/h);cellOf[k]=j*nx+i;tot[j*nx+i]+=this.mS[k]*this.mA[k];}
    let done=0;
    for(let k=0;k<this.nm;k++){const c=cellOf[k];if(c<0||!(tot[c]>0)||dvCell[c]===0)continue;const dv=dvCell[c]*this.mS[k]*this.mA[k]/tot[c];apply(k,dv);done+=dv;}
    return done;}
  // erosion on the sheet: denudation rising steeply with elevation, ero5 km/Myr at 5 km times (e / 5 km)^eroExp
  // (0.5 km/Myr in high ranges, about 20 m/Myr on lowlands at 0.5 km: the magnitudes of Ahnert's relation between
  // denudation and relief). A smooth stand-in for rivers, which the sheet's grid is too coarse to draw (the drainage
  // pass draws them at the end). Eroded rock leaves the system, stripped from the top of the rock columns.
  erode(dt){const {nx,ny,o,act,e}=this,nc=nx*ny,h=o.h,dv=new Float64Array(nc);
    for(let c=0;c<nc;c++){if(!act[c]||!(e[c]>0))continue;const d=Math.min(e[c],o.ero5*Math.pow(e[c]/5,o.eroExp)*dt);dv[c]=-d*h*h;this.eroded[c]+=d;}
    const rem=-this.toMarkers(dv,(k,d)=>{const S=this.mS[k]*this.mA[k],cut=-d/S,col=this.mCol[k];let left=cut;
      while(left>0&&col.length>1){const top=col[col.length-1];if(top[2]>left){top[2]-=left;left=0;}else{left-=top[2];col.pop();}}
      if(left>0&&col.length)col[0][2]-=left;const sum=col.reduce((a,u)=>a+u[2],0)||1;for(const u of col)u[2]/=sum;this.mS[k]+=d/this.mA[k];});
    this.budget.eroded+=rem;}
  // one coupled step: L0 advances (recording its tectonics), then the sheet runs in substeps, then volcanism, erosion
  advance(){const w=this.world,N=w.N,rec=this._rec||(this._rec=new Float32Array(N));
    rec.fill(0);w.tecRec=rec;w.step();w.tecRec=null;const dt=w.dt;
    this.locate();this.Qref=this.refQ();
    // substeps limited by a Courant number of 0.5 on the fastest marker
    let vmax=1;for(let k=0;k<this.nm;k++)if(this.mIn[k]){const s=Math.hypot(this.mUx[k],this.mUy[k]);if(s>vmax)vmax=s;}
    const ns=Math.max(1,Math.ceil(dt*vmax/(0.5*this.o.h)));
    for(let s=0;s<ns;s++){this.plateVelocities();this.deposit();this.solve();this.advect(dt/ns);}
    this.deposit();this.addVolcanism(rec);this.deposit();this.elevation();this.erode(dt);this.deposit();this.elevation();
    this.t+=dt;this.steps++;this.stats(dt,ns);}
  stats(dt,ns){const {nx,ny,o,act,e,S,dS,eps}=this,nc=nx*ny;let land=0,se=0,mx=-1e9,smax=0,dsmax=0,emax=0;
    for(let c=0;c<nc;c++){if(!act[c])continue;if(S[c]>smax)smax=S[c];if(dS[c]>dsmax)dsmax=dS[c];if(eps[c]>emax)emax=eps[c];if(e[c]>0){land++;se+=e[c];if(e[c]>mx)mx=e[c];}}
    this.history.push({t:this.t,landKm2:land*o.h*o.h,meanElevKm:land?se/land:0,maxElevKm:land?mx:0,maxThicknessKm:smax,maxThickeningKm:dsmax,maxStrainRate:emax,substeps:ns||0,cg:this.cgIters||0});}
  balance(){const b=this.budget,now=this.volume(),expect=b.initial+b.arcsHotspots-b.eroded-b.left;return{now,expect,error:(now-expect)/b.initial};}
  // the surface rock of each cell: the commonest top unit among its markers, by volume
  surfaceRock(){const {nx,ny,o}=this,h=o.h,acc=new Float64Array(nx*ny*8);
    for(let k=0;k<this.nm;k++){if(!this.mIn[k])continue;const i=Math.floor((this.mx[k]-this.x0)/h),j=Math.floor((this.my[k]-this.y0)/h),col=this.mCol[k];if(!col.length)continue;acc[(j*nx+i)*8+col[col.length-1][0]]+=this.mS[k]*this.mA[k];}
    for(let c=0;c<nx*ny;c++){let best=0,bv=0;for(let r=1;r<8;r++)if(acc[c*8+r]>bv){bv=acc[c*8+r];best=r;}this.rockTop[c]=best;}return this.rockTop;}
  // map plane -> the L0 elevation there now (the backdrop outside L1's crust)
  l0ElevAt(x,y,tmp){const w=this.world,b=this.proj.inv(x,y,tmp||[0,0,0]),wv=mv(this.Qref,b[0],b[1],b[2],[0,0,0]),out=[0,0,0],nb=w.cube.nb;w.cube.lookupFrac(wv[0],wv[1],wv[2],out);
    const c=out[0],tx=out[1],ty=out[2],cB=nb[4*c+(tx>=0?0:1)],cC=nb[4*c+(ty>=0?2:3)],cD=nb[4*cB+(ty>=0?2:3)],ax=Math.abs(tx),ay=Math.abs(ty),E=w.elev;
    return(1-ax)*(1-ay)*E[c]+ax*(1-ay)*E[cB]+(1-ax)*ay*E[cC]+ax*ay*E[cD];}
}

// ---------- the drainage pass ----------
// With tectonics frozen at the end, the sheet's surface is refined (bilinear, plus relief L0 and the sheet cannot
// resolve, fixed to the planet) and eroded for a few Myr by stream power on a finer grid, only to settle the
// drainage network; each outlet's drainage basin is a natural tile for L2.
function basins(sheet,o){
  o=Object.assign({km:10,myr:5,steps:5,K:0.1,m:0.5,relief:0.15},o);
  const f=Math.max(1,Math.round(sheet.o.h/o.km)),nx=sheet.nx*f,ny=sheet.ny*f,h=sheet.o.h/f,nc=nx*ny,e=new Float64Array(nc),land=new Uint8Array(nc);
  const nz=makeNoise3(hashStr(String(sheet.world.opts.seed)+"|L1basins")|0),b=[0,0,0],H=sheet.o.h,snx=sheet.nx,sny=sheet.ny,se=sheet.e,f0=R_KM/(3*H),oct=4;
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const x=sheet.x0+(i+.5)*h,y=sheet.y0+(j+.5)*h,gx=(x-sheet.x0)/H-0.5,gy=(y-sheet.y0)/H-0.5;
    let i0=Math.floor(gx),j0=Math.floor(gy);const fx=gx-i0,fy=gy-j0;let s=0,w=0;
    for(const [di,dj,ww] of [[0,0,(1-fx)*(1-fy)],[1,0,fx*(1-fy)],[0,1,(1-fx)*fy],[1,1,fx*fy]]){const ii=i0+di,jj=j0+dj;if(ii<0||jj<0||ii>=snx||jj>=sny)continue;const v=se[jj*snx+ii];if(isNaN(v))continue;s+=ww*v;w+=ww;}
    const c=j*nx+i;if(w<0.5){e[c]=-1;continue;}let v=s/w;sheet.proj.inv(x,y,b);const n=nz(b[0],b[1],b[2],f0,oct);v+=(o.relief+o.relief*Math.max(0,v))*2*(n-0.5);e[c]=v;land[c]=v>0?1:0;}
  const G=planarGrid(nx,ny,h),S=scratch(nc),alive=new Uint8Array(nc).fill(1),out=new Uint8Array(nc);
  for(let c=0;c<nc;c++)out[c]=land[c]?0:1;
  const dt=o.myr/o.steps;for(let s=0;s<o.steps;s++){route(G,e,alive,out,S,null);const hN=streamPower(e,S,o.K,o.m,dt);for(let c=0;c<nc;c++)if(!out[c])e[c]=hN[c];}
  route(G,e,alive,out,S,null);
  const basin=new Int32Array(nc).fill(-1),size=new Map();
  for(let q=0;q<S.ns;q++){const l=S.stack[q],r=S.rec[l];basin[l]=r===l?l:basin[r];}
  for(let c=0;c<nc;c++)if(land[c]){const id=basin[c];size.set(id,(size.get(id)||0)+1);}
  const sizes=[...size.values()].sort((a,b)=>b-a);
  return{nx,ny,h,x0:sheet.x0,y0:sheet.y0,e,land,A:Float64Array.from(S.A),rec:Int32Array.from(S.rec),basin,nBasins:size.size,sizes,landCells:land.reduce((a,v)=>a+v,0),S};
}
return {Replay,Sheet,Projection,landmassAt,basins,route,streamPower,planarGrid,scratch,Heap};
})();
// ===== L1-END =====
