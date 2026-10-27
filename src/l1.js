// ===== L1-START =====
// L1: one landmass at the present, refined, with the last few tens of Myr of its tectonic history replayed as a
// thin viscous sheet: where plates converge, its crust shortens and thickens into mountain belts and plateaus at L1
// resolution, carried by markers with L0's rock columns. A short frozen erosion pass then settles its drainage
// basins, the tiles for L2.
const L1 = (() => {
const {R_KM,C,makeLookup,elevOf,makeNoise3,hashStr,getCube}=SIM;
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
// A square grid with 8 neighbours. With jitter > 0, each cell's node (used for distances, hence slopes) is displaced
// within the cell by up to jitter x cell width, by a fixed hash of (x + ix, y + iy) and salt, so the neighbour
// directions differ from cell to cell and rivers are not locked to the lattice's 0 and 45 degree directions.
function jit2(x,y,salt){let h=Math.imul(x,0x27d4eb2d)^Math.imul(y,0x165667b1)^salt;h=Math.imul(h^(h>>>16),0x7feb352d);h=Math.imul(h^(h>>>15),0x846ca68b);h^=h>>>16;return[(h&0xffff)/65536-0.5,(h>>>16)/65536-0.5];}
function planarGrid(nx,ny,dx,jitter=0,ix=0,iy=0,salt=0){const nC=nx*ny,nbr=new Int32Array(8*nC),dist=new Float32Array(8*nC),area=new Float32Array(nC).fill(dx*dx);
  const D=[[1,0],[-1,0],[0,1],[0,-1],[1,1],[-1,1],[1,-1],[-1,-1]],jx=new Float32Array(nC),jy=new Float32Array(nC);
  if(jitter>0)for(let y=0;y<ny;y++)for(let x=0;x<nx;x++){const [u,v]=jit2(x+ix,y+iy,salt);jx[y*nx+x]=2*jitter*u;jy[y*nx+x]=2*jitter*v;}
  for(let y=0;y<ny;y++)for(let x=0;x<nx;x++){const l=y*nx+x;for(let k=0;k<8;k++){const X=x+D[k][0],Y=y+D[k][1],m=X<0||Y<0||X>=nx||Y>=ny?-1:Y*nx+X;
    nbr[8*l+k]=m;dist[8*l+k]=m<0?dx*Math.hypot(D[k][0],D[k][1]):dx*Math.hypot(D[k][0]+jx[m]-jx[l],D[k][1]+jy[m]-jy[l]);}}
  return{nC,nbr,dist,area};}


// ---------- small vector helpers (3x3 rotations row-major; body -> world is Q b) ----------
function mat(A,B){const o=new Array(9);for(let r=0;r<3;r++)for(let c=0;c<3;c++)o[r*3+c]=A[r*3]*B[c]+A[r*3+1]*B[3+c]+A[r*3+2]*B[6+c];return o;}
function matT(A){return[A[0],A[3],A[6],A[1],A[4],A[7],A[2],A[5],A[8]];}
function mv(A,x,y,z,o){o[0]=A[0]*x+A[1]*y+A[2]*z;o[1]=A[3]*x+A[4]*y+A[5]*z;o[2]=A[6]*x+A[7]*y+A[8]*z;return o;}
function mtv(A,x,y,z,o){o[0]=A[0]*x+A[3]*y+A[6]*z;o[1]=A[1]*x+A[4]*y+A[7]*z;o[2]=A[2]*x+A[5]*y+A[8]*z;return o;}
function faceDir(f,a,b,t){let x,y,z;switch(f){case 0:x=1;y=a;z=b;break;case 1:x=-1;y=a;z=b;break;case 2:x=a;y=1;z=b;break;case 3:x=a;y=-1;z=b;break;case 4:x=a;y=b;z=1;break;default:x=a;y=b;z=-1;}
  const L=Math.hypot(x,y,z);t[0]=x/L;t[1]=y/L;t[2]=z/L;return t;}

// ---------- the landmass: picked at the end, its history found by replay ----------
// The object is one landmass at the present. It is grown from the picked cell in four steps:
// 1. Land: the fixed-grid land cells connected to it, corners included (two cells meeting only at a corner are one
//    coast at L0's resolution).
// 2. Buffer: the sea within bufferKm of that land's coast (great-circle distance to the nearest land cell, less a cell). An island with a
//    cell in the buffer joins, and the buffer grows around it, until nothing more joins: islands are chained to the
//    landmass by steps of at most bufferKm (single-linkage clustering of the land at that distance). The buffer is at
//    least the ring of continental crust around the land, so the coast can move.
// 3. Enclosed seas: water the region now surrounds (the pieces of the rest of the sphere other than the largest, the
//    open ocean) joins, with any islands in it. A strait narrower than about twice the buffer is closed by the
//    buffer, so its sea counts as enclosed.
// 4. Crust: of the region's continental crust, only what is connected to the landmass's land through the region's
//    continental crust is kept (a sliver of another continent's shelf that the buffer reaches stays out).
// Cells: the land first (`land` of them), then the rest of the region. The sheet carries the continental crust of all
// of it; the region's oceanic floor stays L0's. Its crust is every parcel that ends up there.
const CORNER=[2,3,1,0];
function landmassAt(world,cell,bufferKm=0){
  const N=world.N,nb=world.cube.nb,E=world.elev,K=world.wK,pos=world.cube.pos;if(E[cell]<=0)return null;
  const nb8=(a,f)=>{for(let k=0;k<4;k++){const b=nb[4*a+k];f(b);f(nb[4*b+CORNER[k]]);}};   // the 4 sides and the 4 corners (NE, SW, NW, SE)
  const inL=new Uint8Array(N),landQ=[],dist=new Float64Array(N).fill(Infinity),src=new Int32Array(N).fill(-1),Q=[];
  // distances are coast to coast: between cell centres, less one cell
  const B=bufferKm>0?(bufferKm+world.cellKm)/R_KM:0,cosB=Math.cos(Math.min(Math.PI,B));
  // add the land component containing c; its cells become buffer sources
  const addLand=c=>{const q=[c];inL[c]=1;for(let h=0;h<q.length;h++){const a=q[h];landQ.push(a);dist[a]=0;src[a]=a;Q.push(a);nb8(a,b=>{if(!inL[b]&&E[b]>0){inL[b]=1;q.push(b);}});}};
  addLand(cell);let islands=0;
  // the buffer, by label-correcting search: each cell keeps its nearest land cell and the great-circle distance to it
  for(let h=0;h<Q.length;h++){const a=Q[h],s=src[a];
    nb8(a,b=>{const c=pos[3*b]*pos[3*s]+pos[3*b+1]*pos[3*s+1]+pos[3*b+2]*pos[3*s+2];if(c<cosB-1e-12)return;
      const d=Math.acos(Math.min(1,c));if(d>=dist[b]-1e-12)return;
      if(E[b]>0&&!inL[b]){islands++;addLand(b);return;}dist[b]=d;src[b]=s;Q.push(b);});}
  const inR=new Uint8Array(N);for(let i=0;i<N;i++)if(dist[i]<Infinity)inR[i]=1;
  for(const a of landQ)for(let k=0;k<4;k++){const b=nb[4*a+k];if(K[b])inR[b]=1;}
  // enclosed seas: the pieces of the rest of the sphere, all but the largest
  const piece=new Int32Array(N).fill(-1),pieces=[];
  for(let i=0;i<N;i++){if(inR[i]||piece[i]>=0)continue;const q=[i];piece[i]=pieces.length;for(let h=0;h<q.length;h++){const a=q[h];for(let k=0;k<4;k++){const b=nb[4*a+k];if(!inR[b]&&piece[b]<0){piece[b]=pieces.length;q.push(b);}}}pieces.push(q);}
  let open=-1;for(let p=0;p<pieces.length;p++)if(open<0||pieces[p].length>pieces[open].length)open=p;
  let seaCells=0;for(let p=0;p<pieces.length;p++)if(p!==open){seaCells+=pieces[p].length;for(const c of pieces[p])inR[c]=1;}
  for(let i=0;i<N;i++)if(inR[i]&&E[i]>0&&!inL[i]){islands++;addLand(i);}
  // the region's continental crust connected to the land through continental crust in the region
  const keep=new Uint8Array(N),q=[];for(const a of landQ){keep[a]=1;q.push(a);}
  for(let h=0;h<q.length;h++){const a=q[h];for(let k=0;k<4;k++){const b=nb[4*a+k];if(inR[b]&&K[b]&&!keep[b]){keep[b]=1;q.push(b);}}}
  const cells=landQ.slice();for(let i=0;i<N;i++)if(inR[i]&&!inL[i]&&(keep[i]||!K[i]))cells.push(i);
  const land=landQ.length,votes=new Map();for(let h=0;h<land;h++){const p=world.owner[cells[h]];votes.set(p,(votes.get(p)||0)+1);}
  let ref=-1,bv=-1;for(const [p,v] of votes)if(v>bv){bv=v;ref=p;}
  return{cells:Int32Array.from(cells),land,ref,bufferKm,islands,seaCells};
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
  finish(){const w=this.world,N=w.N,nb=w.cube.nb,last=this.frames[this.frames.length-1],obj=new Uint8Array(this.idSpace);let nObj=0;
    for(const c of this.pick.cells){const id=last.id[c];if(id>=0&&!obj[id]){obj[id]=1;nObj++;}}
    // plus the landmass's land on oceanic crust (islands built by arcs and hotspots), marked 2: the sheet carries it too
    for(let h=0;h<this.pick.land;h++){const id=last.id[this.pick.cells[h]];if(id<-1&&!obj[-2-id]){obj[-2-id]=2;nObj++;}}
    // and the coastal rim: oceanic cells touching the land (corners included) on the land cell's own plate. L0's coast
    // is the zero of its bilinear surface, half a cell or so out into such a cell, so the sheet's crust has to reach
    // past it, or its edge cuts the coast short in straight lines. Ocean on another plate (across a trench) stays out
    for(let h=0;h<this.pick.land;h++){const a=this.pick.cells[h],pa=w.owner[a];
      for(let k=0;k<4;k++){const b=nb[4*a+k];for(const c of [b,nb[4*b+CORNER[k]]]){const id=last.id[c];if(id<-1&&w.owner[c]===pa&&!obj[-2-id]){obj[-2-id]=2;nObj++;}}}}
    // plus the crust it consumes: a continental parcel that touches the object's crust at some frame and is gone by
    // the end (overrun in a collision, or relabelled by L0) was squeezed into it, so it belongs to the sheet
    const atEnd=new Uint8Array(this.idSpace);for(let i=0;i<N;i++){const v=last.id[i];if(v>=0)atEnd[v]=1;else if(v<-1)atEnd[-2-v]=1;}
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
      for(let i=0;i<N;i++){const v=fr.id[i],id=v<-1?-2-v:v;if(id<0||!obj[id])continue;mtv(Q,pos[3*i],pos[3*i+1],pos[3*i+2],t);proj.fwd(t,xy);
        if(xy[0]<x0)x0=xy[0];if(xy[0]>x1)x1=xy[0];if(xy[1]<y0)y0=xy[1];if(xy[1]>y1)y1=xy[1];}}
    const t0Steps=this.frames[0].steps;
    // L0 is now back at the present (the replay reached it); the sheet's pass restores and replays again
    return{entry:this.entry,endSteps:this.endSteps,t0Steps,window:this.endT-this.frames[0].t,obj,nObj,chain,proj,bbox:[x0,x1,y0,y1],pick:this.pick,idSpace:this.idSpace};}
}
// exp of the 2x2 matrix [[a,b],[c,d]] into E (row-major): e^t (cosh s I + sinh(s)/s N), N = M - t I, s^2 = -det N
function expm2(a,b,c,d,E){const t=(a+d)/2,p=(a-d)/2,s2=p*p+b*c,et=Math.exp(t);let ch,sh;
  if(s2>1e-14){const s=Math.sqrt(s2);ch=Math.cosh(s);sh=Math.sinh(s)/s;}else if(s2<-1e-14){const s=Math.sqrt(-s2);ch=Math.cos(s);sh=Math.sin(s)/s;}else{ch=1+s2/2;sh=1+s2/6;}
  E[0]=et*(ch+sh*p);E[1]=et*sh*b;E[2]=et*sh*c;E[3]=et*(ch-sh*p);return E;}
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
// ---------- dynamic topography from L0's slab record ----------
// Slabs sinking into the mantle pull the surface above them down. L0 records the oceanic area consumed at trenches
// (weighted by slab age) in 10 Myr bins on its fixed, mantle grid. A bin of age a has sunk to depth d(a): 20 km/Myr
// to 660 km, then 10 km/Myr. Its pull spreads over a Gaussian of width 300 km + 0.8 d and weakens as exp(-d/1200 km).
// The deflection is -ampKm times the sum of the spread, weighted slab densities (consumed area per area), less its
// mean over the sphere (return flow lifts the surface away from slabs). This is a calibrated stand-in, not mantle flow: amplitudes of about 1 km over 1,000-3,000 km, as observed over subduction.
// Computed on a cube grid of half L0's resolution, but at least 48 cells a face (about 210 km; the narrowest kernel is
// 300 km); blurs are diffusion steps, applied
// oldest bin first so each bin receives exactly its own total width (Gaussians compose by adding variances).
const DYN={ampKm:0.6,vUpper:20,vLower:10,zTrans:660,sigma0:300,sigmaPerKm:0.8,decayKm:1200};
const dynCache=new Map();
function dynGrid(world){const n=world.n;let g=dynCache.get(n);if(g)return g;
  const n2=Math.max(48,Math.round(n/2)),cube2=getCube(n2),N=world.N,map=new Int32Array(N),pos=world.cube.pos;
  for(let i=0;i<N;i++)map[i]=cube2.lookup(pos[3*i],pos[3*i+1],pos[3*i+2]);
  const area2=Float64Array.from(cube2.area,a=>a*R_KM*R_KM),h2=Math.PI/2/n2*R_KM;
  g={cube2,map,area2,h2};dynCache.set(n,g);return g;}
const dynDepth=a=>{const tU=DYN.zTrans/DYN.vUpper;return a<=tU?DYN.vUpper*a:DYN.zTrans+DYN.vLower*(a-tU);};
// the deflection (km, negative down) on the half-resolution grid at the world's present time
function dynTopo(world){const g=dynGrid(world),{cube2,map,area2,h2}=g,N2=cube2.N,nb=cube2.nb,N=world.N,t=world.t,nB=world.slabT.length;
  if(!world.slab)return{g,F:new Float64Array(N2)};
  const bins=[];for(let q=0;q<nB;q++){const t0=world.slabT[q];if(!(t0>-1e8))continue;const mid=Math.min(t,t0+C.SLAB_DT/2),d=dynDepth(Math.max(0,t-mid));
    bins.push({q,d,sig:DYN.sigma0+DYN.sigmaPerKm*d,k:Math.exp(-d/DYN.decayKm)});}
  bins.sort((a,b)=>b.sig-a.sig);
  let F=new Float64Array(N2),T=new Float64Array(N2);const m=new Float64Array(N2),varStep=h2*h2/2;
  // lazy steps (at most half the neighbours' mean) so the kernel does not checkerboard
  const blur=v=>{let it=v/varStep;while(it>1e-6){const al=Math.min(0.5,it);for(let c=0;c<N2;c++){const s=F[nb[4*c]]+F[nb[4*c+1]]+F[nb[4*c+2]]+F[nb[4*c+3]];T[c]=F[c]+al*(s/4-F[c]);}[F,T]=[T,F];it-=al;}};
  for(let b=0;b<bins.length;b++){const {q,k}=bins[b];m.fill(0);for(let i=0;i<N;i++){const v=world.slab[q*N+i];if(v)m[map[i]]+=v;}
    for(let c=0;c<N2;c++)F[c]+=k*m[c]/area2[c];
    const next=b+1<bins.length?bins[b+1].sig:0;blur(bins[b].sig*bins[b].sig-next*next);}
  // mass conservation: where slabs sink, return flow rises elsewhere, so the deflection averages zero over the sphere
  let sm=0,sa=0;for(let c=0;c<N2;c++){sm+=F[c]*area2[c];sa+=area2[c];}const mean=sm/sa;
  for(let c=0;c<N2;c++)F[c]=-DYN.ampKm*(F[c]-mean);return{g,F};}
// sample a half-resolution field at a world point (bilinear)
function dynAt(D,x,y,z,out){const {cube2}=D.g,nb=cube2.nb;cube2.lookupFrac(x,y,z,out);const c=out[0],tx=out[1],ty=out[2],cB=nb[4*c+(tx>=0?0:1)],cC=nb[4*c+(ty>=0?2:3)],cD=nb[4*cB+(ty>=0?2:3)],ax=Math.abs(tx),ay=Math.abs(ty),F=D.F;
  return(1-ax)*(1-ay)*F[c]+ax*(1-ay)*F[cB]+(1-ax)*ay*F[cC]+ax*ay*F[cD];}
class Sheet{
  constructor(world,plan,o){
    o=Object.assign({h:50,Ar:3,n:3,ellKm:500,U0:50,S0:35,ero5:0.5,eroExp:1.4,flexKm:40,strengthKm:3,cratonStrength:10,markers:6,markerJitter:1,rateMyr:3,dyn:true,picard:2,marginKm:400,cgTol:1e-4,cgMax:200},o);this.o=o;
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
    this.dyn=new Float64Array(nc);this.seed();this.budget.initial=this.volume();this.deposit();this.solve();this.updateDyn();this.elevation();this.stats(0);
    this.ePrev=Float64Array.from(this.e);this.uRate=new Float64Array(nc);
  }
  // markers: m x m per parcel of the object present at t0, one at a random point (a fixed hash) in each of m x m
  // sub-cells of its body cell (stratified, so marker rows don't alias with the grid or bunch into lines under shortening)
  seed(){const w=this.world,N=w.N,n=w.n,nn=n*n,m=this.o.markers,mj=this.o.markerJitter,q=Math.PI/4,t=[0,0,0],b=[0,0,0],xy=[0,0],obj=this.plan.obj,cube=w.cube,pool=w.pool,U=SIM.U;
    this.Qref=this.refQ();const Qr=this.Qref;
    // L0 elevation smoothed (three passes) on the fixed grid, sampled bilinearly at each marker: the inherited surface.
    // Land is averaged only with land and sea only with sea, so every cell keeps its sign and the coast stays where L0's
    // own (bilinear) coast is; averaging across it would sink a cell of coastal lowland into the abyss beside it
    let Es=Float64Array.from(w.elev),tmp=new Float64Array(N);const nb=cube.nb,EL=w.elev;
    for(let p=0;p<3;p++){for(let i=0;i<N;i++){const up=EL[i]>0;let s=Es[i],n=1;for(let k=0;k<4;k++){const j=nb[4*i+k];if(EL[j]>0===up){s+=Es[j];n++;}}tmp[i]=s/n;}[Es,tmp]=[tmp,Es];}
    // craton share of the crust (rock older than 500 Myr), also smoothed on the fixed grid (continental crust with
    // continental crust) and sampled at each marker: strength is 10x on cratons, and taken parcel by parcel its edges
    // would be L0's cell edges, staircases along which the weaker crust between cratons crumples into square ridges
    const crat=c=>{if(c<0)return 0;const cnt=pool.cnt[c],tot=pool.sum(c);let s=0;for(let u=0;u<cnt;u++)if(pool.R[c*U+u]===SIM.ROCK.CRATON&&w.t-pool.B[c*U+u]>500)s+=tot>0?pool.T[c*U+u]/tot:1;return s;};
    let Cs=new Float64Array(N);for(let i=0;i<N;i++)if(w.wK[i])Cs[i]=crat(w.wCol[i]);
    for(let p=0;p<3;p++){for(let i=0;i<N;i++){if(!w.wK[i]){tmp[i]=0;continue;}let s=Cs[i],n=1;for(let k=0;k<4;k++){const j=nb[4*i+k];if(w.wK[j]){s+=Cs[j];n++;}}tmp[i]=s/n;}[Cs,tmp]=[tmp,Cs];}
    const out=[0,0,0],X=[],Y=[],V=[],S0=[],E0=[],PID=[],PL=[],COL=[],STR=[];
    // Every object parcel present now is seeded, continental or not: oceanic parcels marked 2 (land and the coastal
    // rim on oceanic crust), and continental crust that is still oceanic at t0 (accreted during the window).
    // L0 relabels some parcels during the window (resampling copies a neighbour's identity), leaving pinholes in the
    // traced set: a continental parcel with at least three traced neighbours on its plate counts as traced too
    const take=new Map();
    for(const P of w.plates){if(!P.alive||!P.pid)continue;const pid=P.pid,tk=new Uint8Array(N);for(let k=0;k<N;k++)if(pid[k]>=0&&obj[pid[k]]&&P.T[k]>0)tk[k]=1;
      for(let pass=0;pass<2;pass++){const add=[];for(let k=0;k<N;k++){if(tk[k]||!(P.T[k]>0)||!P.K[k])continue;let s=0;for(let d=0;d<4;d++)s+=tk[nb[4*k+d]];if(s>=3)add.push(k);}for(const k of add)tk[k]=1;}
      take.set(P.id,tk);}
    for(const P of w.plates){if(!P.alive||!P.pid)continue;const tk=take.get(P.id);
      for(let k=0;k<N;k++){const id=P.pid[k];if(!tk[k])continue;
        const f=(k/nn)|0,j=((k-f*nn)/n)|0,i=k-f*nn-j*n,area=cube.area[k]*R_KM*R_KM/(m*m);
        // the parcel's rock column, bottom to top, as fractions of its thickness
        const c=P.col[k],cnt=c>=0?pool.cnt[c]:0,tot=c>=0?pool.sum(c):0,col=[];
        for(let u=0;u<cnt;u++)col.push([pool.R[c*U+u],pool.B[c*U+u],tot>0?pool.T[c*U+u]/tot:1]);
        for(let a=0;a<m;a++)for(let bb=0;bb<m;bb++){
          const [ju,jv]=mj>0?jit2(k,a*m+bb,P.id):[0,0];faceDir(f,Math.tan(((i+(a+.5+mj*ju)/m)/n*2-1)*q),Math.tan(((j+(bb+.5+mj*jv)/m)/n*2-1)*q),b);
          mv(P.Q,b[0],b[1],b[2],t);cube.lookupFrac(t[0],t[1],t[2],out);
          const cc=out[0],tx=out[1],ty=out[2],cB=nb[4*cc+(tx>=0?0:1)],cC=nb[4*cc+(ty>=0?2:3)],cD=nb[4*cB+(ty>=0?2:3)],ax=Math.abs(tx),ay=Math.abs(ty);
          const e0=(1-ax)*(1-ay)*Es[cc]+ax*(1-ay)*Es[cB]+(1-ax)*ay*Es[cC]+ax*ay*Es[cD],s0=Math.max(10,C.CONT_BASE+e0/C.ISO);
          const craton=(1-ax)*(1-ay)*Cs[cc]+ax*(1-ay)*Cs[cB]+(1-ax)*ay*Cs[cC]+ax*ay*Cs[cD];
          mtv(Qr,t[0],t[1],t[2],b);this.proj.fwd(b,xy);
          X.push(xy[0]);Y.push(xy[1]);V.push(area);S0.push(s0);E0.push(e0);PID.push(id);PL.push(P.id);COL.push(col);STR.push(1+(this.o.cratonStrength-1)*craton);}}}
    const M=this.nm=X.length;this.mx=Float64Array.from(X);this.my=Float64Array.from(Y);this.mA=Float64Array.from(V);this.mS=Float64Array.from(S0);this.mS0=Float64Array.from(S0);this.mE0=Float64Array.from(E0);
    this.mPid=Int32Array.from(PID);this.mPlate=Int32Array.from(PL);this.mCol=COL;this.mStr=Float64Array.from(STR);this.mIn=new Uint8Array(M).fill(1);this.mUx=new Float64Array(M);this.mUy=new Float64Array(M);
    // finite strain: each marker's 2x2 deformation gradient F (row-major, map plane), identity at t0; det F = area / area at t0
    this.mF=new Float64Array(4*M);for(let k=0;k<M;k++){this.mF[4*k]=1;this.mF[4*k+3]=1;}this.mA0=Float64Array.from(V);
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
  // The velocity gradient: du/dx and dv/dy at cell centres (their sum is the divergence), du/dy and dv/dx at cell
  // corners where all four cells around are crust (zero elsewhere, as at the traction-free edges).
  velGrad(){const {nx,ny,o,u,v,act}=this,h=o.h,nc=nx*ny,nn=(nx+1)*(ny+1);
    const G=this._vg||(this._vg={ux:new Float64Array(nc),vy:new Float64Array(nc),uy:new Float64Array(nn),vx:new Float64Array(nn)});G.uy.fill(0);G.vx.fill(0);
    for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const c=j*nx+i;G.ux[c]=act[c]?(u[j*(nx+1)+i+1]-u[j*(nx+1)+i])/h:0;G.vy[c]=act[c]?(v[(j+1)*nx+i]-v[j*nx+i])/h:0;}
    for(let j=1;j<ny;j++)for(let i=1;i<nx;i++){const a=(j-1)*nx+i-1;if(!(act[a]&&act[a+1]&&act[a+nx]&&act[a+nx+1]))continue;const q=j*(nx+1)+i;
      G.uy[q]=(u[j*(nx+1)+i]-u[(j-1)*(nx+1)+i])/h;G.vx[q]=(v[j*nx+i]-v[j*nx+i-1])/h;}
    return G;}
  // move markers with the sheet (midpoint rule) and strain them by the velocity gradient there: F <- exp(L dt) F.
  // The divergence (trace of L) thins (or, negative, thickens) the crust and grows (shrinks) the marker's area, so
  // det F stays equal to the area ratio; markers leaving the grid leave the budget
  advect(dt){const vel=[0,0],v2=[0,0],{nx,ny,o}=this,h=o.h,X1=this.x0+nx*h,Y1=this.y0+ny*h,G=this.velGrad(),E=[0,0,0,0],F=this.mF;
    for(let k=0;k<this.nm;k++){if(!this.mIn[k])continue;const x=this.mx[k],y=this.my[k];this.velAt(x,y,vel);const xm=x+0.5*dt*vel[0],ym=y+0.5*dt*vel[1];this.velAt(xm,ym,v2);
      const nxk=x+dt*v2[0],nyk=y+dt*v2[1];if(nxk<this.x0||nyk<this.y0||nxk>=X1||nyk>=Y1){this.mIn[k]=0;this.budget.left+=this.mS[k]*this.mA[k];continue;}
      const a=this.sample(G.ux,xm,ym),d=this.sample(G.vy,xm,ym),b=this.sampleN(G.uy,xm,ym),c=this.sampleN(G.vx,xm,ym);expm2(a*dt,b*dt,c*dt,d*dt,E);
      const f0=F[4*k],f1=F[4*k+1],f2=F[4*k+2],f3=F[4*k+3];F[4*k]=E[0]*f0+E[1]*f2;F[4*k+1]=E[0]*f1+E[1]*f3;F[4*k+2]=E[2]*f0+E[3]*f2;F[4*k+3]=E[2]*f1+E[3]*f3;
      const g=Math.exp((a+d)*dt);this.mS[k]/=g;this.mA[k]*=g;this.mx[k]=nxk;this.my[k]=nyk;}}
  // bilinear sample of a corner (node) field
  sampleN(f,x,y){const {nx,ny,o}=this,gx=(x-this.x0)/o.h,gy=(y-this.y0)/o.h,i=Math.max(0,Math.min(nx-1,Math.floor(gx))),j=Math.max(0,Math.min(ny-1,Math.floor(gy))),fx=Math.min(1,Math.max(0,gx-i)),fy=Math.min(1,Math.max(0,gy-j)),W=nx+1;
    return(1-fx)*(1-fy)*f[j*W+i]+fx*(1-fy)*f[j*W+i+1]+(1-fx)*fy*f[(j+1)*W+i]+fx*fy*f[(j+1)*W+i+1];}
  // Finite strain on the grid: per cell, the marker-weighted mean of the left Cauchy-Green tensor B = F F^T and of
  // the rotation, then B's principal stretches. lmin < 1 is horizontal shortening along the angle dir (radians from
  // map east); rot is the rotation about the vertical (counter-clockwise positive). Cached until the sheet moves.
  strainField(){const key=this.steps+":"+this.t;if(this._sf&&this._sf.key===key)return this._sf;const {nx,ny,o}=this,h=o.h,nc=nx*ny,F=this.mF;
    const bxx=new Float64Array(nc),bxy=new Float64Array(nc),byy=new Float64Array(nc),rot=new Float64Array(nc),ws=new Float64Array(nc);
    for(let k=0;k<this.nm;k++){if(!this.mIn[k])continue;const i=Math.floor((this.mx[k]-this.x0)/h),j=Math.floor((this.my[k]-this.y0)/h),c=j*nx+i,w=this.mA0[k];
      const f0=F[4*k],f1=F[4*k+1],f2=F[4*k+2],f3=F[4*k+3];bxx[c]+=w*(f0*f0+f1*f1);bxy[c]+=w*(f0*f2+f1*f3);byy[c]+=w*(f2*f2+f3*f3);rot[c]+=w*Math.atan2(f2-f1,f0+f3);ws[c]+=w;}
    const lmin=new Float64Array(nc).fill(NaN),lmax=new Float64Array(nc).fill(NaN),dir=new Float64Array(nc);
    for(let c=0;c<nc;c++){if(!(ws[c]>0)||!this.act[c])continue;const xx=bxx[c]/ws[c],xy=bxy[c]/ws[c],yy=byy[c]/ws[c],m=(xx+yy)/2,r=Math.hypot((xx-yy)/2,xy);
      lmax[c]=Math.sqrt(m+r);lmin[c]=Math.sqrt(Math.max(0,m-r));dir[c]=0.5*Math.atan2(-2*xy,-(xx-yy));rot[c]/=ws[c];}
    return this._sf={key,lmin,lmax,dir,rot};}
  // bilinear sample of a cell-centred field (inactive cells count as zero)
  sample(f,x,y){const {nx,ny,o}=this,gx=(x-this.x0)/o.h-0.5,gy=(y-this.y0)/o.h-0.5,i=Math.max(0,Math.min(nx-2,Math.floor(gx))),j=Math.max(0,Math.min(ny-2,Math.floor(gy))),fx=Math.min(1,Math.max(0,gx-i)),fy=Math.min(1,Math.max(0,gy-j));
    return(1-fx)*(1-fy)*f[j*nx+i]+fx*(1-fy)*f[j*nx+i+1]+(1-fx)*fy*f[(j+1)*nx+i]+fx*fy*f[(j+1)*nx+i+1];}
  // Elevation: the inherited elevation plus the thickness change, minus a flexural root (Gaussian of width
  // 1.18 x flexural length over the thickness change, the plate yielding beyond strengthKm from local balance)
  elevation(){const {nx,ny,o,S,Sref,eref,act,e,dS}=this,nc=nx*ny,sig=1.18*o.flexKm/o.h;
    for(let c=0;c<nc;c++)dS[c]=act[c]?S[c]-Sref[c]:0;
    const F=this.blur(dS,act,sig),lim=o.strengthKm;
    for(let c=0;c<nc;c++){if(!act[c]){e[c]=NaN;continue;}let rt=ROOT*F[c];const air=ROOT*dS[c];if(rt<air-lim)rt=air-lim;else if(rt>air+lim)rt=air+lim;e[c]=eref[c]+dS[c]-rt+this.dyn[c];}}
  // separable Gaussian blur of f over cells with mask set (normalized by the mask's own blur)
  blur(f,mask,sig){const {nx,ny}=this,nc=nx*ny,R=Math.max(1,Math.ceil(3*sig)),w=[];for(let k=-R;k<=R;k++)w.push(Math.exp(-k*k/(2*sig*sig)));
    const a=new Float64Array(nc),b=new Float64Array(nc),ma=new Float64Array(nc),mb=new Float64Array(nc);
    for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){let s=0,sm=0;for(let k=-R;k<=R;k++){const ii=i+k;if(ii<0||ii>=nx)continue;const c=j*nx+ii;if(!mask[c])continue;s+=w[k+R]*f[c];sm+=w[k+R];}a[j*nx+i]=s;ma[j*nx+i]=sm;}
    for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){let s=0,sm=0;for(let k=-R;k<=R;k++){const jj=j+k;if(jj<0||jj>=ny)continue;const c=jj*nx+i;s+=w[k+R]*a[c];sm+=w[k+R]*ma[c];}b[j*nx+i]=sm>0?s/sm:0;}
    return b;}
  // L0's arcs and hotspots as new crust on the markers beneath them: each marker thickens by what L0 gained there,
  // interpolated bilinearly between L0's cells (a lookup would lay the crust down in blocks of L0's cells)
  addVolcanism(rec){const w=this.world,Qr=this.Qref,b=[0,0,0],t=[0,0,0],out=[0,0,0],nb=w.cube.nb;let done=0;
    for(let k=0;k<this.nm;k++){if(!this.mIn[k])continue;this.proj.inv(this.mx[k],this.my[k],b);mv(Qr,b[0],b[1],b[2],t);w.cube.lookupFrac(t[0],t[1],t[2],out);
      const c=out[0],tx=out[1],ty=out[2],cB=nb[4*c+(tx>=0?0:1)],cC=nb[4*c+(ty>=0?2:3)],cD=nb[4*cB+(ty>=0?2:3)],ax=Math.abs(tx),ay=Math.abs(ty);
      const a=(1-ax)*(1-ay)*rec[c]+ax*(1-ay)*rec[cB]+(1-ax)*ay*rec[cC]+ax*ay*rec[cD];if(!(a>0))continue;
      const col=this.mCol[k],top=col[col.length-1],fr=a/(this.mS[k]+a);for(const u of col)u[2]*=1-fr;
      // a new top unit of arc or hotspot rock (arcs are by far the commoner source on continents)
      if(top&&top[0]===SIM.ROCK.ARC)top[2]+=fr;else col.push([SIM.ROCK.ARC,w.t,fr]);this.mS[k]+=a;done+=a*this.mA[k];}
    this.budget.arcsHotspots+=done;}
  // bilinear sample of a cell-centred field over crust cells only (renormalized; NaN where no crust is near)
  sampleAct(f,x,y){const {nx,ny,o,act}=this,gx=(x-this.x0)/o.h-0.5,gy=(y-this.y0)/o.h-0.5,i0=Math.floor(gx),j0=Math.floor(gy),fx=gx-i0,fy=gy-j0;let s=0,w=0;
    for(let q=0;q<4;q++){const i=i0+(q&1),j=j0+(q>>1);if(i<0||j<0||i>=nx||j>=ny)continue;const c=j*nx+i;if(!act[c])continue;const ww=(q&1?fx:1-fx)*(q>>1?fy:1-fy);s+=ww*f[c];w+=ww;}
    return w>0?s/w:NaN;}
  // erosion on the sheet: denudation rising steeply with elevation, ero5 km/Myr at 5 km times (e / 5 km)^eroExp
  // (0.5 km/Myr in high ranges, about 20 m/Myr on lowlands at 0.5 km: the magnitudes of Ahnert's relation between
  // denudation and relief). A smooth stand-in for rivers, which the sheet's grid is too coarse to draw (the drainage
  // pass draws them at the end). Eroded rock leaves the system, stripped from the top of the rock columns.
  // Each marker loses the depth interpolated (bilinearly) to where it is: a depth per cell shared among the cell's
  // markers would change from cell to cell in steps, and the markers would carry the steps away as blocks.
  erode(dt){const {nx,ny,o,act,e}=this,nc=nx*ny,dc=this._dc||(this._dc=new Float64Array(nc));dc.fill(0);
    for(let c=0;c<nc;c++){if(!act[c]||!(e[c]>0))continue;const d=Math.min(e[c],o.ero5*Math.pow(e[c]/5,o.eroExp)*dt);dc[c]=d;this.eroded[c]+=d;}
    let rem=0;
    for(let k=0;k<this.nm;k++){if(!this.mIn[k])continue;const d0=this.sampleAct(dc,this.mx[k],this.my[k]);if(!(d0>0))continue;
      const d=Math.min(d0,0.5*this.mS[k]),cut=d/this.mS[k],col=this.mCol[k];let left=cut;
      while(left>0&&col.length>1){const top=col[col.length-1];if(top[2]>left){top[2]-=left;left=0;}else{left-=top[2];col.pop();}}
      if(left>0&&col.length)col[0][2]-=left;const sum=col.reduce((a,u)=>a+u[2],0)||1;for(const u of col)u[2]/=sum;
      this.mS[k]-=d;rem+=d*this.mA[k];}
    this.budget.eroded+=rem;}
  // one coupled step: L0 advances (recording its tectonics), then the sheet runs in substeps, then volcanism, erosion
  advance(){const w=this.world,N=w.N,rec=this._rec||(this._rec=new Float32Array(N));
    rec.fill(0);w.tecRec=rec;w.step();w.tecRec=null;const dt=w.dt;
    this.locate();this.Qref=this.refQ();
    // substeps limited by a Courant number of 0.5 on the fastest marker
    let vmax=1;for(let k=0;k<this.nm;k++)if(this.mIn[k]){const s=Math.hypot(this.mUx[k],this.mUy[k]);if(s>vmax)vmax=s;}
    const ns=Math.max(1,Math.ceil(dt*vmax/(0.5*this.o.h)));
    for(let s=0;s<ns;s++){this.plateVelocities();this.deposit();this.solve();this.advect(dt/ns);}
    this.deposit();this.addVolcanism(rec);this.deposit();this.updateDyn();this.elevation();this.tectonicRate(dt);this.erode(dt);this.deposit();this.elevation();
    this.ePrev.set(this.e);this.t+=dt;this.steps++;this.stats(dt,ns);}
  // dynamic topography on the sheet's cells, from L0's slab record now (the mantle frame moves under the map)
  updateDyn(){if(!this.o.dyn)return;const D=dynTopo(this.world),{nx,ny,o}=this,h=o.h,b=[0,0,0],t=[0,0,0],out=[0,0,0],Qr=this.Qref;this.dynD=D;
    for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){this.proj.inv(this.x0+(i+.5)*h,this.y0+(j+.5)*h,b);mv(Qr,b[0],b[1],b[2],t);this.dyn[j*nx+i]=dynAt(D,t[0],t[1],t[2],out);}}
  // The present rate of elevation change from tectonics alone (deformation, volcanism and their isostasy; not
  // erosion): this step's change before erosion, averaged over the last rateMyr. The drainage pass applies it.
  tectonicRate(dt){const {e,ePrev,uRate}=this,a=Math.min(1,dt/this.o.rateMyr);
    for(let c=0;c<e.length;c++){const r=isNaN(e[c])||isNaN(ePrev[c])?0:(e[c]-ePrev[c])/dt;uRate[c]=(1-a)*uRate[c]+a*r;}}
  stats(dt,ns){const {nx,ny,o,act,e,S,dS,eps}=this,nc=nx*ny;let land=0,se=0,mx=-1e9,smax=0,dsmax=0,emax=0;
    for(let c=0;c<nc;c++){if(!act[c])continue;if(S[c]>smax)smax=S[c];if(dS[c]>dsmax)dsmax=dS[c];if(eps[c]>emax)emax=eps[c];if(e[c]>0){land++;se+=e[c];if(e[c]>mx)mx=e[c];}}
    this.history.push({t:this.t,landKm2:land*o.h*o.h,meanElevKm:land?se/land:0,maxElevKm:land?mx:0,maxThicknessKm:smax,maxThickeningKm:dsmax,maxStrainRate:emax,substeps:ns||0,cg:this.cgIters||0});}
  balance(){const b=this.budget,now=this.volume(),expect=b.initial+b.arcsHotspots-b.eroded-b.left;return{now,expect,error:(now-expect)/b.initial};}
  // the surface rock of each cell: the commonest top unit among its markers, by volume
  surfaceRock(){const {nx,ny,o}=this,h=o.h,acc=new Float64Array(nx*ny*8);
    for(let k=0;k<this.nm;k++){if(!this.mIn[k])continue;const i=Math.floor((this.mx[k]-this.x0)/h),j=Math.floor((this.my[k]-this.y0)/h),col=this.mCol[k];if(!col.length)continue;acc[(j*nx+i)*8+col[col.length-1][0]]+=this.mS[k]*this.mA[k];}
    for(let c=0;c<nx*ny;c++){let best=0,bv=0;for(let r=1;r<8;r++)if(acc[c*8+r]>bv){bv=acc[c*8+r];best=r;}this.rockTop[c]=best;}return this.rockTop;}
  // map plane -> the L0 elevation there now (for comparisons with L0)
  l0ElevAt(x,y,tmp){const w=this.world,b=this.proj.inv(x,y,tmp||[0,0,0]),wv=mv(this.Qref,b[0],b[1],b[2],[0,0,0]),out=[0,0,0],nb=w.cube.nb;w.cube.lookupFrac(wv[0],wv[1],wv[2],out);
    const c=out[0],tx=out[1],ty=out[2],cB=nb[4*c+(tx>=0?0:1)],cC=nb[4*c+(ty>=0?2:3)],cD=nb[4*cB+(ty>=0?2:3)],ax=Math.abs(tx),ay=Math.abs(ty),E=w.elev;
    return(1-ax)*(1-ay)*E[c]+ax*(1-ay)*E[cB]+(1-ax)*ay*E[cC]+ax*ay*E[cD];}
}

// ---------- climate for the drainage pass ----------
// A deliberately crude, present-day climate on the fine grid, enough to tell wet from dry:
// - Rain far from relief follows latitude as on Earth: wet equator, dry subtropics, wet westerly belts, dry poles.
// - Winds are zonal: easterly trades below 30 degrees, westerlies to 60, polar easterlies beyond.
// - Moisture blows in from the sea and rains out with distance (1/Lr per km; Lr is long because vegetation recycles
//   rain back into the air) and with uplift (1/He per km climbed),
//   and rain is proportional to the moisture lost, so windward slopes are wet and lee sides lie in a rain shadow.
// - Potential evaporation peaks in the clear-sky subtropics and falls with elevation; land loses water by Fu's form
//   of Budyko's curve.
// Units: m/yr. Returns latitude (radians), rain P, potential evaporation E and runoff R per cell.
const P_LAT=[[0,2.2],[10,1.8],[20,0.6],[28,0.35],[40,0.9],[55,1.0],[70,0.4],[90,0.2]];
// potential evaporation (m/yr) at sea level: highest under the clear skies of the subtropics, lower in the cloudy tropics
const E_LAT=[[0,1.4],[15,1.8],[25,2.0],[35,1.5],[50,0.8],[65,0.4],[90,0.2]];
function lerpTab(T,x){if(x<=T[0][0])return T[0][1];for(let k=1;k<T.length;k++)if(x<=T[k][0]){const [x0,y0]=T[k-1],[x1,y1]=T[k];return y0+(y1-y0)*(x-x0)/(x1-x0);}return T[T.length-1][1];}
// separable Gaussian blur of a grid field (edges renormalized)
function blurGrid(f,nx,ny,sig){const R=Math.max(1,Math.ceil(3*sig)),w=[];for(let k=-R;k<=R;k++)w.push(Math.exp(-k*k/(2*sig*sig)));
  const a=new Float32Array(nx*ny),b=new Float32Array(nx*ny);
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){let s=0,sw=0;for(let k=-R;k<=R;k++){const ii=i+k;if(ii<0||ii>=nx)continue;s+=w[k+R]*f[j*nx+ii];sw+=w[k+R];}a[j*nx+i]=s/sw;}
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){let s=0,sw=0;for(let k=-R;k<=R;k++){const jj=j+k;if(jj<0||jj>=ny)continue;s+=w[k+R]*a[jj*nx+i];sw+=w[k+R];}b[j*nx+i]=s/sw;}
  return b;}
function climate(sheet,nx,ny,h,e,land,o){
  const nc=nx*ny,lat=new Float32Array(nc),wx=new Float32Array(nc),wy=new Float32Array(nc),Qr=sheet.Qref,proj=sheet.proj,b=[0,0,0],p=[0,0,0],eb=[0,0,0],q0=[0,0],q1=[0,0],dl=1e-3;
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const c=j*nx+i,x=sheet.x0+(i+.5)*h,y=sheet.y0+(j+.5)*h;proj.inv(x,y,b);mv(Qr,b[0],b[1],b[2],p);lat[c]=Math.asin(Math.max(-1,Math.min(1,p[2])));
    if(!land[c])continue;const L=Math.hypot(p[0],p[1])||1;mtv(Qr,-p[1]/L,p[0]/L,0,eb);proj.fwd(b,q0);b[0]+=dl*eb[0];b[1]+=dl*eb[1];b[2]+=dl*eb[2];proj.fwd(b,q1);
    const dx=q1[0]-q0[0],dy=q1[1]-q0[1],n=Math.hypot(dx,dy)||1;wx[c]=dx/n;wy[c]=dy/n;}
  // moisture q (1 over the sea), carried downwind: Gauss-Seidel sweeps in the four grid orders, once for a westerly
  // and once for an easterly wind everywhere; the two are blended by latitude (the belts' edges are gradual)
  const qm=new Float32Array(nc),P=new Float32Array(nc),E=new Float32Array(nc),R=new Float32Array(nc);let sgn=1;
  const samp=(f,x,y)=>{const gx=x-0.5,gy=y-0.5,i0=Math.max(0,Math.min(nx-2,Math.floor(gx))),j0=Math.max(0,Math.min(ny-2,Math.floor(gy))),fx=Math.min(1,Math.max(0,gx-i0)),fy=Math.min(1,Math.max(0,gy-j0)),c=j0*nx+i0;
    return(1-fx)*(1-fy)*f[c]+fx*(1-fy)*f[c+1]+(1-fx)*fy*f[c+nx]+fx*fy*f[c+nx+1];};
  // rain responds to relief tens of km across, not to the fine relief: uplift is measured on a smoothed surface
  const eL=blurGrid(Float32Array.from(e,(v,c)=>land[c]?Math.max(0,v):0),nx,ny,o.rainSmoothKm/h);
  const lose=new Float32Array(nc);
  const upd=(i,j)=>{const c=j*nx+i;if(!land[c])return;const x=i+0.5-sgn*wx[c],y=j+0.5-sgn*wy[c],qu=samp(qm,x,y),eu=samp(eL,x,y);
    const k=h/o.Lr+Math.max(0,eL[c]-eu)/o.He;qm[c]=qu*Math.exp(-k);lose[c]=Math.min(o.orogMax/o.Lr,(qu-qm[c])/h);};
  const loseW=new Float32Array(nc);
  for(sgn of [1,-1]){qm.fill(1);lose.fill(0);
    for(let it=0;it<2;it++){for(let j=0;j<ny;j++)for(let i=0;i<nx;i++)upd(i,j);for(let j=0;j<ny;j++)for(let i=nx-1;i>=0;i--)upd(i,j);
      for(let j=ny-1;j>=0;j--)for(let i=0;i<nx;i++)upd(i,j);for(let j=ny-1;j>=0;j--)for(let i=nx-1;i>=0;i--)upd(i,j);}
    if(sgn===1)loseW.set(lose);}
  // westerly weight: 0 in the trades and polar easterlies, 1 in the westerlies, ramps 15 degrees wide at 30 and 60
  const west=a=>Math.min(1,Math.max(0,(a-22.5)/15))*Math.min(1,Math.max(0,(67.5-a)/15));
  for(let c=0;c<nc;c++){const a=Math.abs(lat[c])*180/Math.PI,ww=west(a);P[c]=lerpTab(P_LAT,a)*(land[c]?(ww*loseW[c]+(1-ww)*lose[c])*o.Lr:1);}
  P.set(blurGrid(P,nx,ny,o.rainSmoothKm/h));
  for(let c=0;c<nc;c++){E[c]=lerpTab(E_LAT,Math.abs(lat[c])*180/Math.PI)*Math.max(0.2,1-eL[c]/6);
    // Fu (1981): actual evaporation from land, given rain and potential evaporation
    const r=E[c]/Math.max(1e-6,P[c]),aet=P[c]*(1+r-Math.pow(1+Math.pow(r,o.fuW),1/o.fuW));R[c]=Math.max(0,P[c]-aet);}
  return{lat,P,E,R};
}
// Langbein and Schumm (1958): sediment yield against rain rises from deserts to a peak in semi-arid climates (about
// 0.3 m/yr) and falls as vegetation closes in; scaled to 1 at the peak
function langbeinSchumm(P){const f=x=>Math.pow(x,2.3)/(1+Math.pow(x,3.33));return f(1.273*P/0.3)/f(1.273);}
// ---------- the drainage pass ----------
// With the sheet's present tectonic rates (uplift, subsidence) kept running, its surface is refined (bilinear, plus
// relief the sheet cannot resolve, fixed to the planet) and eroded for a few Myr by stream power on a finer grid, to
// settle the network. Then the water is followed down it:
// - Discharge is runoff (with climate on), and it also weights erosion.
// - Lakes trap sediment: a depression its catchment's sediment can fill within the pass (river erosion plus
//   background denudation) becomes a flat the water crosses. Depressions survive where they subside or are big.
// - A surviving lake overflows only if its inflow beats evaporation minus rain over its full surface; otherwise
//   its basin ends there (drains inland), and it shrinks to the area its inflow can supply (the rest a dry floor).
// - In dry country rivers lose water into their beds and floodplains, about the same per km whatever their size
//   (so small rivers die and big ones cross deserts); where a river runs dry the basin ends at a terminal fan.
// - Sea cut off from the open ocean (and land below sea level) is inland water, under the same rules as lakes.
// - Land none of whose water reaches the coast (dry country where no river forms) is areic: it drains nowhere.
function basins(sheet,o){
  o=Object.assign({km:10,myr:5,steps:5,K:0.1,m:0.5,relief:0.15,jitter:0.4,climate:true,Lr:3000,He:3,orogMax:5,fuW:2.6,rainSmoothKm:30,lakeMinKm:0.002,
    bgDenudation:0.02,rates:true,lossMax:0.3,riverMin:30},o);
  const f=Math.max(1,Math.round(sheet.o.h/o.km)),nx=sheet.nx*f,ny=sheet.ny*f,h=sheet.o.h/f,nc=nx*ny,e=new Float64Array(nc),land=new Uint8Array(nc),U=new Float64Array(nc);
  const nz=makeNoise3(hashStr(String(sheet.world.opts.seed)+"|L1basins"+(o.salt||""))|0),b=[0,0,0],H=sheet.o.h,snx=sheet.nx,sny=sheet.ny,se=sheet.e,f0=R_KM/(3*H),oct=4,ur=o.rates&&sheet.uRate;
  const open=new Uint8Array(nc);
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const x=sheet.x0+(i+.5)*h,y=sheet.y0+(j+.5)*h,gx=(x-sheet.x0)/H-0.5,gy=(y-sheet.y0)/H-0.5;
    let i0=Math.floor(gx),j0=Math.floor(gy);const fx=gx-i0,fy=gy-j0;let s=0,w=0,u=0;
    for(const [di,dj,ww] of [[0,0,(1-fx)*(1-fy)],[1,0,fx*(1-fy)],[0,1,(1-fx)*fy],[1,1,fx*fy]]){const ii=i0+di,jj=j0+dj;if(ii<0||jj<0||ii>=snx||jj>=sny)continue;const v=se[jj*snx+ii];if(isNaN(v))continue;s+=ww*v;w+=ww;if(ur)u+=ww*ur[jj*snx+ii];}
    const c=j*nx+i;if(w<0.5){e[c]=-1;open[c]=1;continue;}let v=s/w;U[c]=u/w;sheet.proj.inv(x,y,b);const n=nz(b[0],b[1],b[2],f0,oct);v+=(o.relief+o.relief*Math.max(0,v))*2*(n-0.5);e[c]=v;land[c]=v>0?1:0;}
  // the open ocean: sea connected to the crust's outside or the grid's edge; other sea is inland water
  {const q=[];for(let c=0;c<nc;c++){const i=c%nx,j=(c/nx)|0;if(open[c]||(!land[c]&&(i===0||j===0||i===nx-1||j===ny-1))){open[c]=1;q.push(c);}}
    for(let k=0;k<q.length;k++){const c=q[k],i=c%nx,j=(c/nx)|0;for(const d of [i>0?-1:0,i<nx-1?1:0,j>0?-nx:0,j<ny-1?nx:0]){const n=c+d;if(d&&!open[n]&&!land[n]){open[n]=1;q.push(n);}}}}
  const inlandSea=new Uint8Array(nc);let nSea=0;for(let c=0;c<nc;c++)if(!land[c]&&!open[c]){inlandSea[c]=1;land[c]=1;nSea++;}
  const G=planarGrid(nx,ny,h,o.jitter,Math.round(sheet.x0/h),Math.round(sheet.y0/h),hashStr(String(sheet.world.opts.seed)+"|L1nodes")|0),S=scratch(nc),alive=new Uint8Array(nc).fill(1),out=new Uint8Array(nc);
  for(let c=0;c<nc;c++)out[c]=land[c]?0:1;
  // climate from the surface as it starts (the network settles on it); erosion weights runoff by its land mean
  const cl=o.climate?climate(sheet,nx,ny,h,e,land,o):null;let wgt=null;
  if(cl){let sr=0,nl=0;for(let c=0;c<nc;c++)if(land[c]){sr+=cl.R[c];nl++;}const mr=sr/(nl||1)||1;wgt=new Float32Array(nc);for(let c=0;c<nc;c++)wgt[c]=cl.R[c]/mr;}
  // erodibility: uniform, or o.Kmul, a multiplier per sheet cell, interpolated in log
  let Kc=o.K;if(o.Kmul){Kc=new Float64Array(nc);const lk=Float64Array.from(o.Kmul,v=>Math.log(v));
    for(let j=0;j<ny;j++)for(let i=0;i<nx;i++)Kc[j*nx+i]=o.K*Math.exp(sheet.sample(lk,sheet.x0+(i+.5)*h,sheet.y0+(j+.5)*h));}
  // the pass: tectonic rates, then erosion, each step; ero records the rock each cell lost (sediment supply)
  const ero=new Float64Array(nc),dt=o.myr/o.steps;
  for(let s=0;s<o.steps;s++){for(let c=0;c<nc;c++)if(!out[c])e[c]+=U[c]*dt;route(G,e,alive,out,S,wgt);const hN=streamPower(e,S,Kc,o.m,dt);
    for(let c=0;c<nc;c++)if(!out[c]){ero[c]+=Math.max(0,e[c]-hN[c]);e[c]=hN[c];}}
  route(G,e,alive,out,S,null);
  const rec=Int32Array.from(S.rec),stack=S.stack.subarray(0,S.ns),A=new Float64Array(nc),Q=new Float64Array(nc),lake=new Uint8Array(nc);let nClosed=0,nDry=0;
  // lakes: cells the filling raised by more than lakeMinKm (or inland sea); each lake's exit is a lake cell draining out
  const fl=S.f,isLake=c=>land[c]&&(fl[c]-e[c]>o.lakeMinKm||inlandSea[c]);
  const lakeArea=new Float64Array(nc),lakeLoss=new Float64Array(nc),lakeVol=new Float64Array(nc),exitOf=new Int32Array(nc).fill(-1),sed=new Float64Array(nc);
  for(let q=0;q<stack.length;q++){const l=stack[q],r=rec[l];if(!isLake(l))continue;exitOf[l]=r!==l&&isLake(r)?exitOf[r]:l;}
  for(let c=0;c<nc;c++)if(exitOf[c]>=0){const x=exitOf[c];lakeArea[x]+=h*h;lakeVol[x]+=Math.max(0,fl[c]-e[c])*h*h;if(cl)lakeLoss[x]+=(cl.E[c]-cl.P[c])*h*h;lake[c]=1;}
  // Down the network from the top: area, discharge (runoff, minus channel losses) and sediment
  const filled=new Uint8Array(nc),qOut=new Float64Array(nc),sedIn=new Float64Array(nc);
  for(let q=stack.length-1;q>=0;q--){const l=stack[q];A[l]+=h*h;
    if(land[l]){if(cl)Q[l]+=cl.R[l]*h*h;sed[l]+=(ero[l]+o.bgDenudation*(cl?langbeinSchumm(cl.P[l]):1)*o.myr)*h*h;}
    let lk=lakeArea[l]>0&&rec[l]!==l;
    if(lk){sedIn[l]=sed[l];if(sed[l]>=lakeVol[l]&&!inlandSea[l]){sed[l]-=lakeVol[l];filled[l]=1;lk=false;}else sed[l]=0;}
    if(lk&&cl&&Q[l]<=lakeLoss[l]){rec[l]=l;nClosed++;lake[l]=2;continue;}
    let qs=lk&&cl?Math.max(0,Q[l]-lakeLoss[l]):Q[l];
    // channel losses in dry country: flood pulses soak into beds and floodplains, about the same per km whatever the
    // river's size (small rivers die, big ones cross deserts), rising from 0 at an aridity index (P/E) of 0.65
    const r=rec[l];
    if(cl&&r!==l&&land[l]&&!lake[l]){const ai=cl.P[l]/Math.max(1e-6,cl.E[l]),k=o.lossMax*Math.max(0,1-ai/0.65);
      if(k>0){const qin=qs;qs=Math.max(0,qs-k*S.rdist[l]);
        if(qin>=o.riverMin&&qs<o.riverMin){rec[l]=l;nDry++;lake[l]=4;continue;}}}
    qOut[l]=qs;if(r!==l){A[r]+=A[l];sed[r]+=sed[l];Q[r]+=qs;}}
  // filled depressions become sediment flats at their spill level
  let nFilled=0;for(let c=0;c<nc;c++){const x=exitOf[c];if(x>=0&&filled[x]){lake[c]=0;e[c]=fl[c];if(c===x)nFilled++;}}
  // labels, now that some basins end inland; closed lakes shrink to the area their inflow can supply
  const basin=new Int32Array(nc).fill(-1),size=new Map(),closedB=new Uint8Array(nc);
  for(let q=0;q<stack.length;q++){const l=stack[q],r=rec[l];basin[l]=r===l?l:basin[r];}
  if(cl){const cells=new Map();for(let c=0;c<nc;c++)if(exitOf[c]>=0&&lake[exitOf[c]]===2){const x=exitOf[c];if(!cells.has(x))cells.set(x,[]);cells.get(x).push(c);}
    for(const [x,cs] of cells){cs.sort((a,b2)=>e[a]-e[b2]);let need=Q[x];
      for(const c of cs){const loss=(cl.E[c]-cl.P[c])*h*h;if(need>0&&loss>0){need-=loss;lake[c]=2;}else lake[c]=3;}lake[x]=lake[x]===3?2:lake[x];}
    // Land drains inland if its basin ends in a lake with no outlet or where a river runs dry (endorheic, 1). Otherwise
    // it is areic (2) if its water never reaches the coast: some cell on its way passes on no water at all (the
    // wash dies in the desert before it joins a river that reaches the sea). Receivers come first in the stack.
    const dead=new Uint8Array(nc);for(let q=0;q<stack.length;q++){const l=stack[q],r=rec[l];if(!land[l])continue;dead[l]=!(qOut[l]>0)||(r!==l&&land[r]&&dead[r])?1:0;}
    for(let c=0;c<nc;c++){if(!land[c])continue;const t=basin[c];if(!out[t]&&lake[t]>=2)closedB[c]=1;else if(dead[c])closedB[c]=2;}}
  for(let c=0;c<nc;c++)if(land[c]){const id=basin[c];size.set(id,(size.get(id)||0)+1);}
  const sizes=[...size.values()].sort((a,b2)=>b2-a);let inland=0,areic=0,nl=0;for(let c=0;c<nc;c++)if(land[c]){nl++;if(closedB[c])inland++;if(closedB[c]===2)areic++;}
  return{nx,ny,h,x0:sheet.x0,y0:sheet.y0,e,land,U,A,Q,qOut,rec,basin,lake,closedB,inlandSea,exitOf,lakeVol,sedIn,K:o.K,m:o.m,myr:o.myr,relief:o.relief,climate:cl,nClosed,nDry,nFilled,nSea,inlandFrac:inland/(nl||1),areicFrac:areic/(nl||1),nBasins:size.size,sizes,landCells:nl,S};
}
// ---------- the L1 -> L2 handoff ----------
// L1 hands L2 its landmass cut into tiles, with what L2 needs to refine each tile on its own:
// - Tiles are drainage basins, so water never crosses a tile edge except at an outlet. Each basin of at least tileKm2
//   is a tile; smaller basins (coastal strips, small inland sinks) are grouped with their neighbours into tiles of
//   about that size. A tile lists its outlets, lakes, and the divides it shares with each neighbouring tile.
// - Fields come at the resolution they were computed: the drainage pass's fine grid (elevation, receivers,
//   discharge, climate, uplift, knickpoint response time) and the sheet's coarse grid (thickness, finite strain,
//   mantle pull, denudation, and the rock stacks of the L0 parcels present).
// - Imprints: what happened to the landmass, each with its class and age (Ma, before the present), oldest first, so
//   L2 can apply them in order (younger structures cut older ones). Structural imprints (rock juxtaposed at a fault,
//   strain fabric) last until erosion strips the rock that carries them; topographic ones decay, at the rates in the
//   lifetimes table.
const LAT=(v)=>[Math.asin(Math.max(-1,Math.min(1,v[2])))*180/Math.PI,Math.atan2(v[1],v[0])*180/Math.PI];
const D8=["E","W","N","S","NE","NW","SE","SW"],M3S=1e6/3.15576e7; // planarGrid's neighbour order; (m/yr x km^2) -> m^3/s
function handoff(sheet,B,o){
  o=Object.assign({tileKm2:1e4,D:0.01,parcels:3},o);
  const w=sheet.world,{nx,ny,h,land,rec,basin,lake,closedB,e,Q}=B,nc=nx*ny,S=B.S,cl=B.climate,minCells=Math.round(o.tileKm2/(h*h)),f=Math.round(sheet.o.h/h);
  const now=w.t,window=sheet.t,r1=v=>Math.round(v*10)/10,r2=v=>Math.round(v*100)/100,r3=v=>Math.round(v*1000)/1000;
  // tiles: major basins, then minor basins grown over their neighbours (first cell order, breadth first) to minCells
  const size=new Map(),first=[];for(let c=0;c<nc;c++)if(land[c]){const b=basin[c];if(!size.has(b)){size.set(b,0);first.push(b);}size.set(b,size.get(b)+1);}
  const tOf=new Map(),tiles=[],minor=b=>size.get(b)<minCells,adj=new Map();
  for(const b of first)if(!minor(b)){tOf.set(b,tiles.length);tiles.push({basins:[b],major:true,n:size.get(b)});}
  // basin adjacency (8 neighbours; count of contacts), for growing groups and placing leftovers
  const link=(a,b)=>{let m=adj.get(a);if(!m)adj.set(a,m=new Map());m.set(b,(m.get(b)||0)+1);};
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const c=j*nx+i;if(!land[c])continue;const a=basin[c];
    for(const d of [i<nx-1?1:0,j<ny-1?nx:0,i<nx-1&&j<ny-1?nx+1:0,i>0&&j<ny-1?nx-1:0]){if(!d)continue;const b=basin[c+d];if(land[c+d]&&b!==a){link(a,b);link(b,a);}}}
  const nbrs=b=>adj.has(b)?[...adj.get(b).keys()].sort((x,y)=>x-y):[];
  for(const b of first){if(!minor(b)||tOf.has(b))continue;const id=tiles.length,T={basins:[b],major:false,n:size.get(b)};tOf.set(b,id);tiles.push(T);
    for(let q=0;q<T.basins.length&&T.n<minCells;q++)for(const x of nbrs(T.basins[q])){if(tOf.has(x)||!minor(x))continue;tOf.set(x,id);T.basins.push(x);T.n+=size.get(x);if(T.n>=minCells)break;}}
  // groups boxed in below half the size join the neighbouring tile they share the longest boundary with (a union of
  // whole basins is still closed to water); isolated islands stay as they are
  for(const T of tiles){if(T.major||T.n>=minCells/2||!T.basins.length)continue;const id=tOf.get(T.basins[0]),con=new Map();
    for(const b of T.basins)if(adj.has(b))for(const [x,cnt] of adj.get(b)){const k=tOf.get(x);if(k!==id)con.set(k,(con.get(k)||0)+cnt);}
    let best=-1,bc=0;for(const [k,cnt] of con)if(cnt>bc||cnt===bc&&k<best){bc=cnt;best=k;}
    if(best>=0){for(const b of T.basins){tOf.set(b,best);tiles[best].basins.push(b);}tiles[best].n+=T.n;T.basins=[];T.n=0;}}
  const keep=[],remap=new Int32Array(tiles.length).fill(-1);tiles.forEach((T,k)=>{if(T.basins.length){remap[k]=keep.length;keep.push(T);}});
  const tileOf=new Int32Array(nc).fill(-1);for(let c=0;c<nc;c++)if(land[c])tileOf[c]=remap[tOf.get(basin[c])];
  // knickpoint response time: the time a change in base level takes to travel up the network from the outlet, at the
  // celerity of stream power with n = 1, K (A_w)^m, A_w the climate-weighted drainage area (as in the pass)
  const Aw=new Float64Array(nc),tau=new Float64Array(nc),stack=S.stack.subarray(0,S.ns);let mr=1;
  if(cl){let s=0,n=0;for(let c=0;c<nc;c++)if(land[c]){s+=cl.R[c];n++;}mr=s/(n||1)||1;}
  for(let q=stack.length-1;q>=0;q--){const l=stack[q];Aw[l]+=h*h*(cl?cl.R[l]/mr:1);const r=rec[l];if(r!==l)Aw[r]+=Aw[l];}
  for(let q=0;q<stack.length;q++){const l=stack[q],r=rec[l];tau[l]=r===l||!land[l]?0:tau[r]+S.rdist[l]/(B.K*Math.pow(Math.max(Aw[l],1e-9),B.m));}
  // the sheet's fields: finite strain, present denudation, and per sheet cell the L0 parcels present (by volume)
  const sf=sheet.strainField(),{o:so}=sheet,snx=sheet.nx,sny=sheet.ny,H=so.h,snc=snx*sny,deno=new Float64Array(snc),arcWin=new Float64Array(snc);
  for(let c=0;c<snc;c++)deno[c]=sheet.act[c]&&sheet.e[c]>0?so.ero5*Math.pow(sheet.e[c]/5,so.eroExp):0;
  const cellP=new Array(snc),t0=sheet.t0;
  for(let k=0;k<sheet.nm;k++){if(!sheet.mIn[k])continue;const i=Math.floor((sheet.mx[k]-sheet.x0)/H),j=Math.floor((sheet.my[k]-sheet.y0)/H),c=j*snx+i,v=sheet.mS[k]*sheet.mA[k],pid=sheet.mPid[k];
    let m=cellP[c];if(!m)m=cellP[c]=new Map();const p=m.get(pid);if(!p)m.set(pid,{v,best:k,bv:v,S:v*sheet.mS[k]});else{p.v+=v;p.S+=v*sheet.mS[k];if(v>p.bv){p.bv=v;p.best=k;}}
    for(const u of sheet.mCol[k])if(u[0]===SIM.ROCK.ARC&&u[1]>=t0-1e-9)arcWin[c]+=u[2]*v/(H*H);}
  const stacksAt=c=>{const m=cellP[c];if(!m)return[];let tot=0;for(const p of m.values())tot+=p.v;
    return [...m.entries()].sort((a,b)=>b[1].v-a[1].v||a[0]-b[0]).slice(0,o.parcels).map(([pid,p])=>{const S0=p.S/p.v;
      return{parcel:pid,share:r3(p.v/tot),thicknessKm:r1(S0),units:sheet.mCol[p.best].map(u=>[SIM.ROCK_NAMES[u[0]],r1(now-u[1]),r2(u[2]*S0)])};});};
  // map plane <-> the present world frame (L0's): world = Qref . inv(x, y)
  const Qr=sheet.Qref,b=[0,0,0],t3=[0,0,0],toLL=(x,y)=>{sheet.proj.inv(x,y,b);mv(Qr,b[0],b[1],b[2],t3);const ll=LAT(t3);return[r2(ll[0]),r2(ll[1])];};
  const fineAt=(x,y)=>{const i=Math.floor((x-B.x0)/h),j=Math.floor((y-B.y0)/h);return i<0||j<0||i>=nx||j>=ny?-1:j*nx+i;};
  // Divides: a tile edge is firm where its local relief (height above the lowest land within one sheet cell) exceeds
  // the peak-to-peak amplitude of the drainage pass's synthetic relief, 2 (relief + relief e); elsewhere contingent,
  // and L2 may move it. How firm is graded, not sharp: re-seeding that relief on the 1 Gyr test landmass, a divide
  // still lies between points 30 km either side of it 77% of the time where firm, 51-53% where contingent, and two
  // random points as far apart straddle some divide 35% of the time (meta.dividePersistence has it by relief).
  const Rw=Math.max(1,Math.round(H/h)),mn1=new Float64Array(nc),rel=new Float64Array(nc),divide=new Uint8Array(nc);
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){let m=Infinity;for(let k=-Rw;k<=Rw;k++){const ii=i+k;if(ii<0||ii>=nx)continue;const c=j*nx+ii;if(land[c]&&e[c]<m)m=e[c];}mn1[j*nx+i]=m;}
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){let m=Infinity;for(let k=-Rw;k<=Rw;k++){const jj=j+k;if(jj<0||jj>=ny)continue;m=Math.min(m,mn1[jj*nx+i]);}const c=j*nx+i;rel[c]=land[c]?e[c]-m:0;}
  const firmAt=(a,c)=>{const r=Math.max(rel[a],rel[c]),ec=Math.max(0,e[a],e[c]);return r>2*B.relief*(1+ec);};
  // per tile: extent, area, outlets, lakes, divides, summaries
  const R=keep.map((T,id)=>({id,kind:"",major:T.major,cells:0,i0:nx,j0:ny,i1:-1,j1:-1,sx:0,sy:0,eMax:-1e9,eSum:0,rain:0,inland:0,areic:0,tauMax:0,tauL:[],shorten:0,arc:0,dynMin:1e9,dynMax:-1e9,coast:0,div:new Map(),outlets:[],lakes:[],imprints:[]}));
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const c=j*nx+i,k=tileOf[c];if(k<0)continue;const T=R[k];T.cells++;T.i0=Math.min(T.i0,i);T.j0=Math.min(T.j0,j);T.i1=Math.max(T.i1,i);T.j1=Math.max(T.j1,j);
    T.sx+=i;T.sy+=j;T.eMax=Math.max(T.eMax,e[c]);T.eSum+=e[c];if(cl)T.rain+=cl.P[c];if(closedB[c])T.inland++;if(closedB[c]===2)T.areic++;T.tauMax=Math.max(T.tauMax,tau[c]);T.tauL.push(tau[c]);
    const sc=Math.floor(j/f)*snx+Math.floor(i/f);if(!isNaN(sf.lmin[sc]))T.shorten=Math.max(T.shorten,1-sf.lmin[sc]);T.arc=Math.max(T.arc,arcWin[sc]);if(sheet.act[sc]){T.dynMin=Math.min(T.dynMin,sheet.dyn[sc]);T.dynMax=Math.max(T.dynMax,sheet.dyn[sc]);}
    for(const d of [i<nx-1?1:0,j<ny-1?nx:0,i>0?-1:0,j>0?-nx:0]){if(!d)continue;const n2=c+d;if(!land[n2]){T.coast+=h;continue;}const k2=tileOf[n2];if(k2===k)continue;
      const p=Math.max(e[c],e[n2]),fm=firmAt(c,n2),dv=T.div.get(k2);divide[c]=Math.max(divide[c],fm?2:1);if(!dv)T.div.set(k2,{len:h,pass:p,firm:fm?h:0});else{dv.len+=h;dv.pass=Math.min(dv.pass,p);if(fm)dv.firm+=h;}}}
  const xKm=i=>B.x0+(i+.5)*h,yKm=j=>B.y0+(j+.5)*h;
  for(const T of R){const K=keep[T.id];
    for(const bId of K.basins){const r=bId,ri=r%nx,rj=(r/nx)|0,kind=!land[r]?"sea":lake[r]===2?"closed lake":lake[r]===4?"river runs dry":S.sink&&S.sink[r]?"enclosed sink":"sea";
      T.outlets.push({kind,xKm:r1(xKm(ri)),yKm:r1(yKm(rj)),latLon:toLL(xKm(ri),yKm(rj)),areaKm2:Math.round(size.get(bId)*h*h),dischargeM3s:r2(Q[r]*M3S)});}
    T.outlets.sort((a,b2)=>b2.areaKm2-a.areaKm2);
    T.kind=T.major?(T.outlets[0].kind==="sea"?"drains to the sea":"drains inland"):T.outlets.some(x=>x.kind==="sea")?"coastal strip":"inland group";}
  // lakes: connected water (open, closed, or the dry floor of a shrunk lake), with their fill time by the sediment
  // reaching them (river erosion and background denudation over the pass)
  {const seen=new Uint8Array(nc);for(let c=0;c<nc;c++){if(seen[c]||!(lake[c]>=1&&lake[c]<=3)||tileOf[c]<0)continue;const q=[c];seen[c]=1;let water=0,playa=0,closed=false,x=-1;
      for(let k=0;k<q.length;k++){const a=q[k];if(lake[a]===3)playa++;else water++;if(lake[a]===2)closed=true;if(B.exitOf[a]>=0)x=B.exitOf[a];const ai=a%nx,aj=(a/nx)|0;
        for(const d of [ai>0?-1:0,ai<nx-1?1:0,aj>0?-nx:0,aj<ny-1?nx:0])if(d&&!seen[a+d]&&lake[a+d]>=1&&lake[a+d]<=3){seen[a+d]=1;q.push(a+d);}}
      const vol=x>=0?B.lakeVol[x]:0,rate=x>=0?B.sedIn[x]/B.myr:0,T=R[tileOf[c]];
      T.lakes.push({kind:closed?"no outlet":"open",inlandSea:!!B.inlandSea[c],waterKm2:water*h*h,dryFloorKm2:playa*h*h,volumeKm3:r1(vol),fillMyr:rate>0?r1(vol/rate):null,xKm:r1(xKm(c%nx)),yKm:r1(yKm((c/nx)|0))});}}
  // imprints, oldest first
  const imprints=[],add=(I,ts)=>{I.id=imprints.length;I.tiles=[...ts].sort((a,b2)=>a-b2);imprints.push(I);for(const k of I.tiles)R[k].imprints.push(I.id);};
  {const hf=w.extractFaults(),xy=[0,0];
    for(const F of hf.faults){const pts=F.trace.map(v=>{mtv(Qr,v[0],v[1],v[2],b);sheet.proj.fwd(b,xy);return[xy[0],xy[1]];}),ts=new Set();
      for(let q=0;q+1<pts.length;q++){const [ax,ay]=pts[q],[bx,by]=pts[q+1],n=Math.max(1,Math.ceil(Math.hypot(bx-ax,by-ay)/(h/2)));
        for(let s=0;s<=n;s++){const c=fineAt(ax+(bx-ax)*s/n,ay+(by-ay)*s/n);if(c>=0&&tileOf[c]>=0)ts.add(tileOf[c]);}}
      if(!ts.size)continue;
      add({class:"fault",startMa:r1(F.startMa),endMa:r1(F.endMa),ongoing:F.active,leaves:["structure","relief"],kind:F.kindName,sense:F.sense,slipAcrossKm:r1(F.slipAcrossKm),slipAlongKm:r1(F.slipAlongKm),
        traceKm:pts.map(p=>[r1(p[0]),r1(p[1])])},ts);}}
  const tsWhere=fn=>R.filter(fn).map(T=>T.id),ongoingStrain=(()=>{let m=0;for(let c=0;c<snc;c++)if(sheet.act[c])m=Math.max(m,sheet.eps[c]);return m>0.05*sheet.e0;})();
  add({class:"crumpling",startMa:r1(window),endMa:0,ongoing:ongoingStrain,leaves:["structure","relief"],note:"finite strain of the sheet over L1's window: coarse fields shortening, stretch, direction, rotation"},tsWhere(T=>T.shorten>0.1));
  add({class:"volcanism",startMa:r1(window),endMa:0,ongoing:true,leaves:["rock","relief"],note:"arc and hotspot rock added during the window: coarse field arcAddedKm, and the top units of the stacks"},tsWhere(T=>T.arc>0.1));
  add({class:"mantle pull",startMa:r1(SIM.C.SLAB_BINS*SIM.C.SLAB_DT),endMa:0,ongoing:true,leaves:["relief"],note:"dynamic topography from slabs subducted over this span; already in the elevations; coarse field mantlePullKm"},tsWhere(()=>true));
  add({class:"drainage",startMa:r1(B.myr),endMa:0,ongoing:true,leaves:["relief","water"],note:"the river network, lakes and inland drainage settled by the drainage pass"},tsWhere(()=>true));
  imprints.sort((a,b2)=>b2.startMa-a.startMa||b2.endMa-a.endMa);const reid=new Map();imprints.forEach((I,k)=>{reid.set(I.id,k);I.id=k;});for(const T of R)T.imprints=T.imprints.map(x=>reid.get(x)).sort((a,b2)=>a-b2);
  // lifetimes of topographic imprints, by erasure mechanism
  const p=so.eroExp,a5=so.ero5/Math.pow(5,p),half=(k,e0)=>r1(Math.pow(e0,1-p)*(Math.pow(2,p-1)-1)/(k*a5*(p-1)));
  const lifetimes=[
    {mechanism:"hillslope blurring",appliesTo:"features narrower than a valley: scarps, small cones, terrace risers",law:"a feature blurs over a width sigma = sqrt(2 D t), so a feature of width w lasts about w^2 / (2 D)",
      params:{D_km2PerMyr:o.D,note:"D in km^2/Myr equals D in m^2/yr; soil-mantled hillslopes range about 0.001-0.05 (literature value, not simulated)"},
      examples:[0.1,0.3,1].map(wk=>({widthKm:wk,lifetimeMyr:r1(wk*wk/(2*o.D))}))},
    {mechanism:"knickpoint sweep",appliesTo:"a change at the outlet (base level, capture, a fault scarp crossing a river): it travels up the network and is gone below the knickpoint",
      law:"celerity K A_w^m (stream power with n = 1); a cell's response time is the integral of distance / celerity from its outlet: fine field responseMyr",
      params:{K:B.K,m:B.m,note:"A_w is drainage area weighted by runoff relative to the land mean, km^2; K is L1's drainage-pass value, not calibrated to absolute rates"},
      examples:[100,1e4,1e6].map(A=>({areaKm2:A,celerityKmPerMyr:r1(B.K*Math.pow(A,B.m))}))},
    {mechanism:"relief decay",appliesTo:"mountains and plateaus as a whole",law:"L1's denudation law, rate = ero5 (e / 5 km)^p of rock; wide relief is isostatically compensated (the surface drops (1 - root share) of what is eroded), narrow relief is not; half-life from elevation e0",
      params:{ero5KmPerMyr:so.ero5,p,rootShare:r2(ROOT)},
      examples:[0.5,1,3].map(e0=>({elevationKm:e0,halfLifeWideMyr:half(1-ROOT,e0),halfLifeNarrowMyr:half(1,e0)}))},
    {mechanism:"lake filling",appliesTo:"lakes and enclosed seas",law:"volume / sediment supply (river erosion plus background denudation in the catchment): per lake, fillMyr",params:{}},
    {mechanism:"stripping",appliesTo:"rock units at the surface",law:"a unit of thickness d at the surface lasts about d / denudation: coarse field denudationKmPerMyr",params:{}}];
  const ref=toLL(0,0),north=toLL(0,100);
  const meta={generator:"L1 landmass (thin viscous sheet and drainage pass)",seed:w.opts.seed,presentMyr:r1(now),windowMyr:r1(window),
    units:{distance:"km",age:"Ma, million years before the present",elevation:"m on the fine grid, km on the coarse grid",discharge:"m^3/s",rain:"mm/yr",rate:"mm/yr = km/Myr"},
    frame:{projection:"azimuthal equidistant, co-moving with the landmass's reference plate",centreLatLon:ref,pointNorthOfCentre:north,
      note:"x east, y north at the centre (km). Latitude and longitude are in L0's present frame (lat = asin z, lon = atan2(y, x))."},
    grids:{fine:{cellKm:h,nx,ny,x0Km:B.x0,y0Km:B.y0,note:"cell (i, j) is centred at x0 + (i + 0.5) h, y0 + (j + 0.5) h"},coarse:{cellKm:H,nx:snx,ny:sny,x0Km:sheet.x0,y0Km:sheet.y0}},
    tiles:`drainage basins of at least ${o.tileKm2} km^2; smaller basins grouped with neighbours into tiles of about that size`,
    order:"imprints are listed oldest first (by start); L2 applies them in that order, so younger structures cut older ones",
    permanence:"structure (rock juxtaposed at faults, strain fabric) lasts until the rock carrying it is eroded; relief decays as in lifetimes",
    receivers:`receiver codes are the direction of steepest descent, ${D8.map((d,k)=>k+" "+d).join(", ")}; -1 is a basin end (lake with no outlet, river that runs dry, enclosed sink)`,
    divides:"a divide is firm where its local relief (height above the lowest land within one coarse cell) exceeds the peak-to-peak amplitude of the synthetic fine relief, 2 x "+B.relief+" x (1 + elevation in km); elsewhere it is contingent and L2 may move it. firmKm is the firm length of each shared divide. Firmness is graded, not sharp: see dividePersistence",
    dividePersistence:{test:"re-seeding the synthetic fine relief (two seeds) on the 1 Gyr test landmass: share of point pairs 30 km either side of a divide that still straddle a divide, against random pairs as far apart",
      firm:0.77,contingent:0.52,byLocalReliefKm:[["0-0.2",0.45,0.25],["0.2-0.4",0.54,0.38],["0.4-0.8",0.67,0.45],["0.8-1.6",0.76,0.5],[">1.6",0.9,0.62]],columns:["relief","divides","random pairs"]},
    caveats:["Which outlet a lowland cell reaches is contingent on relief the sheet cannot resolve (re-seeding that relief moves outlets a median 110-130 km); coasts and outlets to within a few hundred km are what L1 determines, and divides more firmly the more relief they have.",
      "Fault traces are L0's (100 km cells), placed where L0's rigid plates put them; in crumple zones the sheet has moved the crust by up to the shortening.",
      "Parcels in a coarse cell lie side by side (the sheet does not stack crust vertically); share is their share of the cell's volume."]};
  const outTiles=R.map(T=>{const tl=T.tauL.sort((a,b2)=>a-b2);return{id:T.id,kind:T.kind,areaKm2:T.cells*h*h,
    fine:{i0:T.i0,j0:T.j0,nx:T.i1-T.i0+1,ny:T.j1-T.j0+1},centreKm:[r1(xKm(T.sx/T.cells)),r1(yKm(T.sy/T.cells))],centreLatLon:toLL(xKm(T.sx/T.cells),yKm(T.sy/T.cells)),
    maxElevationM:Math.round(T.eMax*1000),meanElevationM:Math.round(T.eSum/T.cells*1000),meanRainMmYr:cl?Math.round(T.rain/T.cells*1000):null,
    inlandShare:r3(T.inland/T.cells),areicShare:r3(T.areic/T.cells),maxShortening:r2(T.shorten),arcAddedKm:r1(T.arc),mantlePullKm:[r2(T.dynMin),r2(T.dynMax)],
    responseMyr:{median:r1(tl[tl.length>>1]||0),max:r1(T.tauMax)},coastKm:T.coast,
    outlets:T.outlets,lakes:T.lakes,divides:[...T.div.entries()].sort((a,b2)=>a[0]-b2[0]).map(([k,d])=>({tile:k,lengthKm:d.len,firmKm:d.firm,passElevationM:Math.round(d.pass*1000)})),imprints:T.imprints};});
  // one tile's fields: the fine window over its extent and the coarse cells under it
  function tile(id){const T=outTiles[id],{i0,j0,nx:tx,ny:ty}=T.fine,n=tx*ty,F={inTile:new Array(n),elevationM:new Array(n),receiver:new Array(n),areaKm2:new Array(n),dischargeM3s:new Array(n),
      lake:new Array(n),drainage:new Array(n),divide:new Array(n),rainMmYr:new Array(n),evaporationMmYr:new Array(n),upliftMmYr:new Array(n),responseMyr:new Array(n)};
    for(let j=0;j<ty;j++)for(let i=0;i<tx;i++){const c=(j0+j)*nx+i0+i,q=j*tx+i,r=rec[c];let code=-1;
      if(r!==c){const dx=r%nx-c%nx,dy=((r/nx)|0)-((c/nx)|0);code=D8.indexOf((dy>0?"N":dy<0?"S":"")+(dx>0?"E":dx<0?"W":""));}
      F.inTile[q]=tileOf[c]===id?1:0;F.elevationM[q]=Math.round(e[c]*1000);F.receiver[q]=land[c]?code:-1;F.areaKm2[q]=Math.round(B.A[c]);F.dischargeM3s[q]=r3((land[c]?B.qOut[c]:Q[c])*M3S);
      F.lake[q]=lake[c];F.drainage[q]=closedB[c];F.divide[q]=divide[c];F.rainMmYr[q]=cl?Math.round(cl.P[c]*1000):null;F.evaporationMmYr[q]=cl?Math.round(cl.E[c]*1000):null;F.upliftMmYr[q]=r3(B.U[c]);F.responseMyr[q]=r1(tau[c]);}
    const si0=Math.floor(i0/f),sj0=Math.floor(j0/f),si1=Math.floor((i0+tx-1)/f),sj1=Math.floor((j0+ty-1)/f),cx=si1-si0+1,cy=sj1-sj0+1,G={i0:si0,j0:sj0,nx:cx,ny:cy,
      crust:[],thicknessKm:[],thickeningKm:[],shortening:[],stretch:[],shorteningAzimuthDeg:[],rotationDeg:[],strainRatePerMyr:[],mantlePullKm:[],denudationKmPerMyr:[],arcAddedKm:[],stacks:[]};
    for(let j=sj0;j<=sj1;j++)for(let i=si0;i<=si1;i++){const c=j*snx+i,a=sheet.act[c];G.crust.push(a);
      G.thicknessKm.push(a?r1(sheet.S[c]):null);G.thickeningKm.push(a?r1(sheet.dS[c]):null);G.shortening.push(isNaN(sf.lmin[c])?null:r3(1-sf.lmin[c]));G.stretch.push(isNaN(sf.lmax[c])?null:r3(sf.lmax[c]-1));
      // azimuth of the shortening direction, degrees clockwise from map north (0-180)
      G.shorteningAzimuthDeg.push(isNaN(sf.lmin[c])?null:Math.round(((90-sf.dir[c]*180/Math.PI)%180+180)%180));G.rotationDeg.push(isNaN(sf.lmin[c])?null:r1(sf.rot[c]*180/Math.PI));
      G.strainRatePerMyr.push(a?r3(sheet.eps[c]):null);G.mantlePullKm.push(a?r2(sheet.dyn[c]):null);G.denudationKmPerMyr.push(a?r3(deno[c]):null);G.arcAddedKm.push(a?r2(arcWin[c]):null);G.stacks.push(stacksAt(c));}
    return{meta:{generator:meta.generator,seed:meta.seed,presentMyr:meta.presentMyr,tile:id,units:meta.units,frame:meta.frame,receivers:meta.receivers,
        stacks:"per coarse cell, the L0 parcels present (largest first, up to "+o.parcels+"): share of the cell's crust, thickness, and units bottom to top as [rock, formed Ma, thickness km]",
        drainage:"drainage: 0 reaches the sea, 1 ends inland (lake with no outlet or a river that runs dry), 2 areic (no water reaches a river); lake: 1 open, 2 no outlet, 3 dry floor of a shrunk lake, 4 where a river runs dry; divide: 1 a contingent edge with another tile (L2 may move it), 2 a firm one"},
      tile:T,imprints:T.imprints.map(k=>imprints[k]),lifetimes,
      fine:Object.assign({cellKm:h,i0,j0,nx:tx,ny:ty,x0Km:B.x0+i0*h,y0Km:B.y0+j0*h},F),coarse:Object.assign({cellKm:H,x0Km:sheet.x0+si0*H,y0Km:sheet.y0+sj0*H},G)};}
  return{meta,lifetimes,imprints,tiles:outTiles,tileOf,divide,responseMyr:tau,tile,summary:()=>({meta,lifetimes,imprints,tiles:outTiles})};
}
return {Replay,Sheet,Projection,landmassAt,basins,handoff,route,streamPower,planarGrid,scratch,Heap,expm2,dynTopo,dynAt,DYN};
})();
// ===== L1-END =====
