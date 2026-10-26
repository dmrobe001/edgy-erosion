// ===== L1-START =====
// L1: one continent, refined from L0 and evolved forward while L0 keeps running, together with whatever crust
// will interact with it during the window. L0 supplies tectonic thickening (arcs, collisions, hotspots, collapse);
// L1 replaces L0's surface processes with stream-power river erosion, derives elevation from crust thickness
// through flexural isostasy, and resolves where its own pieces of crust meet, at its own resolution.
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

// ---------- refined cube-sphere over a set of L0 cells ----------
// Each L0 cell c in the domain gets a slot s; its r x r children have local index s*r*r + jj*r + ii.
// Equal-angle subdivision nests exactly, so every fine cell has one parent and positions are analytic.
// With jitter > 0, each cell's node (used for distances, hence slopes) is displaced within the cell by up to
// jitter x cell width, by a fixed hash of its global index, so the 8 neighbour directions differ from cell to
// cell and rivers are not locked to the lattice's 0 and 45 degree directions. pos stays at the cell centre.
const D8=[[1,0],[-1,0],[0,1],[0,-1],[1,1],[-1,1],[1,-1],[-1,-1]]; // E W N S NE NW SE SW
function faceDir(f,a,b,t){let x,y,z;switch(f){case 0:x=1;y=a;z=b;break;case 1:x=-1;y=a;z=b;break;case 2:x=a;y=1;z=b;break;case 3:x=a;y=-1;z=b;break;case 4:x=a;y=b;z=1;break;default:x=a;y=b;z=-1;}
  const L=Math.hypot(x,y,z);t[0]=x/L;t[1]=y/L;t[2]=z/L;return t;}
function buildGrid(n,slots,slotOf,r,jitter){
  const nS=slots.length,rr=r*r,nC=nS*rr,NR=n*r,nn=n*n,q=Math.PI/4,lk=makeLookup(NR),t=[0,0,0];
  const pos=new Float32Array(3*nC),area=new Float32Array(nC),nbr=new Int32Array(8*nC),dist=new Float32Array(8*nC),gid=new Int32Array(nC);
  const dA=(Math.PI/2/NR)**2*R_KM*R_KM,ang=v=>Math.tan(((v+.5)/NR*2-1)*q);
  const localOf=g=>{const f=(g/(NR*NR))|0,rem=g-f*NR*NR,J=(rem/NR)|0,I=rem-J*NR,s=slotOf[f*nn+((J/r)|0)*n+((I/r)|0)];return s<0?-1:s*rr+(J%r)*r+(I%r);};
  for(let s=0;s<nS;s++){const c=slots[s],f=(c/nn)|0,j=((c-f*nn)/n)|0,i=c-f*nn-j*n;
    for(let jj=0;jj<r;jj++)for(let ii=0;ii<r;ii++){const l=s*rr+jj*r+ii,I=i*r+ii,J=j*r+jj,a=ang(I),b=ang(J);gid[l]=f*NR*NR+J*NR+I;
      faceDir(f,a,b,t);pos[3*l]=t[0];pos[3*l+1]=t[1];pos[3*l+2]=t[2];area[l]=(1+a*a)*(1+b*b)/Math.pow(1+a*a+b*b,1.5)*dA;
      for(let k=0;k<8;k++){const I2=I+D8[k][0],J2=J+D8[k][1];let m;
        if(I2>=0&&J2>=0&&I2<NR&&J2<NR){const s2=slotOf[f*nn+((J2/r)|0)*n+((I2/r)|0)];m=s2<0?-1:s2*rr+(J2%r)*r+(I2%r);}
        else{faceDir(f,ang(I2),ang(J2),t);m=localOf(lk(t[0],t[1],t[2]));if(m===l)m=-1;}
        nbr[8*l+k]=m;}}}
  let node=pos;
  if(jitter>0){node=new Float32Array(3*nC);
    for(let l=0;l<nC;l++){const g=gid[l],f=(g/(NR*NR))|0,rem=g-f*NR*NR,J=(rem/NR)|0,I=rem-J*NR;
      let h=Math.imul(g^(g>>>16),0x7feb352d);h=Math.imul(h^(h>>>15),0x846ca68b);h^=h>>>16;const u=(h&0xffff)/65536-0.5,v=(h>>>16)/65536-0.5;
      faceDir(f,ang(I+2*jitter*u),ang(J+2*jitter*v),t);node[3*l]=t[0];node[3*l+1]=t[1];node[3*l+2]=t[2];}}
  for(let l=0;l<nC;l++)for(let k=0;k<8;k++){const m=nbr[8*l+k];if(m<0)continue;
    dist[8*l+k]=Math.hypot(node[3*l]-node[3*m],node[3*l+1]-node[3*m+1],node[3*l+2]-node[3*m+2])*R_KM;}
  return{n,r,rr,nS,nC,NR,pos,node,area,nbr,dist,gid,lookup:lk,localOf,jitter};
}

// 3x3 rotations, row-major
function mat(A,B){const o=new Array(9);for(let r=0;r<3;r++)for(let c=0;c<3;c++)o[r*3+c]=A[r*3]*B[c]+A[r*3+1]*B[3+c]+A[r*3+2]*B[6+c];return o;}
function matT(A){return[A[0],A[3],A[6],A[1],A[4],A[7],A[2],A[5],A[8]];}
// Pieces (one per plate at t0, each in that plate's frame) are concatenated into one graph; cell indices of
// piece p are offset by its first slot x r^2. There are no static links between pieces.
function concatGrids(list){if(list.length===1)return Object.assign({},list[0],{parts:list,offsets:[0]});const g0=list[0];let nC=0,nS=0;for(const g of list){nC+=g.nC;nS+=g.nS;}
  const out={n:g0.n,r:g0.r,rr:g0.rr,nS,nC,NR:g0.NR,pos:new Float32Array(3*nC),node:new Float32Array(3*nC),area:new Float32Array(nC),nbr:new Int32Array(8*nC),dist:new Float32Array(8*nC),gid:new Int32Array(nC),parts:list,offsets:[],jitter:g0.jitter};
  let o=0;for(const g of list){out.offsets.push(o);out.pos.set(g.pos,3*o);out.node.set(g.node,3*o);out.area.set(g.area,o);out.gid.set(g.gid,o);out.dist.set(g.dist,8*o);
    for(let q=0;q<8*g.nC;q++)out.nbr[8*o+q]=g.nbr[q]<0?-1:g.nbr[q]+o;o+=g.nC;}
  return out;}

// ---------- continent selection ----------
// The connected continental crust of the clicked cell's plate, in that plate's body frame, plus every body
// cell of the same plate within `margin` L0 cells of it (the ocean margin is base level for the rivers).
function select(world,cell,margin){return selectBody(world,world.owner[cell],world.bodyIdx[cell],margin);}
function selectBody(world,plate,k0,margin){
  const P=world.plates[plate],N=world.N,nb=world.cube.nb;
  if(!P||!P.alive||!P.K[k0]||P.T[k0]<=0)return null;
  const d=new Int8Array(N).fill(-1),q=[k0];d[k0]=0;
  for(let h=0;h<q.length;h++){const a=q[h];for(let k=0;k<4;k++){const b=nb[4*a+k];if(d[b]<0&&P.K[b]&&P.T[b]>0){d[b]=0;q.push(b);}}}
  const nCont=q.length;
  for(let h=0;h<q.length;h++){const a=q[h];if(d[a]>=margin)continue;for(let k=0;k<4;k++){const b=nb[4*a+k];if(d[b]<0&&P.T[b]>0){d[b]=d[a]+1;q.push(b);}}}
  const slots=Int32Array.from(q).sort();
  return{plate:P.id,k0,slots,nCont,dist:d};
}

// ---------- which crust matters: a lookahead and a backward pass ----------
// L0 tags every crust parcel present at t0 with an identity (world.trackParcels), then runs ahead through the
// window from a snapshot, recording for every step which parcel sits on each fixed cell and whether it is
// continental. Crust is relevant at a step if it lies within reachKm of the target (the picked continent's
// parcels, wherever they are, including every fragment it breaks into), through continental crust contiguous at
// that step. So a landmass that mashes against the target matters while it is attached and within reach, and
// drops out once it has separated. (Chaining relevance back through time, "within reach of whatever will be
// relevant next step", compounds the reach every step and swallows whole landmasses.) Each parcel is later
// brought in a warm-up before its first relevant time. Then the snapshot is restored, so the real run replays
// the same future, identities included.
class Lookahead{
  constructor(world,cell,window,o){
    o=Object.assign({reachKm:1200,margin:3},o);this.o=o;
    const sel=select(world,cell,o.margin);if(!sel)throw new Error("not a continental cell");
    world.trackParcels();
    const N=world.N,P=world.plates[sel.plate];this.world=world;this.N=N;this.t0=world.t;this.window=window;this.ids=world.plates.length*N;
    this.target=new Uint8Array(this.ids);for(const k of sel.slots)if(sel.dist[k]===0&&P.pid[k]>=0)this.target[P.pid[k]]=1;
    this.snap=world.snapshot();this.frames=[];this.record();this.done=false;}
  record(){const w=this.world,N=w.N,id=new Int32Array(N),K=new Uint8Array(N);
    for(let i=0;i<N;i++){const P=w.plates[w.owner[i]];id[i]=P.pid?P.pid[w.bodyIdx[i]]:-1;K[i]=w.wK[i];}
    this.frames.push({t:w.t-this.t0,id,K});}
  // run up to maxMs of L0 steps; returns true when the window is covered
  advance(maxMs){const w=this.world,t1=Date.now();
    while(!this.done&&Date.now()-t1<maxMs){w.step();this.record();if(w.t-this.t0>=this.window-1e-9)this.done=true;}
    return this.done;}
  // restore L0 to t0; return each relevant parcel's identity with the first and last time (Myr after t0) it matters
  finish(){const w=this.world;w.restore(this.snap);
    const N=this.N,nb=w.cube.nb,depth=Math.max(1,Math.round(this.o.reachKm/w.cellKm)),ids=this.ids;
    const relNow=new Uint8Array(ids),first=new Float32Array(ids).fill(NaN),last=new Float32Array(ids).fill(NaN);
    const dist=new Int16Array(N),q=new Int32Array(N),touched=[],present=new Uint8Array(ids),presList=[],end=this.frames[this.frames.length-1].t;
    for(let s=this.frames.length-1;s>=0;s--){const {t,id,K}=this.frames[s];let qh=0,qt=0;dist.fill(-1);
      for(let i=0;i<N;i++){const p=id[i];if(K[i]&&p>=0&&this.target[p]){dist[i]=0;q[qt++]=i;}}
      while(qh<qt){const a=q[qh++];if(dist[a]>=depth)continue;for(let k=0;k<4;k++){const b=nb[4*a+k];if(dist[b]<0&&K[b]){dist[b]=dist[a]+1;q[qt++]=b;}}}
      for(const p of touched)relNow[p]=0;touched.length=0;
      // a parcel last relevant now that is absent from L0 at the next step was removed by L0 (consumed or overwritten),
      // not left behind: keep it to the end, and let L1 resolve its crust itself
      for(let h=0;h<qt;h++){const p=id[q[h]];if(p<0||relNow[p])continue;relNow[p]=1;touched.push(p);if(isNaN(last[p]))last[p]=s<this.frames.length-1&&!present[p]?end:t;first[p]=t;}
      for(const p of presList)present[p]=0;presList.length=0;for(let i=0;i<N;i++){const p=id[i];if(p>=0&&!present[p]){present[p]=1;presList.push(p);}}}
    const pids=[];for(let p=0;p<ids;p++)if(!isNaN(first[p])){pids.push(p);if(last[p]>=end-1e-3)last[p]=Infinity;} // relevant to the end: never dropped (Float32 times, so a loose tolerance)
    return{N,t0:this.t0,window:this.window,pids:Int32Array.from(pids),first:Float32Array.from(pids,p=>first[p]),last:Float32Array.from(pids,p=>last[p]),
      target:Uint8Array.from(pids,p=>this.target[p]),reachKm:this.o.reachKm};}
}

// undirected edges (a[i], b[i], length d[i]) as CSR adjacency {off, adj, dist, count}
function csr(n,a,b,d){const off=new Int32Array(n+1);for(let q=0;q<a.length;q++){off[a[q]+1]++;off[b[q]+1]++;}for(let l=0;l<n;l++)off[l+1]+=off[l];
  const adj=new Int32Array(off[n]),dist=d?new Float32Array(off[n]):null,fill=off.slice(0,n);
  for(let q=0;q<a.length;q++){adj[fill[a[q]]]=b[q];if(d)dist[fill[a[q]]]=d[q];fill[a[q]]++;adj[fill[b[q]]]=a[q];if(d)dist[fill[b[q]]]=d[q];fill[b[q]]++;}
  return{off,adj,dist,count:a.length};}

// ---------- a spatial index in the fixed (world) frame ----------
// Fine cells of every piece, at their current world positions, hashed by the global fine-grid cell they fall in.
class CellIndex{
  constructor(NR,cap){this.NR=NR;this.lk=makeLookup(NR);let P=1;while(P<2*cap)P<<=1;this.mask=P-1;this.key=new Int32Array(P);this.head=new Int32Array(P);this.next=new Int32Array(cap);this.g=new Int32Array(cap).fill(-1);this.around=new Int32Array(9);}
  clear(){this.key.fill(-1);this.g.fill(-1);}
  insert(l,x,y,z){const g=this.lk(x,y,z);let h=(Math.imul(g,0x9E3779B1)>>>0)&this.mask;
    while(this.key[h]!==-1&&this.key[h]!==g)h=(h+1)&this.mask;if(this.key[h]===-1){this.key[h]=g;this.head[h]=-1;}
    this.next[l]=this.head[h];this.head[h]=l;this.g[l]=g;return g;}
  first(g){let h=(Math.imul(g,0x9E3779B1)>>>0)&this.mask;for(;;){const k=this.key[h];if(k===-1)return -1;if(k===g)return this.head[h];h=(h+1)&this.mask;}}
  // the global cell g and its in-face neighbours (across a cube edge they are skipped); returns how many
  neighbours(g){const NR=this.NR,F=NR*NR,f=(g/F)|0,rem=g-f*F,J=(rem/NR)|0,I=rem-J*NR;let n=0;
    for(let dj=-1;dj<=1;dj++){const J2=J+dj;if(J2<0||J2>=NR)continue;for(let di=-1;di<=1;di++){const I2=I+di;if(I2<0||I2>=NR)continue;this.around[n++]=f*F+J2*NR+I2;}}
    return n;}
}

class Region{
  // o.cell: the picked fixed-grid cell. o.relevance: from Lookahead.finish(), which also switched on L0's parcel
  // identities. Every parcel that is ever relevant gets a slot; slots are grouped into pieces by the plate that
  // held them at t0 (each piece in that plate's frame), and each slot is active from (first relevant - warm-up)
  // to its last relevant time. Ocean margins around relevant crust share the activity of what they border.
  constructor(world,o){
    o=Object.assign({r:8,margin:3,K:0.1,m:0.5,flexKm:40,flowKm2PerMyr:300,strengthKm:3,jitter:0.4,window:100,warmupMyr:15},o);this.o=o;this.world=world;
    const rel=o.relevance;if(!rel||!world.tracking)throw new Error("Region needs a lookahead (Lookahead.finish())");
    const n=world.n,N=world.N,cube=world.cube,nb=cube.nb,t0=world.t,r=o.r,rr=r*r;
    this.t0=t0;this.t=0;this.steps=0;
    const sel=select(world,o.cell,o.margin);if(!sel)throw new Error("not a continental cell");
    // group relevant parcels by their plate at t0; the picked continent's plate comes first
    const groups=new Map();const gget=p=>{let g=groups.get(p);if(!g)groups.set(p,g=new Map());return g;};gget(sel.plate);
    for(let q=0;q<rel.pids.length;q++){const pid=rel.pids[q],p0=(pid/N)|0,k0=pid-p0*N;
      gget(p0).set(k0,{spawn:Math.max(0,rel.first[q]-o.warmupMyr),drop:rel.last[q],target:rel.target[q]});}
    this.pieces=[];let nS=0;const recs=[];
    for(const [p0,g] of groups){const P=world.plates[p0];if(!P||!P.alive||!g.size)continue;
      // ocean margin: m rings of the plate's own crust around relevant parcels, active whenever a neighbour is
      for(let ring=0;ring<o.margin;ring++){const add=[];
        for(const [k,v] of g)for(let d=0;d<4;d++){const b=nb[4*k+d];if(P.T[b]<=0)continue;const u=g.get(b);
          if(!u)add.push([b,v]);else if(u.margin&&(v.spawn<u.spawn||v.drop>u.drop)){u.spawn=Math.min(u.spawn,v.spawn);u.drop=Math.max(u.drop,v.drop);}}
        for(const [b,v] of add){const u=g.get(b);if(!u)g.set(b,{spawn:v.spawn,drop:v.drop,target:0,margin:1});else{u.spawn=Math.min(u.spawn,v.spawn);u.drop=Math.max(u.drop,v.drop);}}}
      const slots=Int32Array.from([...g.keys()]).sort(),slotOf=new Int32Array(N).fill(-1);for(let s=0;s<slots.length;s++)slotOf[slots[s]]=s;
      this.pieces.push({plate0:p0,s0:nS,nS:slots.length,slots,slotOf,frames:new Map([[p0,[1,0,0,0,1,0,0,0,1]]])});
      for(const k of slots)recs.push(g.get(k));nS+=slots.length;}
    this.pieceOfPlate0=new Int32Array(world.plates.length).fill(-1);this.pieces.forEach((pc,i)=>this.pieceOfPlate0[pc.plate0]=i);
    const slots=this.slots=new Int32Array(nS),pieceOf=this.pieceOf=new Uint8Array(nS);
    this.slotSpawn=new Float32Array(nS);this.slotDrop=new Float32Array(nS);this.slotTarget=new Uint8Array(nS);this.slotMargin=new Uint8Array(nS);
    this.pieces.forEach((pc,pi)=>{for(let s=0;s<pc.nS;s++){const S=pc.s0+s,v=recs[S];slots[S]=pc.slots[s];pieceOf[S]=pi;this.slotSpawn[S]=v.spawn;this.slotDrop[S]=v.drop;this.slotTarget[S]=v.target;this.slotMargin[S]=v.margin?1:0;}});
    this.slotCont=this.slotTarget; // the picked continent's own parcels (for the globe overlay and the map centre)
    const so=s=>{const pc=this.pieces[pieceOf[s]];return c=>{const v=pc.slotOf[c];return v<0?-1:v+pc.s0;};};
    this.slotPlate=new Int32Array(nS);this.slotBody=Int32Array.from(slots);this.slotFound=new Uint8Array(nS);this.slotActive=new Uint8Array(nS);
    for(let s=0;s<nS;s++)this.slotPlate[s]=this.pieces[pieceOf[s]].plate0;
    this.slotNb=new Int32Array(4*nS);this.quad=new Int32Array(12*nS);
    for(let s=0;s<nS;s++){const f=so(s),c=slots[s];for(let k=0;k<4;k++)this.slotNb[4*s+k]=f(nb[4*c+k]);
      // bilinear stencils: for each slot and quadrant, the three other parents (as slots, -1 outside the domain)
      for(let qd=0;qd<4;qd++){const cB=nb[4*c+(qd&1?1:0)],cC=nb[4*c+(qd&2?3:2)],cD=nb[4*cB+(qd&2?3:2)];
        this.quad[12*s+3*qd]=f(cB);this.quad[12*s+3*qd+1]=f(cC);this.quad[12*s+3*qd+2]=f(cD);}}
    const G=this.G=concatGrids(this.pieces.map(pc=>buildGrid(n,pc.slots,pc.slotOf,r,o.jitter))),nC=G.nC;
    this.subQ=new Uint8Array(rr);this.subW=new Float32Array(2*rr);
    for(let jj=0;jj<r;jj++)for(let ii=0;ii<r;ii++){const tx=(ii+.5)/r-.5,ty=(jj+.5)/r-.5,sb=jj*r+ii;this.subQ[sb]=(tx<0?1:0)+(ty<0?2:0);this.subW[2*sb]=Math.abs(tx);this.subW[2*sb+1]=Math.abs(ty);}
    // flexure and crustal flow run on a coarser nested grid (b x b fine cells per mid cell), then are interpolated back
    const sigma=1.18*o.flexKm,cellKm=world.cellKm;let rf=1;for(let d=r;d>=1;d--)if(r%d===0&&cellKm/d>=sigma/2.5){rf=d;break;}
    this.rf=rf;this.Gm=concatGrids(this.pieces.map(pc=>buildGrid(n,pc.slots,pc.slotOf,rf,0)));this.buildMidMaps();this.flexPasses=Math.max(1,Math.round(sigma*sigma/(0.4*(cellKm/rf)**2)));
    this.T=new Float64Array(nC);this.e=new Float64Array(nC);this.cont=new Uint8Array(nC);this.alive=new Uint8Array(nC).fill(1);this.live=new Uint8Array(nC);this.eOcean=new Float32Array(nC);
    this.T0=new Float64Array(nC);this.e0=new Float32Array(nC);this.eroded=new Float32Array(nC);this.uplifted=new Float32Array(nC);this.outlet=new Uint8Array(nC);this.root=new Float64Array(nC);
    this.W3=new Float32Array(3*nC);this.Wm=new Float32Array(3*this.Gm.nC);
    this.S=scratch(nC);this.S.seaEdge=new Uint8Array(nC);this.rate=new Float64Array(nS);this.rec0=new Float32Array(N);this.precip=null;this.K=o.K;
    this.index=new CellIndex(G.NR,nC);this.midIndex=new CellIndex(this.Gm.NR,this.Gm.nC);
    this.nz=makeNoise3(hashStr(String(world.opts.seed)+"|L1")|0);
    this.budget={initial:0,activated:0,dropped:0,tectonicApplied:0,tectonicL0:0,unappliedOcean:0,eroded:0,capped:0,consumed:0,consumedOutside:0,refilled:0,flowResidual:0};
    this.events=[];this.history=[];
    this.updateFrames();this.locate();this.activate(0);this.budget.initial=this.budget.activated;this.budget.activated=0;this.locate();
    this.place();this.flexure();this.classify();this.seaEdges();route(this.effGraph(),this.e,this.live,this.outlet,this.S,this.precip,this.X);this.stats(0,0);
    // map orientation: the picked continent's centroid, and world north at t0, in the first piece's frame
    let cx=0,cy=0,cz=0;for(let s=0;s<this.pieces[0].nS;s++){if(!this.slotTarget[s])continue;const c=slots[s];cx+=cube.pos[3*c];cy+=cube.pos[3*c+1];cz+=cube.pos[3*c+2];}
    const cl=Math.hypot(cx,cy,cz)||1,Q=world.plates[this.pieces[0].plate0].Q;this.center=[cx/cl,cy/cl,cz/cl];this.northB=[Q[6],Q[7],Q[8]];
  }
  // bilinear interpolation of a per-slot field onto fine cells (parents outside the domain, or masked out, are skipped);
  // with `cells`, only those fine cells are written
  interp(field,out,mask,cells){const {G,quad,subQ,subW}=this,rr=G.rr,nC=G.nC;out=out||new Float64Array(nC);
    const one=l=>{const s=(l/rr)|0,sb=l-s*rr,qd=subQ[sb],wx=subW[2*sb],wy=subW[2*sb+1];
      const sB=quad[12*s+3*qd],sC=quad[12*s+3*qd+1],sD=quad[12*s+3*qd+2];
      let v=0,w=0,ww;
      ww=(1-wx)*(1-wy);if(!mask||mask[s]){v+=ww*field[s];w+=ww;}
      ww=wx*(1-wy);if(sB>=0&&(!mask||mask[sB])){v+=ww*field[sB];w+=ww;}
      ww=(1-wx)*wy;if(sC>=0&&(!mask||mask[sC])){v+=ww*field[sC];w+=ww;}
      ww=wx*wy;if(sD>=0&&(!mask||mask[sD])){v+=ww*field[sD];w+=ww;}
      out[l]=w>0?v/w:field[s];};
    if(cells)for(const l of cells)one(l);else for(let l=0;l<nC;l++)one(l);
    return out;}
  volume(){const {T,cont,live,G}=this;let v=0;for(let l=0;l<G.nC;l++)if(cont[l]&&live[l])v+=T[l]*G.area[l];return v;}
  // ---- the coarser nested grid used for long-wavelength processes (flexure, lower-crustal flow) ----
  // Per fine cell: its mid cell, and a bilinear stencil (three neighbouring mid cells + weights) from a table.
  buildMidMaps(){const G=this.G,Gm=this.Gm,nC=G.nC,r=G.r,rr=G.rr,rf=Gm.r,b=r/rf,nb=Gm.nbr,mC=Gm.nC;
    this.mid=new Int32Array(nC);this.midQ=new Uint8Array(nC);this.midS=new Uint16Array(nC);
    this.mq=new Int32Array(12*mC);for(let m=0;m<mC;m++)for(let qd=0;qd<4;qd++){const sx=qd&1,sy=qd&2;
      this.mq[12*m+3*qd]=nb[8*m+(sx?1:0)];this.mq[12*m+3*qd+1]=nb[8*m+(sy?3:2)];this.mq[12*m+3*qd+2]=nb[8*m+(sx?(sy?7:5):(sy?6:4))];}
    this.mw=new Float32Array(2*b*b);
    for(let y=0;y<b;y++)for(let x=0;x<b;x++){const tx=(x+.5)/b-.5,ty=(y+.5)/b-.5;this.mw[2*(y*b+x)]=Math.abs(tx);this.mw[2*(y*b+x)+1]=Math.abs(ty);}
    for(let l=0;l<nC;l++){const s=(l/rr)|0,sb=l-s*rr,jj=(sb/r)|0,ii=sb-jj*r,x=ii%b,y=jj%b;
      this.mid[l]=s*rf*rf+((jj/b)|0)*rf+((ii/b)|0);this.midQ[l]=((x+.5)/b-.5<0?1:0)+((y+.5)/b-.5<0?2:0);this.midS[l]=y*b+x;}}
  // area-weighted mean of fine values src over cells with use[l] set; W gets the contributing area
  toMid(src,use,X,W){const {G,mid}=this,nC=G.nC,area=G.area;X.fill(0);W.fill(0);
    for(let l=0;l<nC;l++){if(!use[l])continue;const m=mid[l],a=area[l];X[m]+=a*src[l];W[m]+=a;}
    for(let m=0;m<X.length;m++)if(W[m]>0)X[m]/=W[m];}
  // bilinear interpolation of a mid-grid field back to every fine cell with use[l] set (mid cells with W = 0, or on
  // another plate, are skipped)
  fromMid(A,W,use,out){const {G,mid,midQ,midS,mq,mw}=this,nC=G.nC,mp=this.midPlate;
    for(let l=0;l<nC;l++){if(!use[l])continue;const m=mid[l],o=12*m+3*midQ[l],wx=mw[2*midS[l]],wy=mw[2*midS[l]+1],mB=mq[o],mC=mq[o+1],mD=mq[o+2],p=mp[m];
      let v=0,w=0,ww=(1-wx)*(1-wy);if(W[m]>0){v+=ww*A[m];w+=ww;}
      ww=wx*(1-wy);if(mB>=0&&W[mB]>0&&mp[mB]===p){v+=ww*A[mB];w+=ww;}
      ww=(1-wx)*wy;if(mC>=0&&W[mC]>0&&mp[mC]===p){v+=ww*A[mC];w+=ww;}
      ww=wx*wy;if(mD>=0&&W[mD]>0&&mp[mD]===p){v+=ww*A[mD];w+=ww;}
      out[l]=w>0?v/w:0;}
    return out;}
  fineBuf(name){return this[name]||(this[name]=new Float64Array(this.G.nC));}
  midBuf(name){return this[name]||(this[name]=new Float64Array(this.Gm.nC));}
  // Flexural isostasy, approximated, applied to changes since each cell became active (the inherited state is
  // taken as in equilibrium): e = e0 + dT - (1 - ISO) * F[dT], where F is a Gaussian of width 1.18 x flexural
  // length (the half-power point of the plate's response 1/(1+(k*alpha)^4)). Long wavelengths follow L0's Airy
  // isostasy; loads narrower than the flexural length are held up by the plate's strength, up to a limit.
  // The smoothing crosses contacts between pieces, so a load on one side bends the other.
  flexure(){
    const {G,T,T0,e0,cont,live,eOcean,root,e}=this,nC=G.nC,X=this.midBuf("_X"),W=this.midBuf("_W"),Y=this.midBuf("_Y"),mC=X.length;
    const dT=this.fineBuf("_dT");for(let l=0;l<nC;l++)dT[l]=cont[l]?T[l]-T0[l]:0;
    this.toMid(dT,live,X,W);
    let A=X,B=Y;
    const {off,adj}=this.Am;
    for(let p=0;p<this.flexPasses;p++){for(let m=0;m<mC;m++){if(!(W[m]>0)){B[m]=0;continue;}let v=A[m],w=1;
        for(let q=off[m];q<off[m+1];q++){const j=adj[q];if(W[j]>0){v+=A[j];w++;}}B[m]=v/w;}[A,B]=[B,A];}
    // the plate yields rather than hold more than strengthKm of topography away from local (Airy) balance
    const lim=this.o.strengthKm;
    const F=this.fromMid(A,W,live,this.fineBuf("_F"));
    for(let l=0;l<nC;l++){if(!live[l])continue;let rt=ROOT*F[l];
      const air=ROOT*dT[l];if(rt<air-lim)rt=air-lim;else if(rt>air+lim)rt=air+lim;root[l]=rt;e[l]=cont[l]?e0[l]+dT[l]-rt:eOcean[l];}
  }
  // Lower-crustal flow: thickened crust spreads sideways, on the mid grid only (so it relaxes collision piles
  // without filling valleys), across contacts too. It diffuses the thickness change since activation, not the
  // thickness itself: inherited contrasts (thin rifted margins) are equilibrium, as in flexure(); diffusing total
  // thickness pumped crust out to the margins and raised a ring of islands. Channel flux grows with the cube of
  // channel thickness. Fluxes are symmetric; the residue of interpolating back to fine cells is budget.flowResidual.
  crustFlow(dt){
    const kappa=this.o.flowKm2PerMyr;if(!(kappa>0))return;
    const {G,T,T0,cont,live}=this,nC=G.nC,X=this.midBuf("_FX"),W=this.midBuf("_FW"),V=this.midBuf("_FV"),Tm=this.midBuf("_FT"),mC=X.length,dx=this.world.cellKm/this.rf;
    const use=this._cu||(this._cu=new Uint8Array(nC));for(let l=0;l<nC;l++)use[l]=live[l]&&cont[l]?1:0;
    const dT=this.fineBuf("_FdT");for(let l=0;l<nC;l++)dT[l]=T[l]-T0[l];
    this.toMid(T,use,Tm,W);this.toMid(dT,use,X,W);const X0=this.midBuf("_FX0");X0.set(X);
    const cap=8*mC,eI=this._eI&&this._eI.length>=cap?this._eI:(this._eI=new Int32Array(cap)),eJ=this._eJ&&this._eJ.length>=cap?this._eJ:(this._eJ=new Int32Array(cap)),eC=this._eC&&this._eC.length>=cap?this._eC:(this._eC=new Float64Array(cap));
    let ne=0,cmax=0;
    const {off,adj}=this.Am;
    for(let m=0;m<mC;m++){if(!(W[m]>0))continue;for(let q=off[m];q<off[m+1];q++){const j=adj[q];if(j<=m||!(W[j]>0)||ne>=cap)continue;
      const t=(Tm[m]+Tm[j])/(2*C.CONT_BASE),c=kappa*t*t*t;eI[ne]=m;eJ[ne]=j;eC[ne]=c*Math.min(W[m],W[j]);ne++;if(c>cmax)cmax=c;}}
    const D=cmax*dt/(dx*dx),nsub=Math.max(1,Math.ceil(D/0.2)),h=dt/nsub/(dx*dx);
    for(let it=0;it<nsub;it++){V.fill(0);
      for(let q=0;q<ne;q++){const m=eI[q],j=eJ[q],f=h*eC[q]*(X[j]-X[m]);V[m]+=f;V[j]-=f;}
      for(let m=0;m<mC;m++)if(W[m]>0)X[m]+=V[m]/W[m];}
    for(let m=0;m<mC;m++)X[m]-=X0[m];
    const add=this.fromMid(X,W,use,this.fineBuf("_F"));let res=0;for(let l=0;l<nC;l++){if(!use[l])continue;T[l]+=add[l];res+=add[l]*G.area[l];}
    this.budget.flowResidual+=res;
  }
  // Outlets are sea: ocean cells and drowned continental cells.
  classify(){const {G,cont,live,e,outlet}=this,nC=G.nC;for(let l=0;l<nC;l++)outlet[l]=live[l]&&(!cont[l]||e[l]<=0)?1:0;}

  // ---- following parcels and their frames ----
  // Each piece keeps a frame for every plate it has parcels on: for parcels on plate X, world = Q(X) M b, with b
  // the parcel's position in the piece's own (t0) body frame. A split copies the parent's M to the new plate;
  // a weld of L into W composes M with the fixed rotation between their frames.
  updateFrames(){const plates=this.world.plates,pcs=this.pieces;
    for(let id=this.framesSeen===undefined?plates.length:this.framesSeen;id<plates.length;id++){const C=plates[id];if(!C.alive||C.parent===undefined)continue;
      for(const pc of pcs){const M=pc.frames.get(C.parent);if(M&&!pc.frames.has(id))pc.frames.set(id,M.slice());}}
    this.framesSeen=plates.length;
    for(const pc of pcs)for(const [p,M] of [...pc.frames]){let P=plates[p],Mn=M;
      for(let g=0;g<16&&!P.alive&&P.mergedInto!==undefined;g++){const W=plates[P.mergedInto];Mn=mat(matT(W.Q),mat(P.Q,Mn));P=W;}
      if(P.id!==p){pc.frames.delete(p);if(P.alive&&!pc.frames.has(P.id))pc.frames.set(P.id,Mn);}}}
  // Where is each parcel now? Scan the plates for its identity. A parcel L0 no longer has (consumed, or its cell
  // overwritten) keeps riding the plate it was last on: in L1 its crust lives on until L1 itself resolves it.
  locate(){const w=this.world,N=w.N,plates=w.plates,{slotPlate,slotBody,slotFound,pieces,pieceOfPlate0}=this,np0=pieceOfPlate0.length;slotFound.fill(0);
    for(const P of plates){if(!P.alive||!P.pid)continue;const pid=P.pid,T=P.T;
      for(let k=0;k<N;k++){const id=pid[k];if(id<0||T[k]<=0)continue;const p0=(id/N)|0;if(p0>=np0)continue;const pi=pieceOfPlate0[p0];if(pi<0)continue;
        const pc=pieces[pi],s0=pc.slotOf[id-p0*N];if(s0<0)continue;const s=s0+pc.s0;if(slotFound[s]||!pc.frames.has(P.id))continue;
        slotFound[s]=1;slotPlate[s]=P.id;slotBody[s]=k;}}
    for(let s=0;s<slotPlate.length;s++){if(slotFound[s])continue;let P=plates[slotPlate[s]];
      for(let g=0;g<16&&!P.alive&&P.mergedInto!==undefined;g++)P=plates[P.mergedInto];slotPlate[s]=P.id;}
    // which rigid parts exist: (piece, plate) with active parcels
    const on=pieces.map(()=>new Map());for(let s=0;s<slotPlate.length;s++)if(this.slotActive[s]){const m=on[this.pieceOf[s]],p=slotPlate[s];m.set(p,(m.get(p)||0)+1);}
    this.piecePlates=on;}
  // L0's current elevation of each located parcel, smoothed (three passes over located neighbours)
  slotElevation(){const w=this.world,nS=this.slots.length,raw=new Float64Array(nS),kf=new Float64Array(nS),ok=new Uint8Array(nS);
    for(let s=0;s<nS;s++){if(!this.slotFound[s])continue;const P=w.plates[this.slotPlate[s]],k=this.slotBody[s];raw[s]=elevOf(P.K[k],P.T[k],P.B[k],P.S[k],w.t);kf[s]=P.K[k];ok[s]=1;}
    let es=raw,tmp=new Float64Array(nS);
    for(let pass=0;pass<3;pass++){for(let s=0;s<nS;s++){if(!ok[s]){tmp[s]=es[s];continue;}let v=es[s],c=1;for(let k=0;k<4;k++){const u=this.slotNb[4*s+k];if(u>=0&&ok[u]){v+=es[u];c++;}}tmp[s]=v/c;}[es,tmp]=[tmp,es];}
    return{es,kf,ok};}
  // Bring in slots whose time has come (and that L0 still has), starting them from L0's state now; drop slots
  // whose relevance has ended. Each newly active cell's inherited state counts as equilibrium (T0, e0).
  activate(t){const nS=this.slots.length,rr=this.G.rr,{slotActive,slotSpawn,slotDrop,slotFound}=this,G=this.G,b=this.budget;
    const on=[],off=[];
    for(let s=0;s<nS;s++){if(!slotActive[s]&&slotSpawn[s]<=t+1e-9&&slotDrop[s]>=t&&slotFound[s])on.push(s);else if(slotActive[s]&&slotDrop[s]<t-1e-9)off.push(s);}
    for(const s of off){slotActive[s]=0;for(let l=s*rr;l<(s+1)*rr;l++){if(this.alive[l]&&this.cont[l])b.dropped+=this.T[l]*G.area[l];this.live[l]=0;}}
    if(off.length)this.events.push({t,type:"drop",slots:off.length});
    if(!on.length)return;
    const {es,kf,ok}=this.slotElevation(),cells=[];for(const s of on)for(let l=s*rr;l<(s+1)*rr;l++)cells.push(l);
    const ei=this.interp(es,this.fineBuf("_ae"),ok,cells),ki=this.interp(kf,this.fineBuf("_ak"),ok,cells);
    // Relief L0 cannot resolve: fractal noise from 1.5 L0 cells down to 2 L1 cells (50 m in lowlands, more in
    // highlands), fixed in each piece's frame. It is elevation only (no thickness), and gives erosion something to
    // organize; on a bilinear surface D8 routing makes parallel streams.
    const r=G.r,f0=R_KM/(1.5*this.world.cellKm),oct=Math.max(1,Math.ceil(Math.log2(1.5*r/2))+1),p=G.pos;
    for(const l of cells){this.alive[l]=1;this.eroded[l]=0;this.uplifted[l]=0;
      if(ki[l]>=0.5||ei[l]>0){const amp=0.05+0.05*Math.max(0,ei[l]),v=this.nz(p[3*l],p[3*l+1],p[3*l+2],f0,oct);
        this.cont[l]=1;this.T[l]=this.T0[l]=C.CONT_BASE+ei[l]/C.ISO;this.e0[l]=ei[l]+amp*2*(v-0.5)*2;b.activated+=this.T[l]*G.area[l];}
      else{this.cont[l]=0;this.eOcean[l]=ei[l];}}
    for(const s of on)slotActive[s]=1;
    for(const l of cells)this.live[l]=1;
    if(t>0)this.events.push({t,type:"spawn",slots:on.length});}
  // rigid transform (piece frame -> world) for slot s
  partF(s){const pc=this.pieces[this.pieceOf[s]],p=this.slotPlate[s],M=pc.frames.get(p);return M?mat(this.world.plates[p].Q,M):null;}

  // ---- the world-space layer ----
  // Every parcel of L1 crust spans about one L0 cell, so the layer works parcel-first: each active parcel's centre
  // is placed in the world and its L0 cell found (one lookup per parcel). Parcels whose L0 neighbourhood holds
  // another rigid part are "near"; parcels whose neighbourhood L0 shows on another plate are "exposed". Only
  // their fine cells are placed, indexed and searched:
  //  - overlap (closer than 0.6 cell widths, different rigid parts): continental beats oceanic; otherwise the part
  //    that L0 shows on top there wins (majority over the overlap), else the thicker crust. The loser is consumed.
  //  - crust outside L1 that L0 shows on top of a cell (bilinear weight over half) consumes it too.
  //  - contact (within 1.5 cell widths, different parts): a temporary link, rebuilt every step, so rivers,
  //    flexure and crustal flow cross collision zones before (and after) any weld.
  place(){
    const w=this.world,{G,live,alive,cont,T,W3,slotPlate,pieceOf,index}=this,nC=G.nC,rr=G.rr,pos=G.pos,nS=this.slots.length,cube=w.cube,cpos=cube.pos,nb=cube.nb,own=w.owner,N=w.N;
    const partOfSlot=s=>pieceOf[s]*65536+slotPlate[s],part=l=>partOfSlot((l/rr)|0);
    // parcels: transform, world L0 cell, and per L0 cell which part is there (or "mixed")
    const Fs=this.Fs=new Array(nS),cw=this.slotCell&&this.slotCell.length===nS?this.slotCell:(this.slotCell=new Int32Array(nS));
    const first=this._cF&&this._cF.length===N?this._cF:(this._cF=new Int32Array(N)),mixed=this._cM&&this._cM.length===N?this._cM:(this._cM=new Uint8Array(N));first.fill(-1);mixed.fill(0);
    for(let s=0;s<nS;s++){cw[s]=-1;if(!this.slotActive[s])continue;const F=this.partF(s);if(!F)continue;Fs[s]=F;const k=this.slots[s],x=cpos[3*k],y=cpos[3*k+1],z=cpos[3*k+2];
      const c=cube.lookup(F[0]*x+F[1]*y+F[2]*z,F[3]*x+F[4]*y+F[5]*z,F[6]*x+F[7]*y+F[8]*z);cw[s]=c;const p=partOfSlot(s);if(first[c]<0)first[c]=p;else if(first[c]!==p)mixed[c]=1;}
    // the parcels in each L0 cell, for drawing
    {const off=new Int32Array(N+1);for(let s=0;s<nS;s++)if(cw[s]>=0)off[cw[s]+1]++;for(let c=0;c<N;c++)off[c+1]+=off[c];const list=new Int32Array(off[N]),fill=off.slice(0,N);
      for(let s=0;s<nS;s++)if(cw[s]>=0)list[fill[cw[s]]++]=s;this.cellSlots={off,list};}
    const ring=c=>[c,nb[4*c],nb[4*c+1],nb[4*c+2],nb[4*c+3],nb[4*nb[4*c]+2],nb[4*nb[4*c]+3],nb[4*nb[4*c+1]+2],nb[4*nb[4*c+1]+3]];
    const near=this._near&&this._near.length===nS?this._near:(this._near=new Uint8Array(nS)),exposed=this._exp&&this._exp.length===nS?this._exp:(this._exp=new Uint8Array(nS));
    for(let s=0;s<nS;s++){near[s]=0;exposed[s]=0;if(cw[s]<0)continue;const p=partOfSlot(s),sp=slotPlate[s];
      for(const c of ring(cw[s])){if(mixed[c]||(first[c]>=0&&first[c]!==p))near[s]=1;if(own[c]!==sp)exposed[s]=1;}}
    // fine cells of near or exposed parcels: placed and indexed
    index.clear();const placed=this._pl&&this._pl.length===nC?this._pl:(this._pl=new Uint8Array(nC));placed.fill(0);
    for(let s=0;s<nS;s++){if(!Fs[s]){for(let l=s*rr;l<(s+1)*rr;l++)live[l]=0;continue;}if(!(near[s]||exposed[s]))continue;const F=Fs[s];
      for(let l=s*rr;l<(s+1)*rr;l++){if(!live[l])continue;const x=pos[3*l],y=pos[3*l+1],z=pos[3*l+2];
        W3[3*l]=F[0]*x+F[1]*y+F[2]*z;W3[3*l+1]=F[3]*x+F[4]*y+F[5]*z;W3[3*l+2]=F[6]*x+F[7]*y+F[8]*z;placed[l]=1;index.insert(l,W3[3*l],W3[3*l+1],W3[3*l+2]);}}
    // 1. overlaps and contacts between L1 parts, in one search; overlaps decided per pair of parts by L0's view
    const pairs=[],votes=new Map(),lk=cube.lookup,cand=[],deadAt=this.deadAt||(this.deadAt=new Float32Array(nC));
    for(let a=0;a<nC;a++){if(!placed[a]||!live[a]||!near[(a/rr)|0])continue;const pa=part(a),h=Math.sqrt(G.area[a])/R_KM,ho=0.36*h*h,hc=2.25*h*h,na=index.neighbours(index.g[a]);
      for(let q=0;q<na;q++)for(let b=index.first(index.around[q]);b>=0;b=index.next[b]){if(b<=a||!live[b])continue;const pb=part(b);if(pb===pa)continue;
        const dx=W3[3*a]-W3[3*b],dy=W3[3*a+1]-W3[3*b+1],dz=W3[3*a+2]-W3[3*b+2],d2=dx*dx+dy*dy+dz*dz;if(d2>hc)continue;
        cand.push(a,b,d2);if(d2>ho)continue;
        pairs.push(a,b);if(cont[a]===cont[b]){const o=own[lk(W3[3*a],W3[3*a+1],W3[3*a+2])],key=pa<pb?pa+"|"+pb:pb+"|"+pa,v=votes.get(key)||[0,0];
          if(o===slotPlate[(a/rr)|0])v[pa<pb?0:1]++;else if(o===slotPlate[(b/rr)|0])v[pa<pb?1:0]++;votes.set(key,v);}}}
    let consumed=0;
    for(let q=0;q<pairs.length;q+=2){const a=pairs[q],b=pairs[q+1];if(!live[a]||!live[b])continue;let lose;
      if(cont[a]!==cont[b])lose=cont[a]?b:a;
      else{const pa=part(a),pb=part(b),v=votes.get(pa<pb?pa+"|"+pb:pb+"|"+pa)||[0,0],lo=pa<pb?a:b,hi=pa<pb?b:a;
        lose=v[0]>v[1]?hi:v[1]>v[0]?lo:(T[a]>=T[b]?b:a);}
      if(cont[lose])consumed+=T[lose]*G.area[lose];alive[lose]=0;live[lose]=0;deadAt[lose]=this.t;}
    this.budget.consumed+=consumed;
    // 2. crust outside L1 that L0 puts on top of a cell
    const out=[0,0,0],plates=w.plates,bi=w.bodyIdx;let outside=0;
    const notL1=i=>{const P=plates[own[i]],id=P.pid?P.pid[bi[i]]:-1;if(id<0)return true;const p0=(id/N)|0,pi=p0<this.pieceOfPlate0.length?this.pieceOfPlate0[p0]:-1;if(pi<0)return true;
      const pc=this.pieces[pi],s=pc.slotOf[id-p0*N];return s<0||!this.slotActive[s+pc.s0];};
    // bilinear weight, at a world position, of L0 crust on another plate that L1 does not simulate
    const outsideWeight=(X,Y,Z,sp)=>{cube.lookupFrac(X,Y,Z,out);
      const c=out[0],tx=out[1],ty=out[2],cB=nb[4*c+(tx>=0?0:1)],cC=nb[4*c+(ty>=0?2:3)],cD=nb[4*cB+(ty>=0?2:3)],wx=Math.abs(tx),wy=Math.abs(ty);
      if(own[c]===sp&&own[cB]===sp&&own[cC]===sp&&own[cD]===sp)return 0;
      let wt=0;if(own[c]!==sp&&notL1(c))wt+=(1-wx)*(1-wy);if(own[cB]!==sp&&notL1(cB))wt+=wx*(1-wy);if(own[cC]!==sp&&notL1(cC))wt+=(1-wx)*wy;if(own[cD]!==sp&&notL1(cD))wt+=wx*wy;return wt;};
    for(let l=0;l<nC;l++){if(!placed[l]||!live[l]||!exposed[(l/rr)|0])continue;
      if(outsideWeight(W3[3*l],W3[3*l+1],W3[3*l+2],slotPlate[(l/rr)|0])>0.5){if(cont[l])outside+=T[l]*G.area[l];alive[l]=0;live[l]=0;deadAt[l]=this.t;}}
    this.budget.consumedOutside+=outside;
    // 2b. refill holes, as L0 does: a consumed cell of an active parcel that has stayed consumed for 5 Myr (gaps that
    // open for a step at a moving collision front are left alone) comes back when neither rule above would consume
    // it again (no live cell of another part within overlap distance, no outside crust on top by the same bilinear
    // test) and at least three neighbours on the same part are alive. It takes their mean state (the inherited
    // reference too), and the volume is booked as refilled.
    {const nbr=G.nbr,T0=this.T0,e0=this.e0,eO=this.eOcean;let refill=0;const back=[];
      for(let l=0;l<nC;l++){if(alive[l]||this.t-deadAt[l]<5)continue;const s=(l/rr)|0;if(!this.slotActive[s]||!Fs[s])continue;
        let n=0,nc=0,sT=0,sT0=0,se0=0,so=0;for(let k=0;k<8;k++){const j=nbr[8*l+k];if(j<0||!live[j]||slotPlate[(j/rr)|0]!==slotPlate[s])continue;n++;
          if(cont[j]){nc++;sT+=T[j];sT0+=T0[j];se0+=e0[j];}else so+=eO[j];}
        if(n<3)continue;
        const F=Fs[s],x=pos[3*l],y=pos[3*l+1],z=pos[3*l+2],X=F[0]*x+F[1]*y+F[2]*z,Y=F[3]*x+F[4]*y+F[5]*z,Z=F[6]*x+F[7]*y+F[8]*z;
        if(exposed[s]&&outsideWeight(X,Y,Z,slotPlate[s])>0.5)continue;
        if(near[s]){const pl=partOfSlot(s),h=0.6*Math.sqrt(G.area[l])/R_KM;
          let blocked=false;const na=index.neighbours(index.lk(X,Y,Z));for(let q=0;q<na&&!blocked;q++)for(let b=index.first(index.around[q]);b>=0;b=index.next[b]){if(!live[b]||part(b)===pl)continue;
            const dx=W3[3*b]-X,dy=W3[3*b+1]-Y,dz=W3[3*b+2]-Z;if(dx*dx+dy*dy+dz*dz<h*h){blocked=true;break;}}
          if(blocked)continue;}
        back.push([l,nc*2>=n,nc?sT/nc:0,nc?sT0/nc:0,nc?se0/nc:0,n-nc?so/(n-nc):0]);}
      for(const [l,isC,t,t0,ee0,eo] of back){alive[l]=1;live[l]=1;this.eroded[l]=0;this.uplifted[l]=0;
        if(isC){cont[l]=1;T[l]=t;T0[l]=t0;e0[l]=ee0;refill+=t*G.area[l];}else{cont[l]=0;eO[l]=eo;}}
      this.budget.refilled+=refill;}
    // contacts: the surviving candidate pairs
    {const pa=[],pb=[],pd=[];for(let q=0;q<cand.length;q+=3){const a=cand[q],b=cand[q+1];if(live[a]&&live[b]){pa.push(a);pb.push(b);pd.push(Math.max(1e-3,Math.sqrt(cand[q+2])*R_KM));}}
      this.X=csr(nC,pa,pb,pd);}
    // 3. the same on the mid grid (flexure and crustal flow cross contacts), and its adjacency for this step
    const Gm=this.Gm,mC=Gm.nC,rf2=this.rf*this.rf,Wm=this.Wm,mlive=this._ml||(this._ml=new Uint8Array(mC)),mp=this.midPlate||(this.midPlate=new Int32Array(mC));
    mlive.fill(0);for(let l=0;l<nC;l++)if(live[l])mlive[this.mid[l]]=1;
    const mi=this.midIndex;mi.clear();const mnear=this._mn||(this._mn=new Uint8Array(mC));
    for(let m=0;m<mC;m++){const s=(m/rf2)|0;mp[m]=partOfSlot(s);mnear[m]=0;if(!mlive[m]||!near[s])continue;const F=Fs[s];
      const x=Gm.pos[3*m],y=Gm.pos[3*m+1],z=Gm.pos[3*m+2];Wm[3*m]=F[0]*x+F[1]*y+F[2]*z;Wm[3*m+1]=F[3]*x+F[4]*y+F[5]*z;Wm[3*m+2]=F[6]*x+F[7]*y+F[8]*z;mnear[m]=1;mi.insert(m,Wm[3*m],Wm[3*m+1],Wm[3*m+2]);}
    this.Xm=this.contacts(Gm,Wm,mnear,mi,m=>mp[m],1.5);
    {const nbm=Gm.nbr,X=this.Xm,ma=[],mb=[];for(let m=0;m<mC;m++){if(!mlive[m])continue;for(let k=0;k<4;k++){const j=nbm[8*m+k];if(j>m&&mlive[j]&&mp[j]===mp[m]){ma.push(m);mb.push(j);}}
        for(let q=X.off[m];q<X.off[m+1];q++){const j=X.adj[q];if(j>m){ma.push(m);mb.push(j);}}}
      this.Am=csr(mC,ma,mb,null);}
  }
  // contact links between cells of different parts within `reach` cell widths, as CSR (both directions)
  contacts(G,P3,live,index,part,reach,near){const nC=G.nC,pa=[],pb=[],pd=[];
    for(let a=0;a<nC;a++){if(!live[a]||(near&&!near[a]))continue;const ka=part(a),h=reach*Math.sqrt(G.area[a])/R_KM,h2=h*h,na=index.neighbours(index.g[a]);
      for(let q=0;q<na;q++)for(let b=index.first(index.around[q]);b>=0;b=index.next[b]){if(b<=a||!live[b]||part(b)===ka)continue;
        const dx=P3[3*a]-P3[3*b],dy=P3[3*a+1]-P3[3*b+1],dz=P3[3*a+2]-P3[3*b+2],d2=dx*dx+dy*dy+dz*dz;if(d2>h2)continue;
        pa.push(a);pb.push(b);pd.push(Math.max(1e-3,Math.sqrt(d2)*R_KM));}}
    return csr(nC,pa,pb,pd);}
  // static links that are still valid: both ends live and on the same plate (a rift breaks them)
  effGraph(){const {G,live,slotPlate}=this,nC=G.nC,rr=G.rr,nbr=G.nbr,eff=this._eff||(this._eff=new Int32Array(8*nC));
    for(let l=0;l<nC;l++){const p=slotPlate[(l/rr)|0];for(let k=0;k<8;k++){const j=nbr[8*l+k];eff[8*l+k]=j>=0&&live[l]&&live[j]&&slotPlate[(j/rr)|0]===p?j:-1;}}
    return{nC,nbr:eff,dist:G.dist,area:G.area};}
  // Land cells with a missing neighbour (a rift, the edge of the domain, consumed crust) border open sea when
  // L0 shows sea at the missing neighbour's position and no L1 cell is there; they then drain to the sea.
  // Where L0 shows land outside L1, the edge stays a drainage divide.
  seaEdges(){const w=this.world,{G,live,cont,e,outlet,W3,index}=this,nC=G.nC,rr=G.rr,nbr=G.nbr,sea=this.S.seaEdge,lk=w.cube.lookup,el=w.elev,t=[0,0,0],sp=this.slotPlate;
    sea.fill(0);const NR=G.NR,q=Math.PI/4,ang=v=>Math.tan(((v+.5)/NR*2-1)*q);let F=null,Fs=-1;
    for(let l=0;l<nC;l++){if(!live[l]||!cont[l]||outlet[l])continue;const s=(l/rr)|0,p=sp[s];let miss=false;
      for(let k=0;k<8&&!miss;k++){const j=nbr[8*l+k];if(j<0||!live[j]||sp[(j/rr)|0]!==p)miss=true;}
      if(!miss)continue;if(Fs!==s){F=this.partF(s);Fs=s;}if(!F)continue;
      const g=G.gid[l],f=(g/(NR*NR))|0,rem=g-f*NR*NR,J=(rem/NR)|0,I=rem-J*NR,h=0.75*Math.sqrt(G.area[l])/R_KM;
      for(let k=0;k<8;k++){const j=nbr[8*l+k];if(j>=0&&live[j]&&sp[(j/rr)|0]===p)continue;
        faceDir(f,ang(I+D8[k][0]),ang(J+D8[k][1]),t);const x=F[0]*t[0]+F[1]*t[1]+F[2]*t[2],y=F[3]*t[0]+F[4]*t[1]+F[5]*t[2],z=F[6]*t[0]+F[7]*t[1]+F[8]*t[2];
        if(el[lk(x,y,z)]>0)continue;
        // occupied by some live L1 cell? then it is a contact, not open sea
        let occ=false;const gg=index.lk(x,y,z),na=index.neighbours(gg);
        for(let a=0;a<na&&!occ;a++)for(let b=index.first(index.around[a]);b>=0;b=index.next[b]){if(!live[b])continue;const dx=W3[3*b]-x,dy=W3[3*b+1]-y,dz=W3[3*b+2]-z;if(dx*dx+dy*dy+dz*dz<h*h){occ=true;break;}}
        if(!occ){sea[l]=1;break;}}}}

  // tectonic thickening recorded by L0 this step, mapped onto slots by parcel identity, then spread twice
  // (volume-conserving) over neighbouring parcels: against L0's resampling noise, and so that input concentrated in
  // one L0 cell becomes a rounded bell rather than a square pyramid from bilinear interpolation
  collectForcing(dt){
    const w=this.world,N=w.N,rec=this.rec0,{slotNb,rate,slotActive,pieces,pieceOfPlate0}=this,nS=this.slots.length,np0=pieceOfPlate0.length;
    const raw=this._raw||(this._raw=new Float64Array(nS));raw.fill(0);
    const own=w.owner,bi=w.bodyIdx,areaC=w.cube.area,kmA=R_KM*R_KM,plates=w.plates;let tot=0;
    for(let i=0;i<N;i++){const v=rec[i];if(v===0)continue;const P=plates[own[i]],id=P.pid?P.pid[bi[i]]:-1;if(id<0)continue;const p0=(id/N)|0;if(p0>=np0)continue;
      const pi=pieceOfPlate0[p0];if(pi<0)continue;const pc=pieces[pi],s0=pc.slotOf[id-p0*N];if(s0<0)continue;const s=s0+pc.s0;if(!slotActive[s])continue;raw[s]+=v;tot+=v*areaC[i]*kmA;}
    this.budget.tectonicL0+=tot;
    let src=raw;
    for(let pass=0;pass<2;pass++){const dst=pass===0?(this._raw2||(this._raw2=new Float64Array(nS))):rate;dst.fill(0);
      for(let s=0;s<nS;s++){if(!slotActive[s]||src[s]===0)continue;let c=0;for(let k=0;k<4;k++){const u=slotNb[4*s+k];if(u>=0&&slotActive[u])c++;}
        dst[s]+=src[s]*(1-0.125*c);for(let k=0;k<4;k++){const u=slotNb[4*s+k];if(u>=0&&slotActive[u])dst[u]+=src[s]*0.125;}}
      src=dst;}
    for(let s=0;s<nS;s++)rate[s]/=dt;
  }
  // one coupled step: L0 advances (recording its tectonics), then L1 follows its parcels, brings crust in or
  // drops it, resolves where crust meets, applies L0's tectonics, and erodes
  advance(){
    const w=this.world;this.rec0.fill(0);w.tecRec=this.rec0;w.step();w.tecRec=null;const dt=w.dt;
    this.updateFrames();this.locate();this.activate(this.t+dt);this.place();this.collectForcing(dt);this.step(dt);this.t+=dt;this.steps++;
  }
  step(dt){
    const {G,T,cont,live,eroded,uplifted,S,budget}=this,nC=G.nC,area=G.area;
    const u=this.interp(this.rate,this._u||(this._u=new Float64Array(nC)),this.slotActive);
    let up=0;
    for(let l=0;l<nC;l++){if(!live[l]||u[l]===0)continue;const d=u[l]*dt;
      if(cont[l]){T[l]+=d;uplifted[l]+=d;up+=d*area[l];}else if(d>0)budget.unappliedOcean+=d*area[l];}
    budget.tectonicApplied+=up;
    this.crustFlow(dt);
    let cap=0;for(let l=0;l<nC;l++)if(live[l]&&cont[l]&&T[l]>C.CAP){cap+=(T[l]-C.CAP)*area[l];T[l]=C.CAP;}
    budget.capped+=cap;
    this.flexure();this.classify();this.seaEdges();
    route(this.effGraph(),this.e,live,this.outlet,S,this.precip,this.X);
    const hN=streamPower(this.e,S,this.K,this.o.m,dt),e=this.e,outlet=this.outlet;
    let ero=0;
    for(let l=0;l<nC;l++){if(!live[l]||!cont[l]||outlet[l])continue;const d=e[l]-hN[l];if(d<=0)continue;
      T[l]-=d;eroded[l]+=d;ero+=d*area[l];e[l]=hN[l];}
    budget.eroded+=ero;
    this.stats(ero/dt,up/dt);
  }
  stats(ero,up){const {G,e,cont,live}=this,nC=G.nC;let nl=0,se=0,mx=-1e9,al=0,cells=0;
    for(let l=0;l<nC;l++){if(!live[l])continue;cells++;if(!cont[l]||e[l]<=0)continue;nl++;se+=e[l]*G.area[l];al+=G.area[l];if(e[l]>mx)mx=e[l];}
    this.history.push({t:this.t,cells,landKm2:al,meanElevKm:al?se/al:0,maxElevKm:nl?mx:0,erosionKm3PerMyr:ero,upliftKm3PerMyr:up,erosionRateKmPerMyr:al?ero/al:0,upliftRateKmPerMyr:al?up/al:0,contacts:this.X?this.X.count:0});}
  // crust volume: what is there now, against what the bookkeeping says should be there
  balance(){const b=this.budget,now=this.volume(),expect=b.initial+b.activated-b.dropped+b.tectonicApplied+b.flowResidual+b.refilled-b.eroded-b.capped-b.consumed-b.consumedOutside;
    return{now,expect,error:(now-expect)/b.initial};}
  // coarse-grained L1 thickness against L0's own thickness, over parcels that are continental in both
  divergence(){const w=this.world,{G,T,cont,live}=this,rr=G.rr,nS=this.slots.length;let n=0,s1=0,s2=0,e2=0;
    for(let s=0;s<nS;s++){if(!this.slotActive[s]||!this.slotFound[s])continue;const P=w.plates[this.slotPlate[s]],k=this.slotBody[s];if(!P.K[k])continue;
      let v=0,a=0,ok=true;for(let l=s*rr;l<(s+1)*rr;l++){if(!live[l]||!cont[l]){ok=false;break;}v+=T[l]*G.area[l];a+=G.area[l];}if(!ok)continue;
      const d=v/a-P.T[k],de=C.ISO*d;n++;s1+=d;s2+=d*d;e2+=de*de;}
    return{parcels:n,meanThicknessDiffKm:n?s1/n:0,rmsThicknessDiffKm:n?Math.sqrt(s2/n):0,rmsElevationDiffKm:n?Math.sqrt(e2/n):0};}
  // The map is drawn in the frame of the plate holding most of the first piece's active parcels. Every piece's
  // parcels on each plate (a "part") move rigidly; for each part, A maps a map-frame direction into the piece's
  // own frame. F0 maps the map frame to the world.
  mapFrame(){const plates=this.world.plates,pc=this.pieces[0],on=this.piecePlates?this.piecePlates[0]:new Map([[pc.plate0,1]]);
    let X0=-1,bn=-1;for(const [q,c] of on)if(c>bn&&pc.frames.has(q)){bn=c;X0=q;}if(X0<0)X0=[...pc.frames.keys()][0];
    return mat(plates[X0].Q,pc.frames.get(X0));}
  viewParts(){const plates=this.world.plates,F0=this.mapFrame(),out=[];
    (this.piecePlates||[]).forEach((on,p)=>{const pc=this.pieces[p];for(const [q] of on){const M=pc.frames.get(q);if(!M)continue;out.push({piece:p,plate:q,A:mat(matT(mat(plates[q].Q,M)),F0)});}});
    return out;}
  // world direction -> the live fine cell there, or -1: try the parcels in that L0 cell and its neighbours, each
  // through its own rigid transform and its piece's grid
  cellAtWorld(x,y,z){const w=this.world,cube=w.cube,nb=cube.nb,{off,list}=this.cellSlots,G=this.G,grids=G.parts,offs=G.offsets,rr=G.rr,c=cube.lookup(x,y,z);
    const cs=[c,nb[4*c],nb[4*c+1],nb[4*c+2],nb[4*c+3]];
    for(const cc of cs)for(let q=off[cc];q<off[cc+1];q++){const s=list[q],F=this.Fs[s];if(!F)continue;const pi=this.pieceOf[s],g=grids[pi];
      const X=F[0]*x+F[3]*y+F[6]*z,Y=F[1]*x+F[4]*y+F[7]*z,Z=F[2]*x+F[5]*y+F[8]*z,m=g.localOf(g.lookup(X,Y,Z));if(m<0)continue;const l=m+offs[pi],sl=(l/rr)|0;
      if(this.live[l]&&this.pieceOf[sl]===pi&&this.slotPlate[sl]===this.slotPlate[s])return l;}
    return -1;}
  // the slot holding the crust on fixed-grid cell i (by parcel identity), or -1
  slotAtFixed(i){const w=this.world,N=w.N,P=w.plates[w.owner[i]],id=P.pid?P.pid[w.bodyIdx[i]]:-1;if(id<0)return -1;const p0=(id/N)|0;
    const pi=p0<this.pieceOfPlate0.length?this.pieceOfPlate0[p0]:-1;if(pi<0)return -1;const pc=this.pieces[pi],s=pc.slotOf[id-p0*N];return s<0?-1:s+pc.s0;}
  basins(){const {S,G}=this,out=this._basin||(this._basin=new Int32Array(G.nC).fill(-1));out.fill(-1);
    for(let h=0;h<S.ns;h++){const l=S.stack[h],r=S.rec[l];out[l]=r===l?l:out[r];}return out;}
}
return {Region,Lookahead,CellIndex,select,selectBody,buildGrid,route,streamPower,planarGrid,scratch,Heap};
})();
// ===== L1-END =====
