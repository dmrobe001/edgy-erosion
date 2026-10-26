// ===== L1-START =====
// L1: one continent, refined from L0 and evolved forward in its own body frame while L0 keeps running.
// L0 supplies tectonic thickening (arcs, collisions, hotspots, collapse); L1 replaces L0's surface processes
// with stream-power river erosion, and derives elevation from crust thickness through flexural isostasy.
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
// Route flow over a graph with 8 neighbour slots per cell (-1 = none). Outlets are fixed base level.
// Depressions are filled with an epsilon gradient (Barnes et al. 2014, Priority-Flood+epsilon) so every
// cell drains to an outlet; receivers are steepest descent on the filled surface; area accumulates down the tree.
function route(G,e,alive,outlet,S,precip){
  const {nC,nbr,dist,area}=G,{f,closed,rec,rdist,stack,ndon,don,A,heap,pit}=S;
  heap.n=0;let ph=0,pt=0;
  for(let l=0;l<nC;l++){closed[l]=0;if(!alive[l]){closed[l]=1;rec[l]=-1;}else if(outlet[l]){closed[l]=1;f[l]=e[l];}}
  for(let l=0;l<nC;l++){if(!alive[l]||!outlet[l])continue;
    for(let k=0;k<8;k++){const j=nbr[8*l+k];if(j>=0&&alive[j]&&!outlet[j]){heap.push(l,f[l]);break;}}}
  // Zhou et al. (2016) variant: a cell strictly higher than a finalized neighbour is final at its own height,
  // so whole hillslopes are swept uphill with a plain queue; only cells that still have lower, unresolved
  // neighbours wait in the heap. Same filled surface as the plain algorithm, far fewer heap operations.
  const up=S.up||(S.up=new Int32Array(nC));
  const flood=()=>{while(ph<pt||heap.n>0){const c=ph<pt?pit[ph++]:heap.pop(),fc=f[c];
    for(let k=0;k<8;k++){const j=nbr[8*c+k];if(j<0||closed[j])continue;closed[j]=1;
      if(e[j]<=fc){f[j]=fc+EPS;pit[pt++]=j;continue;}
      f[j]=e[j];let uh=0,ut=0;up[ut++]=j;
      while(uh<ut){const a=up[uh++],fa=f[a];let defer=false;
        for(let q=0;q<8;q++){const b=nbr[8*a+q];if(b<0||closed[b])continue;
          if(e[b]>fa){closed[b]=1;f[b]=e[b];up[ut++]=b;}else defer=true;}
        if(defer)heap.push(a,fa);}}}};
  flood();
  // land the flood never reached is cut off from the sea (e.g. enclosed by consumed parcels): each such piece
  // becomes an inland basin draining to its lowest cell
  const sink=S.sink||(S.sink=new Uint8Array(nC));sink.fill(0);S.sinks=0;
  for(;;){let lo=-1;for(let l=0;l<nC;l++)if(!closed[l]&&(lo<0||e[l]<e[lo]))lo=l;if(lo<0)break;
    closed[lo]=1;sink[lo]=1;f[lo]=e[lo];heap.push(lo,f[lo]);S.sinks++;flood();}
  for(let l=0;l<nC;l++){if(!alive[l])continue;
    if(outlet[l]||sink[l]){rec[l]=l;rdist[l]=0;continue;}
    let best=l,bs=0,bd=0;
    for(let k=0;k<8;k++){const j=nbr[8*l+k];if(j<0||!alive[j])continue;const d=dist[8*l+k],s=(f[l]-f[j])/d;if(s>bs){bs=s;best=j;bd=d;}}
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
// A cell lower than its receiver's new height (a filled depression) is left alone. Returns hN, the new heights.
function streamPower(e,S,K,m,dt){
  const {rec,rdist,stack,A,hN,ns}=S;
  for(let h=0;h<ns;h++){const l=stack[h],r=rec[l];
    if(r===l){hN[l]=e[l]>0?e[l]:0;continue;}
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
// Link edge cells of piece i (neighbour slot -1) to the cell of piece j one step beyond, if any. Rji maps piece j's
// frame into piece i's, Rij the reverse. Distances are measured between nodes in piece i's frame.
function stitchGrid(G,i,j,Rij,Rji,use){const Gi=G.parts[i],Gj=G.parts[j],oi=G.offsets[i],oj=G.offsets[j],NR=Gi.NR,q=Math.PI/4,t=[0,0,0];let links=0;
  const ang=v=>Math.tan(((v+.5)/NR*2-1)*q);
  for(let l=0;l<Gi.nC;l++){const L=l+oi;if(use&&!use(L))continue;let hasEdge=false;for(let k=0;k<8;k++)if(G.nbr[8*L+k]<0){hasEdge=true;break;}if(!hasEdge)continue;
    const g=Gi.gid[l],f=(g/(NR*NR))|0,rem=g-f*NR*NR,J=(rem/NR)|0,I=rem-J*NR;
    for(let k=0;k<8;k++){if(G.nbr[8*L+k]>=0)continue;faceDir(f,ang(I+D8[k][0]),ang(J+D8[k][1]),t);
      const x=Rij[0]*t[0]+Rij[1]*t[1]+Rij[2]*t[2],y=Rij[3]*t[0]+Rij[4]*t[1]+Rij[5]*t[2],z=Rij[6]*t[0]+Rij[7]*t[1]+Rij[8]*t[2];
      const m=Gj.localOf(Gj.lookup(x,y,z));if(m<0||(use&&!use(m+oj)))continue;const M=m+oj,nd=G.node;
      const a=nd[3*M],b=nd[3*M+1],c=nd[3*M+2],X=Rji[0]*a+Rji[1]*b+Rji[2]*c,Y=Rji[3]*a+Rji[4]*b+Rji[5]*c,Z=Rji[6]*a+Rji[7]*b+Rji[8]*c;
      G.nbr[8*L+k]=M;G.dist[8*L+k]=Math.hypot(nd[3*L]-X,nd[3*L+1]-Y,nd[3*L+2]-Z)*R_KM;links++;}}
  return links;}
// Pieces (one per continent, each in its own plate's frame) are concatenated into one graph; cell indices of
// piece p are offset by its first slot x r^2. There are no links between pieces until they are stitched.
function concatGrids(list){if(list.length===1)return list[0];const g0=list[0];let nC=0,nS=0;for(const g of list){nC+=g.nC;nS+=g.nS;}
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

// ---------- looking ahead for sutures ----------
// Run L0 forward through the window from a snapshot, following the picked continent's plate lineage (splits add
// the new plate; a weld with a plate outside the lineage makes it a partner). Each step, remember where outside
// continental crust touches the lineage's continental crust, by plate and body cell. When a partner welds on,
// those last contacts identify its continent: each is traced back in time to where that crust was at t0, through
// splits (a split keeps body indices) and welds (L0's merge log says which absorbed cell filled each surviving
// cell). Then restore the snapshot, so the real run replays the same future.
class Lookahead{
  constructor(world,cell,window,margin){
    const sel=select(world,cell,margin||3);if(!sel)throw new Error("not a continental cell");
    this.world=world;this.t0=world.t;this.window=window;this.snap=world.snapshot();this.seen0=world.plates.length;world.mergeLog=[];this.splits=[];
    this.lineage=new Set([sel.plate]);this.deadSeen=new Set(world.plates.filter(P=>!P.alive).map(P=>P.id));this.platesSeen=world.plates.length;
    this.contacts=new Map();this.found=[];this.done=false;}
  // run up to maxMs of L0 steps; returns true when the window is covered
  advance(maxMs){const w=this.world,t1=Date.now();
    while(!this.done&&Date.now()-t1<maxMs){w.step();const plates=w.plates,lin=this.lineage;
      for(let id=this.platesSeen;id<plates.length;id++){if(!plates[id].alive)continue;this.splits.push({id,parent:plates[id].parent,t:w.t});if(lin.has(plates[id].parent))lin.add(id);}
      this.platesSeen=plates.length;
      for(const P of plates){if(P.alive||P.mergedInto===undefined||this.deadSeen.has(P.id))continue;this.deadSeen.add(P.id);
        const a=P.id,b=P.mergedInto,partner=lin.has(a)&&!lin.has(b)?b:lin.has(b)&&!lin.has(a)?a:-1;
        if(partner>=0){const c=this.contacts.get(partner);this.found.push({plate:partner,t:w.t-this.t0,contacts:c?c.ks:[],tc:c?c.t:w.t});lin.add(a);lin.add(b);}}
      // contacts this step: continental cells of an outside plate next to continental cells of the lineage
      const now=new Map(),O=w.owner,K=w.wK,nb=w.cube.nb;
      for(let i=0;i<w.N;i++){if(!K[i]||!lin.has(O[i]))continue;for(let k=0;k<4;k++){const j=nb[4*i+k];const q=O[j];if(!K[j]||lin.has(q))continue;
        let arr=now.get(q);if(!arr)now.set(q,arr=[]);if(arr.length<40)arr.push(w.bodyIdx[j]);}}
      for(const [q,arr] of now)this.contacts.set(q,{t:w.t,ks:arr});
      if(w.t-this.t0>=this.window-1e-9)this.done=true;}
    return this.done;}
  // restore L0 to t0 and return the partners as [{plate, k, t, nCont}], each a continent present at t0
  // where the crust at (plate q, body cell k) at time t was at t0, or null (crust created after t0)
  trace(q,k,t){const ev=[...this.merges.filter(m=>m.t<=t).map(m=>({m,t:m.t})),...this.splits.filter(x=>x.t<=t).map(x=>({x,t:x.t}))].sort((a,b)=>b.t-a.t);
    for(const e of ev){if(e.m&&e.m.into===q&&e.m.map[k]>=0){k=e.m.map[k];q=e.m.from;}else if(e.x&&e.x.id===q)q=e.x.parent;}
    return q<this.seen0?{q,k}:null;}
  finish(){const w=this.world;this.merges=w.mergeLog||[];w.restore(this.snap);delete w.mergeLog;const out=[];
    for(const f of this.found)for(const k0 of f.contacts){const o=this.trace(f.plate,k0,f.tc);if(!o)continue;
      if(out.some(x=>x.plate===o.q&&x.sel.dist[o.k]===0))continue;
      const sel=selectBody(w,o.q,o.k,1);if(!sel||sel.nCont<8)continue;out.push({plate:o.q,k:o.k,t:f.t,nCont:sel.nCont,sel});}
    return out.map(({sel,...x})=>x);}
}

class Region{
  // o.cell: the picked fixed-grid cell. o.partners: continents that will weld onto it within the window, as
  // [{plate, k}] (a plate and body cell at t0; see lookahead()). Each continent is a piece in its own plate's frame.
  constructor(world,o){
    o=Object.assign({r:8,margin:3,K:0.1,m:0.5,flexKm:40,flowKm2PerMyr:300,strengthKm:3,jitter:0.4,window:100,partners:[]},o);this.o=o;this.world=world;
    const sel=select(world,o.cell,o.margin);if(!sel)throw new Error("not a continental cell");
    const n=world.n,N=world.N,cube=world.cube,nb=cube.nb,t0=world.t,r=o.r,rr=r*r;
    // continents on the same plate share a frame, so they are one piece (the union of their slots)
    const sels=[sel];for(const pt of o.partners){if(sels.some(x=>x.plate===pt.plate&&x.dist[pt.k]>=0))continue;const ps=selectBody(world,pt.plate,pt.k,o.margin);if(!ps)continue;
      const same=sels.find(x=>x.plate===ps.plate);if(!same){sels.push(ps);continue;}
      const d=same.dist;for(let k=0;k<N;k++)if(ps.dist[k]>=0&&(d[k]<0||ps.dist[k]<d[k]))d[k]=ps.dist[k];
      same.slots=Int32Array.from(d.reduce((a,v,k)=>{if(v>=0)a.push(k);return a;},[]));same.nCont+=ps.nCont;}
    this.t0=t0;this.t=0;this.steps=0;this.plate0=sel.plate;this.Q0=Float64Array.from(world.plates[sel.plate].Q);
    // pieces: slot range, the plate whose body frame holds the piece, and the frame's fixed offset M (world = Q(plate) M b)
    this.pieces=[];let nS=0;for(const ps of sels){const slotOf=new Int32Array(N).fill(-1);for(let s=0;s<ps.slots.length;s++)slotOf[ps.slots[s]]=s;
      this.pieces.push({plate0:ps.plate,s0:nS,nS:ps.slots.length,sel:ps,slotOf,frames:new Map([[ps.plate,[1,0,0,0,1,0,0,0,1]]]),nCont:ps.nCont});nS+=ps.slots.length;}
    this.framesSeen=world.plates.length;
    const slots=this.slots=new Int32Array(nS),pieceOf=this.pieceOf=new Uint8Array(nS);this.slotCont=new Uint8Array(nS);
    for(const [pi,pc] of this.pieces.entries())for(let s=0;s<pc.nS;s++){const k=pc.sel.slots[s];slots[pc.s0+s]=k;pieceOf[pc.s0+s]=pi;this.slotCont[pc.s0+s]=pc.sel.dist[k]===0?1:0;}
    this.nCont=sel.nCont;this.slotOf=this.pieces[0].slotOf;
    const so=s=>{const pc=this.pieces[pieceOf[s]];return c=>{const v=pc.slotOf[c];return v<0?-1:v+pc.s0;};};
    // slot lineage: which plate and body cell holds each parcel now
    this.platesSeen=world.plates.length;this.slotPlate=new Int32Array(nS);this.slotBody=Int32Array.from(slots);this.slotAlive=new Uint8Array(nS).fill(1);
    for(let s=0;s<nS;s++)this.slotPlate[s]=this.pieces[pieceOf[s]].plate0;
    this.slotNb=new Int32Array(4*nS);this.quad=new Int32Array(12*nS);
    for(let s=0;s<nS;s++){const f=so(s),c=slots[s];for(let k=0;k<4;k++)this.slotNb[4*s+k]=f(nb[4*c+k]);
      // bilinear stencils: for each slot and quadrant, the three other parents (as slots, -1 outside the domain)
      for(let qd=0;qd<4;qd++){const cB=nb[4*c+(qd&1?1:0)],cC=nb[4*c+(qd&2?3:2)],cD=nb[4*cB+(qd&2?3:2)];
        this.quad[12*s+3*qd]=f(cB);this.quad[12*s+3*qd+1]=f(cC);this.quad[12*s+3*qd+2]=f(cD);}}
    const G=this.G=concatGrids(this.pieces.map(pc=>buildGrid(n,pc.sel.slots,pc.slotOf,r,o.jitter))),nC=G.nC;
    this.subQ=new Uint8Array(rr);this.subW=new Float32Array(2*rr);
    for(let jj=0;jj<r;jj++)for(let ii=0;ii<r;ii++){const tx=(ii+.5)/r-.5,ty=(jj+.5)/r-.5,sb=jj*r+ii;this.subQ[sb]=(tx<0?1:0)+(ty<0?2:0);this.subW[2*sb]=Math.abs(tx);this.subW[2*sb+1]=Math.abs(ty);}
    // flexure and crustal flow run on a coarser nested grid (b x b fine cells per mid cell), then are interpolated back
    const sigma=1.18*o.flexKm,cellKm=world.cellKm;let rf=1;for(let d=r;d>=1;d--)if(r%d===0&&cellKm/d>=sigma/2.5){rf=d;break;}
    this.rf=rf;this.Gm=concatGrids(this.pieces.map(pc=>buildGrid(n,pc.sel.slots,pc.slotOf,rf,0)));this.buildMidMaps();this.flexPasses=Math.max(1,Math.round(sigma*sigma/(0.4*(cellKm/rf)**2)));
    // initial state from L0: elevation smoothed over the domain (three passes), interpolated, inverted through Airy isostasy
    const e0=new Float64Array(nS),kf=new Float64Array(nS);
    for(let s=0;s<nS;s++){const k=slots[s],P=world.plates[this.pieces[pieceOf[s]].plate0];e0[s]=elevOf(P.K[k],P.T[k],P.B[k],P.S[k],t0);kf[s]=P.K[k];}
    let es=e0,tmp=new Float64Array(nS);
    for(let pass=0;pass<3;pass++){for(let s=0;s<nS;s++){let v=es[s],w=1;for(let k=0;k<4;k++){const u=this.slotNb[4*s+k];if(u>=0){v+=es[u];w++;}}tmp[s]=v/w;}[es,tmp]=[tmp,es];}
    this.T=new Float64Array(nC);this.e=new Float64Array(nC);this.cont=new Uint8Array(nC);this.alive=new Uint8Array(nC).fill(1);this.eOcean=new Float32Array(nC);
    this.T0=new Float64Array(nC);this.e0=new Float32Array(nC);this.eroded=new Float32Array(nC);this.uplifted=new Float32Array(nC);this.outlet=new Uint8Array(nC);this.root=new Float64Array(nC);
    const ei=this.interp(es),ki=this.interp(kf);
    // Relief L0 cannot resolve: fractal noise from 1.5 L0 cells down to 2 L1 cells (50 m in lowlands, more in highlands),
    // fixed in the body frame so the same continent looks alike at every L1 resolution. It is elevation only (no
    // thickness), and gives erosion something to organize; on a bilinear surface D8 routing makes parallel streams.
    const nz=makeNoise3(hashStr(String(world.opts.seed)+"|L1")|0),f0=R_KM/(1.5*world.cellKm),oct=Math.max(1,Math.ceil(Math.log2(1.5*r/2))+1);
    for(let l=0;l<nC;l++){
      if(ki[l]>=0.5||ei[l]>0){const p=G.pos,amp=0.05+0.05*Math.max(0,ei[l]),v=nz(p[3*l],p[3*l+1],p[3*l+2],f0,oct);
        this.cont[l]=1;this.T[l]=this.T0[l]=C.CONT_BASE+ei[l]/C.ISO;this.e0[l]=ei[l]+amp*2*(v-0.5)*2;}else{this.eOcean[l]=ei[l];}}
    this.S=scratch(nC);this.rate=new Float64Array(nS);this.rec0=new Float32Array(N);this.slotAt=new Map();this.touched=[];
    this.precip=null;this.K=o.K;this.welds=[];
    this.budget={initial:this.volume(),tectonicApplied:0,tectonicL0:0,unappliedOcean:0,eroded:0,capped:0,consumed:0,rifted:0,flowResidual:0};
    this.history=[];this.flexure();this.classify();route(G,this.e,this.alive,this.outlet,this.S,this.precip);this.stats(0,0);
    // map orientation: the picked continent's centroid, and world north at t0, in its body coordinates
    let cx=0,cy=0,cz=0;for(let s=0;s<this.pieces[0].nS;s++){if(!this.slotCont[s])continue;const c=slots[s];cx+=cube.pos[3*c];cy+=cube.pos[3*c+1];cz+=cube.pos[3*c+2];}
    const cl=Math.hypot(cx,cy,cz);this.center=[cx/cl,cy/cl,cz/cl];const Q=this.Q0;this.northB=[Q[6],Q[7],Q[8]];
  }
  // bilinear interpolation of a per-slot field onto fine cells (parents outside the domain, or dead, are skipped)
  interp(field,out,mask){const {G,quad,subQ,subW}=this,rr=G.rr,nC=G.nC;out=out||new Float64Array(nC);
    for(let l=0;l<nC;l++){const s=(l/rr)|0,sb=l-s*rr,qd=subQ[sb],wx=subW[2*sb],wy=subW[2*sb+1];
      const sB=quad[12*s+3*qd],sC=quad[12*s+3*qd+1],sD=quad[12*s+3*qd+2];
      let v=0,w=0,ww;
      ww=(1-wx)*(1-wy);if(!mask||mask[s]){v+=ww*field[s];w+=ww;}
      ww=wx*(1-wy);if(sB>=0&&(!mask||mask[sB])){v+=ww*field[sB];w+=ww;}
      ww=(1-wx)*wy;if(sC>=0&&(!mask||mask[sC])){v+=ww*field[sC];w+=ww;}
      ww=wx*wy;if(sD>=0&&(!mask||mask[sD])){v+=ww*field[sD];w+=ww;}
      out[l]=w>0?v/w:0;}
    return out;}
  volume(){const {T,cont,alive,G}=this;let v=0;for(let l=0;l<G.nC;l++)if(cont[l]&&alive[l])v+=T[l]*G.area[l];return v;}
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
  // bilinear interpolation of a mid-grid field back to every fine cell with use[l] set (mid cells with W = 0 are skipped)
  fromMid(A,W,use,out){const {G,mid,midQ,midS,mq,mw}=this,nC=G.nC;
    for(let l=0;l<nC;l++){if(!use[l])continue;const m=mid[l],o=12*m+3*midQ[l],wx=mw[2*midS[l]],wy=mw[2*midS[l]+1],mB=mq[o],mC=mq[o+1],mD=mq[o+2];
      let v=0,w=0,ww=(1-wx)*(1-wy);if(W[m]>0){v+=ww*A[m];w+=ww;}
      ww=wx*(1-wy);if(mB>=0&&W[mB]>0){v+=ww*A[mB];w+=ww;}
      ww=(1-wx)*wy;if(mC>=0&&W[mC]>0){v+=ww*A[mC];w+=ww;}
      ww=wx*wy;if(mD>=0&&W[mD]>0){v+=ww*A[mD];w+=ww;}
      out[l]=w>0?v/w:0;}
    return out;}
  fineBuf(name){return this[name]||(this[name]=new Float64Array(this.G.nC));}
  midBuf(name){return this[name]||(this[name]=new Float64Array(this.Gm.nC));}
  // Flexural isostasy, approximated, applied to changes since the start (the inherited state is taken as in
  // equilibrium): e = e0 + dT - (1 - ISO) * F[dT], where F is a Gaussian of width 1.18 x flexural length
  // (the half-power point of the plate's response 1/(1+(k*alpha)^4)). Long wavelengths follow L0's Airy
  // isostasy; loads narrower than the flexural length are held up by the plate's strength, up to a limit.
  flexure(){
    const {G,T,T0,e0,cont,alive,eOcean,root,e}=this,nC=G.nC,X=this.midBuf("_X"),W=this.midBuf("_W"),Y=this.midBuf("_Y"),nb=this.Gm.nbr,mC=X.length;
    const dT=this.fineBuf("_dT");for(let l=0;l<nC;l++)dT[l]=cont[l]?T[l]-T0[l]:0;
    this.toMid(dT,alive,X,W);
    let A=X,B=Y;
    for(let p=0;p<this.flexPasses;p++){for(let m=0;m<mC;m++){if(!(W[m]>0)){B[m]=0;continue;}let v=A[m],w=1;
        for(let k=0;k<4;k++){const j=nb[8*m+k];if(j>=0&&W[j]>0){v+=A[j];w++;}}B[m]=v/w;}[A,B]=[B,A];}
    // the plate yields rather than hold more than strengthKm of topography away from local (Airy) balance
    const lim=this.o.strengthKm;
    const F=this.fromMid(A,W,alive,this.fineBuf("_F"));
    for(let l=0;l<nC;l++){if(!alive[l])continue;let rt=ROOT*F[l];
      const air=ROOT*dT[l];if(rt<air-lim)rt=air-lim;else if(rt>air+lim)rt=air+lim;root[l]=rt;e[l]=cont[l]?e0[l]+dT[l]-rt:eOcean[l];}
  }
  // Lower-crustal flow: thickened crust spreads sideways. At L0's scale its thickness diffusion stands in for this;
  // here it acts only at long wavelengths (on the mid grid), so it relaxes collision piles without filling valleys.
  // It diffuses the thickness change since the start, not the thickness itself: inherited contrasts (thin rifted
  // margins against a thick interior) are part of the equilibrium state, as in flexure(). Diffusing total
  // thickness pumped interior crust out to the margins and raised a ring of islands along the shelf edge.
  // Fluxes between mid cells are symmetric (volume-conserving); interpolating the change back to fine cells is
  // not exactly so, and that residue is booked as budget.flowResidual.
  crustFlow(dt){
    const kappa=this.o.flowKm2PerMyr;if(!(kappa>0))return;
    const {G,T,T0,cont,alive}=this,nC=G.nC,X=this.midBuf("_FX"),W=this.midBuf("_FW"),V=this.midBuf("_FV"),Tm=this.midBuf("_FT"),nb=this.Gm.nbr,mC=X.length,dx=this.world.cellKm/this.rf;
    const use=this._cu||(this._cu=new Uint8Array(nC));for(let l=0;l<nC;l++)use[l]=alive[l]&&cont[l]?1:0;
    const dT=this.fineBuf("_FdT");for(let l=0;l<nC;l++)dT[l]=T[l]-T0[l];
    this.toMid(T,use,Tm,W);this.toMid(dT,use,X,W);const X0=this.midBuf("_FX0");X0.set(X);
    // channel flux grows with the cube of channel thickness (Poiseuille flow), so hot thick crust spreads much faster;
    // edge coefficients are fixed for the step (thickness barely changes within it), then explicit substeps
    const eI=this._eI||(this._eI=new Int32Array(4*mC)),eJ=this._eJ||(this._eJ=new Int32Array(4*mC)),eC=this._eC||(this._eC=new Float64Array(4*mC));
    let ne=0,cmax=0;
    for(let m=0;m<mC;m++){if(!(W[m]>0))continue;for(let k=0;k<4;k++){const j=nb[8*m+k];if(j<=m||!(W[j]>0))continue;
      const t=(Tm[m]+Tm[j])/(2*C.CONT_BASE),c=kappa*t*t*t;eI[ne]=m;eJ[ne]=j;eC[ne]=c*Math.min(W[m],W[j]);ne++;
      if(c>cmax)cmax=c;}}
    const D=cmax*dt/(dx*dx),nsub=Math.max(1,Math.ceil(D/0.2)),h=dt/nsub/(dx*dx);
    for(let it=0;it<nsub;it++){V.fill(0);
      for(let q=0;q<ne;q++){const m=eI[q],j=eJ[q],f=h*eC[q]*(X[j]-X[m]);V[m]+=f;V[j]-=f;}
      for(let m=0;m<mC;m++)if(W[m]>0)X[m]+=V[m]/W[m];}
    for(let m=0;m<mC;m++)X[m]-=X0[m];
    const add=this.fromMid(X,W,use,this.fineBuf("_F"));let res=0;for(let l=0;l<nC;l++){if(!use[l])continue;T[l]+=add[l];res+=add[l]*G.area[l];}
    this.budget.flowResidual+=res;
  }
  // Outlets are sea: ocean cells and drowned continental cells. A land edge of the domain (consumed parcels, or crust
  // beyond the margin) is a closed boundary, a drainage divide, so rivers there drain back into the domain.
  classify(){const {G,cont,alive,e,outlet}=this,nC=G.nC;for(let l=0;l<nC;l++)outlet[l]=alive[l]&&(!cont[l]||e[l]<=0)?1:0;}
  // Follow each parcel through welds (mapped into the survivor's frame) and splits. A split is caught in the step
  // it happens: the new plate holds exactly the parcels that moved, at their old body indices. (Searching older
  // descendants later is wrong: their grids refill those indices with unrelated crust.) A parcel whose body cell
  // is empty afterwards has been consumed.
  updateLineage(){
    const w=this.world,plates=w.plates,pos=w.cube.pos,lk=w.cube.lookup,{slotPlate,slotBody,slotAlive}=this,nS=this.slots.length;
    for(let s=0;s<nS;s++){if(!slotAlive[s])continue;let P=plates[slotPlate[s]],k=slotBody[s];
      for(let g=0;g<16&&!P.alive&&P.mergedInto!==undefined;g++){const Q=P.Q,bx=pos[3*k],by=pos[3*k+1],bz=pos[3*k+2];
        const x=Q[0]*bx+Q[1]*by+Q[2]*bz,y=Q[3]*bx+Q[4]*by+Q[5]*bz,z=Q[6]*bx+Q[7]*by+Q[8]*bz,Wn=plates[P.mergedInto],QW=Wn.Q;
        k=lk(QW[0]*x+QW[3]*y+QW[6]*z,QW[1]*x+QW[4]*y+QW[7]*z,QW[2]*x+QW[5]*y+QW[8]*z);P=Wn;}
      slotPlate[s]=P.id;slotBody[s]=k;}
    for(let id=this.platesSeen;id<plates.length;id++){const Cp=plates[id];if(!Cp.alive||Cp.parent===undefined)continue;
      for(let s=0;s<nS;s++)if(slotAlive[s]&&slotPlate[s]===Cp.parent&&Cp.T[slotBody[s]]>0)slotPlate[s]=id;}
    this.platesSeen=plates.length;
    const died=[];
    for(let s=0;s<nS;s++){if(!slotAlive[s])continue;const P=plates[slotPlate[s]];if(!P.alive||P.T[slotBody[s]]<=0){slotAlive[s]=0;died.push(s);}}
    if(died.length){const {G,T,cont,alive}=this,rr=G.rr;
      for(const s of died)for(let l=s*rr;l<(s+1)*rr;l++){if(alive[l]&&cont[l])this.budget.consumed+=T[l]*G.area[l];alive[l]=0;}}
  }
  // tectonic thickening recorded by L0 this step, mapped onto slots, smoothed once (volume-conserving) against resampling noise
  collectForcing(dt){
    const w=this.world,N=w.N,rec=this.rec0,{slotPlate,slotBody,slotAlive,slotNb,rate}=this,nS=this.slots.length;
    for(const [arr,k] of this.touched)arr[k]=-1;this.touched.length=0;
    for(let s=0;s<nS;s++){if(!slotAlive[s])continue;const p=slotPlate[s];let arr=this.slotAt.get(p);if(!arr){arr=new Int32Array(N).fill(-1);this.slotAt.set(p,arr);}
      arr[slotBody[s]]=s;this.touched.push([arr,slotBody[s]]);}
    const raw=this._raw||(this._raw=new Float64Array(nS));raw.fill(0);
    const own=w.owner,bi=w.bodyIdx,areaC=w.cube.area,kmA=R_KM*R_KM;let tot=0;
    for(let i=0;i<N;i++){const v=rec[i];if(v===0)continue;const arr=this.slotAt.get(own[i]);if(!arr)continue;const s=arr[bi[i]];if(s<0)continue;raw[s]+=v;tot+=v*areaC[i]*kmA;}
    this.budget.tectonicL0+=tot;
    rate.fill(0);
    for(let s=0;s<nS;s++){if(!slotAlive[s]||raw[s]===0)continue;let c=0;for(let k=0;k<4;k++){const u=slotNb[4*s+k];if(u>=0&&slotAlive[u])c++;}
      rate[s]+=raw[s]*(1-0.125*c);for(let k=0;k<4;k++){const u=slotNb[4*s+k];if(u>=0&&slotAlive[u])rate[u]+=raw[s]*0.125;}}
    for(let s=0;s<nS;s++)rate[s]/=dt;
  }
  // one coupled step: L0 advances (recording its tectonics), then L1 applies that forcing and erodes
  advance(){
    const w=this.world;this.rec0.fill(0);w.tecRec=this.rec0;w.step();w.tecRec=null;const dt=w.dt;
    this.updateLineage();this.updateFrames();this.collectForcing(dt);this.step(dt);this.t+=dt;this.steps++;
  }
  // Each piece keeps a frame for every plate it has parcels on: for parcels on plate X, world = Q(X) M b, with b
  // the parcel's position in the piece's own (t0) body frame. A split copies the parent's M to the new plate;
  // a weld of L into W composes M with the fixed rotation between their frames. Two pieces have welded once
  // they have parcels on the same plate; they are then stitched, using that plate's frames.
  updateFrames(){const plates=this.world.plates,pcs=this.pieces;
    for(let id=this.framesSeen;id<plates.length;id++){const C=plates[id];if(!C.alive||C.parent===undefined)continue;
      for(const pc of pcs){const M=pc.frames.get(C.parent);if(M&&!pc.frames.has(id))pc.frames.set(id,M.slice());}}
    this.framesSeen=plates.length;
    for(const pc of pcs)for(const [p,M] of [...pc.frames]){let P=plates[p],Mn=M;
      for(let g=0;g<16&&!P.alive&&P.mergedInto!==undefined;g++){const W=plates[P.mergedInto];Mn=mat(matT(W.Q),mat(P.Q,Mn));P=W;}
      if(P.id!==p){pc.frames.delete(p);if(P.alive&&!pc.frames.has(P.id))pc.frames.set(P.id,Mn);}}
    // where each piece's live parcels are now
    const on=pcs.map(()=>new Map());for(let s=0;s<this.slots.length;s++)if(this.slotAlive[s]){const m=on[this.pieceOf[s]],p=this.slotPlate[s];m.set(p,(m.get(p)||0)+1);}
    this.piecePlates=on;
    for(let i=0;i<pcs.length;i++)for(let j=i+1;j<pcs.length;j++){if(this.welds.some(x=>x.i===i&&x.j===j))continue;
      for(const [p] of on[i])if(on[j].has(p)&&pcs[i].frames.has(p)&&pcs[j].frames.has(p)){this.welds.push({i,j,t:this.t,plate:p,links:this.stitch(i,j,p)});break;}}}
  // Link the outward-facing edge cells of pieces i and j that are on plate p, on the fine grid and the coarse
  // (flexure) grid, to the cell of the other piece one step beyond each edge. Gaps (consumed crust) stay unlinked.
  stitch(i,j,p){const Mi=this.pieces[i].frames.get(p),Mj=this.pieces[j].frames.get(p),Rji=mat(matT(Mi),Mj),Rij=matT(Rji);let links=0;
    const sp=this.slotPlate,sa=this.slotAlive;
    for(const G of [this.G,this.Gm]){const per=G.r*G.r,onP=l=>{const s=(l/per)|0;return sa[s]&&sp[s]===p;};
      links+=stitchGrid(G,i,j,Rij,Rji,onP)+stitchGrid(G,j,i,Rji,Rij,onP);}
    this.buildMidMaps();return links;}
  // The map is drawn in the frame of the plate holding most of the first piece. Every piece's parcels on each
  // plate (a "part") move rigidly; for each part, A maps a map-frame direction into the piece's own frame.
  viewParts(){const plates=this.world.plates,pcs=this.pieces,on=this.piecePlates||pcs.map(pc=>new Map([[pc.plate0,1]]));
    let X0=-1,bn=-1;for(const [q,c] of on[0])if(c>bn&&pcs[0].frames.has(q)){bn=c;X0=q;}if(X0<0)X0=[...pcs[0].frames.keys()][0];
    const F0=mat(plates[X0].Q,pcs[0].frames.get(X0)),out=[];
    pcs.forEach((pc,p)=>{for(const [q] of on[p]){const M=pc.frames.get(q);if(!M)continue;out.push({piece:p,plate:q,A:mat(matT(mat(plates[q].Q,M)),F0)});}});
    return out;}
  welded(){return this.welds.length>=this.pieces.length-1;}
  step(dt){
    const {G,T,cont,alive,eroded,uplifted,S,budget}=this,nC=G.nC,nbr=G.nbr,rr=G.rr,area=G.area;
    const u=this.interp(this.rate,this._u||(this._u=new Float64Array(nC)),this.slotAlive);
    let up=0;
    for(let l=0;l<nC;l++){if(!alive[l]||u[l]===0)continue;const d=u[l]*dt;
      if(cont[l]){T[l]+=d;uplifted[l]+=d;up+=d*area[l];}else if(d>0)budget.unappliedOcean+=d*area[l];}
    budget.tectonicApplied+=up;
    // parcels now on different plates have rifted apart: the cells along the break become sea. Single parcels that
    // flip plate are L0 resampling noise near a boundary, so the break is judged on a majority-filtered plate map.
    const sp=this.coherentPlates();
    for(let l=0;l<nC;l++){if(!alive[l]||!cont[l])continue;const p=sp[(l/rr)|0];
      for(let k=0;k<8;k++){const j=nbr[8*l+k];if(j>=0&&alive[j]&&sp[(j/rr)|0]!==p){budget.rifted+=T[l]*area[l];cont[l]=0;this.eOcean[l]=-0.3;break;}}}
    this.crustFlow(dt);
    let cap=0;for(let l=0;l<nC;l++)if(alive[l]&&cont[l]&&T[l]>C.CAP){cap+=(T[l]-C.CAP)*area[l];T[l]=C.CAP;}
    budget.capped+=cap;
    this.flexure();this.classify();
    route(G,this.e,alive,this.outlet,S,this.precip);
    const hN=streamPower(this.e,S,this.K,this.o.m,dt),e=this.e,outlet=this.outlet;
    let ero=0;
    for(let l=0;l<nC;l++){if(!alive[l]||!cont[l]||outlet[l])continue;const d=e[l]-hN[l];if(d<=0)continue;
      T[l]-=d;eroded[l]+=d;ero+=d*area[l];e[l]=hN[l];}
    budget.eroded+=ero;
    this.stats(ero/dt,up/dt);
  }
  coherentPlates(){const {slotPlate,slotAlive,slotNb}=this,nS=slotPlate.length,out=this._cp||(this._cp=new Int32Array(nS));
    for(let s=0;s<nS;s++){out[s]=slotPlate[s];if(!slotAlive[s])continue;let same=0,n=0,alt=-1,altc=0;
      for(let k=0;k<4;k++){const u=slotNb[4*s+k];if(u<0||!slotAlive[u])continue;n++;if(slotPlate[u]===slotPlate[s])same++;
        else{let c=0;for(let q=0;q<4;q++){const v=slotNb[4*s+q];if(v>=0&&slotAlive[v]&&slotPlate[v]===slotPlate[u])c++;}if(c>altc){altc=c;alt=slotPlate[u];}}}
      if(n>=2&&altc>same&&altc*2>n)out[s]=alt;}
    return out;}
  stats(ero,up){const {G,e,cont,alive}=this,nC=G.nC;let nl=0,se=0,mx=-1e9,al=0;
    for(let l=0;l<nC;l++){if(!alive[l]||!cont[l]||e[l]<=0)continue;nl++;se+=e[l]*G.area[l];al+=G.area[l];if(e[l]>mx)mx=e[l];}
    this.history.push({t:this.t,landKm2:al,meanElevKm:al?se/al:0,maxElevKm:nl?mx:0,erosionKm3PerMyr:ero,upliftKm3PerMyr:up,erosionRateKmPerMyr:al?ero/al:0,upliftRateKmPerMyr:al?up/al:0});}
  // crust volume: what is there now, against what the bookkeeping says should be there
  balance(){const b=this.budget,now=this.volume(),expect=b.initial+b.tectonicApplied+b.flowResidual-b.eroded-b.capped-b.consumed-b.rifted;return{now,expect,error:(now-expect)/b.initial};}
  // coarse-grained L1 thickness against L0's own thickness, over parcels that are continental in both
  divergence(){const w=this.world,{G,T,cont,alive}=this,rr=G.rr,nS=this.slots.length;let n=0,s1=0,s2=0,e1=0,e2=0;
    for(let s=0;s<nS;s++){if(!this.slotAlive[s])continue;const P=w.plates[this.slotPlate[s]],k=this.slotBody[s];if(!P.K[k])continue;
      let v=0,a=0,ok=true;for(let l=s*rr;l<(s+1)*rr;l++){if(!alive[l]||!cont[l]){ok=false;break;}v+=T[l]*G.area[l];a+=G.area[l];}if(!ok)continue;
      const d=v/a-P.T[k],de=C.ISO*d;n++;s1+=d;s2+=d*d;e1+=de;e2+=de*de;}
    return{parcels:n,meanThicknessDiffKm:n?s1/n:0,rmsThicknessDiffKm:n?Math.sqrt(s2/n):0,rmsElevationDiffKm:n?Math.sqrt(e2/n):0};}
  // body-frame direction -> fine cell (or -1)
  // map-frame direction -> fine cell (or -1). Without parts, the map frame is the first piece's own frame.
  // With parts (from viewParts), a cell counts only if its parcel is on that part's plate now.
  localAt(x,y,z,parts){const G=this.G,grids=G.parts||[G],offs=G.offsets||[0],rr=G.rr;
    if(!parts)return grids[0].localOf(grids[0].lookup(x,y,z));let first=-1;
    for(const pt of parts){const R=pt.A,g=grids[pt.piece],X=R[0]*x+R[1]*y+R[2]*z,Y=R[3]*x+R[4]*y+R[5]*z,Z=R[6]*x+R[7]*y+R[8]*z;
      const m=g.localOf(g.lookup(X,Y,Z));if(m<0)continue;const l=m+offs[pt.piece],s=(l/rr)|0;
      if(this.slotPlate[s]===pt.plate&&this.alive[l])return l;if(first<0)first=l;}
    return first;}
  // more than one rigid part (several pieces, or a piece split by a rift): the map changes as they move
  moving(){return (this.piecePlates||[]).reduce((a,m)=>a+m.size,0)>1;}
  basins(){const {S,G}=this,out=this._basin||(this._basin=new Int32Array(G.nC).fill(-1));out.fill(-1);
    for(let h=0;h<S.ns;h++){const l=S.stack[h],r=S.rec[l];out[l]=r===l?l:out[r];}return out;}
}
return {Region,Lookahead,select,selectBody,buildGrid,route,streamPower,planarGrid,scratch,Heap};
})();
// ===== L1-END =====
