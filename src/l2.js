// ===== L2-START =====
// L2: drainage basins selected on the refined L1 map, at 500 m. So far: the surface L2 starts from, as a function of
// position. Any point of it can be evaluated without the rest, which is what burning in basin boundaries in strips
// needs (README, section 6, L2).
// The module is a function of SIM and L1 so that its source can be handed to Web Workers (see L2work below).
function L2module(SIM,L1){
const {hashStr}=SIM;

// 2D gradient noise (Perlin's improved noise, quintic fade), 16 gradient directions, lattice spacing 1; octave k
// draws its gradients from its own hash, so octaves are independent, and is turned by its own angle (k times the
// golden angle): a lattice's axes show faintly in its noise, and octaves sharing them add up to a plaid
const G16=Array.from({length:16},(_,k)=>[Math.cos(2*Math.PI*(k+0.5)/16),Math.sin(2*Math.PI*(k+0.5)/16)]);
function makeGrad2(seed){
  const h=(i,j,k)=>{let n=(Math.imul(i,374761393)+Math.imul(j,668265263)+Math.imul(k,1440662683)+seed)|0;n=Math.imul(n^(n>>>13),1274126177);n^=n>>>16;return G16[n&15];};
  const f=t=>t*t*t*(t*(t*6-15)+10);
  const GA=Math.PI*(3-Math.sqrt(5)),rc=[],rs=[];for(let k=0;k<64;k++){rc.push(Math.cos(k*GA));rs.push(Math.sin(k*GA));}
  return(x0,y0,k)=>{const q=k&63,x=x0*rc[q]-y0*rs[q],y=x0*rs[q]+y0*rc[q],i=Math.floor(x),j=Math.floor(y),u=x-i,v=y-j,a=h(i,j,k),b=h(i+1,j,k),c=h(i,j+1,k),d=h(i+1,j+1,k);
    const n00=a[0]*u+a[1]*v,n10=b[0]*(u-1)+b[1]*v,n01=c[0]*u+c[1]*(v-1),n11=d[0]*(u-1)+d[1]*(v-1),fu=f(u),fv=f(v);
    return n00+(n10-n00)*fu+(n01-n00)*fv+(n00-n10-n01+n11)*fu*fv;};}

// separable Gaussian blur on a grid, over the cells where mask is set (normalized by the mask's own blur); cells with
// no masked cell within reach get 0
function gblur(f,mask,nx,ny,sig){const R=Math.max(1,Math.ceil(3*sig)),w=[];for(let k=-R;k<=R;k++)w.push(Math.exp(-k*k/(2*sig*sig)));
  const a=new Float64Array(nx*ny),am=new Float64Array(nx*ny),b=new Float64Array(nx*ny);
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){let s=0,sm=0;for(let k=-R;k<=R;k++){const ii=i+k;if(ii<0||ii>=nx)continue;const c=j*nx+ii;if(!mask[c])continue;s+=w[k+R]*f[c];sm+=w[k+R];}a[j*nx+i]=s;am[j*nx+i]=sm;}
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){let s=0,sm=0;for(let k=-R;k<=R;k++){const jj=j+k;if(jj<0||jj>=ny)continue;const c=jj*nx+i;s+=w[k+R]*a[c];sm+=w[k+R]*am[c];}b[j*nx+i]=sm>1e-9?s/sm:0;}
  return b;}

// ---------- the relief L2 starts from ----------
// z(x, y) on the L1 sheet's map (km, fixed to the reference plate), the sum of three parts:
// 1. Base: L1's surface (the drainage pass's 10 km grid, rivers and lakes included), Catmull-Rom bicubic, so it
//    is smooth. L1's values are cell means, and an interpolant's cell mean is not the value at the cell centre, so
//    the values it passes through are corrected first (a few iterations of e' += e - M e', M the cell mean of the
//    interpolant, a small separable filter): averaged over each L1 cell, the base is L1's surface.
// 2. Added relief in the band L1 cannot resolve: gradient-noise octaves of lattice spacing lamMaxKm, half that, ...
//    down to lamMinKm (16, 8, 4, 2, 1 km; L1's grid resolves 20 km and up). Each octave's RMS is L1's local roughness
//    r1 continued down the spectrum: r1 (lambda / refKm)^H. r1 is the RMS, over a Gaussian window of windowKm, of
//    L1's land surface less its blur at one cell (relief at about 20-60 km), on land only, so the coast's step from
//    land to L1's flat sea does not count as roughness. H = 0.5: L1's own spectrum falls by about 0.7 per octave
//    (section 5.4), as does Earth's topography at these scales (Huang and Turcotte 1989).
//    Where L1 shortened the crust, the relief is stretched along strike (across the shortening axis) by
//    a = 1 + aniso x shortening, at most anisoMax: ridges and valleys of fold-thrust belts run along strike, and with
//    them the divides (trellis drainage). The direction and amount come from a lattice of nodes latticeKm apart, each
//    with its own rotated, stretched copy of the noise, blended with smooth weights; where the nodes agree (no
//    shortening) the copies are one field, otherwise they are independent and the blend is renormalized to keep the
//    RMS. Noise of lattice spacing 16 km still has power at wavelengths L1 resolves, so the added relief's own L1
//    cell means (8 x 8 samples of the octaves of 4 km and up, per cell, computed when first needed) are interpolated
//    as the base is and subtracted: averaged over L1's cells, the added relief is about zero, and the band is L2's.
// 3. L1's trunk rivers burnt in: along every L1 river draining at least burnKm2, the surface is at most the river's
//    level (L1's filled surface, linear along each 10 km reach) plus (d - floor) x g beyond a floor floorKm either
//    side, d the distance to the river: a narrow trench with a floor about an L2 cell wide (so a 500 m grid always
//    has cells on it, wherever the river runs), so the added relief cannot dam L1's rivers, but does not flatten
//    their valleys. Its sides rise at the
//    steeper of L1's valley side (the steepest slope of L1's surface beside the reach) and burnRise RMS of the added
//    relief over burnWidthKm. L1 knows its rivers' paths only to a cell, and its reaches run in its grid's 0 and 45
//    degree directions, so the river is drawn where the reach is in coordinates warped by a smooth noise (RMS
//    warpKm, lattice warpLamKm): it wanders within L1's cell instead of running as a canal. A river entering the sea
//    is held at sea level, not taken down to L1's sea floor. Small streams are L2's own.
// Octaves shorter than minKm (a pixel, when drawing) are left out. Deterministic in the world seed.
// What the relief needs from L1, as plain data (typed arrays and numbers), so that a worker can build the same
// function: the seed, the drainage pass's surface and rivers, and the shortening lattice (from the sheet).
function reliefInputs(sheet,B,o){
  o=Object.assign({H:0.5,lamMaxKm:16,lamMinKm:1,refKm:40,windowKm:30,aniso:5,anisoMax:4,latticeKm:20,burnKm2:1000,burnRise:3,burnWidthKm:1,burnMaxKm:15,floorKm:0.6,warpKm:1.5,warpLamKm:16},o);
  const {nx,ny,h}=B,seed=hashStr(String(sheet.world.opts.seed)+"|L2relief")|0;
  // shortening, as an axial vector field (s cos 2t, s sin 2t) on the sheet's grid, t the shortening axis, sampled on a
  // lattice latticeKm apart: per node the strike's direction and the stretch along it
  const sf=sheet.strainField(),snc=sheet.nx*sheet.ny,qc=new Float64Array(snc),qs=new Float64Array(snc);
  for(let c=0;c<snc;c++){if(isNaN(sf.lmin[c]))continue;const s=Math.max(0,1-sf.lmin[c]);qc[c]=s*Math.cos(2*sf.dir[c]);qs[c]=s*Math.sin(2*sf.dir[c]);}
  const Lb=o.latticeKm,lnx=Math.ceil(nx*h/Lb)+2,lny=Math.ceil(ny*h/Lb)+2,cs=new Float64Array(lnx*lny),sn=new Float64Array(lnx*lny),an=new Float64Array(lnx*lny).fill(1);
  for(let J=0;J<lny;J++)for(let I=0;I<lnx;I++){const x=B.x0+I*Lb,y=B.y0+J*Lb,a=sheet.sample(qc,x,y),b=sheet.sample(qs,x,y),s=Math.hypot(a,b),k=J*lnx+I;
    const st=Math.min(o.anisoMax,1+o.aniso*s);if(st<1.02){cs[k]=1;sn[k]=0;continue;}
    const strike=0.5*Math.atan2(b,a)+Math.PI/2;cs[k]=Math.cos(strike);sn[k]=Math.sin(strike);an[k]=st;}
  return{o,seed,B:{nx,ny,h,x0:B.x0,y0:B.y0,e:B.e,land:B.land,rec:B.rec,A:B.A,f:B.S.f},lat:{lnx,lny,cs,sn,an}};}
function relief(sheet,B,o){return reliefFrom(reliefInputs(sheet,B,o));}
function reliefFrom(D){
  const o=D.o,B=D.B,{nx,ny,h,e,land}=B,nc=nx*ny,seed=D.seed,noise=makeGrad2(seed);
  // RMS of one noise octave (sampled once on a fine, irregular set of points), to scale octaves to a given RMS
  let ss=0;const NS=40000;for(let k=0;k<NS;k++){const v=noise(k*0.6180339887%97.3,k*0.7548776662%89.1,1000);ss+=v*v;}const nrm=1/Math.sqrt(ss/NS);
  // 1. L1 roughness r1, land only
  const lm=Uint8Array.from(land,(v,c)=>v&&e[c]>0?1:0),hp=new Float64Array(nc),bl=gblur(e,lm,nx,ny,1);
  for(let c=0;c<nc;c++)hp[c]=lm[c]?(e[c]-bl[c])**2:0;
  const r1=gblur(hp,lm,nx,ny,o.windowKm/h);for(let c=0;c<nc;c++)r1[c]=Math.sqrt(r1[c]);
  // 2. the shortening lattice (reliefInputs)
  const Lb=o.latticeKm,lx0=B.x0,ly0=B.y0,{lnx,lny,cs,sn,an}=D.lat;
  // 3. trunk rivers: segments from each cell's centre to its receiver's, bucketed by the cells within burnMaxKm
  const f=B.f,rec=B.rec,A=B.A,gside=new Float64Array(nc),segs=[];
  // (octave amplitudes first: the trench's sides depend on them)
  const nOct=Math.round(Math.log2(o.lamMaxKm/o.lamMinKm))+1,amp=[];for(let k=0;k<nOct;k++)amp.push(nrm*Math.pow(o.lamMaxKm/2**k/o.refKm,o.H));
  const sigma=Math.sqrt(amp.reduce((a,v)=>a+(v/nrm)**2,0)); // RMS of the added relief per unit r1 (octaves independent)
  for(let c=0;c<nc;c++){if(!land[c])continue;const i=c%nx,j=(c/nx)|0;let g=o.burnRise*sigma*r1[c]/o.burnWidthKm;
    for(let dj=-1;dj<=1;dj++)for(let di=-1;di<=1;di++){const I=i+di,J=j+dj;if((di||dj)&&I>=0&&J>=0&&I<nx&&J<ny){const d=e[J*nx+I]-e[c];if(d>0)g=Math.max(g,d/(h*Math.hypot(di,dj)));}}gside[c]=g;}
  for(let c=0;c<nc;c++){const r=rec[c];if(!land[c]||r===c||r<0||!(A[c]>=o.burnKm2))continue;segs.push(c);}
  const reach=Math.ceil(o.burnMaxKm/h),cnt=new Int32Array(nc+1);
  const cover=(c,fn)=>{const r=rec[c],i0=Math.min(c%nx,r%nx)-reach,i1=Math.max(c%nx,r%nx)+reach,j0=Math.min((c/nx)|0,(r/nx)|0)-reach,j1=Math.max((c/nx)|0,(r/nx)|0)+reach;
    for(let J=Math.max(0,j0);J<=Math.min(ny-1,j1);J++)for(let I=Math.max(0,i0);I<=Math.min(nx-1,i1);I++)fn(J*nx+I);};
  for(const c of segs)cover(c,q=>cnt[q+1]++);for(let q=0;q<nc;q++)cnt[q+1]+=cnt[q];
  const bucket=new Int32Array(cnt[nc]),fill=cnt.slice(0,nc);for(const c of segs)cover(c,q=>{bucket[fill[q]++]=c;});
  const cx=c=>B.x0+(c%nx+.5)*h,cy=c=>B.y0+(((c/nx)|0)+.5)*h;
  // samplers on the 10 km grid
  const at=(F,i,j)=>F[Math.max(0,Math.min(ny-1,j))*nx+Math.max(0,Math.min(nx-1,i))];
  const cr=(p0,p1,p2,p3,t)=>p1+0.5*t*(p2-p0+t*(2*p0-5*p1+4*p2-p3+t*(3*(p1-p2)+p3-p0)));
  const cubic=(F,x,y)=>{const gx=(x-B.x0)/h-0.5,gy=(y-B.y0)/h-0.5,i=Math.floor(gx),j=Math.floor(gy),tx=gx-i,ty=gy-j,r=[0,0,0,0];
    for(let m=0;m<4;m++)r[m]=cr(at(F,i-1,j-1+m),at(F,i,j-1+m),at(F,i+1,j-1+m),at(F,i+2,j-1+m),tx);return cr(r[0],r[1],r[2],r[3],ty);};
  // the cell mean of the 1D interpolant of a unit impulse at offsets -2..2: a separable 5-point filter
  const ker=[-2,-1,0,1,2].map(k=>{let s=0;const NK=400;for(let q=0;q<NK;q++){const u=(q+.5)/NK-0.5,i=Math.floor(u),t=u-i,p=d=>i+d===k?1:0;s+=cr(p(-1),p(0),p(1),p(2),t);}return s/NK;});
  const conv=(F,dx_,dy_)=>{const o2=new Float64Array(nc);for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){let s=0;for(let k=-2;k<=2;k++)s+=ker[k+2]*at(F,i+k*dx_,j+k*dy_);o2[j*nx+i]=s;}return o2;};
  const eb=Float64Array.from(e);for(let it=0;it<4;it++){const m=conv(conv(eb,1,0),0,1);for(let c=0;c<nc;c++)eb[c]+=e[c]-m[c];}
  const base=(x,y)=>cubic(eb,x,y);
  function bil(F,x,y){const gx=(x-B.x0)/h-0.5,gy=(y-B.y0)/h-0.5,i=Math.floor(gx),j=Math.floor(gy),tx=gx-i,ty=gy-j;
    return(1-tx)*(1-ty)*at(F,i,j)+tx*(1-ty)*at(F,i+1,j)+(1-tx)*ty*at(F,i,j+1)+tx*ty*at(F,i+1,j+1);}
  // the added relief at (x, y) in units of r1: octaves of lattice spacing lamMaxKm / 2^k, in each node's frame
  function band(x,y,minKm){let kMax=nOct;while(kMax>0&&o.lamMaxKm/2**(kMax-1)<minKm)kMax--;if(!kMax)return 0;
    const gx=(x-lx0)/Lb,gy=(y-ly0)/Lb,I=Math.floor(gx),J=Math.floor(gy),fx=gx-I,fy=gy-J,sx=fx*fx*(3-2*fx),sy=fy*fy*(3-2*fy);
    let tot=0,w2=0,iso=0;
    for(let q=0;q<4;q++){const ii=Math.max(0,Math.min(lnx-1,I+(q&1))),jj=Math.max(0,Math.min(lny-1,J+(q>>1))),k=jj*lnx+ii,w=(q&1?sx:1-sx)*(q>>1?sy:1-sy);if(w<=0)continue;
      if(an[k]===1){iso+=w;continue;}
      const u=(x*cs[k]+y*sn[k])/an[k],v=-x*sn[k]+y*cs[k];let s=0;for(let o2=0;o2<kMax;o2++){const L=o.lamMaxKm/2**o2;s+=amp[o2]*noise(u/L,v/L,o2);}tot+=w*s;w2+=w*w;}
    if(iso>0){let s=0;for(let o2=0;o2<kMax;o2++){const L=o.lamMaxKm/2**o2;s+=amp[o2]*noise(x/L,y/L,o2);}tot+=iso*s;w2+=iso*iso;}
    return tot/Math.sqrt(w2);}
  // the burnt-in trunk rivers: the lowest trench over the reaches near (x, y), or Infinity
  const lev=c=>land[c]?f[c]:0,wa=o.warpKm*nrm,wl=o.warpLamKm;
  function trench(x0,y0){const x=x0+wa*noise(x0/wl,y0/wl,50),y=y0+wa*noise(x0/wl,y0/wl,51),i=Math.floor((x-B.x0)/h),j=Math.floor((y-B.y0)/h);if(i<0||j<0||i>=nx||j>=ny)return Infinity;const q=j*nx+i;let best=Infinity;
    for(let p=cnt[q];p<cnt[q+1];p++){const c=bucket[p],r=rec[c],ax=cx(c),ay=cy(c),bx=cx(r),by=cy(r),dx=bx-ax,dy=by-ay,L2=dx*dx+dy*dy;
      const t=Math.max(0,Math.min(1,((x-ax)*dx+(y-ay)*dy)/L2)),d=Math.hypot(x-ax-t*dx,y-ay-t*dy);if(d>o.burnMaxKm)continue;
      const z=lev(c)+t*(lev(r)-lev(c))+Math.max(0,d-o.floorKm)*gside[c];if(z<best)best=z;}
    return best;}
  // the added relief's L1 cell means, when first needed; then the added relief less their interpolation
  const cm=new Float64Array(nc).fill(NaN),NM=8;
  const cellMean=(i,j)=>{i=Math.max(0,Math.min(nx-1,i));j=Math.max(0,Math.min(ny-1,j));const c=j*nx+i;if(!isNaN(cm[c]))return cm[c];let s=0;
    for(let b=0;b<NM;b++)for(let a=0;a<NM;a++){const x=B.x0+(i+(a+.5)/NM)*h,y=B.y0+(j+(b+.5)/NM)*h;s+=bil(r1,x,y)*band(x,y,4);}return cm[c]=s/(NM*NM);};
  function added(x,y,minKm=0){const gx=(x-B.x0)/h-0.5,gy=(y-B.y0)/h-0.5,i=Math.floor(gx),j=Math.floor(gy),tx=gx-i,ty=gy-j,r=[0,0,0,0];
    for(let m=0;m<4;m++)r[m]=cr(cellMean(i-1,j-1+m),cellMean(i,j-1+m),cellMean(i+1,j-1+m),cellMean(i+2,j-1+m),tx);
    return bil(r1,x,y)*band(x,y,minKm)-cr(r[0],r[1],r[2],r[3],ty);}
  // the surface; minKm leaves out octaves shorter than that (default: all)
  function z(x,y,minKm=0){return Math.min(base(x,y)+added(x,y,minKm),trench(x,y));}
  // the same on a grid: nx2 x ny2 cells of dx km from (x0, y0) (cell centres), octaves down to 2 dx
  function grid(x0,y0,dx,nx2,ny2,minKm=2*dx){const out=new Float32Array(nx2*ny2);for(let j=0;j<ny2;j++)for(let i=0;i<nx2;i++)out[j*nx2+i]=z(x0+(i+.5)*dx,y0+(j+.5)*dx,minKm);return out;}
  return{o,z,grid,base,band,added,trench,r1,roughAt:(x,y)=>bil(r1,x,y),strikeAt:(x,y)=>{const I=Math.round((x-lx0)/Lb),J=Math.round((y-ly0)/Lb);if(I<0||J<0||I>=lnx||J>=lny)return null;const k=J*lnx+I;return{stretch:an[k],strike:Math.atan2(sn[k],cs[k])};},
    nSegments:segs.length,octaves:Array.from({length:nOct},(_,k)=>o.lamMaxKm/2**k)};
}
// ---------- basin boundaries at 500 m ----------
// Which L1 basin a point of L2's surface belongs to: the basin (at L1) of the place where its water leaves the land, at
// the sea, at an inland sink (an L1 lake with no outlet, a dry lake floor or a river's dry end; or land below sea
// level), or at the edge of the window routed. Depressions are filled (priority flood, as at L1), so every point
// drains somewhere. On a window: flood from those outlets, lowest first; each cell takes the label of the cell it was
// reached from, so it carries the label of the outlet its water reaches (Barnes et al. 2014).
// L1 sea cells take the basin of the nearest L1 land cell, so a coast that L2 moves seaward keeps its basins.
function l1Labels(B){const {nx,ny,land,basin}=B,nc=nx*ny,lab=new Int32Array(nc).fill(-1),q=[];
  for(let c=0;c<nc;c++)if(land[c]){lab[c]=basin[c];q.push(c);}
  for(let k=0;k<q.length;k++){const c=q[k],i=c%nx,j=(c/nx)|0;for(const [di,dj] of [[1,0],[-1,0],[0,1],[0,-1]]){const I=i+di,J=j+dj;if(I<0||J<0||I>=nx||J>=ny)continue;const n=J*nx+I;if(lab[n]<0){lab[n]=lab[c];q.push(n);}}}
  return lab;}
// labels on a window of W x H cells with heights z: seeds (lab0 >= 0) are outlets; returns each cell's label, and
// the seed (outlet cell) it came from
function flood(z,W,H,lab0){const n=W*H,src=new Int32Array(n).fill(-1),key=new Float64Array(n),done=new Uint8Array(n),hp=new L1.Heap(9*n+8);
  for(let q=0;q<n;q++)if(lab0[q]>=0){src[q]=q;key[q]=z[q];hp.push(q,z[q]);}
  while(hp.n){const q=hp.pop();if(done[q])continue;done[q]=1;const i=q%W,j=(q/W)|0,kq=key[q];
    for(let dj=-1;dj<=1;dj++)for(let di=-1;di<=1;di++){if(!di&&!dj)continue;const I=i+di,J=j+dj;if(I<0||J<0||I>=W||J>=H)continue;const p=J*W+I;
      if(done[p]||lab0[p]>=0)continue;const kp=Math.max(z[p],kq);if(src[p]<0||kp<key[p]){src[p]=src[q];key[p]=kp;hp.push(p,kp);}}}
  const lab=new Int32Array(n);for(let q=0;q<n;q++)lab[q]=src[q]>=0?lab0[src[q]]:-1;
  return{lab,src};}
// A flood with given seeds: seed[q] set for outlets, each with its starting level key0[q], label lab0[q] and origin
// org0[q]; every other point takes the label and origin of the point it was reached from, at the level
// max(its height, that point's level). Returns lab, org and the level key.
function floodSeeded(z,W,H,seed,key0,lab0,org0){const n=W*H,key=new Float64Array(n).fill(Infinity),lab=new Int32Array(n).fill(-1),org=new Int32Array(n).fill(-1),done=new Uint8Array(n),hp=new L1.Heap(9*n+8);
  for(let q=0;q<n;q++)if(seed[q]){key[q]=key0[q];lab[q]=lab0[q];org[q]=org0[q];hp.push(q,key[q]);}
  while(hp.n){const q=hp.pop();if(done[q])continue;done[q]=1;const i=q%W,j=(q/W)|0,kq=key[q];
    for(let dj=-1;dj<=1;dj++)for(let di=-1;di<=1;di++){if(!di&&!dj)continue;const I=i+di,J=j+dj;if(I<0||J<0||I>=W||J>=H)continue;const p=J*W+I;
      if(done[p]||seed[p])continue;const kp=Math.max(z[p],kq);if(kp<key[p]){key[p]=kp;lab[p]=lab[q];org[p]=org[q];hp.push(p,kp);}}}
  return{lab,org,key};}
// Pass A on one tile of W x W points (a pure function of its heights and outlets, so it can run anywhere): flood from
// every edge point and every outlet point, each its own seed, and return the graph between seeds: for each pair whose
// catchments touch, the lowest level at which water crosses between them (the greater of two touching points'
// levels, least over all touching pairs). Minimax ("lowest highest point") distances between seeds are the same on
// this graph as on the points (Barnes et al. 2014, 2016).
function tileGraph(z,W,isOut){const n=W*W,seed=new Uint8Array(n),key0=new Float64Array(n),ids=new Int32Array(n);
  const per=[];for(let i=0;i<W;i++)per.push(i);for(let j=1;j<W;j++)per.push(j*W+W-1);for(let i=W-2;i>=0;i--)per.push((W-1)*W+i);for(let j=W-2;j>=1;j--)per.push(j*W);
  for(const q of per)seed[q]=1;for(let q=0;q<n;q++)if(isOut[q])seed[q]=1;
  for(let q=0;q<n;q++){key0[q]=z[q];ids[q]=q;}
  const {org,key}=floodSeeded(z,W,W,seed,key0,ids,ids),best=new Map();
  for(let j=0;j<W;j++)for(let i=0;i<W;i++){const a=j*W+i,sa=org[a];
    for(const [di,dj] of [[1,0],[0,1],[1,1],[-1,1]]){const I=i+di,J=j+dj;if(I<0||J<0||I>=W||J>=W)continue;const b=J*W+I,sb=org[b];if(sa===sb||isOut[sa]&&isOut[sb])continue;
      const w=Math.max(key[a],key[b]),k=sa<sb?sa*n+sb:sb*n+sa,o=best.get(k);if(o===undefined||w<o)best.set(k,w);}}
  const eu=new Int32Array(best.size),ev=new Int32Array(best.size),ew=new Float64Array(best.size);let m=0;for(const [k,w] of best){eu[m]=Math.floor(k/n);ev[m]=k%n;ew[m]=w;m++;}
  return{per:Int32Array.from(per),eu,ev,ew};}
// Pass A for one tile from scratch: its heights (the relief, evaluated), outlets (sea, land below sea level, L1 sinks;
// outside the grid) with their L1 labels, then its graph. Shared by the page's main thread and its workers, which
// build it from plain data (reliefInputs, and the drainage pass's land, lakes and basins).
function tileMaker(R,B,o){const {nx,ny,h,x0,y0,land,lake}=B,f=Math.round(h/o.dx),dx=h/f,W=o.T*f,lab1=l1Labels(B),sink=c=>land[c]&&lake&&lake[c]>=2;
  function make(TI,TJ){const z=new Float32Array(W*W),out=new Uint8Array(W*W),lab=new Int32Array(W*W).fill(-1),cell=new Int32Array(W*W).fill(-1);
    for(let bj=0;bj<o.T;bj++)for(let bi=0;bi<o.T;bi++){const I=TI*o.T+bi,J=TJ*o.T+bj,inG=I<nx&&J<ny,c=inG?J*nx+I:-1;
      for(let v=0;v<f;v++)for(let u=0;u<f;u++){const q=(bj*f+v)*W+bi*f+u;if(!inG){z[q]=-1;out[q]=1;continue;}
        z[q]=R.z(x0+I*h+(u+.5)*dx,y0+J*h+(v+.5)*dx,R.o.lamMinKm);cell[q]=c;lab[q]=lab1[c];if(z[q]<=0||sink(c))out[q]=1;}}
    return Object.assign({z,out,lab,cell,TI,TJ},tileGraph(z,W,out));}
  return{make,lab1,f,dx,W,sink};}
// The selection's outline at L2's resolution, stitched from tiles. Tiles are T x T L1 cells on a fixed lattice;
// each is evaluated and reduced to its graph (pass A) once. For a selection, a domain of tiles is routed:
//   - stitch: the tiles' graphs joined across their shared edges (touching edge points, crossing level the greater
//     of the two heights) and flooded as one graph, from the outlets and from the domain's outer edge (labelled by
//     L1); this gives the exact level, label and origin arriving at every tile's edge points;
//   - pass B: each tile flooded again from its outlets and its edge points at those levels: the same labels as
//     routing the whole domain at once (ties aside).
// A label is a guess if its water reached the domain's outer edge where L1's label is not sure (an L1 cell within a
// cell of another basin). The domain starts as the tiles of the L1 cells within `band` cells of the selection's outer
// edge and grows by the tiles beyond the edge points that guesses came from, until no point in those L1 cells is a
// guess (or maxTiles). Every L1 cell of the domain with no guess gets L2's labels; the rest keep L1's.
function boundary(R,B,o){
  o=Object.assign({dx:0.5,T:4,band:2,sureR:1,grow:0,maxTiles:4000},o);const stats=[];
  const {nx,ny,h,land}=B,M=tileMaker(R,B,o),{f,dx,W,lab1,sink}=M,graphs=new Map(),tnx=Math.ceil(nx/o.T),tny=Math.ceil(ny/o.T);
  // L1 cells whose label is sure: every L1 cell within sureR cells shares it
  const sure=new Uint8Array(nx*ny);for(let c=0;c<nx*ny;c++){const i=c%nx,j=(c/nx)|0,d=o.sureR;let ok=1;for(let dj=-d;dj<=d&&ok;dj++)for(let di=-d;di<=d;di++){const I=i+di,J=j+dj;if(I>=0&&J>=0&&I<nx&&J<ny&&lab1[J*nx+I]!==lab1[c]){ok=0;break;}}sure[c]=ok;}
  // pass A, once per tile (here, or by a worker: putGraph)
  const graph=k=>{let g=graphs.get(k);if(!g){g=M.make(Math.floor(k/tny),k%tny);graphs.set(k,g);}return g;};
  const putGraph=(k,g)=>{if(!graphs.has(k))graphs.set(k,g);};
  const tileOfCell=c=>Math.floor((c%nx)/o.T)*tny+Math.floor(((c/nx)|0)/o.T);
  // the selection as a flag per L1 cell
  const selFlags=sel=>{const inSel=new Set(sel.basins),fl=new Uint8Array(nx*ny);for(let c=0;c<nx*ny;c++)if(lab1[c]>=0&&inSel.has(lab1[c]))fl[c]=1;return fl;};
  // stitch and pass B over a domain (a set of tile keys); returns per tile the labels and whether each is a guess,
  // and the domain-edge points the guesses came from
  // final: pass B over every tile and the edge relabelling; otherwise (while the domain grows) pass B only over the
  // tiles in `only`, enough to find the guesses
  function solve(dom,final,only){const T0=Date.now(),keys=[...dom],gi=new Map();keys.forEach((k,t)=>gi.set(k,t));
    // global node per (tile, seed point): edge points, then the outlets the tile's graph mentions
    const off=[],nodeOf=[];let N=0;for(const k of keys){const g=graphs.get(k),m=new Map();for(const q of g.per)m.set(q,N++);for(let e=0;e<g.eu.length;e++)for(const q of [g.eu[e],g.ev[e]])if(!m.has(q))m.set(q,N++);nodeOf.push(m);}
    const nz=new Float64Array(N),seed=new Uint8Array(N),nlab=new Int32Array(N).fill(-1),guessy=new Uint8Array(N),ntile=new Int32Array(N),npt=new Int32Array(N),U=[],V=[],Wt=[];
    keys.forEach((k,t)=>{const g=graphs.get(k);for(const [q,id] of nodeOf[t]){nz[id]=g.z[q];ntile[id]=t;npt[id]=q;if(g.out[q]){seed[id]=1;nlab[id]=g.lab[q];}}
      for(let e=0;e<g.eu.length;e++){U.push(nodeOf[t].get(g.eu[e]));V.push(nodeOf[t].get(g.ev[e]));Wt.push(g.ew[e]);}
      // across the tile's edges: touching edge points of a neighbouring tile in the domain, or the domain's outer edge
      const {TI,TJ}=g;for(const q of g.per){const u=q%W,v=(q/W)|0,X=TI*W+u,Y=TJ*W+v,a=nodeOf[t].get(q);let outer=false;
        for(let dj=-1;dj<=1;dj++)for(let di=-1;di<=1;di++){if(!di&&!dj)continue;const u2=u+di,v2=v+dj;if(u2>=0&&v2>=0&&u2<W&&v2<W)continue;
          const X2=X+di,Y2=Y+dj,T2I=Math.floor(X2/W),T2J=Math.floor(Y2/W),k2=T2I*tny+T2J;
          if(X2<0||Y2<0||T2I>=tnx||T2J>=tny||!gi.has(k2)){outer=true;continue;}
          if(k2<k)continue;const g2=graphs.get(k2),q2=(Y2-T2J*W)*W+(X2-T2I*W),b=nodeOf[gi.get(k2)].get(q2);U.push(a);V.push(b);Wt.push(Math.max(g.z[q],g2.z[q2]));}
        if(outer&&!seed[a]){seed[a]=1;nlab[a]=g.lab[q];const c=g.cell[q];guessy[a]=c>=0&&!sure[c]&&!sink(c)?1:0;}}});
    const T1=Date.now();
    // the stitched graph, flooded
    const deg=new Int32Array(N+1);for(let e=0;e<U.length;e++){deg[U[e]+1]++;deg[V[e]+1]++;}for(let i=0;i<N;i++)deg[i+1]+=deg[i];
    const adj=new Int32Array(deg[N]),aw=new Float64Array(deg[N]),fp=deg.slice(0,N);for(let e=0;e<U.length;e++){adj[fp[U[e]]]=V[e];aw[fp[U[e]]++]=Wt[e];adj[fp[V[e]]]=U[e];aw[fp[V[e]]++]=Wt[e];}
    const K=new Float64Array(N).fill(Infinity),L=new Int32Array(N).fill(-1),O=new Int32Array(N).fill(-1),dn=new Uint8Array(N),hp=new L1.Heap(deg[N]+N+8);
    for(let i=0;i<N;i++)if(seed[i]){K[i]=nz[i];L[i]=nlab[i];O[i]=i;hp.push(i,K[i]);}
    while(hp.n){const a=hp.pop();if(dn[a])continue;dn[a]=1;for(let p=deg[a];p<deg[a+1];p++){const b=adj[p];if(dn[b]||seed[b])continue;const k=Math.max(aw[p],K[a]);if(k<K[b]){K[b]=k;L[b]=L[a];O[b]=O[a];hp.push(b,k);}}}
    const T2=Date.now();
    // pass B per tile: edge points seeded at the stitched levels
    const res=new Map();
    const passB=(k,t)=>{const g=graphs.get(k),sd=new Uint8Array(W*W),k0=new Float64Array(W*W),l0=new Int32Array(W*W),o0=new Int32Array(W*W).fill(-1);
      for(let q=0;q<W*W;q++)if(g.out[q]){sd[q]=1;k0[q]=g.z[q];l0[q]=g.lab[q];}
      for(const q of g.per){const id=nodeOf[t].get(q);sd[q]=1;k0[q]=K[id];l0[q]=L[id];o0[q]=O[id];}
      const r=floodSeeded(g.z,W,W,sd,k0,l0,o0),gs=new Uint8Array(W*W);for(let q=0;q<W*W;q++){const og=r.org[q];if(og>=0&&guessy[og])gs[q]=1;}
      res.set(k,{lab:r.lab,org:r.org,key:r.key,guess:gs});};
    keys.forEach((k,t)=>{if(final||only.has(k))passB(k,t);});const T3=Date.now();let iters=0,reflooded=0;
    // The stitched levels are exact, but at an edge point reached at the same level from several sides, the graph
    // chose a side by its own order. Routing the whole domain gives such a point the label of its lowest neighbour
    // (where its water goes), so edge points take it, across seams too; tiles whose edge labels changed are flooded
    // again, until none change.
    for(let it=0;final&&it<50;it++){const dirty=new Set();
      keys.forEach((k,t)=>{const g=graphs.get(k),r=res.get(k);for(const q of g.per){const id=nodeOf[t].get(q);if(seed[id])continue;const u=q%W,v=(q/W)|0,X=g.TI*W+u,Y=g.TJ*W+v;
        let bk=Infinity,bl=-1,bo=-1;for(let dj=-1;dj<=1;dj++)for(let di=-1;di<=1;di++){if(!di&&!dj)continue;const X2=X+di,Y2=Y+dj,T2I=Math.floor(X2/W),T2J=Math.floor(Y2/W),k2=T2I*tny+T2J;
          if(X2<0||Y2<0||T2I>=tnx||T2J>=tny)continue;const r2=res.get(k2);if(!r2)continue;const q2=(Y2-T2J*W)*W+X2-T2I*W;if(r2.key[q2]<bk){bk=r2.key[q2];bl=r2.lab[q2];bo=r2.org[q2];}}
        if(bk<K[id]&&bl!==L[id]){L[id]=bl;O[id]=bo;dirty.add(t);}}});
      iters++;if(!dirty.size)break;reflooded+=dirty.size;for(const t of dirty)passB(keys[t],t);}
    stats.push({tiles:keys.length,nodes:N,edges:U.length,build:T1-T0,graph:T2-T1,passB:T3-T2,relabel:Date.now()-T3,iters,reflooded});
    return{res,guessy,ntile,npt,keys};}
  // the domain-edge point a guess came from, and the tiles beyond it (and within `grow` tiles of those: growing faster
  // took fewer rounds but a larger domain, and cost more)
  function beyond(S,og,add){const g=graphs.get(S.keys[S.ntile[og]]),q=S.npt[og],u=q%W,v=(q/W)|0,X=g.TI*W+u,Y=g.TJ*W+v,G=o.grow;
    for(let dj=-1;dj<=1;dj++)for(let di=-1;di<=1;di++){const T2I=Math.floor((X+di)/W),T2J=Math.floor((Y+dj)/W);
      for(let b=-G;b<=G;b++)for(let a=-G;a<=G;a++){const I=T2I+a,J=T2J+b;if(I<0||J<0||I>=tnx||J>=tny)continue;add(I*tny+J);}}}
  // The selection refined, as a generator: yields after each tile evaluated, so the page can spread the work over
  // frames; returns the result. Tiles evaluated before (for any selection) are reused.
  // Yields {need: tile keys} whenever tiles must be evaluated (pass A); whoever drives it may supply them with putGraph
  // before resuming (workers), or leave them to be evaluated here.
  function* refineSteps(sel){const t0=Date.now(),fl=selFlags(sel),selC=c=>fl[c]===1;
    // core: L1 cells within `band` of an edge cell (selected with an unselected neighbour on land, or the reverse)
    const core=new Set();for(let c=0;c<nx*ny;c++){const i=c%nx,j=(c/nx)|0;let on=false;
      for(let dj=-1;dj<=1&&!on;dj++)for(let di=-1;di<=1;di++){const I=i+di,J=j+dj;if(I<0||J<0||I>=nx||J>=ny)continue;if(fl[J*nx+I]!==fl[c]&&(land[c]||land[J*nx+I])){on=true;break;}}
      if(on){const d=o.band;for(let dj=-d;dj<=d;dj++)for(let di=-d;di<=d;di++){const I=i+di,J=j+dj;if(I>=0&&J>=0&&I<nx&&J<ny)core.add(J*nx+I);}}}
    const dom=new Set();for(const c of core)dom.add(tileOfCell(c));const coreTiles=new Set(dom);let S=null,rounds=0,evaluated=0;
    const fill=function*(){const need=[...dom].filter(k=>!graphs.has(k));if(need.length){evaluated+=need.length;yield{need,evaluated,domain:dom.size,rounds};for(const k of need)graph(k);}};
    for(;;){rounds++;yield*fill();
      S=solve(dom,false,coreTiles);const add=new Set();
      for(const c of core){const k=tileOfCell(c),r=S.res.get(k),bi=(c%nx)-Math.floor(k/tny)*o.T,bj=((c/nx)|0)-(k%tny)*o.T;
        for(let v=0;v<f;v++)for(let u=0;u<f;u++){const q=(bj*f+v)*W+bi*f+u;if(r.guess[q])beyond(S,r.org[q],kk=>{if(!dom.has(kk))add.add(kk);});}}
      if(!add.size||dom.size>=o.maxTiles)break;for(const k of add){if(dom.size>=o.maxTiles)break;dom.add(k);}}
    S=solve(dom,true);
    // labels per L1 cell: the core always, and any other L1 cell of the domain with no guess
    const cells=new Map();let guessCells=0;
    for(const k of dom){const r=S.res.get(k),TI=Math.floor(k/tny),TJ=k%tny;for(let bj=0;bj<o.T;bj++)for(let bi=0;bi<o.T;bi++){const I=TI*o.T+bi,J=TJ*o.T+bj;if(I>=nx||J>=ny)continue;const c=J*nx+I,Lc=new Int32Array(f*f);let g=0;
      for(let v=0;v<f;v++)for(let u=0;u<f;u++){const q=(bj*f+v)*W+bi*f+u;Lc[v*f+u]=r.lab[q];g+=r.guess[q];}
      if(core.has(c)){cells.set(c,Lc);if(g)guessCells++;}else if(!g)cells.set(c,Lc);}}
    const inSel=new Set(sel.basins);
    const inside=(x,y)=>{const i=Math.floor((x-B.x0)/h),j=Math.floor((y-B.y0)/h);if(i<0||j<0||i>=nx||j>=ny)return false;const c=j*nx+i,Lc=cells.get(c);
      if(!Lc)return selC(c);const u=Math.min(f-1,Math.floor((x-B.x0-i*h)/dx)),v=Math.min(f-1,Math.floor((y-B.y0-j*h)/dx));return inSel.has(Lc[v*f+u]);};
    // land area: L1 cells refined count their land points in the selection; the rest count whole if selected
    let area=0,areaL1=0;for(let c=0;c<nx*ny;c++){if(!land[c])continue;if(selC(c))areaL1+=h*h;const Lc=cells.get(c);if(!Lc){if(selC(c))area+=h*h;continue;}const g=graphs.get(tileOfCell(c)),bi=(c%nx)-Math.floor(tileOfCell(c)/tny)*o.T,bj=((c/nx)|0)-(tileOfCell(c)%tny)*o.T;let n=0;for(let v=0;v<f;v++)for(let u=0;u<f;u++)if(inSel.has(Lc[v*f+u])&&g.z[(bj*f+v)*W+bi*f+u]>0)n++;area+=n*dx*dx;}
    return{inside,cells,solved:S,tileKey:(x,y)=>{const X=Math.floor((x-B.x0)/dx),Y=Math.floor((y-B.y0)/dx),TI=Math.floor(X/W),TJ=Math.floor(Y/W),k=TI*tny+TJ,r=S.res.get(k);if(!r)return null;const q=(Y-TJ*W)*W+X-TI*W;return{lab:r.lab[q],key:r.key[q],u:X-TI*W,v:Y-TJ*W};},nTiles:dom.size,evaluated,rounds,nCells:cells.size,coreCells:core.size,guessCells,areaKm2:area,areaL1Km2:areaL1,ms:Date.now()-t0};}
  function refine(sel){const it=refineSteps(sel);let r;while(!(r=it.next()).done);return r.value;}
  // the inputs a worker needs to build the same tiles: the relief's inputs (passed in), the grid and the outlets
  const workerInputs=()=>({B:{nx,ny,h,x0:B.x0,y0:B.y0,land:B.land,lake:B.lake,basin:B.basin},o});
  // the same labelling routed over a whole box at once (the reference the stitching is tested against): nX x nY
  // points of dx from (x0, y0), the box's edge an outlet with L1's labels
  function full(x0,y0,nX,nY){const z=R.grid(x0,y0,dx,nX,nY,R.o.lamMinKm),l0=new Int32Array(nX*nY).fill(-1);
    for(let j=0;j<nY;j++)for(let i=0;i<nX;i++){const q=j*nX+i,x=x0+(i+.5)*dx,y=y0+(j+.5)*dx,ci=Math.floor((x-B.x0)/h),cj=Math.floor((y-B.y0)/h);
      if(ci<0||cj<0||ci>=nx||cj>=ny){l0[q]=-2;continue;}const c=cj*nx+ci;if(i===0||j===0||i===nX-1||j===nY-1||z[q]<=0||sink(c))l0[q]=lab1[c];}
    return{z,lab:flood(z,nX,nY,l0).lab};}
  return{o,dx,f,lab1,tny,refine,refineSteps,putGraph,ensure:graph,workerInputs,full,graphs,stats};}
// ---------- the L2 pass ----------
// Erodibility from rock, per sheet cell (as tests/exp_erodibility.js grounded it): a factor by the rock at the surface
// (sediments weak, cratonic basement strong), young sediments weaker still (up to 10x when fresh, e-folding over
// 30 Myr), log-averaged over the cell's markers by volume and scaled to a geometric mean of 1 on land.
const KTYPE=[1,0.7,0.7,1.5,0.5,0.7,3,3]; // none, ocean basalt, hotspot basalt, arc, craton, orogen, shelf, clastic
function erodibility(sheet){const {nx,ny,o}=sheet,h=o.h,nc=nx*ny,s=new Float64Array(nc),ws=new Float64Array(nc),t=sheet.world.t,R=SIM.ROCK;
  for(let k=0;k<sheet.nm;k++){if(!sheet.mIn[k])continue;const col=sheet.mCol[k];if(!col.length)continue;const u=col[col.length-1],i=Math.floor((sheet.mx[k]-sheet.x0)/h),j=Math.floor((sheet.my[k]-sheet.y0)/h),c=j*nx+i;
    const young=u[0]===R.SHELF||u[0]===R.CLASTIC?1+9*Math.exp(-(t-u[1])/30):1,v=sheet.mS[k]*sheet.mA[k];s[c]+=v*Math.log((KTYPE[u[0]]||1)*young);ws[c]+=v;}
  let g=0,n=0;for(let c=0;c<nc;c++)if(ws[c]>0&&sheet.act[c]&&sheet.e[c]>0){g+=s[c]/ws[c];n++;}g/=n||1;
  const lnK=new Float64Array(nc);for(let c=0;c<nc;c++)lnK[c]=ws[c]>0?s[c]/ws[c]-g:0;
  return{lnK,nx,ny,h,x0:sheet.x0,y0:sheet.y0};}
// node positions jittered within their cells by a fixed hash of the global point index (as L1's drainage pass), so
// steepest descent is not locked to the grid's 0 and 45 degree directions
function jit(gx,gy,salt){let h=Math.imul(gx,0x27d4eb2d)^Math.imul(gy,0x165667b1)^salt;h=Math.imul(h^(h>>>16),0x7feb352d);h=Math.imul(h^(h>>>15),0x846ca68b);h^=h>>>16;return[(h&0xffff)/65536-0.5,(h>>>16)/65536-0.5];}
// The L2 pass on a selection whose outline is burnt in (bd: its boundary, ref: refine's result): landscape evolution at
// 500 m for `myr` at present rates, a generator yielding after each step.
// - Domain: the selected land points. Outlets: sea (land below sea level included) at sea level, and L1's inland sinks
//   (lakes with no outlet, dry lake floors, rivers' dry ends) at their own height, as in the burn-in. The rest of the
//   outline is a water divide held fixed.
// - Each step: uplift at L1's present tectonic rate; routing (priority flood with an epsilon gradient, then steepest
//   descent between jittered nodes); discharge, drainage area weighted by L1's runoff over its land mean; incision by
//   stream power, implicit, n = 1: dh/dt = U - K k_rock A_w^m S (Braun and Willett 2013), K and m as L1's pass. A point
//   lower than its receiver's new height (in a filled depression) is left alone.
// - No hillslope diffusion: at 0.001-0.01 m^2/yr it changes a 500 m cell by D t / dx^2 = 0.4-4% of its relief in
//   10 Myr; hillslopes are L3's.
// - Depressions, as at L1: one its catchment's sediment can fill within the run (the run's incision plus background
//   denudation, bgKmMyr) becomes a flat the water crosses; the others stay lakes.
// - Then L1's scale is restored: the run's net change, averaged over each L1 cell, is interpolated (as the relief's
//   base) and removed, so L2 changes only the band it owns. Its size before removal is reported.
// Arrays hold only the points that matter (selected land and the outlets beside it), indexed through the box.
function* runPass(R,B,bd,sel,ref,kin,o){
  o=Object.assign({myr:2,dt:0.25,K:0.1,m:0.5,jitter:0.4,eps:1e-6,salt:0,channelKm2:1,bgKmMyr:0.02},o);
  const t0=Date.now(),{nx,ny,h,land,lake}=B,f=bd.f,dx=bd.dx,lab1=bd.lab1,inS=new Set(sel.basins),cA=dx*dx;
  // the box: L1 cells holding selected points, and a cell around
  let i0=1e9,i1=-1,j0=1e9,j1=-1;
  for(let c=0;c<nx*ny;c++){const L=ref.cells.get(c);let has=!L&&inS.has(lab1[c]);if(L)for(let q=0;q<L.length;q++)if(inS.has(L[q])){has=true;break;}
    if(!has)continue;const i=c%nx,j=(c/nx)|0;if(i<i0)i0=i;if(i>i1)i1=i;if(j<j0)j0=j;if(j>j1)j1=j;}
  i0=Math.max(0,i0-1);j0=Math.max(0,j0-1);i1=Math.min(nx-1,i1+1);j1=Math.min(ny-1,j1+1);
  const W=(i1-i0+1)*f,H=(j1-j0+1)*f,gx0=i0*f,gy0=j0*f,X0=B.x0+i0*h,Y0=B.y0+j0*h;
  const gxOf=q=>q%W,gyOf=q=>(q/W)|0,px=q=>X0+(q%W+.5)*dx,py=q=>Y0+(((q/W)|0)+.5)*dx,cellOf=q=>(j0+((((q/W)|0)/f)|0))*nx+i0+(((q%W)/f)|0);
  const sink=c=>land[c]&&lake&&lake[c]>=2,T=bd.o.T,TW=T*f;
  const zAt=q=>{const GX=gx0+q%W,GY=gy0+((q/W)|0),g=bd.graphs.get(Math.floor(GX/TW)*bd.tny+Math.floor(GY/TW));return g?g.z[(GY%TW)*TW+GX%TW]:R.z(px(q),py(q),R.o.lamMinKm);};
  // selected land points, then the outlets beside them; idx maps a box point to its index here (or -1)
  const idx=new Int32Array(W*H).fill(-1),pts=[],kd=[],zz=[];
  // (yielding every 32 rows: setting up a large selection takes seconds)
  for(let q=0;q<W*H;q++){if(q%(32*W)===0)yield{step:0,steps:0,setup:q/(W*H)};if(!ref.inside(px(q),py(q)))continue;const zq=zAt(q),c=cellOf(q);if(zq<=0||sink(c))continue;idx[q]=pts.length;pts.push(q);kd.push(1);zz.push(zq);}
  const nAct0=pts.length;
  for(let k=0;k<nAct0;k++){const q=pts[k],i=q%W,j=(q/W)|0;for(let dj=-1;dj<=1;dj++)for(let di=-1;di<=1;di++){const I=i+di,J=j+dj;if(I<0||J<0||I>=W||J>=H)continue;const p=J*W+I;if(idx[p]>=0||idx[p]===-2)continue;
    const zp=zAt(p),c=cellOf(p);if(zp<=0||sink(c)){idx[p]=pts.length;pts.push(p);kd.push(2);zz.push(zp);}else idx[p]=-2;}}
  for(let q=0;q<W*H;q++)if(idx[q]===-2)idx[q]=-1;
  const n=pts.length,kind=Uint8Array.from(kd),z=Float32Array.from(zz),z0=Float32Array.from(zz),base=k=>z[k]<=0?0:z[k];
  // neighbours (8, as indices here, -1 none) and the jittered distances to them
  const D8=[[1,0],[-1,0],[0,1],[0,-1],[1,1],[-1,1],[1,-1],[-1,-1]],nb=new Int32Array(8*n),nd=new Float32Array(8*n),jx=new Float32Array(n),jy=new Float32Array(n);
  for(let k=0;k<n;k++){const [a,b]=jit(gx0+pts[k]%W,gy0+((pts[k]/W)|0),o.salt);jx[k]=2*o.jitter*a;jy[k]=2*o.jitter*b;}
  for(let k=0;k<n;k++){const q=pts[k],i=q%W,j=(q/W)|0;for(let d=0;d<8;d++){const [di,dj]=D8[d],I=i+di,J=j+dj,m=I<0||J<0||I>=W||J>=H?-1:idx[J*W+I];nb[8*k+d]=m;if(m>=0)nd[8*k+d]=dx*Math.hypot(di+jx[m]-jx[k],dj+jy[m]-jy[k]);}}
  yield{step:0,steps:0,setup:true};
  // inputs per selected point: uplift (L1's present tectonic rate), runoff weight (L1's climate over its land mean),
  // erodibility (rock; log-interpolated from the sheet's cells)
  const bil=(F,x,y)=>{const gx=(x-B.x0)/h-0.5,gy=(y-B.y0)/h-0.5,i=Math.floor(gx),j=Math.floor(gy),tx=gx-i,ty=gy-j,at=(I,J)=>F[Math.max(0,Math.min(ny-1,J))*nx+Math.max(0,Math.min(nx-1,I))];
    return(1-tx)*(1-ty)*at(i,j)+tx*(1-ty)*at(i+1,j)+(1-tx)*ty*at(i,j+1)+tx*ty*at(i+1,j+1);};
  const cl=B.climate;let mr=1;if(cl){let s2=0,k2=0;for(let c=0;c<nx*ny;c++)if(land[c]){s2+=cl.R[c];k2++;}mr=s2/(k2||1)||1;}
  const kAt=(x,y)=>{const gx=(x-kin.x0)/kin.h-0.5,gy=(y-kin.y0)/kin.h-0.5,i=Math.floor(gx),j=Math.floor(gy),tx=gx-i,ty=gy-j,at=(I,J)=>kin.lnK[Math.max(0,Math.min(kin.ny-1,J))*kin.nx+Math.max(0,Math.min(kin.nx-1,I))];
    return Math.exp((1-tx)*(1-ty)*at(i,j)+tx*(1-ty)*at(i+1,j)+(1-tx)*ty*at(i,j+1)+tx*ty*at(i+1,j+1));};
  const U=new Float32Array(n),wgt=new Float32Array(n).fill(1),Kp=new Float32Array(n),inc=new Float32Array(n);
  for(let k=0;k<nAct0;k++){const x=px(pts[k]),y=py(pts[k]);U[k]=B.U?bil(B.U,x,y):0;if(cl)wgt[k]=Math.max(0,bil(cl.R,x,y))/mr;Kp[k]=o.K*kAt(x,y);}
  const fz=new Float64Array(n),rec=new Int32Array(n),rlen=new Float32Array(n),stack=new Int32Array(n),A=new Float64Array(n),done=new Uint8Array(n),pushed=new Uint8Array(n),hp=new L1.Heap(8*n+8);
  let ns=0,nAct=nAct0,firstRec=null;
  // routing: priority flood from the outlets (epsilon gradient: every point drains), receivers by steepest descent on
  // the filled surface between jittered nodes, the stack (outlets first) and discharge-weighted area
  function route(){done.fill(0);pushed.fill(0);hp.n=0;
    for(let k=0;k<n;k++){rec[k]=-1;if(kind[k]===2){fz[k]=base(k);pushed[k]=1;hp.push(k,fz[k]);}}
    let left=nAct;
    for(;;){while(hp.n){const k=hp.pop();if(done[k])continue;done[k]=1;if(kind[k]===1)left--;
        for(let d=0;d<8;d++){const m=nb[8*k+d];if(m<0||kind[m]!==1||pushed[m])continue;fz[m]=Math.max(z[m],fz[k]+o.eps);pushed[m]=1;hp.push(m,fz[m]);}}
      // points no outlet reaches (cut off by the fixed divide): their lowest point becomes a sink
      if(left<=0)break;let lo=-1;for(let k=0;k<n;k++)if(kind[k]===1&&!done[k]&&(lo<0||z[k]<z[lo]))lo=k;if(lo<0)break;kind[lo]=2;nAct--;left--;fz[lo]=z[lo];pushed[lo]=1;hp.push(lo,fz[lo]);}
    for(let k=0;k<n;k++){if(kind[k]!==1)continue;let best=-1,bs=0,bl=0;
      for(let d=0;d<8;d++){const m=nb[8*k+d];if(m<0||!kind[m])continue;const fm=kind[m]===2?base(m):fz[m];if(fm>=fz[k])continue;const L=nd[8*k+d],sl=(fz[k]-fm)/L;if(sl>bs){bs=sl;best=m;bl=L;}}
      rec[k]=best;rlen[k]=bl;}
    const cnt=new Int32Array(n+1);for(let k=0;k<n;k++)if(kind[k]===1&&rec[k]>=0)cnt[rec[k]+1]++;for(let k=0;k<n;k++)cnt[k+1]+=cnt[k];
    const don=new Int32Array(cnt[n]),fp=cnt.slice(0,n);for(let k=0;k<n;k++)if(kind[k]===1&&rec[k]>=0)don[fp[rec[k]]++]=k;
    ns=0;for(let k=0;k<n;k++)if(kind[k]===2)stack[ns++]=k;for(let s2=0;s2<ns;s2++){const k=stack[s2];for(let p=cnt[k];p<cnt[k+1];p++)stack[ns++]=don[p];}
    A.fill(0);for(let s2=ns-1;s2>=0;s2--){const k=stack[s2];if(kind[k]!==1)continue;A[k]+=cA*wgt[k];const r=rec[k];if(r>=0)A[r]+=A[k];}}
  const steps=Math.max(1,Math.round(o.myr/o.dt)),dt=o.myr/steps;
  for(let s2=0;s2<steps;s2++){
    for(let k=0;k<n;k++)if(kind[k]===1)z[k]+=U[k]*dt;
    route();if(!firstRec)firstRec=Int32Array.from(rec);
    // implicit stream power, receivers first
    for(let s3=0;s3<ns;s3++){const k=stack[s3];if(kind[k]!==1)continue;const r=rec[k];if(r<0)continue;const zr=kind[r]===2?base(r):z[r];if(z[k]<=zr)continue;
      const F=Kp[k]*Math.pow(A[k],o.m)*dt/rlen[k],zn=(z[k]+F*zr)/(1+F);inc[k]+=z[k]-zn;z[k]=zn;}
    yield{step:s2+1,steps};}
  route();let changed=0;for(let k=0;k<n;k++)if(kind[k]===1&&rec[k]!==firstRec[k])changed++;
  // depressions: connected flooded points; filled if their catchment's sediment (incision and background denudation
  // over the run) reaches their volume
  const lid=new Int32Array(n).fill(-1),vol=[],sed=[],lakes=[];
  for(let k=0;k<n;k++){if(kind[k]!==1||lid[k]>=0||!(fz[k]>z[k]+1e-6))continue;const id=vol.length,q=[k];lid[k]=id;let v=0;
    for(let a=0;a<q.length;a++){const u=q[a];v+=(fz[u]-z[u])*cA;for(let d=0;d<8;d++){const m=nb[8*u+d];if(m>=0&&kind[m]===1&&lid[m]<0&&fz[m]>z[m]+1e-6){lid[m]=id;q.push(m);}}}vol.push(v);sed.push(0);lakes.push(q);}
  const qs=new Float64Array(n);for(let s3=ns-1;s3>=0;s3--){const k=stack[s3];if(kind[k]!==1)continue;qs[k]+=(inc[k]+o.bgKmMyr*o.myr)*cA;
    if(lid[k]>=0){sed[lid[k]]+=qs[k];continue;}const r=rec[k];if(r>=0&&kind[r]===1)qs[r]+=qs[k];}
  let filledN=0,filledKm2=0,lakeKm2=0,lakeN=0;const lakeFlag=new Uint8Array(n);
  for(let id=0;id<vol.length;id++){if(sed[id]>=vol[id]){filledN++;for(const u of lakes[id]){z[u]=fz[u];filledKm2+=cA;}}else{lakeN++;for(const u of lakes[id]){lakeFlag[u]=1;lakeKm2+=cA;}}}
  // L1's scale: the net change averaged per L1 cell (over the cell's selected points), interpolated, removed
  const cw=i1-i0+1,ch=j1-j0+1,cm=new Float64Array(cw*ch),cn=new Float64Array(cw*ch);
  for(let k=0;k<n;k++){if(kind[k]!==1)continue;const q=pts[k],ci=((q%W)/f)|0,cj=(((q/W)|0)/f)|0;cm[cj*cw+ci]+=z[k]-z0[k];cn[cj*cw+ci]++;}
  let drift2=0,dn=0;for(let c=0;c<cw*ch;c++)if(cn[c]){cm[c]/=cn[c];drift2+=cm[c]*cm[c]*cn[c];dn+=cn[c];}
  const crS=(p0,p1,p2,p3,t)=>p1+0.5*t*(p2-p0+t*(2*p0-5*p1+4*p2-p3+t*(3*(p1-p2)+p3-p0))),cat=(I,J)=>cm[Math.max(0,Math.min(ch-1,J))*cw+Math.max(0,Math.min(cw-1,I))];
  for(let k=0;k<n;k++){if(kind[k]!==1)continue;const q=pts[k],gx=(q%W+.5)/f-0.5,gy=(((q/W)|0)+.5)/f-0.5,I=Math.floor(gx),J=Math.floor(gy),tx=gx-I,ty=gy-J,r4=[0,0,0,0];
    for(let m2=0;m2<4;m2++)r4[m2]=crS(cat(I-1,J-1+m2),cat(I,J-1+m2),cat(I+1,J-1+m2),cat(I+2,J-1+m2),tx);z[k]-=crS(r4[0],r4[1],r4[2],r4[3],ty);}
  route();
  // after: change, Strahler order on the channels (weighted area at least channelKm2: a channel's order rises by one
  // only where two of the highest order among its tributaries meet), main-stem length (longest flow path upstream)
  const dz=new Float32Array(n);let er=0,maxInc=0;for(let k=0;k<n;k++)if(kind[k]===1){dz[k]=z[k]-z0[k];er+=dz[k];if(-dz[k]>maxInc)maxInc=-dz[k];}
  const order=new Uint8Array(n),len=new Float32Array(n),mo=new Uint8Array(n),mc=new Uint8Array(n);
  for(let s3=ns-1;s3>=0;s3--){const k=stack[s3];if(kind[k]!==1)continue;const r=rec[k];
    if(A[k]>=o.channelKm2){order[k]=!mo[k]?1:mc[k]>=2?mo[k]+1:mo[k];if(r>=0&&kind[r]===1){if(order[k]>mo[r]){mo[r]=order[k];mc[r]=1;}else if(order[k]===mo[r])mc[r]++;}}
    if(r>=0&&kind[r]===1){const L=len[k]+rlen[k];if(L>len[r])len[r]=L;}}
  // Hack's law: main-stem length against (weighted) drainage area on the channels, least squares in logs
  let sx=0,sy=0,sxx=0,sxy=0,nh=0;for(let k=0;k<n;k++){if(kind[k]!==1||A[k]<o.channelKm2||!(len[k]>0))continue;const a=Math.log(A[k]),l=Math.log(len[k]);sx+=a;sy+=l;sxx+=a*a;sxy+=a*l;nh++;}
  const hack=nh>2?(sxy-sx*sy/nh)/(sxx-sx*sx/nh):NaN;let maxOrder=0;for(let k=0;k<n;k++)if(order[k]>maxOrder)maxOrder=order[k];
  return{W,H,dx,X0,Y0,pts:Int32Array.from(pts),idx,kind,z,z0,dz,rec,A,order,len,lake:lakeFlag,steps,dt,myr:o.myr,nAct,o,
    stats:{hack,maxOrder,meanChangeM:1000*er/nAct,maxIncisionM:1000*maxInc,l1DriftM:1000*Math.sqrt(drift2/(dn||1)),receiversChangedFrac:changed/nAct,depressions:vol.length,filled:filledN,filledKm2,lakes:lakeN,lakeKm2,activeKm2:nAct*cA,ms:Date.now()-t0},
    at:(x,y)=>{const i=Math.floor((x-X0)/dx),j=Math.floor((y-Y0)/dx);return i<0||j<0||i>=W||j>=H?-1:idx[j*W+i];},
    // the surface for drawing: bilinear between the selected points around (x, y) (drawn point by point, the
    // hillshade lit only the points' edges)
    zAt:(x,y)=>{const gx=(x-X0)/dx-0.5,gy=(y-Y0)/dx-0.5,i=Math.floor(gx),j=Math.floor(gy),tx=gx-i,ty=gy-j;let s=0,w=0;
      for(let c=0;c<4;c++){const I=i+(c&1),J=j+(c>>1);if(I<0||J<0||I>=W||J>=H)continue;const k=idx[J*W+I];if(k<0||kind[k]!==1)continue;const ww=(c&1?tx:1-tx)*(c>>1?ty:1-ty);s+=ww*z[k];w+=ww;}
      return w>0?s/w:NaN;}};}
function pass(R,B,bd,sel,ref,kin,o){const it=runPass(R,B,bd,sel,ref,kin,o);let r;while(!(r=it.next()).done);return r.value;}

return {relief,reliefInputs,reliefFrom,boundary,tileMaker,flood,floodSeeded,tileGraph,l1Labels,makeGrad2,gblur,erodibility,runPass,pass};
}
const L2=L2module(SIM,L1);
// ===== L2-END =====
