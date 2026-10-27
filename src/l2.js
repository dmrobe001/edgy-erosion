// ===== L2-START =====
// L2: drainage basins selected on the refined L1 map, at 500 m. So far: the surface L2 starts from, as a function of
// position. Any point of it can be evaluated without the rest, which is what burning in basin boundaries in strips
// needs (README, section 6, L2).
const L2 = (() => {
const {R_KM,hashStr}=SIM;

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
function relief(sheet,B,o){
  o=Object.assign({H:0.5,lamMaxKm:16,lamMinKm:1,refKm:40,windowKm:30,aniso:5,anisoMax:4,latticeKm:20,burnKm2:1000,burnRise:3,burnWidthKm:1,burnMaxKm:15,floorKm:0.6,warpKm:1.5,warpLamKm:16},o);
  const {nx,ny,h,e,land}=B,nc=nx*ny,seed=hashStr(String(sheet.world.opts.seed)+"|L2relief")|0,noise=makeGrad2(seed);
  // RMS of one noise octave (sampled once on a fine, irregular set of points), to scale octaves to a given RMS
  let ss=0;const NS=40000;for(let k=0;k<NS;k++){const v=noise(k*0.6180339887%97.3,k*0.7548776662%89.1,1000);ss+=v*v;}const nrm=1/Math.sqrt(ss/NS);
  // 1. L1 roughness r1, land only
  const lm=Uint8Array.from(land,(v,c)=>v&&e[c]>0?1:0),hp=new Float64Array(nc),bl=gblur(e,lm,nx,ny,1);
  for(let c=0;c<nc;c++)hp[c]=lm[c]?(e[c]-bl[c])**2:0;
  const r1=gblur(hp,lm,nx,ny,o.windowKm/h);for(let c=0;c<nc;c++)r1[c]=Math.sqrt(r1[c]);
  // 2. shortening, as an axial vector field (s cos 2t, s sin 2t) on the sheet's grid, t the shortening axis
  const sf=sheet.strainField(),snc=sheet.nx*sheet.ny,qc=new Float64Array(snc),qs=new Float64Array(snc);
  for(let c=0;c<snc;c++){if(isNaN(sf.lmin[c]))continue;const s=Math.max(0,1-sf.lmin[c]);qc[c]=s*Math.cos(2*sf.dir[c]);qs[c]=s*Math.sin(2*sf.dir[c]);}
  const Lb=o.latticeKm,lx0=B.x0,ly0=B.y0,lnx=Math.ceil(nx*h/Lb)+2,lny=Math.ceil(ny*h/Lb)+2,cs=new Float64Array(lnx*lny),sn=new Float64Array(lnx*lny),an=new Float64Array(lnx*lny).fill(1);
  for(let J=0;J<lny;J++)for(let I=0;I<lnx;I++){const x=lx0+I*Lb,y=ly0+J*Lb,a=sheet.sample(qc,x,y),b=sheet.sample(qs,x,y),s=Math.hypot(a,b),k=J*lnx+I;
    const st=Math.min(o.anisoMax,1+o.aniso*s);if(st<1.02){cs[k]=1;sn[k]=0;continue;}
    const strike=0.5*Math.atan2(b,a)+Math.PI/2;cs[k]=Math.cos(strike);sn[k]=Math.sin(strike);an[k]=st;}
  // 3. trunk rivers: segments from each cell's centre to its receiver's, bucketed by the cells within burnMaxKm
  const f=B.S.f,rec=B.rec,A=B.A,gside=new Float64Array(nc),segs=[];
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
// The selection's outline at L2's resolution, by strips: labels are computed on a fixed lattice of tiles (T x T L1
// cells), each routed on itself plus a halo of k L1 cells, the halo's outer edge an outlet with L1's own labels. A
// tile's labels depend only on its position, never on the selection, so two selections sharing a boundary compute
// the same outline there. The halo is wide enough when no core cell takes its label from a guess: a core cell's
// water reaches a real outlet (the sea, a sink), or the halo's edge where L1's label is sure (the L1 cells around all
// share it); if it reaches the edge near an L1 divide, the halo doubles, up to kMax. (Two halos agreeing
// with each other proved nothing: both sent the same water to the same uncertain edge.) The first halo is how far
// L2's relief can move a divide: the added relief's RMS over L1's regional
// slope at the cell (a contour on a slope s moves by about sigma / s), in L1 cells, from kMin to kMax, the most
// over the tile. The regional slope is that of L1's surface blurred over slopeKm: the steepest drop to a neighbour
// is usually the wall of a small valley L1 carved, and in flat country it hid how far water wanders at 500 m.
// Tiles are routed for the L1 cells within max(band, k / 2) cells of the selection's outer edge; elsewhere a point
// keeps its L1 cell's label.
function boundary(R,B,o){
  o=Object.assign({dx:0.5,T:4,kMin:2,kMax:24,band:2,slopeKm:20},o);
  const {nx,ny,h,e,land,lake}=B,f=Math.round(h/o.dx),dx=h/f,lab1=l1Labels(B),blocks=new Map(),tiles=new Map();
  // the relief's RMS per unit of L1 roughness (all octaves, as drawn; an upper bound once the band is high-passed)
  const sig=Math.sqrt(R.octaves.reduce((a,L)=>a+Math.pow(L/R.o.refKm,2*R.o.H),0));
  const lm=Uint8Array.from(land),eS=gblur(e,lm,nx,ny,o.slopeKm/h),kc=new Uint8Array(nx*ny);
  for(let c=0;c<nx*ny;c++){const i=c%nx,j=(c/nx)|0,at=(I,J)=>eS[Math.max(0,Math.min(ny-1,J))*nx+Math.max(0,Math.min(nx-1,I))];
    const s=Math.hypot(at(i+1,j)-at(i-1,j),at(i,j+1)-at(i,j-1))/(2*h),w=sig*R.r1[c]/Math.max(s,1e-6);kc[c]=Math.max(o.kMin,Math.min(o.kMax,Math.ceil(w/h)));}
  const kOf=c=>kc[c];
  // the surface at L2's resolution in one L1 cell, f x f values, when first needed
  const block=c=>{let b=blocks.get(c);if(b)return b;b=new Float32Array(f*f);const i=c%nx,j=(c/nx)|0;
    for(let v=0;v<f;v++)for(let u=0;u<f;u++)b[v*f+u]=R.z(B.x0+i*h+(u+.5)*dx,B.y0+j*h+(v+.5)*dx,R.o.lamMinKm);blocks.set(c,b);return b;};
  const sink=c=>land[c]&&lake&&lake[c]>=2;
  // L1 cells whose label is sure: every L1 cell within a cell shares it
  const sure=new Uint8Array(nx*ny);for(let c=0;c<nx*ny;c++){const i=c%nx,j=(c/nx)|0;let ok=1;for(let dj=-1;dj<=1&&ok;dj++)for(let di=-1;di<=1;di++){const I=i+di,J=j+dj;if(I>=0&&J>=0&&I<nx&&J<ny&&lab1[J*nx+I]!==lab1[c]){ok=0;break;}}sure[c]=ok;}
  // a tile's labels: per core L1 cell, f x f labels
  function tile(TI,TJ){const key=TI+","+TJ;let t=tiles.get(key);if(t)return t;const i0=TI*o.T,j0=TJ*o.T;let k=o.kMin;
    for(let j=j0;j<j0+o.T;j++)for(let i=i0;i<i0+o.T;i++)if(i<nx&&j<ny)k=Math.max(k,kOf(j*nx+i));
    // double the halo while a core cell's label is a guess (or to kMax)
    let r=route(i0,j0,k),tries=1;while(r.guesses&&k<o.kMax){k=Math.min(o.kMax,2*k);r=route(i0,j0,k);tries++;}
    t={k,tries,guesses:r.guesses,cells:r.out};tiles.set(key,t);return t;}
  // the core cells' labels of the tile at (i0, j0), routed with a halo of k cells
  function route(i0,j0,k){
    const a0=i0-k,b0=j0-k,n1=o.T+2*k,W=n1*f,z=new Float32Array(W*W).fill(1e9),l0=new Int32Array(W*W).fill(-1);
    for(let bj=0;bj<n1;bj++)for(let bi=0;bi<n1;bi++){const I=a0+bi,J=b0+bj;const inG=I>=0&&J>=0&&I<nx&&J<ny,c=inG?J*nx+I:-1,bl=inG?block(c):null,edge=bi===0||bj===0||bi===n1-1||bj===n1-1;
      for(let v=0;v<f;v++)for(let u=0;u<f;u++){const q=(bj*f+v)*W+bi*f+u;if(!inG){z[q]=-1;l0[q]=-2;continue;}z[q]=bl[v*f+u];
        const ring=edge&&(bi===0&&u===0||bj===0&&v===0||bi===n1-1&&u===f-1||bj===n1-1&&v===f-1);if(ring||z[q]<=0||sink(c))l0[q]=lab1[c];}}
    const {lab,src}=flood(z,W,W,l0),out=new Map();let guesses=0;
    for(let bj=k;bj<k+o.T;bj++)for(let bi=k;bi<k+o.T;bi++){const I=a0+bi,J=b0+bj;if(I<0||J<0||I>=nx||J>=ny)continue;const L=new Int32Array(f*f);
      for(let v=0;v<f;v++)for(let u=0;u<f;u++){const q=(bj*f+v)*W+bi*f+u;L[v*f+u]=lab[q];
        // a guess: the outlet is on the halo's edge (not sea or sink) in an L1 cell near a divide
        const sq=src[q];if(sq>=0&&z[sq]>0){const sb=(sq/W/f)|0,sa=((sq%W)/f)|0,cc=(b0+sb)*nx+a0+sa;if(cc>=0&&cc<nx*ny&&!sink(cc)&&!sure[cc])guesses++;}}
      out.set(J*nx+I,L);}
    return{out,guesses};}
  // The selection refined. Edge cells: selected with an unselected neighbour, or the reverse; then every L1 cell
  // within max(band, k / 2) of one, k its tile's halo, gets L2 labels from its tile.
  const tKey=c=>Math.floor((c%nx)/o.T)+","+Math.floor(((c/nx)|0)/o.T);
  const tileAt=c=>tile(Math.floor((c%nx)/o.T),Math.floor(((c/nx)|0)/o.T));
  function edgeCells(selC){const edge=[];for(let c=0;c<nx*ny;c++){const i=c%nx,j=(c/nx)|0;let on=false;
      for(let dj=-1;dj<=1&&!on;dj++)for(let di=-1;di<=1;di++){const I=i+di,J=j+dj;if(I<0||J<0||I>=nx||J>=ny)continue;if(selC(J*nx+I)!==selC(c)&&(land[c]||land[J*nx+I])){on=true;break;}}if(on)edge.push(c);}
    return edge;}
  function nearCells(selC,edge){const near=new Set();
    for(const c of edge){const d=Math.max(o.band,Math.ceil(Math.max(kc[c],tileAt(c).k)/2)),i=c%nx,j=(c/nx)|0;for(let dj=-d;dj<=d;dj++)for(let di=-d;di<=d;di++){const I=i+di,J=j+dj;if(I>=0&&J>=0&&I<nx&&J<ny)near.add(J*nx+I);}}
    return near;}
  // tiles a refinement still needs routed, [TI, TJ] each: the edge cells' first (the band's width depends on their
  // halos), then the band's; empty when refine(sel) can run at once. Lets the page route a tile per frame.
  // the selection as a flag per L1 cell (a set lookup per neighbour made the scan of a large grid slow)
  const selFlags=sel=>{const inSel=new Set(sel.basins),fl=new Uint8Array(nx*ny);for(let c=0;c<nx*ny;c++)if(lab1[c]>=0&&inSel.has(lab1[c]))fl[c]=1;return c=>fl[c]===1;};
  function pending(sel){const selC=selFlags(sel),edge=edgeCells(selC),need=new Map();
    for(const c of edge){const k=tKey(c);if(!tiles.has(k))need.set(k,k.split(",").map(Number));}if(need.size)return[...need.values()];
    for(const c of nearCells(selC,edge)){const k=tKey(c);if(!tiles.has(k))need.set(k,k.split(",").map(Number));}return[...need.values()];}
  function refine(sel){const t0=Date.now(),inSel=new Set(sel.basins),selC=selFlags(sel),near=nearCells(selC,edgeCells(selC));
    const cells=new Map(),seen=new Set();for(const c of near){seen.add(tKey(c));cells.set(c,tileAt(c).cells.get(c));}
    const inside=(x,y)=>{const i=Math.floor((x-B.x0)/h),j=Math.floor((y-B.y0)/h);if(i<0||j<0||i>=nx||j>=ny)return false;const c=j*nx+i,L=cells.get(c);
      if(!L)return selC(c);const u=Math.min(f-1,Math.floor((x-B.x0-i*h)/dx)),v=Math.min(f-1,Math.floor((y-B.y0-j*h)/dx));return inSel.has(L[v*f+u]);};
    // land area: L1 cells refined count their land points in the selection; the rest count whole if selected
    let area=0,areaL1=0;for(let c=0;c<nx*ny;c++){if(!land[c])continue;if(selC(c))areaL1+=h*h;const L=cells.get(c);if(!L){if(selC(c))area+=h*h;continue;}const bl=block(c);let n=0;for(let q=0;q<L.length;q++)if(inSel.has(L[q])&&bl[q]>0)n++;area+=n*dx*dx;}
    let kMaxUsed=0,unsettled=0;for(const key of seen){const t=tiles.get(key);kMaxUsed=Math.max(kMaxUsed,t.k);if(t.guesses)unsettled++;}
    return{inside,cells,nTiles:seen.size,nCells:near.size,areaKm2:area,areaL1Km2:areaL1,kMaxUsed,unsettled,ms:Date.now()-t0};}
  // the same labelling routed over a whole box at once (the reference the strips are tested against): nX x nY cells
  // of dx from (x0, y0), the box's edge an outlet with L1's labels
  function full(x0,y0,nX,nY){const z=R.grid(x0,y0,dx,nX,nY,R.o.lamMinKm),l0=new Int32Array(nX*nY).fill(-1);
    for(let j=0;j<nY;j++)for(let i=0;i<nX;i++){const q=j*nX+i,x=x0+(i+.5)*dx,y=y0+(j+.5)*dx,ci=Math.floor((x-B.x0)/h),cj=Math.floor((y-B.y0)/h);
      if(ci<0||cj<0||ci>=nx||cj>=ny){l0[q]=-2;continue;}const c=cj*nx+ci;if(i===0||j===0||i===nX-1||j===nY-1||z[q]<=0||sink(c))l0[q]=lab1[c];}
    return{z,lab:flood(z,nX,nY,l0).lab};}
  return{o,dx,f,lab1,tile,pending,refine,full,kOf,blocks,tiles};}
return {relief,boundary,flood,l1Labels,makeGrad2,gblur};
})();
// ===== L2-END =====
