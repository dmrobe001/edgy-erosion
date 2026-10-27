// ---------- rendering ----------
const {ROCK,R_KM}=SIM;
function hexToRgb(h){const n=parseInt(h.slice(1),16);return[(n>>16)&255,(n>>8)&255,n&255];}
function buildLut(stops,n,lo,hi){const lut=new Uint8Array(n*3),s=stops.map(([v,c])=>[v,hexToRgb(c)]);
  for(let k=0;k<n;k++){const v=lo+(hi-lo)*k/(n-1);let a=0;while(a<s.length-2&&v>s[a+1][0])a++;
    const [v0,c0]=s[a],[v1,c1]=s[a+1];let t=(v-v0)/(v1-v0);t=Math.max(0,Math.min(1,t));for(let c=0;c<3;c++)lut[k*3+c]=Math.round(c0[c]+(c1[c]-c0[c])*t);}
  return lut;}
const ELEV_STOPS=[[-8,"#06132a"],[-6,"#0b2347"],[-4,"#14396b"],[-2,"#235a94"],[-0.6,"#3f86bd"],[-0.01,"#7fbad6"],[0,"#5f8f4e"],[0.4,"#7ea258"],[1.2,"#b6a76a"],[2.5,"#9a7b52"],[4,"#8a7f78"],[5.5,"#e8e6e2"],[8,"#ffffff"]];
const AGE_STOPS=[[0,"#c8102e"],[20,"#f26b21"],[45,"#f6d23c"],[75,"#5fb85a"],[110,"#2a9dbf"],[150,"#2c4fa3"],[200,"#4b2a83"]];
const THICK_STOPS=[[0,"#1d1140"],[8,"#30307a"],[25,"#20718b"],[35,"#3aa67a"],[50,"#b7d24b"],[65,"#f4b740"],[80,"#f0f0d0"]];
const ROCKS=[null,
  {name:"Ocean-ridge basalt",c:"#5a6b86"},{name:"Hotspot and flood basalt",c:"#c8433a"},{name:"Arc volcanics",c:"#ec9468"},
  {name:"Cratonic basement",c:"#d781a8"},{name:"Collisional metamorphics",c:"#8c6cc0"},{name:"Marine shelf sediments",c:"#ead26e"},{name:"Clastic sediments",c:"#c6a06a"}];
const ROCK_RGB=ROCKS.map(r=>r?hexToRgb(r.c):[0,0,0]);
const LUT_E=buildLut(ELEV_STOPS,1024,-8,8),LUT_A=buildLut(AGE_STOPS,256,0,200),LUT_T=buildLut(THICK_STOPS,256,0,80);
function gradCss(stops,lo,hi){return "linear-gradient(90deg,"+stops.map(([v,c])=>`${c} ${((v-lo)/(hi-lo)*100).toFixed(1)}%`).join(",")+")";}
function hsl(h,s,l){s/=100;l/=100;const k=n=>(n+h/30)%12,a=s*Math.min(l,1-l),f=n=>l-a*Math.max(-1,Math.min(k(n)-3,Math.min(9-k(n),1)));return[Math.round(f(0)*255),Math.round(f(8)*255),Math.round(f(4)*255)];}
function getVar(n){return getComputedStyle(document.documentElement).getPropertyValue(n).trim()||"#888888";}
const $=id=>document.getElementById(id);

const canvas=$("map"),wrap=$("mapwrap");
let world,view="elev",proj="globe",running=false,msEma=0,rateEma=0,lastFrame=performance.now(),loggedCount=0;
let cellRGB=null,colorsDirty=true,lightW=null,handoff=null,downloadsP=null,selCell=-1,l1=null,l1Time=0,geo=null,geoDirty=true;
const fineNoise=SIM.makeNoise3(918273);
let yaw=0.6,pitch=0.35,camDirty=true,pixIdx=null,pixLimb=null,pixDir=null,imgData=null,pixW=0,pixH=0,Mw=null;
// zoom on the globe (1 = whole disc) and on the map (1 = whole map, centred on mC in map units); run target (Myr)
let gZoom=1,mZoom=1,mC=[0,0],goode=null,goodeSteps=-1,runStart=0,runTarget=0,selMask=null,selMaskKey="";

function prepLight(cube){ // per-neighbour weights for hillshading under a fixed sun
  if(cube._lw)return cube._lw;const N=cube.N,p=cube.pos,nb=cube.nb,L=[-0.45,0.55,0.70];const cr=Math.PI/2/cube.n,lw=new Float32Array(4*N);
  for(let i=0;i<N;i++)for(let k=0;k<4;k++){const j=nb[4*i+k];const dx=p[3*j]-p[3*i],dy=p[3*j+1]-p[3*i+1],dz=p[3*j+2]-p[3*i+2],d2=dx*dx+dy*dy+dz*dz||1;
    lw[4*i+k]=(dx*L[0]+dy*L[1]+dz*L[2])/d2*cr;}
  return cube._lw=lw;}

function computeColors(){
  const w=world,N=w.N,E=w.elev,K=w.wK,O=w.owner,T=w.wT,B=w.wB,R=w.wR,RA=w.wRA,tr=w.trench,nb=w.cube.nb,lw=lightW,pos=w.cube.pos,pl=w.plates,t=w.t;
  if(!cellRGB||cellRGB.length!==N*3)cellRGB=new Uint8Array(N*3);
  const pc=view==="plates"?pl.map(p=>[hsl(p.hue,48,64),hsl(p.hue,42,36)]):null;
  const conv=hexToRgb(getVar("--conv")),divc=hexToRgb(getVar("--div")),tran=hexToRgb(getVar("--trans"));
    for(let i=0;i<N;i++){
    let e=E[i];const land=e>0;let g=0;for(let k=0;k<4;k++)g+=(E[nb[4*i+k]]-e)*lw[4*i+k];
    let s=land?1-0.35*g:1-0.1*g;s=s<0.6?0.6:s>1.25?1.25:s;
    let r,gg,b;
    if(view==="elev"){if(!K[i]&&tr[i]>0.05)e-=2.2*Math.min(1,tr[i]);let k=Math.round((e+8)/16*1023);k=k<0?0:k>1023?1023:k;r=LUT_E[k*3];gg=LUT_E[k*3+1];b=LUT_E[k*3+2];}
    else if(view==="rock"){const c=ROCK_RGB[R[i]]||ROCK_RGB[1];const age=t-RA[i];const f=1-0.32*Math.min(1,age/1500);r=c[0]*f;gg=c[1]*f;b=c[2]*f;s=1+(s-1)*0.6;}
    else{const p=O[i],c=pc[p][K[i]?0:1];r=c[0];gg=c[1];b=c[2];s=1+(s-1)*0.4;
      for(let k=0;k<4;k++){const j=nb[4*i+k],q=O[j];if(q===p)continue;
        const P=pl[p],Q=pl[q],dx=Q.W[0]-P.W[0],dy=Q.W[1]-P.W[1],dz=Q.W[2]-P.W[2],x=pos[3*i],y=pos[3*i+1],z=pos[3*i+2];
        const vx=dy*z-dz*y,vy=dz*x-dx*z,vz=dx*y-dy*x;let ex=pos[3*j]-x,ey=pos[3*j+1]-y,ez=pos[3*j+2]-z;const el=Math.hypot(ex,ey,ez)||1;
        const rel=(vx*ex+vy*ey+vz*ez)/el*R_KM/10;const bc=rel<-0.7?conv:rel>0.7?divc:tran;r=bc[0];gg=bc[1];b=bc[2];s=1;break;}}
    cellRGB[3*i]=Math.min(255,r*s);cellRGB[3*i+1]=Math.min(255,gg*s);cellRGB[3*i+2]=Math.min(255,b*s);
  }
  colorsDirty=false;
}

function camMatrix(){ // columns: world directions of camera right, up, toward-viewer
  const cp=Math.cos(pitch),sp=Math.sin(pitch),cy=Math.cos(yaw),sy=Math.sin(yaw);
  const f=v=>{const y1=v[1]*cp+v[2]*sp,z1=-v[1]*sp+v[2]*cp;const wx=z1,wy=v[0],wz=y1;return[wx*cy-wy*sy,wx*sy+wy*cy,wz];};
  const a=f([1,0,0]),b=f([0,1,0]),c=f([0,0,1]);return[a,b,c];}

// Interrupted Goode homolosine (equal area): sinusoidal within PHI1 of the equator, Mollweide beyond, shifted to meet.
// Map units: x is longitude offset at the equator (radians, -pi..pi), y in the same units. The cuts are placed for this
// planet: the map edge on the meridian crossing the least land, then one more cut in the north and two in the south,
// each on the meridian crossing the least land in its hemisphere (counting land within 4 degrees), 50 degrees apart.
const G_PHI1=0.7109889,G_DY=0.0528035,G_YMAX=Math.SQRT2-G_DY,G_ASP=2*Math.PI/(2*G_YMAX);
function goodeLat(phi){const a=Math.abs(phi);if(a<=G_PHI1)return[phi,Math.cos(phi)];
  let t=a,tg=Math.PI*Math.sin(a);for(let k=0;k<30;k++){const f=2*t+Math.sin(2*t)-tg,d=2+2*Math.cos(2*t);if(Math.abs(d)<1e-12)break;t-=f/d;if(Math.abs(f)<1e-12)break;}
  t=Math.min(t,Math.PI/2);return[Math.sign(phi)*(Math.SQRT2*Math.sin(t)-G_DY),2*Math.SQRT2/Math.PI*Math.cos(t)];}
function goodeLatInv(y){const a=Math.abs(y);if(a<=G_PHI1)return[y,Math.cos(y)];if(a>G_YMAX)return null;
  const t=Math.asin(Math.min(1,(a+G_DY)/Math.SQRT2)),phi=Math.asin(Math.max(-1,Math.min(1,(2*t+Math.sin(2*t))/Math.PI)));return[Math.sign(y)*phi,2*Math.SQRT2/Math.PI*Math.cos(t)];}
function goodeLayout(){const lk=world.cube.lookup,E=world.elev,nl=[new Float64Array(360),new Float64Array(360)];
  for(let k=0;k<360;k++){const lon=(k+0.5-180)*Math.PI/180,c=Math.cos(lon),sn=Math.sin(lon);
    for(let j=0;j<90;j++){const lat=(j+0.5)*Math.PI/180,cl=Math.cos(lat),sl=Math.sin(lat);if(E[lk(cl*c,cl*sn,sl)]>0)nl[0][k]++;if(E[lk(cl*c,cl*sn,-sl)]>0)nl[1][k]++;}}
  const cost=nl.map(a=>Float64Array.from(a,(v,k)=>{let s=0;for(let d=-4;d<=4;d++)s+=a[(k+d+360)%360];return s;})),far=(k,cuts)=>cuts.every(c=>Math.min(Math.abs(k-c),360-Math.abs(k-c))>=50);
  let edge=0;for(let k=1;k<360;k++)if(cost[0][k]+cost[1][k]<cost[0][edge]+cost[1][edge])edge=k;
  const pick=(h,n)=>{const cuts=[edge];for(let q=0;q<n;q++){let b=-1;for(let k=0;k<360;k++)if(far(k,cuts)&&(b<0||cost[h][k]<cost[h][b]))b=k;if(b>=0)cuts.push(b);}return cuts.slice(1);};
  const lc=(edge+0.5-180)*Math.PI/180+Math.PI,rel=k=>{let r=(k+0.5-180)*Math.PI/180-lc;while(r<=-Math.PI)r+=2*Math.PI;while(r>Math.PI)r-=2*Math.PI;return r;};
  const lobes=cuts=>{const xs=[-Math.PI,...cuts.map(rel).sort((a,b)=>a-b),Math.PI],L=[];for(let q=0;q+1<xs.length;q++)L.push({a:xs[q],b:xs[q+1],m:(xs[q]+xs[q+1])/2});return L;};
  goode={lc,edgeLonDeg:(edge+0.5-180),north:lobes(pick(0,1)),south:lobes(pick(1,2))};goodeSteps=world.steps;camDirty=true;}
function goodeFwd(lat,lon){let r=lon-goode.lc;while(r<=-Math.PI)r+=2*Math.PI;while(r>Math.PI)r-=2*Math.PI;const L=lat>=0?goode.north:goode.south;
  let q=0;while(q<L.length-1&&r>=L[q].b)q++;const [y,k]=goodeLat(lat);return[L[q].m+(r-L[q].m)*k,y,(lat>=0?0:10)+q];}
function goodeInv(x,y,out){const li=goodeLatInv(y);if(!li||x<-Math.PI||x>Math.PI)return false;const L=y>=0?goode.north:goode.south;
  let q=0;while(q<L.length-1&&x>=L[q].b)q++;const lb=L[q],r=lb.m+(x-lb.m)/li[1];if(r<lb.a||r>lb.b)return false;
  const lon=goode.lc+r,cl=Math.cos(li[0]);out[0]=cl*Math.cos(lon);out[1]=cl*Math.sin(lon);out[2]=Math.sin(li[0]);return true;}
// map units of a canvas pixel, and back
const mapXY=(px,py)=>[mC[0]+((px+0.5)/pixW-0.5)*2*Math.PI/mZoom,mC[1]+(0.5-(py+0.5)/pixH)*2*G_YMAX/mZoom];
function clampMap(){const hx=Math.PI*(1-1/mZoom),hy=G_YMAX*(1-1/mZoom);mC[0]=Math.max(-hx,Math.min(hx,mC[0]));mC[1]=Math.max(-hy,Math.min(hy,mC[1]));}
function sizeCanvas(){
  const dpr=window.devicePixelRatio||1,cw=wrap.clientWidth;
  if(proj==="flat"&&(!goode||(!running&&goodeSteps!==world.steps)))goodeLayout();
  let W,H;if(proj==="globe"){W=H=Math.min(520,Math.round(cw*dpr));}else{W=Math.min(1024,Math.round(cw*dpr));W-=W%2;H=Math.round(W/G_ASP);}
  if(W!==pixW||H!==pixH||!pixIdx){pixW=W;pixH=H;canvas.width=W;canvas.height=H;imgData=canvas.getContext("2d").createImageData(W,H);pixIdx=new Int32Array(W*H);pixLimb=new Uint8Array(W*H);camDirty=true;}
}
function buildPixIdx(){
  const lk=world.cube.lookup,W=pixW,H=pixH;
  if(proj==="globe"){
    Mw=camMatrix();const [a,b,c]=Mw;
    for(let py=0;py<H;py++){const sy=(1-(py+0.5)/H*2)/gZoom;for(let px=0;px<W;px++){const sx=((px+0.5)/W*2-1)/gZoom,r2=sx*sx+sy*sy,q=py*W+px;
      if(r2>1){pixIdx[q]=-1;continue;}const sz=Math.sqrt(1-r2);
      pixIdx[q]=lk(a[0]*sx+b[0]*sy+c[0]*sz,a[1]*sx+b[1]*sy+c[1]*sz,a[2]*sx+b[2]*sy+c[2]*sz);pixLimb[q]=Math.round(255*(0.5+0.5*Math.sqrt(sz)));}}
  } else {
    const v=[0,0,0];for(let py=0;py<H;py++)for(let px=0;px<W;px++){const q=py*W+px,[x,y]=mapXY(px,py);
      if(goodeInv(x,y,v)){pixIdx[q]=lk(v[0],v[1],v[2]);pixLimb[q]=255;}else pixIdx[q]=-1;}
  }
  camDirty=false;geoDirty=true;
}
function buildGeo(){ // per pixel: four cells and bilinear weights on the cube grid, plus stand-in fine relief noise
  const cube=world.cube,nb=cube.nb,W=pixW,H=pixH,P=W*H,out=[0,0,0];
  if(!geo||geo.n.length!==P)geo={c:new Int32Array(4*P),wx:new Float32Array(P),wy:new Float32Array(P),n:new Float32Array(P)};
  const put=(q,x,y,z)=>{cube.lookupFrac(x,y,z,out);const c=out[0],tx=out[1],ty=out[2];
    const cB=nb[4*c+(tx>=0?0:1)],cC=nb[4*c+(ty>=0?2:3)],cD=nb[4*cB+(ty>=0?2:3)];
    geo.c[4*q]=c;geo.c[4*q+1]=cB;geo.c[4*q+2]=cC;geo.c[4*q+3]=cD;geo.wx[q]=Math.abs(tx);geo.wy[q]=Math.abs(ty);
    const v=fineNoise(x,y,z,22,5),r=1-Math.abs(2*v-1);geo.n[q]=r*r;};
  if(proj==="globe"){const [a,b,c]=Mw;
    for(let py=0;py<H;py++){const sy=(1-(py+0.5)/H*2)/gZoom;for(let px=0;px<W;px++){const sx=((px+0.5)/W*2-1)/gZoom,r2=sx*sx+sy*sy,q=py*W+px;if(r2>1)continue;const sz=Math.sqrt(1-r2);
      put(q,a[0]*sx+b[0]*sy+c[0]*sz,a[1]*sx+b[1]*sy+c[1]*sz,a[2]*sx+b[2]*sy+c[2]*sz);}}
  } else {
    const v=[0,0,0];for(let py=0;py<H;py++)for(let px=0;px<W;px++){const [x,y]=mapXY(px,py);if(goodeInv(x,y,v))put(py*W+px,v[0],v[1],v[2]);}
  }
  geoDirty=false;
}
function ensureL1(){const now=performance.now();
  if(!l1||(l1.steps!==world.steps&&(!running||now-l1Time>2500))){l1=handoff&&handoff.l1&&handoff.l1.steps===world.steps?handoff.l1:world.buildL1();l1Time=now;}}
function paintL1(){
  const d=imgData.data,n=pixW*pixH,{Es,domain,sOff,sId,sTh,keys,rockOf}=l1,K=world.wK;
  const cs=[0,0,0,0],ws=[0,0,0,0],ptr=[0,0,0,0],lo=[0,0,0,0];
  for(let q=0;q<n;q++){const o=4*q;if(pixIdx[q]<0){d[o+3]=0;continue;}
    const wx=geo.wx[q],wy=geo.wy[q];ws[0]=(1-wx)*(1-wy);ws[1]=wx*(1-wy);ws[2]=(1-wx)*wy;ws[3]=wx*wy;
    for(let k=0;k<4;k++)cs[k]=geo.c[4*q+k];
    const e=ws[0]*Es[cs[0]]+ws[1]*Es[cs[1]]+ws[2]*Es[cs[2]]+ws[3]*Es[cs[3]];
    const nn=geo.n[q],A=0.06+0.22*Math.max(0,e),eL=e+A*(2*nn-1);
    let r,g,b,s=1;
    if(eL<=0){let k=Math.round((eL+8)/16*1023);k=k<0?0:k>1023?1023:k;r=LUT_E[k*3];g=LUT_E[k*3+1];b=LUT_E[k*3+2];s=0.92+0.12*nn;}
    else{
      // which domain governs this pixel: the continental cell with the largest weight
      let dom=-1,bw=0;for(let k=0;k<4;k++)if(K[cs[k]]&&ws[k]>bw){bw=ws[k];dom=domain[cs[k]];}
      let rock=ROCK.SHELF;
      if(dom>=0){let wt=0;for(let k=0;k<4;k++){if(K[cs[k]]&&domain[cs[k]]===dom){ptr[k]=sOff[cs[k]+1]-1;lo[k]=sOff[cs[k]];wt+=ws[k];}else{ptr[k]=-1;lo[k]=0;}}
        const cut=Math.max(0,e+0.6*A-eL);let acc=0,last=-1;
        for(let guard=0;guard<48;guard++){let bk=-Infinity,f=-1;
          for(let k=0;k<4;k++)if(ptr[k]>=lo[k]&&ptr[k]>=0){const fk=keys[sId[ptr[k]]];if(fk>bk){bk=fk;f=sId[ptr[k]];}}
          if(f<0)break;let th=0;for(let k=0;k<4;k++)if(ptr[k]>=lo[k]&&ptr[k]>=0&&sId[ptr[k]]===f){th+=ws[k]*sTh[ptr[k]];ptr[k]--;}
          acc+=th/wt;last=f;if(acc>=cut)break;}
        if(last>=0)rock=rockOf[last];}
      const c=ROCK_RGB[rock];r=c[0];g=c[1];b=c[2];s=0.7+0.5*nn;}
    const l=pixLimb[q]*s;d[o]=Math.min(255,r*l/255);d[o+1]=Math.min(255,g*l/255);d[o+2]=Math.min(255,b*l/255);d[o+3]=255;}
}
function project(v){ // world unit vector -> canvas px, or null if hidden
  if(proj==="globe"){const [a,b,c]=Mw;const x=v[0]*a[0]+v[1]*a[1]+v[2]*a[2],y=v[0]*b[0]+v[1]*b[1]+v[2]*b[2],z=v[0]*c[0]+v[1]*c[1]+v[2]*c[2];
    if(z<0.05)return null;return[(x*gZoom+1)/2*pixW,(1-y*gZoom)/2*pixH];}
  const lat=Math.asin(Math.max(-1,Math.min(1,v[2]))),lon=Math.atan2(v[1],v[0]),g=goodeFwd(lat,lon);
  return[((g[0]-mC[0])*mZoom/(2*Math.PI)+0.5)*pixW,(0.5-(g[1]-mC[1])*mZoom/(2*G_YMAX))*pixH,g[2]];
}
function render(){
  sizeCanvas();if(camDirty)buildPixIdx();
  const d=imgData.data,n=pixW*pixH;
  if(view==="l1"){ensureL1();if(geoDirty)buildGeo();paintL1();}
  else{if(colorsDirty)computeColors();
    for(let q=0;q<n;q++){const i=pixIdx[q],o=4*q;if(i<0){d[o+3]=0;continue;}const l=pixLimb[q];
      d[o]=cellRGB[3*i]*l>>8;d[o+1]=cellRGB[3*i+1]*l>>8;d[o+2]=cellRGB[3*i+2]*l>>8;d[o+3]=255;}}
  const mask=selGlobeMask();
  // the selected landmass: everything else dimmed, its submerged rim slightly less, its land outlined in white
  if(mask){const W=pixW,H=pixH,inL=q=>pixIdx[q]>=0&&mask[pixIdx[q]]===2;
    for(let q=0;q<n;q++){const i=pixIdx[q];if(i<0)continue;const m=mask[i],o=4*q;
      if(m===2){const x=q%W,y=(q/W)|0;if(x>0&&!inL(q-1)||x<W-1&&!inL(q+1)||y>0&&!inL(q-W)||y<H-1&&!inL(q+W)){d[o]=d[o+1]=d[o+2]=255;}continue;}
      const f=m===1?0.7:0.45;d[o]*=f;d[o+1]*=f;d[o+2]*=f;}}
  const ctx=canvas.getContext("2d");ctx.clearRect(0,0,pixW,pixH);ctx.putImageData(imgData,0,0);
  const dpr=pixW/wrap.clientWidth;
  const marks=$("marks").checked;
  if(marks){for(const h of world.hotspots){const p=project(h.v);if(!p)continue;const rad=Math.max(4,5*dpr);
    ctx.lineWidth=Math.max(1.5,1.6*dpr);ctx.strokeStyle="rgba(0,0,0,.55)";ctx.beginPath();ctx.arc(p[0],p[1],rad+1,0,7);ctx.stroke();
    ctx.strokeStyle="#ffcf4a";ctx.beginPath();ctx.arc(p[0],p[1],rad,0,7);ctx.stroke();}}
  if(marks){
    for(const P of world.plates){if(!P.alive||P.area<world.N*0.006)continue;const c=P.c;
      const v=[P.W[1]*c[2]-P.W[2]*c[1],P.W[2]*c[0]-P.W[0]*c[2],P.W[0]*c[1]-P.W[1]*c[0]];const k=28;
      const e=[c[0]+v[0]*k,c[1]+v[1]*k,c[2]+v[2]*k];const el=Math.hypot(...e);e[0]/=el;e[1]/=el;e[2]/=el;
      const a=project(c),b=project(e);if(!a||!b)continue;if(proj==="flat"&&a[2]!==b[2])continue;drawArrow(ctx,a[0],a[1],b[0],b[1],dpr);}
  }
  drawTraces(ctx,dpr);
}
function drawTraces(ctx,dpr){
  if(!handoff||!$("traces").checked)return;
  const cols=[null,getVar("--conv"),getVar("--div"),getVar("--trans")];
  for(const f of handoff.faults){
    const pts=f.trace.map(project);
    for(const [col,wd] of [["rgba(0,0,0,.65)",4*dpr],[cols[f.kind],2*dpr]]){
      ctx.strokeStyle=col;ctx.lineWidth=wd;ctx.lineCap="round";ctx.lineJoin="round";ctx.setLineDash(f.active||col.startsWith("rgba")?[]:[5*dpr,4*dpr]);
      ctx.beginPath();let pen=false;
      for(let k=0;k<pts.length;k++){const p=pts[k];if(!p){pen=false;continue;}
        if(pen&&proj==="flat"&&p[2]!==pts[k-1][2])pen=false;
        if(pen)ctx.lineTo(p[0],p[1]);else{ctx.moveTo(p[0],p[1]);pen=true;}}
      ctx.stroke();}
  }
  ctx.setLineDash([]);
}
function drawArrow(ctx,x0,y0,x1,y1,dpr){
  const a=Math.atan2(y1-y0,x1-x0),hl=6*dpr;
  for(const [col,w] of [["rgba(0,0,0,.6)",4*dpr],["#ffffff",1.8*dpr]]){
    ctx.strokeStyle=col;ctx.fillStyle=col;ctx.lineWidth=w;ctx.lineCap="round";
    ctx.beginPath();ctx.moveTo(x0,y0);ctx.lineTo(x1,y1);ctx.stroke();
    ctx.beginPath();ctx.moveTo(x1+Math.cos(a)*hl*0.4,y1+Math.sin(a)*hl*0.4);ctx.lineTo(x1-Math.cos(a-0.5)*hl,y1-Math.sin(a-0.5)*hl);ctx.lineTo(x1-Math.cos(a+0.5)*hl,y1-Math.sin(a+0.5)*hl);ctx.closePath();ctx.fill();}
  ctx.fillStyle="#fff";ctx.beginPath();ctx.arc(x0,y0,2.2*dpr,0,7);ctx.fill();
}

// ---------- ui ----------
const EMPTY_LOG='<li class="empty" style="border:0">No events yet.</li>';
function build(){
  world=new SIM.World({n:+$("grid").value,seed:$("seed").value||"genesis",nPlates:+$("nPlates").value,mode:$("mode").value,nnr:true,contTarget:+$("contT").value/100});
  actions=[];l0LastInspect=l1LastInspect=null;l1Reset();snaps=[];snapTick();lightW=prepLight(world.cube);colorsDirty=true;handoff=null;l1=null;geoDirty=true;selCell=-1;selMask=null;goode=null;runStart=runTarget=0;showColumn();renderSummary();camDirty=true;loggedCount=0;$("log").innerHTML=EMPTY_LOG;msEma=0;rateEma=0;render();updateUi();logAction({a:"build",settings:{seed:world.opts.seed,resolution:world.n,initialPlates:world.opts.nPlates,start:world.opts.mode,continentalTarget:world.contTarget}});
}
function fmtT(t){return t>=1000?(t/1000).toFixed(2):Math.round(t).toString();}
function updateUi(){
  const t=world.t;$("clock").innerHTML=t>=1000?fmtT(t)+"<small>billion years</small>":fmtT(t)+"<small>million years</small>";
  $("sLand").textContent=(world.landFrac*100).toFixed(1)+"%";$("sDxDt").textContent=`${Math.round(world.cellKm)} km / ${+world.dt.toFixed(2)} Myr`;
  $("sMs").textContent=msEma?msEma.toFixed(1)+" ms":"–";
  const ev=world.events,last=world._lastN||0;
  if(ev.length&&ev[ev.length-1].n>last){
    const ul=$("log");if(loggedCount===0)ul.innerHTML="";
    for(const e of ev){if(e.n<=last)continue;const li=document.createElement("li");li.className=e.type;
      li.innerHTML=`<span class="t num">${fmtT(e.t)}${e.t>=1000?" Gyr":" Myr"}</span>${e.text}<span class="w num">${e.where}</span>`;ul.prepend(li);loggedCount++;}
    while(ul.children.length>80)ul.lastChild.remove();world._lastN=ev[ev.length-1].n;
  }
}
function updateLegend(){
  const L=$("legend");
  if(view==="elev")L.innerHTML=`<span class="num">−8 km</span><span class="bar" style="background:${gradCss(ELEV_STOPS,-8,8)}"></span><span class="num">+8 km</span>`;
  else if(view==="rock"||view==="l1")L.innerHTML=ROCKS.slice(1).map(r=>`<span><i class="swatch" style="background:${r.c}"></i>${r.name}</span>`).join("");
  else L.innerHTML=`<span><i class="swatch" style="background:var(--conv)"></i>Convergent</span><span><i class="swatch" style="background:var(--div)"></i>Divergent</span><span><i class="swatch" style="background:var(--trans)"></i>Transform</span>`;
}
function frame(now){
  let need=false;
  if(running){const n=+$("speed").value,t0=performance.now();let k=0;for(;k<n&&world.t<runTarget-1e-9;k++){world.step();snapTick();}
    if(k){const ms=(performance.now()-t0)/k;msEma=msEma?msEma*0.9+ms*0.1:ms;
      const fdt=(now-lastFrame)/1000;if(fdt>0&&fdt<1){const rate=k*world.dt/fdt;rateEma=rateEma?rateEma*0.9+rate*0.1:rate;}
      colorsDirty=true;need=true;if(handoff)staleHandoff();updateUi();if(selCell>=0)showColumn();l1WorldMoved();}
    if(world.t>=runTarget-1e-9)stopRun();}
  if(l1Busy()){l1Tick();lastFrame=now;requestAnimationFrame(frame);return;}
  if(need)render();lastFrame=now;requestAnimationFrame(frame);
}
// Run advances L0 by the Gyr slider from where it starts, then stops; Run again continues for that long
function stopRun(){running=false;logAction({a:"stop"});$("play").textContent="Run";render();l1SelectionMoved();}
$("play").onclick=()=>{if(l1Busy())return;if(running){stopRun();return;}
  if(world.t>=runTarget-1e-9){runStart=world.t;runTarget=runStart+1000*+$("gyr").value;}logAction({a:"run",untilMyr:runTarget});running=true;$("play").textContent="Pause";};
$("gyr").oninput=e=>{$("gyrVal").textContent=e.target.value;if(runTarget>runStart)runTarget=runStart+1000*+e.target.value;};
$("stepBtn").onclick=()=>{if(l1Busy())return;if(running)stopRun();l1WorldMoved();const t0=performance.now();world.step();snapTick();msEma=performance.now()-t0;colorsDirty=true;if(handoff)staleHandoff();render();updateUi();if(selCell>=0)showColumn();l1SelectionMoved();logAction({a:"step"});};
$("speed").oninput=e=>$("speedVal").textContent=e.target.value;
$("nPlates").oninput=e=>$("nPlatesVal").textContent=e.target.value;
$("contT").oninput=e=>$("contTVal").textContent=e.target.value+"%";
$("contT").onchange=build;
$("reset").onclick=build;
$("randSeed").onclick=()=>{const w=["basalt","craton","ophiolite","gneiss","tethys","rodinia","zircon","karst","obsidian","pelagic","moho","flysch","horst","graben","olivine","andesite"];
  $("seed").value=w[Math.floor(Math.random()*w.length)]+"-"+Math.floor(Math.random()*1000);build();};
$("seed").addEventListener("keydown",e=>{if(e.key==="Enter")build();});
["mode","grid"].forEach(id=>$(id).onchange=build);
["marks","traces"].forEach(id=>$(id).onchange=()=>render());
document.querySelectorAll("#views button").forEach(b=>b.onclick=()=>{view=b.dataset.v;
  document.querySelectorAll("#views button").forEach(x=>x.setAttribute("aria-pressed",x===b?"true":"false"));colorsDirty=true;updateLegend();render();});
document.querySelectorAll("#proj button").forEach(b=>b.onclick=()=>{proj=b.dataset.p;
  document.querySelectorAll("#proj button").forEach(x=>x.setAttribute("aria-pressed",x===b?"true":"false"));
  wrap.className="mapwrap "+proj;pixIdx=null;render();});
// Pointer gestures for a map element: one-finger drag, two-finger pinch or the wheel to zoom (about a point), and a tap
// (press and release without moving) to select
function gestures(el,{drag,zoom,tap}){const pts=new Map();let moved=false,sx=0,sy=0,lastD=0;
  const mid=()=>{const a=[...pts.values()];return[(a[0][0]+a[1][0])/2,(a[0][1]+a[1][1])/2];},dist=()=>{const a=[...pts.values()];return Math.hypot(a[0][0]-a[1][0],a[0][1]-a[1][1]);};
  el.addEventListener("pointerdown",e=>{try{el.setPointerCapture(e.pointerId);}catch(_){}pts.set(e.pointerId,[e.clientX,e.clientY]);
    if(pts.size===1){moved=false;sx=e.clientX;sy=e.clientY;}else{moved=true;lastD=dist();}});
  el.addEventListener("pointermove",e=>{const p=pts.get(e.pointerId);if(!p)return;const dx=e.clientX-p[0],dy=e.clientY-p[1];pts.set(e.pointerId,[e.clientX,e.clientY]);
    if(pts.size>=2){const d=dist();if(lastD>0&&d>0){const m=mid();zoom(d/lastD,m[0],m[1]);}lastD=d;return;}
    if(Math.abs(e.clientX-sx)+Math.abs(e.clientY-sy)>5)moved=true;if(moved)drag(dx,dy);});
  const up=e=>{const had=pts.delete(e.pointerId);if(had&&pts.size===0&&!moved&&e.type==="pointerup")tap(e);if(pts.size<2)lastD=0;};
  el.addEventListener("pointerup",up);el.addEventListener("pointercancel",up);
  el.addEventListener("wheel",e=>{e.preventDefault();zoom(Math.exp(-e.deltaY*(e.deltaMode===1?0.05:0.0015)),e.clientX,e.clientY);},{passive:false});}
let dragging=false;
gestures(wrap,{
  drag:(dx,dy)=>{if(proj==="globe"){const k=3/(wrap.clientWidth*gZoom);yaw-=dx*k;pitch=Math.max(-1.45,Math.min(1.45,pitch+dy*k));}
    else{const r=canvas.getBoundingClientRect();mC[0]-=dx/r.width*2*Math.PI/mZoom;mC[1]+=dy/r.height*2*G_YMAX/mZoom;clampMap();}camDirty=true;render();},
  zoom:(f,cx,cy)=>{if(proj==="globe")gZoom=Math.max(1,Math.min(8,gZoom*f));
    else{const r=canvas.getBoundingClientRect(),px=(cx-r.left)/r.width*pixW,py=(cy-r.top)/r.height*pixH,[x0,y0]=mapXY(px,py);mZoom=Math.max(1,Math.min(12,mZoom*f));
      const [x1,y1]=mapXY(px,py);mC[0]+=x0-x1;mC[1]+=y0-y1;clampMap();}
    camDirty=true;render();},
  tap:e=>inspect(e)});
function inspect(e){
  const r=canvas.getBoundingClientRect(),px=Math.floor((e.clientX-r.left)/r.width*pixW),py=Math.floor((e.clientY-r.top)/r.height*pixH);
  if(px<0||py<0||px>=pixW||py>=pixH)return;const i=pixIdx[py*pixW+px];if(i<0)return;
  selCell=i;selMask=null;showColumn();l1Select(i);logAction({a:"select",cell:i});render();
  const w=world,e2=w.elev[i],box=$("inspect");box.style.display="block";
  box.textContent=`${cellLoc(i)}: ${e2>=0?"elevation":"depth"} ${Math.abs(e2).toFixed(1)} km`;l0LastInspect={cell:i,text:box.textContent};
  clearTimeout(box._t);box._t=setTimeout(()=>box.style.display="none",3000);
}
function cellLoc(i){const p=world.cube.pos,lat=Math.asin(p[3*i+2])*180/Math.PI,lon=Math.atan2(p[3*i+1],p[3*i])*180/Math.PI;
  return `${Math.abs(lat).toFixed(0)}°${lat>=0?"N":"S"} ${Math.abs(lon).toFixed(0)}°${lon>=0?"E":"W"}`;}
function fmtAgo(a){a=Math.max(0,a);return a>=1000?(a/1000).toFixed(2)+" Gyr":Math.round(a)+" Myr";}
function fmtKm(v){const a=Math.abs(v);return a>=100?Math.round(a).toLocaleString("en-US"):a.toFixed(a>=10?0:1);}
const COMPASS=["north","northeast","east","southeast","south","southwest","west","northwest"];
function showColumn(){
  const body=$("colBody"),panel=$("colpanel");
  if(selCell<0||!world){panel.hidden=true;return;}
  const w=world,i=selCell,P=w.plates[w.owner[i]],e=w.elev[i],v=w.velocityAt(i);
  const units=w.pool.units(w.wCol[i]).reverse();
  const crust=w.wK[i]?"Continental crust":`Oceanic crust, ${fmtAgo(w.t-w.wB[i])} old`;
  let html=`<div class="colhead"><b>${cellLoc(i)}</b><span>${crust}, ${w.wT[i].toFixed(0)} km thick</span><span>${e>=0?"Elevation":"Depth"} ${Math.abs(e).toFixed(1)} km</span><span>Plate ${P.id}, moving ${v.speed.toFixed(1)} cm/yr toward the ${COMPASS[Math.round(v.bearing/45)%8]}</span></div><div class="strat">`;
  for(const u of units){const h=Math.round(8+46*Math.log10(1+u.th)/Math.log10(76));const R=ROCKS[u.rock]||ROCKS[1];
    html+=`<div class="unit"><div class="blk" style="height:${h}px;background:${R.c}"></div><div class="lbl">${R.name}<br><span class="num">${u.th<0.1?u.th.toFixed(2):u.th.toFixed(1)} km, formed ${fmtAgo(w.t-u.born)} ago</span></div></div>`;}
  html+=`</div>`;
  const eps=w.faultEpisodes(i);
  if(eps.length){html+=`<ul class="eps">`;
    for(const ep of eps.reverse()){const k=ep.kind==="strike-slip"?`Strike-slip (${ep.sense})`:ep.kind[0].toUpperCase()+ep.kind.slice(1);
      html+=`<li><b>${k}</b>, ${ep.active?"active now":`from ${fmtAgo(w.t-ep.t0)} to ${fmtAgo(w.t-ep.t1)} ago`}: <span class="num">${fmtKm(ep.across)} km ${ep.across>=0?"of convergence":"of extension"}, ${fmtKm(ep.along)} km sideways</span></li>`;}
    html+=`</ul>`;}
  else html+=`<p class="muted" style="margin:10px 0 0;font-size:13.5px">No fault slip recorded in this crust.</p>`;
  body.innerHTML=html;panel.hidden=false;
}
function staleHandoff(){handoff=null;$("download").hidden=true;$("tracesLbl").hidden=true;$("summary").innerHTML=`<div class="muted">The world has moved on since the last snapshot. Extract again to refresh the handoff.</div>`;}
function renderSummary(){
  const S=$("summary");if(!handoff){S.innerHTML="";$("download").hidden=true;$("tracesLbl").hidden=true;return;}
  const f=handoff.faults,by=[0,0,0,0],act=[0,0,0,0],len=[0,0,0,0];for(const x of f){by[x.kind]++;if(x.active)act[x.kind]++;len[x.kind]+=x.lengthKm;}
  const row=(k,name)=>`<div><span class="swatch" style="background:${["",'var(--conv)','var(--div)','var(--trans)'][k]}"></span>${by[k]} ${name} faults, ${act[k]} active, ${Math.round(len[k]).toLocaleString("en-US")} km of trace</div>`;
  S.innerHTML=`<div class="num"><b>Snapshot at ${fmtT(handoff.t)} ${handoff.t>=1000?"billion":"million"} years</b></div>${row(1,"thrust")}${row(2,"normal")}${row(3,"strike-slip")}
    <div>${world.N.toLocaleString("en-US")} rock columns, averaging ${world.meanUnits.toFixed(1)} units on continents</div>
    <div>${handoff.l1.forms.filter(f=>f.cells>=20).length} formations spanning 20 or more columns (${handoff.l1.forms.length} in all), in ${handoff.l1.nDomains} domains</div><div>${world.events.length} catalogued events, ${world.plates.filter(p=>p.alive).length} plate motions</div>
    <div class="muted" id="dlNote"></div>`;
  $("tracesLbl").hidden=false;
  if(downloadsP)downloadsP.then(d=>{$("download").hidden=!d;if(!d)$("dlNote").textContent="File downloads aren’t available in this view.";});
}
$("extract").onclick=()=>{if(running){running=false;$("play").textContent="Run";}
  handoff=world.extractFaults();handoff.l1=world.buildL1();logAction({a:"extract"});l1=handoff.l1;l1Time=performance.now();renderSummary();render();};
$("download").onclick=async()=>{const d=await downloadsP;if(!d||!handoff)return;
  const data=JSON.stringify(world.exportHandoff(handoff));
  const name=`l0-handoff-${($("seed").value||"genesis").replace(/[^a-z0-9-]/gi,"_")}-${Math.round(handoff.t)}myr.json`;
  try{await d.save({filename:name,data});$("dlNote").textContent="Saved.";}
  catch(err){$("dlNote").textContent=err&&err.code==="declined"?"Download cancelled.":"The download didn’t go through.";}};
downloadsP=(window.claude&&typeof window.claude.use==="function")?window.claude.use("downloads").catch(()=>null):Promise.resolve(null);
window.addEventListener("resize",()=>{pixIdx=null;render();l1Pix=null;if(l1St&&(l1St.sh||l1St.pv))l1Render();});
// ---------- L1 landmass ----------
// L0 keeps its whole state every 10 Myr (the last four), so a landmass picked at the present can be replayed
let snaps=[],l1St=null,l1View="terrain",l1Zoom=1,l1Pan=[0,0],l1Pix=null,l1PW=0,l1PH=0,l1Img=null,l1Geo=null,l1Ms=0,l1Mask=null;
const SNAP_EVERY=10,SNAP_KEEP=4;
function snapTick(){if(!world)return;const last=snaps[snaps.length-1];
  if(!last||world.t-last.t>=SNAP_EVERY-1e-9){snaps.push({t:world.t,steps:world.steps,snap:world.snapshot()});while(snaps.length>SNAP_KEEP)snaps.shift();}}
const l1Busy=()=>!!(l1St&&/^(replay|forward|run|basins)$/.test(l1St.phase));
const THICK_DIV=[[-15,"#2c5aa0"],[-5,"#9db8de"],[0,"#f2f0ea"],[5,"#e6a07a"],[15,"#b8322a"],[30,"#5e0f12"]];
const STRAIN_STOPS=[[-3,"#f3f1ea"],[-2,"#d9e2c0"],[-1.25,"#9cc38c"],[-0.5,"#3f9b8a"],[0,"#265c8f"],[0.5,"#2b1f5c"]];
const FIN_STOPS=[[-0.5,"#2c5aa0"],[-0.2,"#9db8de"],[0,"#f2f0ea"],[0.15,"#e9b08a"],[0.35,"#c4452f"],[0.6,"#5e0f12"]];
const DYN_STOPS=[[-1,"#2c5aa0"],[-0.4,"#9db8de"],[0,"#f2f0ea"],[0.3,"#e9b08a"],[0.6,"#b8322a"]];
const LUT_DYN=buildLut(DYN_STOPS,256,-1,0.6),RAIN_STOPS=[[-1.3,"#96703c"],[-0.7,"#dcc88c"],[-0.3,"#b4dc96"],[0.1,"#3ca078"],[0.5,"#145a8c"]];
const LUT_RAIN=buildLut(RAIN_STOPS,256,-1.3,0.5),LUT_FIN=buildLut(FIN_STOPS,256,-0.5,0.6),LUT_DS=buildLut(THICK_DIV,256,-15,30),LUT_EPS=buildLut(STRAIN_STOPS,256,-3,0.5);
function l1Reset(){l1St=null;l1Mask=null;l1Pix=null;$("l1Run").disabled=true;
  $("l1stage").hidden=true;$("l1controls").hidden=true;$("l1stats").innerHTML="";$("l1budget").innerHTML="";$("l1legend").innerHTML="";
  $("l1Info").innerHTML="";l1TileSel=-1;l1TileInfo();}
// the newest saved state at least a window before the present
const l1Entry=win=>{for(let k=snaps.length-1;k>=0;k--)if(snaps[k].t<=world.t-win-1e-9)return snaps[k];return null;};
// A finished run belongs to its moment of L0: once L0 moves on it is marked stale and kept until the selection changes.
function l1WorldMoved(){if(l1St&&l1St.phase==="done"&&!l1St.stale){l1St.stale=true;l1Status();}}
// the selected point's landmass, highlighted on L0's map (land strongly, its submerged rim faintly)
function selGlobeMask(){if(!world||selCell<0||l1Busy())return null;const key=selCell+":"+world.steps+":"+l1Buf();if(selMaskKey===key)return selMask;selMaskKey=key;
  const pick=L1.landmassAt(world,selCell,l1Buf());if(!pick)return selMask=null;const m=new Uint8Array(world.N);for(let h=0;h<pick.cells.length;h++)m[pick.cells[h]]=h<pick.land?2:1;return selMask=m;}
// Selection drives the L1 panel: a landmass shows at once, as L0 has it now (refining replays from a window earlier).
// A finished run stays while the selection is on its landmass.
function l1Select(i,force){if(l1Busy()||!world)return;
  const S=l1St;if(!force&&S&&S.phase==="done"&&!S.stale&&S.pick.cells.indexOf(i)>=0&&S.pick.cells.indexOf(i)<S.pick.land){l1Status();return;}
  const pick=L1.landmassAt(world,i,l1Buf());
  if(!pick){l1Reset();$("l1Info").innerHTML=`<div class="muted">The selected point is under the sea.</div>`;return;}
  l1Reset();l1St={phase:"picked",pick,cell:i,pv:l1MakePreview(pick)};l1Zoom=1;l1Pan=[0,0];l1View="terrain";
  $("l1stage").hidden=false;$("l1controls").hidden=false;l1SetViews();l1Legend();l1Status();l1Render();}
// after L0 moves (a step, or a run stopping): an unrefined selection follows L0
function l1SelectionMoved(){if(selCell>=0&&(!l1St||l1St.phase==="picked"))l1Select(selCell,true);else l1Status();}
// The preview frame: azimuthal equidistant about the landmass's centroid now, north up, in L0's fixed frame
function l1MakePreview(pick){const pos=world.cube.pos;let cx=0,cy=0,cz=0;for(let h=0;h<pick.land;h++){const c=pick.cells[h];cx+=pos[3*c];cy+=pos[3*c+1];cz+=pos[3*c+2];}
  const L=Math.hypot(cx,cy,cz)||1,c=[cx/L,cy/L,cz/L];let nr=[-c[2]*c[0],-c[2]*c[1],1-c[2]*c[2]],ln=Math.hypot(...nr);if(ln<1e-9){nr=[1,0,0];ln=1;}nr=nr.map(v=>v/ln);
  const ea=[nr[1]*c[2]-nr[2]*c[1],nr[2]*c[0]-nr[0]*c[2],nr[0]*c[1]-nr[1]*c[0]],proj=new L1.Projection(c,ea,nr),xy=[0,0],b=[0,0,0];
  let x0=Infinity,x1=-Infinity,y0=Infinity,y1=-Infinity;for(const k of pick.cells){b[0]=pos[3*k];b[1]=pos[3*k+1];b[2]=pos[3*k+2];proj.fwd(b,xy);x0=Math.min(x0,xy[0]);x1=Math.max(x1,xy[0]);y0=Math.min(y0,xy[1]);y1=Math.max(y1,xy[1]);}
  const h=50,mg=400;x0=Math.floor((x0-mg)/h)*h;y0=Math.floor((y0-mg)/h)*h;
  const m=new Uint8Array(world.N);for(let q=0;q<pick.cells.length;q++)m[pick.cells[q]]=q<pick.land?2:1;
  return{proj,mask:m,x0,y0,nx:Math.ceil((x1+mg-x0)/h),ny:Math.ceil((y1+mg-y0)/h),o:{h}};}
// L0's surface at a world direction (bilinear on the cube grid)
function l0Elev(x,y,z){const w=world,out=[0,0,0],nb=w.cube.nb;w.cube.lookupFrac(x,y,z,out);
  const c=out[0],tx=out[1],ty=out[2],cB=nb[4*c+(tx>=0?0:1)],cC=nb[4*c+(ty>=0?2:3)],cD=nb[4*cB+(ty>=0?2:3)],ax=Math.abs(tx),ay=Math.abs(ty),E=w.elev;
  return(1-ax)*(1-ay)*E[c]+ax*(1-ay)*E[cB]+(1-ax)*ay*E[cC]+ax*ay*E[cD];}
// status lines under the L1 button; the button refines the current selection when a saved state is old enough
function l1Status(){const S=l1St,el=$("l1Info");if(!world||l1Busy())return;const win=+$("l1Win").value,entry=l1Entry(win);
  const pick=S&&S.pick;$("l1Run").disabled=!(pick&&entry&&selCell>=0);if(!pick)return;const a=pick.land*world.cellKm*world.cellKm;
  let html=`<div><b>Landmass</b>: ${(a/1e6).toFixed(1)} M km² of land (${pick.land.toLocaleString("en-US")} L0 cells), mostly on plate ${pick.ref}${pick.islands?`, with ${pick.islands} island${pick.islands>1?"s":""} within ${pick.bufferKm} km`:""}; with the buffer and enclosed seas, ${((pick.cells.length-pick.land)*world.cellKm*world.cellKm/1e6).toFixed(1)} M km² of sea and shelf.</div>`;
  if(S.phase==="picked")html+=`<div class="muted">${entry?`Showing L0 now; refining replays from ${fmtT(entry.t)} Myr.`:`No saved L0 state ${win} Myr back yet: run L0 longer or shorten the window.`}</div>`;
  else if(S.phase==="done")html+=`<div class="muted">Refined over ${fmtT(S.plan.window)} Myr from L0 at ${fmtT(S.entry.t)} Myr${S.stale?"; L0 has moved on since":""}.</div>`;
  el.innerHTML=html;}
// the buffer changes which landmass is selected: the L0 highlight and an unrefined preview follow it
function l1Buf(){return +$("l1Buf").value;}
$("l1Buf").oninput=e=>{$("l1BufVal").textContent=(+e.target.value).toLocaleString("en-US")+" km";if(l1Busy())return;selMaskKey="";if(selCell>=0&&(!l1St||l1St.phase==="picked"))l1Select(selCell,true);render();};
$("l1Win").oninput=e=>{$("l1WinVal").textContent=e.target.value+" Myr";l1Status();};
$("l1Ell").oninput=e=>$("l1EllVal").textContent=(+e.target.value).toLocaleString("en-US")+" km";
$("l1Run").onclick=()=>{if(l1Busy()||selCell<0)return;const win=+$("l1Win").value,entry=l1Entry(win);if(!entry)return;
  const pick=l1St&&l1St.phase==="picked"?l1St.pick:L1.landmassAt(world,selCell,l1Buf());if(!pick)return;
  if(running)stopRun();$("l1Run").disabled=true;
  const o={h:+$("l1H").value,Ar:+$("l1Ar").value,ellKm:+$("l1Ell").value};let rp;
  try{rp=new L1.Replay(world,entry,win,pick);}catch(err){$("l1Info").insertAdjacentHTML("beforeend",`<div class="muted">Could not start: ${err.message}</div>`);return;}
  l1TileSel=-1;l1TileInfo();$("l1stats").innerHTML="";$("l1budget").innerHTML="";
  l1St={phase:"replay",pick,cell:selCell,rp,entry,win,o,t0Wall:performance.now()};logAction({a:"refine",cell:selCell,bufferKm:pick.bufferKm,entrySteps:entry.steps,windowMyr:win,sheet:o});l1Mask=null;$("exportBtn").disabled=true;
  $("l1Info").innerHTML=`<div class="muted" id="l1Prog">Replaying L0 to trace the landmass’s crust…</div>`;};
function l1Prog(t){const el=$("l1Prog");if(el)el.textContent=t;}
// one slice of work per animation frame: replay, replay again to the window's start, the sheet a step at a time, the basins
function l1Tick(){const S=l1St,w=world;colorsDirty=true;
  if(S.phase==="replay"){const done=S.rp.advance(120);l1Prog(`Pass 1, tracing the crust: L0 at ${fmtT(w.t)} Myr (present ${fmtT(S.rp.endT)})…`);
    if(!done)return;S.plan=S.rp.finish();S.rp=null;w.restore(S.entry.snap);w.trackParcels();S.phase="forward";return;}
  if(S.phase==="forward"){const t1=performance.now();while(w.steps<S.plan.t0Steps&&performance.now()-t1<120)w.step();l1Prog(`Pass 2: L0 back at ${fmtT(w.t)} Myr, heading for the window’s start…`);
    if(w.steps<S.plan.t0Steps)return;
    try{S.sh=new L1.Sheet(w,S.plan,S.o);}catch(err){l1Prog(`Could not build the sheet: ${err.message}`);w.untrackParcels();l1St=null;$("exportBtn").disabled=false;return;}
    S.phase="run";l1Zoom=1;l1Pan=[0,0];l1Pix=null;l1View="terrain";l1SetViews();$("l1stage").hidden=false;$("l1controls").hidden=false;
    $("l1Info").insertAdjacentHTML("beforeend",`<div class="muted">Sheet: ${S.sh.nx}×${S.sh.ny} cells of ${S.o.h} km, ${S.sh.nm.toLocaleString("en-US")} markers carrying ${S.plan.nObj.toLocaleString("en-US")} L0 parcels.</div>`);
    l1Legend();l1Render();l1Stats();return;}
  if(S.phase==="run"){const t0=performance.now();S.sh.advance();const ms=performance.now()-t0;l1Ms=l1Ms?l1Ms*0.7+ms*0.3:ms;
    l1Prog(`Pass 2, the sheet: ${fmtT(S.sh.t)} of ${fmtT(S.plan.window)} Myr…`);updateUi();if(selCell>=0)showColumn();l1Pix=null;l1Render();l1Stats();
    if(w.steps>=S.plan.endSteps)S.phase="basins";return;}
  if(S.phase==="basins"){l1Prog("Tectonics frozen: 5 Myr of erosion on a 10 km grid to settle the drainage basins…");
    // the frame draws the message first; the pass itself runs on the next one
    if(!S.basinsDue){S.basinsDue=true;return;}
    const t0=performance.now();S.bs=L1.basins(S.sh);S.bsMs=performance.now()-t0;S.ho=L1.handoff(S.sh,S.bs);S.hoMs=performance.now()-t0-S.bsMs;w.untrackParcels();l1TileSel=-1;l1TileInfo();
    S.phase="done";$("exportBtn").disabled=false;l1Mask=null;l1View="basins";l1SetViews();l1Pix=null;l1Legend();l1Render();l1Stats();l1Status();updateUi();}}
function l1SetViews(){const has=!!(l1St&&l1St.bs),sheet=!!(l1St&&l1St.sh);document.querySelectorAll("#l1views button").forEach(x=>{if(x.dataset.v!=="terrain")x.disabled=!sheet;if(x.dataset.v==="basins"||x.dataset.v==="rain")x.disabled=!has;if(x.dataset.v==="tiles")x.disabled=!(l1St&&l1St.ho);x.setAttribute("aria-pressed",x.dataset.v===l1View?"true":"false");});}
// canvas geometry: the sheet's grid, zoomed and panned; per pixel its map-plane point
function l1Size(){const sh=l1St.sh||l1St.pv,wrap=$("l1wrap"),X0=sh.x0,X1=sh.x0+sh.nx*sh.o.h,Y0=sh.y0,Y1=sh.y0+sh.ny*sh.o.h,asp=Math.min(1.3,Math.max(0.45,(Y1-Y0)/(X1-X0)));
  wrap.style.height=Math.round(wrap.clientWidth*asp)+"px";const dpr=Math.min(2,window.devicePixelRatio||1),W=Math.min(1100,Math.round(wrap.clientWidth*dpr)),H=Math.round(W*asp);
  if(W!==l1PW||H!==l1PH||!l1Img){l1PW=W;l1PH=H;const cv=$("l1map");cv.width=W;cv.height=H;l1Img=cv.getContext("2d").createImageData(W,H);l1Pix=null;}
  const s0=Math.max((X1-X0)/W,(Y1-Y0)/H),bx=(X0+X1)/2,by=(Y0+Y1)/2;l1Geo={s:s0/l1Zoom,s0,bx,by,cx:bx+l1Pan[0],cy:by+l1Pan[1]};}
const l1XY=(px,py)=>[l1Geo.cx+(px+0.5-l1PW/2)*l1Geo.s,l1Geo.cy-(py+0.5-l1PH/2)*l1Geo.s];
// per pixel: the field values to draw (recomputed when the sheet moves or the view changes)
function l1Fields(){if(!l1St.sh)return l1PreviewFields();const S=l1St,sh=S.sh,bs=S.bs,W=l1PW,H=l1PH,n=W*H,h=sh.o.h,{nx,ny}=sh;
  const P=l1Pix={cell:new Int32Array(n),fine:new Int32Array(n),e:new Float32Array(n),back:new Float32Array(n)};
  for(let py=0;py<H;py++)for(let px=0;px<W;px++){const q=py*W+px,[x,y]=l1XY(px,py),i=Math.floor((x-sh.x0)/h),j=Math.floor((y-sh.y0)/h);
    const c=i>=0&&j>=0&&i<nx&&j<ny?j*nx+i:-1;P.cell[q]=c>=0&&sh.act[c]?c:-1;P.fine[q]=-1;P.back[q]=NaN;
    if(bs&&P.cell[q]>=0){const fi=Math.floor((x-bs.x0)/bs.h),fj=Math.floor((y-bs.y0)/bs.h);if(fi>=0&&fj>=0&&fi<bs.nx&&fj<bs.ny)P.fine[q]=fj*bs.nx+fi;}
    // land: the drainage pass's 10 km surface, bilinear between its land cells (drawn cell by cell, the hillshade lit
    // only the cells' edges, a crosshatch); land only, so the coast stays where the pass has it
    if(P.fine[q]>=0&&bs.e[P.fine[q]]>0){const gx=(x-bs.x0)/bs.h-0.5,gy=(y-bs.y0)/bs.h-0.5,i0=Math.floor(gx),j0=Math.floor(gy),fx=gx-i0,fy=gy-j0;let s=0,ws=0;
      for(let k=0;k<4;k++){const ii=i0+(k&1),jj=j0+(k>>1);if(ii<0||jj<0||ii>=bs.nx||jj>=bs.ny)continue;const v=bs.e[jj*bs.nx+ii];if(!(v>0))continue;const ww=(k&1?fx:1-fx)*(k>>1?fy:1-fy);s+=ww*v;ws+=ww;}
      P.e[q]=ws>0?s/ws:bs.e[P.fine[q]];}
    else{// the sheet's surface, bilinear over active neighbours; at sea it fades into the plain ocean outside the crust
      // over the edge cell, so the crust's 50 km cells don't show as a staircase in the sea
      const gx=(x-sh.x0)/h-0.5,gy=(y-sh.y0)/h-0.5,i0=Math.floor(gx),j0=Math.floor(gy),fx=gx-i0,fy=gy-j0;let s=0,ws=0;
      for(const [di,dj,ww] of [[0,0,(1-fx)*(1-fy)],[1,0,fx*(1-fy)],[0,1,(1-fx)*fy],[1,1,fx*fy]]){const ii=i0+di,jj=j0+dj;if(ii<0||jj<0||ii>=nx||jj>=ny)continue;const cc=jj*nx+ii;if(!sh.act[cc])continue;s+=ww*sh.e[cc];ws+=ww;}
      const v=ws>0?s/ws:P.cell[q]>=0?sh.e[P.cell[q]]:NaN;P.e[q]=P.fine[q]>=0||v<=0&&ws<1?Math.min(-0.01,ws*v+(1-ws)*L1_OUT_KM):P.cell[q]>=0?v:NaN;}}}
// before refining: L0's present surface in the preview frame, the landmass (and its submerged rim) as crust, the rest as backdrop
function l1PreviewFields(){const pv=l1St.pv,W=l1PW,H=l1PH,n=W*H,b=[0,0,0],lk=world.cube.lookup;
  const P=l1Pix={cell:new Int32Array(n).fill(-1),fine:new Int32Array(n).fill(-1),e:new Float32Array(n),back:new Float32Array(n)};
  for(let py=0;py<H;py++)for(let px=0;px<W;px++){const q=py*W+px,[x,y]=l1XY(px,py);pv.proj.inv(x,y,b);const v=l0Elev(b[0],b[1],b[2]);
    if(pv.mask[lk(b[0],b[1],b[2])]){P.e[q]=v;P.back[q]=NaN;}else{P.e[q]=NaN;P.back[q]=v;}}}
const L1_OUT_KM=-5,L1_OUT=(()=>{const k=Math.round((L1_OUT_KM+8)/16*1023);return[LUT_E[k*3],LUT_E[k*3+1],LUT_E[k*3+2]];})();
function l1Render(){if(!l1St||!(l1St.sh||l1St.pv))return;l1Size();if(!l1Pix)l1Fields();
  const S=l1St,sh=S.sh,bs=S.bs,W=l1PW,H=l1PH,d=l1Img.data,P=l1Pix,E=P.e,ex=10/l1Geo.s,rivers=$("l1Rivers").checked;
  const rock=l1View==="rock"?sh.surfaceRock():null,sf=l1View==="finite"?sh.strainField():null,fAmin=2e4;
  const bcol=new Map(),bsize=bs?l1BasinSizes():null;
  for(let py=0;py<H;py++)for(let px=0;px<W;px++){const q=py*W+px,o=4*q,c=P.cell[q],f=P.fine[q];let r,g,b;
    // outside L1's crust: before refining, L0's own surface, land greyed; after, plain ocean (L0's coast there is
    // not L1's, and drawing both side by side shows two shorelines)
    if(isNaN(E[q])){
      const hb=P.back[q];if(isNaN(hb)){r=L1_OUT[0];g=L1_OUT[1];b=L1_OUT[2];}else{let k=Math.round((hb+8)/16*1023);k=k<0?0:k>1023?1023:k;r=LUT_E[k*3];g=LUT_E[k*3+1];b=LUT_E[k*3+2];
        if(hb>0){const m=(r+g+b)/3;r=0.35*r+0.65*m;g=0.35*g+0.65*m;b=0.35*b+0.65*m;}r*=0.72;g*=0.72;b*=0.72;}
      d[o]=r;d[o+1]=g;d[o+2]=b;d[o+3]=255;continue;}
    const h=E[q],land=h>0;
    // hillshade from the drawn surface (vertical exaggeration 10x)
    const eL=px>0&&!isNaN(E[q-1])?E[q-1]:h,eR=px<W-1&&!isNaN(E[q+1])?E[q+1]:h,eU=py>0&&!isNaN(E[q-W])?E[q-W]:h,eD=py<H-1&&!isNaN(E[q+W])?E[q+W]:h;
    const gx=(eR-eL)*ex/2,gy=(eU-eD)*ex/2,nl=1/Math.sqrt(1+gx*gx+gy*gy);let sh_=land?Math.max(0.35,Math.min(1.25,(0.55*-gx+0.55*gy+0.63)*nl/0.63)):1;
    if(l1View==="terrain"||!land&&l1View!=="thick"&&l1View!=="strain"){let k=Math.round((h+8)/16*1023);k=k<0?0:k>1023?1023:k;r=LUT_E[k*3]*sh_;g=LUT_E[k*3+1]*sh_;b=LUT_E[k*3+2]*sh_;}
    else if(l1View==="thick"){const v=c>=0?sh.dS[c]:0;let k=Math.round((v+15)/45*255);k=k<0?0:k>255?255:k;const s2=0.8+0.2*sh_;r=LUT_DS[k*3]*s2;g=LUT_DS[k*3+1]*s2;b=LUT_DS[k*3+2]*s2;}
    // finite strain: shortening (1 - least stretch) in reds; extension in blues only where there is next to no
    // shortening (crust squeezed out sideways stretches as much as it shortens, and counts as shortened)
    else if(l1View==="finite"){let v=0;if(c>=0&&!isNaN(sf.lmin[c])){const sh1=1-sf.lmin[c],ex=sf.lmax[c]-1;v=sh1>=0.05?sh1:ex>0.05?-ex:0;}
      let k=Math.round((v+0.5)/1.1*255);k=k<0?0:k>255?255:k;const s2=0.8+0.2*sh_;r=LUT_FIN[k*3]*s2;g=LUT_FIN[k*3+1]*s2;b=LUT_FIN[k*3+2]*s2;}
    else if(l1View==="strain"){const v=c>=0?Math.log10(Math.max(1e-3,sh.eps[c])):-3;let k=Math.round((v+3)/3.5*255);k=k<0?0:k>255?255:k;const s2=0.8+0.2*sh_;r=LUT_EPS[k*3]*s2;g=LUT_EPS[k*3+1]*s2;b=LUT_EPS[k*3+2]*s2;}
    else if(l1View==="rock"){const R=ROCK_RGB[c>=0?rock[c]:0]||[128,128,128];r=R[0]*sh_;g=R[1]*sh_;b=R[2]*sh_;}
    else if(l1View==="dyn"){const v=c>=0?sh.dyn[c]:0;let k=Math.round((v+1)/1.6*255);k=k<0?0:k>255?255:k;const s2=0.8+0.2*sh_;r=LUT_DYN[k*3]*s2;g=LUT_DYN[k*3+1]*s2;b=LUT_DYN[k*3+2]*s2;}
    else if(l1View==="rain"){if(f>=0&&bs.land[f]){let k=Math.round((Math.log10(Math.max(0.05,bs.climate.P[f]))+1.3)/1.8*255);k=k<0?0:k>255?255:k;r=LUT_RAIN[k*3]*sh_;g=LUT_RAIN[k*3+1]*sh_;b=LUT_RAIN[k*3+2]*sh_;}
      else{let k=Math.round((h+8)/16*1023);k=k<0?0:k>1023?1023:k;r=LUT_E[k*3];g=LUT_E[k*3+1];b=LUT_E[k*3+2];}}
    else if(l1View==="basins"){if(f>=0&&bs.land[f]){const id=bs.closedB[f]===2?-2:bs.basin[f];let col=bcol.get(id);
        // basins draining inland in muted earth tones, those reaching the sea in clear colours
        if(!col){const big=bsize.get(id)>=l1BasinMin(),u=(Math.imul(id+1,2654435761)>>>0)/4294967296;col=bs.closedB[f]===2?[226,208,164]:bs.closedB[f]?hsl(25+30*u,35,52+14*((u*13)%1)):big?hsl(u*360,45+20*((u*7)%1),58+12*((u*13)%1)):[205,200,188];bcol.set(id,col);}
        r=col[0]*sh_;g=col[1]*sh_;b=col[2]*sh_;}
      else{let k=Math.round((h+8)/16*1023);k=k<0?0:k>1023?1023:k;r=LUT_E[k*3];g=LUT_E[k*3+1];b=LUT_E[k*3+2];}}
    // L2 tiles: one colour per tile (the selected one bright, the rest faded when one is selected); edges with another
    // tile dark where the divide is firm (set by relief the sheet resolved), pale where it is contingent
    else if(l1View==="tiles"){const ho=S.ho;if(f>=0&&bs.land[f]&&ho.tileOf[f]>=0){const k=ho.tileOf[f],u=(Math.imul(k+7,2654435761)>>>0)/4294967296;let col=hsl(u*360,40+20*((u*7)%1),60+12*((u*13)%1));
        if(l1TileSel>=0&&k!==l1TileSel){const m=(col[0]+col[1]+col[2])/3;col=[0.3*col[0]+0.7*m,0.3*col[1]+0.7*m,0.3*col[2]+0.7*m];}
        r=col[0]*sh_;g=col[1]*sh_;b=col[2]*sh_;const dv=ho.divide[f];if(dv===2){r=45;g=38;b=34;}else if(dv===1){r=0.45*r+140;g=0.45*g+140;b=0.45*b+140;}}
      else{let k=Math.round((h+8)/16*1023);k=k<0?0:k>1023?1023:k;r=LUT_E[k*3];g=LUT_E[k*3+1];b=LUT_E[k*3+2];}}
    // rivers carrying what 20,000 km² of land at 0.3 m/yr of runoff would (the drainage pass draws them; the sheet
    // has none), and lakes: open (blue), closed (pale blue), dry lake floors (white)
    if(rivers&&land&&f>=0&&(l1View==="terrain"||l1View==="basins"||l1View==="rain"||l1View==="tiles")){const A=bs.climate?bs.Q[f]/0.3:bs.A[f],a=A>=fAmin?Math.min(0.9,0.4+0.25*Math.log10(A/fAmin)):0;
      if(a>0){r=r*(1-a)+40*a;g=g*(1-a)+95*a;b=b*(1-a)+170*a;}}
    if(land&&f>=0&&bs.lake[f]&&l1View!=="thick"&&l1View!=="strain"&&l1View!=="finite"&&l1View!=="dyn"){const L=[null,[52,104,190],[120,196,214],[238,234,222],[214,178,112]][bs.lake[f]];r=L[0];g=L[1];b=L[2];}
    d[o]=r;d[o+1]=g;d[o+2]=b;d[o+3]=255;}
  const ctx=$("l1map").getContext("2d");ctx.putImageData(l1Img,0,0);
  // ticks along the shortening axis where crust has shortened by 10% or more, on a lattice of cells about 26 px apart
  if(sf){const {nx,ny}=sh,h=sh.o.h,step=Math.max(1,Math.round(26*l1Geo.s/h)),len=0.42*step*h/l1Geo.s,dpr=l1PW/$("l1wrap").clientWidth;
    ctx.lineCap="round";for(const [col,wd] of [["rgba(255,255,255,.75)",3.2],["rgba(20,20,20,.9)",1.4]]){ctx.strokeStyle=col;ctx.lineWidth=wd*dpr;ctx.beginPath();
      for(let j=Math.floor(step/2);j<ny;j+=step)for(let i=Math.floor(step/2);i<nx;i+=step){const c=j*nx+i;if(isNaN(sf.lmin[c])||1-sf.lmin[c]<0.1)continue;
        const x=sh.x0+(i+.5)*h,y=sh.y0+(j+.5)*h,px=(x-l1Geo.cx)/l1Geo.s+l1PW/2,py=-(y-l1Geo.cy)/l1Geo.s+l1PH/2,dx=Math.cos(sf.dir[c])*len/2,dy=-Math.sin(sf.dir[c])*len/2;
        if(px<0||py<0||px>l1PW||py>l1PH)continue;ctx.moveTo(px-dx,py-dy);ctx.lineTo(px+dx,py+dy);}ctx.stroke();}}}
// an axis's bearing, e.g. N40°E (angle in radians from map east, counter-clockwise; axes have no sense)
function l1Bearing(a){let az=(90-a*180/Math.PI)%180;if(az<0)az+=180;return az<=90?`N${Math.round(az)}°E`:`N${Math.round(180-az)}°W`;}
// basins worth a colour of their own: at least 10,000 km² (100 cells of 10 km)
const l1BasinMin=()=>Math.round(1e4/(l1St.bs.h*l1St.bs.h));
function l1BasinSizes(){const bs=l1St.bs;if(bs._size)return bs._size;const m=new Map();for(let f=0;f<bs.land.length;f++)if(bs.land[f]){const id=bs.basin[f];m.set(id,(m.get(id)||0)+1);}return bs._size=m;}
function l1Legend(){const L=$("l1legend"),bar=(st,lo,hi,a,b)=>`<span class="num">${a}</span><span class="bar" style="background:${gradCss(st,lo,hi)}"></span><span class="num">${b}</span>`;
  if(l1View==="terrain")L.innerHTML=bar(ELEV_STOPS,-8,8,"−8 km","+8 km");
  else if(l1View==="thick")L.innerHTML=bar(THICK_DIV,-15,30,"−15","+30 km");
  else if(l1View==="finite")L.innerHTML=bar(FIN_STOPS,-0.5,0.6,"50% ext.","60% short.");
  else if(l1View==="strain")L.innerHTML=bar(STRAIN_STOPS,-3,0.5,"0.001","3 /Myr");
  else if(l1View==="rock")L.innerHTML=ROCKS.slice(1).map(r=>`<span><i class="swatch" style="background:${r.c}"></i>${r.name}</span>`).join("");
  else if(l1View==="dyn")L.innerHTML=bar(DYN_STOPS,-1,0.6,"−1 km","+0.6 km");
  else if(l1View==="rain")L.innerHTML=bar(RAIN_STOPS,-1.3,0.5,"50 mm","3 m/yr");
  else if(l1View==="tiles")L.innerHTML=`<span><i class="swatch" style="background:#2d2622"></i>firm divide</span><span><i class="swatch" style="background:#d8d4cc"></i>contingent divide</span>`;
  else L.innerHTML=`<span><i class="swatch" style="background:#3468be"></i>lake</span><span><i class="swatch" style="background:#78c4d6"></i>no outlet</span><span><i class="swatch" style="background:#eeeade"></i>dry lake floor</span><span><i class="swatch" style="background:#e2d0a4"></i>areic</span>`;}
function l1MaxShort(sh){const sf=sh.strainField();let m=0;for(let c=0;c<sf.lmin.length;c++)if(!isNaN(sf.lmin[c])&&1-sf.lmin[c]>m)m=1-sf.lmin[c];return Math.round(100*m);}
function l1Stats(){const S=l1St;if(!S||!S.sh)return;const sh=S.sh,h=sh.history[sh.history.length-1],h0=sh.history[0],b=sh.budget,bal=sh.balance(),M=v=>`${(v/1e6).toFixed(2)} M km³`;
  $("l1stats").innerHTML=`<div>L1 time<b class="num">${sh.t.toFixed(1)} / ${S.plan.window.toFixed(1)} Myr</b></div><div>Land<b class="num">${(h.landKm2/1e6).toFixed(1)} M km²</b></div>
    <div>Mean / max elevation<b class="num">${h.meanElevKm.toFixed(2)} / ${h.maxElevKm.toFixed(1)} km</b></div><div>At start<b class="num">${h0.meanElevKm.toFixed(2)} / ${h0.maxElevKm.toFixed(1)} km</b></div>
    <div>Most thickening<b class="num">${h.maxThickeningKm.toFixed(1)} km</b></div><div>Peak strain rate<b class="num">${h.maxStrainRate.toFixed(2)} /Myr</b></div><div>Most shortening<b class="num">${l1MaxShort(sh)}%</b></div>
    <div>Sheet dx / dt<b class="num">${sh.o.h} km / ${+(sh.world.dt/Math.max(1,h.substeps||1)).toFixed(2)} Myr</b></div>${S.bs?`<div>Drainage pass dx / dt<b class="num">${S.bs.h} km / ${+(S.bs.myr/(S.bs.steps||5)).toFixed(2)} Myr</b></div>`:""}
    <div>Compute per step<b class="num">${l1Ms?Math.round(l1Ms)+" ms":"–"}</b></div><div>Markers<b class="num">${sh.nm.toLocaleString("en-US")}</b></div>`;
  let html=`<div><b>Crust budget</b>: started ${M(b.initial)}; ${M(b.arcsHotspots)} added by arcs and hotspots, ${M(b.eroded)} eroded, ${M(b.left)} carried off the grid; now ${M(bal.now)}, off by ${(Math.abs(bal.error)*100).toExponential(1)}%.</div>`;
  if(S.bs){const bs=S.bs,sz=l1BasinSizes(),min=l1BasinMin(),big=[...sz.values()].filter(v=>v>=min).sort((a,b2)=>b2-a),cell=bs.h*bs.h,inBig=big.reduce((a,v)=>a+v,0);
    html+=`<div><b>Drainage basins</b>: ${(100*bs.inlandFrac).toFixed(1)}% of the land drains inland: ${(100*(bs.inlandFrac-bs.areicFrac)).toFixed(1)}% to ${bs.nClosed} lakes with no outlet and ${bs.nDry} rivers that run dry, and ${(100*bs.areicFrac).toFixed(1)}% is desert whose water never reaches a river at all${bs.nSea?` (${(bs.nSea*bs.h*bs.h/1e3).toFixed(0)} thousand km² of sea cut off from the ocean)`:""}. ${big.length} of at least 10,000 km² cover ${(100*inBig/bs.landCells).toFixed(0)}% of the land; the largest are ${big.slice(0,5).map(v=>(v*cell/1e6).toFixed(2)).join(", ")} M km². ${(bs.nBasins-big.length).toLocaleString("en-US")} smaller ones drain the rest, mostly short coastal streams. (${bs.nx}×${bs.ny} cells of ${bs.h} km, ${Math.round(S.bsMs)} ms.)</div>`;}
  $("l1budget").innerHTML=html;}
document.querySelectorAll("#l1views button").forEach(b=>b.onclick=()=>{if(b.disabled)return;l1View=b.dataset.v;l1SetViews();l1Legend();l1Render();});
$("l1Rivers").onchange=()=>l1Render();
// the L1 map: drag to pan, wheel or pinch to zoom about the pointer, tap to inspect
gestures($("l1wrap"),{
  drag:(dx,dy)=>{if(!l1Geo||l1Zoom<=1)return;const k=l1Geo.s*l1PW/$("l1wrap").clientWidth;l1Pan[0]-=dx*k;l1Pan[1]+=dy*k;l1Pix=null;l1Render();},
  zoom:(f,cx,cy)=>{if(!l1St||!l1Geo)return;const r=$("l1map").getBoundingClientRect(),px=(cx-r.left)/r.width*l1PW,py=(cy-r.top)/r.height*l1PH,[x,y]=l1XY(px,py);
    l1Zoom=Math.max(1,Math.min(32,l1Zoom*f));if(l1Zoom===1)l1Pan=[0,0];
    else{const s1=l1Geo.s0/l1Zoom;l1Pan=[x-(px+0.5-l1PW/2)*s1-l1Geo.bx,y+(py+0.5-l1PH/2)*s1-l1Geo.by];}
    l1Pix=null;l1Render();},
  tap:e=>l1Inspect(e)});
function l1Inspect(e){if(!l1St||!l1Pix)return;const cv=$("l1map"),r=cv.getBoundingClientRect(),px=Math.floor((e.clientX-r.left)/r.width*l1PW),py=Math.floor((e.clientY-r.top)/r.height*l1PH);
  if(px<0||py<0||px>=l1PW||py>=l1PH)return;const q=py*l1PW+px,c=l1Pix.cell[q],f=l1Pix.fine[q],sh=l1St.sh,bs=l1St.bs,box=$("l1inspect"),h=l1Pix.e[q];
  if(isNaN(h)){const hb=l1Pix.back[q];if(isNaN(hb)){box.textContent="Outside the landmass’s crust: not modelled at L1";return;}box.textContent=`Not the landmass’s crust; L0 ${hb>=0?"elevation":"depth"} ${Math.abs(hb).toFixed(2)} km`;}
  else{let t=`${h>=0?"Elevation":"Depth"} ${Math.abs(h).toFixed(2)} km`;
    if(c>=0&&sh){const rk=sh.surfaceRock()[c];t+=`, crust ${sh.S[c].toFixed(1)} km (${sh.dS[c]>=0?"+":"−"}${Math.abs(sh.dS[c]).toFixed(1)} in the window), strain rate ${sh.eps[c].toFixed(3)}/Myr, mantle pull ${Math.round(sh.dyn[c]*1000)} m${rk?", "+ROCKS[rk].name.toLowerCase()+" at the surface":""}`;const sf=sh.strainField();if(!isNaN(sf.lmin[c])){const s1=1-sf.lmin[c],e1=sf.lmax[c]-1;t+=`; ${s1>=0.02?`shortened ${Math.round(100*s1)}% along ${l1Bearing(sf.dir[c])}`:"little shortened"}${e1>=0.02?`, stretched ${Math.round(100*e1)}% across`:""}, rotated ${Math.abs(Math.round(sf.rot[c]*180/Math.PI))}° ${sf.rot[c]>=0?"anticlockwise":"clockwise"}`;}}
    if(f>=0&&bs&&bs.land[f]){const sz=l1BasinSizes().get(bs.basin[f])*bs.h*bs.h;t+=`; drains ${Math.round(bs.A[f]).toLocaleString("en-US")} km², in a basin of ${Math.round(sz).toLocaleString("en-US")} km² that ${["reaches the sea","drains inland","never reaches a river (desert)"][bs.closedB[f]]}`;if(bs.climate)t+=`; rain ${Math.round(bs.climate.P[f]*1000)} mm/yr, runoff ${Math.round(bs.climate.R[f]*1000)} mm/yr`;if(bs.lake[f])t+=["","; a lake that overflows","; a lake with no outlet","; a dry lake floor","; where a river runs dry"][bs.lake[f]];if(bs.inlandSea[f])t+=" (a sea cut off from the ocean)";}
    if(l1View==="tiles"&&f>=0&&bs&&bs.land[f]&&l1St.ho){const k=l1St.ho.tileOf[f];if(k>=0&&k!==l1TileSel){l1TileSel=k;l1TileInfo();l1Render();}}
    box.textContent=t;const [xk,yk]=l1XY(px,py);l1LastInspect={xKm:+xk.toFixed(1),yKm:+yk.toFixed(1),text:t};}
  box.style.display="block";clearTimeout(box._t);box._t=setTimeout(()=>box.style.display="none",5000);}
// the L1 -> L2 handoff: the selected tile's summary, and downloads of the landmass summary and of one tile's fields
let l1TileSel=-1;
function l1TileInfo(){const box=$("l1tile"),dl=$("l1dl"),ho=l1St&&l1St.ho;if(!ho){box.hidden=true;dl.hidden=true;return;}
  dl.hidden=false;$("l1DlTile").disabled=l1TileSel<0;box.hidden=false;
  const nF=ho.tiles.length,firm=ho.tiles.reduce((a,T)=>a+T.divides.reduce((s,d)=>s+d.firmKm,0),0),all=ho.tiles.reduce((a,T)=>a+T.divides.reduce((s,d)=>s+d.lengthKm,0),0);
  let html=`<div><b>L2 handoff</b>: ${nF} tiles, ${ho.imprints.length} imprints (oldest first: ${ho.imprints.slice(0,3).map(I=>(I.kind?I.kind+" ":"")+I.class+" from "+I.startMa+" Ma").join(", ")}…); ${all?Math.round(100*firm/all):0}% of the divides between tiles are firm. (${Math.round(l1St.hoMs)} ms.)</div>`;
  if(l1TileSel>=0){const T=ho.tiles[l1TileSel],o=T.outlets,sea=o.filter(x=>x.kind==="sea"),inl=o.length-sea.length,fd=T.divides.reduce((s,d)=>s+d.firmKm,0),ad=T.divides.reduce((s,d)=>s+d.lengthKm,0),fmt=v=>v.toLocaleString("en-US");
    html+=`<div><b>Tile ${T.id}</b>, ${T.kind}: ${fmt(T.areaKm2)} km², ${T.meanElevationM} m mean elevation (highest ${fmt(T.maxElevationM)} m)${T.meanRainMmYr!=null?`, ${T.meanRainMmYr} mm/yr of rain`:""}. ${sea.length?`${sea.length} outlet${sea.length>1?"s":""} to the sea (largest ${fmt(Math.round(sea[0].dischargeM3s))} m³/s)`:"No outlet to the sea"}${inl?`, ${inl} ending inland`:""}${T.lakes.length?`; ${T.lakes.length} lake${T.lakes.length>1?"s":""}`:""}. Divides with ${T.divides.length} tile${T.divides.length===1?"":"s"}, ${ad?Math.round(100*fd/ad):0}% firm. A change at the outlet takes a median ${T.responseMyr.median} Myr to reach a cell (at most ${T.responseMyr.max}). ${T.maxShortening>=0.05?`Shortened up to ${Math.round(100*T.maxShortening)}% in the window. `:""}${T.imprints.length} imprints: ${T.imprints.map(k=>{const I=ho.imprints[k];return(I.kind?I.kind+" ":"")+I.class+" "+I.startMa+"–"+I.endMa+" Ma";}).join("; ")}.</div>`;}
  box.innerHTML=html;}
async function l1Save(kind){const d=await downloadsP,note=$("l1DlNote"),ho=l1St&&l1St.ho;if(!ho)return;if(!d){note.textContent="File downloads aren’t available in this view.";return;}
  const base=`l1-${($("seed").value||"genesis").replace(/[^a-z0-9-]/gi,"_")}-${Math.round(ho.meta.presentMyr)}myr`,data=JSON.stringify(kind==="tile"?ho.tile(l1TileSel):ho.summary());
  try{await d.save({filename:kind==="tile"?`${base}-tile${l1TileSel}.json`:`${base}-landmass.json`,data});note.textContent="Saved.";}
  catch(err){note.textContent=err&&err.code==="declined"?"Download cancelled.":"The download didn’t go through.";}}
$("l1DlTile").onclick=()=>{if(l1TileSel>=0)l1Save("tile");};$("l1DlAll").onclick=()=>l1Save("all");
// ---------- session export ----------
// A session is reproducible from a short log: L0 is deterministic given its settings and step count, and L1 given the
// selected cell, the saved state it started from (L0 at a step count), its window and its settings. The export records
// those, check values to confirm a rerun matches, the viewports (what is on screen, how it is zoomed and coloured), the
// statistics shown or removed from the page, and the actions taken. tools/session.js reruns it.
let actions=[],l0LastInspect=null,l1LastInspect=null;
function logAction(o){if(!world)return;actions.push(Object.assign({a:o.a,steps:world.steps,t:+world.t.toFixed(4)},o));if(actions.length>500)actions.splice(1,1);}
function l0Hash(w){let a=0;for(let i=0;i<w.N;i++)a=(a*31+Math.round(w.elev[i]*1e4)+w.owner[i])|0;return a;}
const PROJ_DATA=p=>p?{center:p.c,east:p.ea,north:p.nr}:null;
function sessionData(){const w=world,S=l1St,p=w.cube.pos,ll=i=>[+(Math.asin(p[3*i+2])*180/Math.PI).toFixed(3),+(Math.atan2(p[3*i+1],p[3*i])*180/Math.PI).toFixed(3)];
  const l1Extent=l1Geo&&S?{xKm:[l1Geo.cx-l1PW/2*l1Geo.s,l1Geo.cx+l1PW/2*l1Geo.s],yKm:[l1Geo.cy-l1PH/2*l1Geo.s,l1Geo.cy+l1PH/2*l1Geo.s],kmPerPx:l1Geo.s}:null;
  const sh=S&&S.sh,bs=S&&S.bs,ho=S&&S.ho;
  return{format:"edgy-erosion session",version:2,build:typeof BUILD_ID!=="undefined"?BUILD_ID:null,savedAt:new Date().toISOString(),
    settings:{seed:w.opts.seed,resolution:w.n,initialPlates:w.opts.nPlates,start:w.opts.mode,continentalTarget:w.contTarget,netRotationRemoved:true,
      stepsPerFrame:+$("speed").value,runGyr:+$("gyr").value,l1:{bufferKm:l1Buf(),windowMyr:+$("l1Win").value,sheetKm:+$("l1H").value,argand:+$("l1Ar").value,couplingKm:+$("l1Ell").value}},
    reproduce:{L0:{steps:w.steps,t:w.t},selectedCell:selCell,faultHandoffAtSteps:handoff?w.steps:null,
      L1:S&&S.entry?{phase:S.phase,cell:S.cell,bufferKm:S.pick.bufferKm,presentSteps:S.plan?S.plan.endSteps:S.rp?S.rp.endSteps:null,entrySteps:S.entry.steps,entryT:S.entry.t,windowMyr:S.win,sheet:S.o,stale:!!S.stale}:null},
    checks:{L0:{hash:l0Hash(w),platesAlive:w.plates.filter(q=>q.alive).length,landFraction:w.landFrac},
      L1:sh?{sheetSteps:sh.steps,markers:sh.nm,t:sh.t,budgetError:sh.balance().error,last:sh.history[sh.history.length-1],
        basins:bs?{nBasins:bs.nBasins,landCells:bs.landCells,inlandFrac:bs.inlandFrac,nClosed:bs.nClosed,nDry:bs.nDry}:null,tiles:ho?ho.tiles.length:null}:null},
    viewports:{
      L0:{colorMode:view,projection:proj,marks:$("marks").checked,faultTraces:!!handoff&&$("traces").checked,canvasPx:[pixW,pixH],
        globe:{yaw,pitch,zoom:gZoom},map:{zoom:mZoom,centre:mC,goode:goode&&{centreLonDeg:goode.lc*180/Math.PI,edgeLonDeg:goode.edgeLonDeg,
          northLobesDeg:goode.north.map(l=>[l.a,l.b].map(v=>+(v*180/Math.PI).toFixed(2))),southLobesDeg:goode.south.map(l=>[l.a,l.b].map(v=>+(v*180/Math.PI).toFixed(2)))}},
        selected:selCell>=0?{cell:selCell,latLon:ll(selCell),elevationKm:w.elev[selCell]}:null,lastInspected:l0LastInspect},
      L1:S?{colorMode:l1View,rivers:$("l1Rivers").checked,zoom:l1Zoom,pan:l1Pan,canvasPx:[l1PW,l1PH],visible:l1Extent,
        frame:PROJ_DATA(sh?S.plan.proj:S.pv&&S.pv.proj),frameNote:sh?"the sheet's map: azimuthal equidistant co-moving with the reference plate (km)":"preview: azimuthal equidistant about the landmass's centroid in L0's fixed frame (km)",
        lastInspected:l1LastInspect}:null,
      L2:ho?{tilesShown:l1View==="tiles",selectedTile:l1TileSel,tile:l1TileSel>=0?ho.tiles[l1TileSel]:null}:null},
    stats:{L0:{t:w.t,steps:w.steps,dt:w.dt,platesAlive:w.plates.filter(q=>q.alive).length,landFraction:w.landFrac,continentalFraction:w.contFrac,meanSeafloorAgeMyr:w.oceanAge,
        meanUnitsPerContinentalColumn:w.meanUnits,netRotationDegPerMyr:w.netRot,msPerStep:msEma,simulatedMyrPerSecond:rateEma,savedStates:snaps.map(x=>({t:x.t,steps:x.steps}))},
      L1:sh?{history:sh.history,budget:sh.budget,msPerStep:l1Ms,basinsMs:S.bsMs,handoffMs:S.hoMs}:null},
    actions,events:w.events,
    landmass:S&&S.pick?{cells:Array.from(S.pick.cells),land:S.pick.land,ref:S.pick.ref}:null,
    handoffSummary:ho?ho.summary():null};}
$("exportBtn").onclick=async()=>{const note=$("exportNote");if(l1Busy()||!world)return;
  const d=await downloadsP;if(!d){note.textContent="File downloads aren’t available in this view.";return;}
  const data=JSON.stringify(sessionData()),name=`session-${String(world.opts.seed).replace(/[^a-z0-9-]/gi,"_")}-${Math.round(world.t)}myr.json`;
  try{await d.save({filename:name,data});note.textContent=`Saved ${name}.`;}
  catch(err){note.textContent=err&&err.code==="declined"?"Download cancelled.":"The export didn’t go through.";}};
$("infoBtn").onclick=()=>{const el=$("info"),open=el.hidden;el.hidden=!open;$("infoBtn").setAttribute("aria-expanded",open?"true":"false");};
build();updateLegend();requestAnimationFrame(frame);
</script>
</body>
</html>
