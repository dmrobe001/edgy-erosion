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
  const fkRGB=[[0,0,0],conv,divc,tran];
  for(let i=0;i<N;i++){
    let e=E[i];const land=e>0;let g=0;for(let k=0;k<4;k++)g+=(E[nb[4*i+k]]-e)*lw[4*i+k];
    let s=land?1-0.35*g:1-0.1*g;s=s<0.6?0.6:s>1.25?1.25:s;
    let r,gg,b;
    if(view==="elev"){if(!K[i]&&tr[i]>0.05)e-=2.2*Math.min(1,tr[i]);let k=Math.round((e+8)/16*1023);k=k<0?0:k>1023?1023:k;r=LUT_E[k*3];gg=LUT_E[k*3+1];b=LUT_E[k*3+2];}
    else if(view==="rock"){const c=ROCK_RGB[R[i]]||ROCK_RGB[1];const age=t-RA[i];const f=1-0.32*Math.min(1,age/1500);r=c[0]*f;gg=c[1]*f;b=c[2]*f;s=1+(s-1)*0.6;}
    else if(view==="faults"){const base=land?120+Math.min(60,e*12):38+Math.max(-20,e*2);r=gg=b=base;s=1+(s-1)*0.5;
      const fi=w.faultCellInfo(i);if(fi){const a=Math.min(1,Math.log10(1+fi.slip)/3.2)*(fi.active?1:0.55),c=fkRGB[fi.kind];r=r*(1-a)+c[0]*a;gg=gg*(1-a)+c[1]*a;b=b*(1-a)+c[2]*a;}}
    else if(view==="age"){if(K[i]){const v=land?150:110;r=gg=b=v;}else{let k=Math.round((t-B[i])/200*255);if(k>255)k=255;if(k<0)k=0;r=LUT_A[k*3];gg=LUT_A[k*3+1];b=LUT_A[k*3+2];s=1+(s-1)*0.4;}}
    else if(view==="thick"){let k=Math.round(T[i]/80*255);k=k<0?0:k>255?255:k;r=LUT_T[k*3];gg=LUT_T[k*3+1];b=LUT_T[k*3+2];s=1+(s-1)*0.5;}
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

function sizeCanvas(){
  const dpr=window.devicePixelRatio||1,cw=wrap.clientWidth;
  let W,H;if(proj==="globe"){W=H=Math.min(520,Math.round(cw*dpr));}else{W=Math.min(1024,Math.round(cw*dpr));W-=W%2;H=W/2;}
  if(W!==pixW||H!==pixH||!pixIdx){pixW=W;pixH=H;canvas.width=W;canvas.height=H;imgData=canvas.getContext("2d").createImageData(W,H);pixIdx=new Int32Array(W*H);pixLimb=new Uint8Array(W*H);camDirty=true;}
}
function buildPixIdx(){
  const lk=world.cube.lookup,W=pixW,H=pixH;
  if(proj==="globe"){
    Mw=camMatrix();const [a,b,c]=Mw;
    for(let py=0;py<H;py++){const sy=1-(py+0.5)/H*2;for(let px=0;px<W;px++){const sx=(px+0.5)/W*2-1,r2=sx*sx+sy*sy,q=py*W+px;
      if(r2>1){pixIdx[q]=-1;continue;}const sz=Math.sqrt(1-r2);
      pixIdx[q]=lk(a[0]*sx+b[0]*sy+c[0]*sz,a[1]*sx+b[1]*sy+c[1]*sz,a[2]*sx+b[2]*sy+c[2]*sz);pixLimb[q]=Math.round(255*(0.5+0.5*Math.sqrt(sz)));}}
  } else {
    for(let py=0;py<H;py++){const lat=Math.PI/2-(py+0.5)/H*Math.PI,cl=Math.cos(lat),sl=Math.sin(lat);
      for(let px=0;px<W;px++){const lon=(px+0.5)/W*2*Math.PI-Math.PI;pixIdx[py*W+px]=lk(cl*Math.cos(lon),cl*Math.sin(lon),sl);pixLimb[py*W+px]=255;}}
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
    for(let py=0;py<H;py++){const sy=1-(py+0.5)/H*2;for(let px=0;px<W;px++){const sx=(px+0.5)/W*2-1,r2=sx*sx+sy*sy,q=py*W+px;if(r2>1)continue;const sz=Math.sqrt(1-r2);
      put(q,a[0]*sx+b[0]*sy+c[0]*sz,a[1]*sx+b[1]*sy+c[1]*sz,a[2]*sx+b[2]*sy+c[2]*sz);}}
  } else {
    for(let py=0;py<H;py++){const lat=Math.PI/2-(py+0.5)/H*Math.PI,cl=Math.cos(lat),sl=Math.sin(lat);
      for(let px=0;px<W;px++){const lon=(px+0.5)/W*2*Math.PI-Math.PI;put(py*W+px,cl*Math.cos(lon),cl*Math.sin(lon),sl);}}
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
    if(z<0.05)return null;return[(x+1)/2*pixW,(1-y)/2*pixH];}
  const lat=Math.asin(Math.max(-1,Math.min(1,v[2]))),lon=Math.atan2(v[1],v[0]);return[(lon+Math.PI)/(2*Math.PI)*pixW,(Math.PI/2-lat)/Math.PI*pixH,lon];
}
function render(){
  sizeCanvas();if(camDirty)buildPixIdx();
  const d=imgData.data,n=pixW*pixH;
  if(view==="l1"){ensureL1();if(geoDirty)buildGeo();paintL1();}
  else{if(colorsDirty)computeColors();
    for(let q=0;q<n;q++){const i=pixIdx[q],o=4*q;if(i<0){d[o+3]=0;continue;}const l=pixLimb[q];
      d[o]=cellRGB[3*i]*l>>8;d[o+1]=cellRGB[3*i+1]*l>>8;d[o+2]=cellRGB[3*i+2]*l>>8;d[o+3]=255;}}
  const ctx=canvas.getContext("2d");ctx.clearRect(0,0,pixW,pixH);ctx.putImageData(imgData,0,0);
  const dpr=pixW/wrap.clientWidth;
  if($("hots").checked){for(const h of world.hotspots){const p=project(h.v);if(!p)continue;const rad=Math.max(4,5*dpr);
    ctx.lineWidth=Math.max(1.5,1.6*dpr);ctx.strokeStyle="rgba(0,0,0,.55)";ctx.beginPath();ctx.arc(p[0],p[1],rad+1,0,7);ctx.stroke();
    ctx.strokeStyle="#ffcf4a";ctx.beginPath();ctx.arc(p[0],p[1],rad,0,7);ctx.stroke();}}
  if($("arrows").checked){
    for(const P of world.plates){if(!P.alive||P.area<world.N*0.006)continue;const c=P.c;
      const v=[P.W[1]*c[2]-P.W[2]*c[1],P.W[2]*c[0]-P.W[0]*c[2],P.W[0]*c[1]-P.W[1]*c[0]];const k=28;
      const e=[c[0]+v[0]*k,c[1]+v[1]*k,c[2]+v[2]*k];const el=Math.hypot(...e);e[0]/=el;e[1]/=el;e[2]/=el;
      const a=project(c),b=project(e);if(!a||!b)continue;if(proj==="flat"&&Math.abs(a[2]-b[2])>Math.PI)continue;drawArrow(ctx,a[0],a[1],b[0],b[1],dpr);}
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
        if(pen&&proj==="flat"&&Math.abs(p[2]-pts[k-1][2])>Math.PI)pen=false;
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
const EMPTY_LOG='<li class="empty" style="border:0">Events appear here as the world evolves: rifts, collisions, sutures, consumed plates, hotspot islands.</li>';
function build(){
  world=new SIM.World({n:+$("grid").value,seed:$("seed").value||"genesis",nPlates:+$("nPlates").value,mode:$("mode").value,nnr:$("nnr").checked,contTarget:+$("contT").value/100});
  lightW=prepLight(world.cube);colorsDirty=true;handoff=null;l1=null;geoDirty=true;selCell=-1;showColumn();renderSummary();camDirty=true;loggedCount=0;$("log").innerHTML=EMPTY_LOG;msEma=0;rateEma=0;render();updateUi();
}
function fmtT(t){return t>=1000?(t/1000).toFixed(2):Math.round(t).toString();}
function updateUi(){
  const t=world.t;$("clock").innerHTML=t>=1000?fmtT(t)+"<small>billion years</small>":fmtT(t)+"<small>million years</small>";
  $("sPlates").textContent=world.plates.filter(p=>p.alive).length;
  $("sLand").textContent=(world.landFrac*100).toFixed(1)+"%";
  $("sCont").textContent=(world.contFrac*100).toFixed(1)+"%";
  $("sAge").textContent=Math.round(world.oceanAge)+" Myr";
  $("sUnits").textContent=world.meanUnits.toFixed(1);
  $("sNet").textContent=(world.netRot||0).toFixed(2)+"°/Myr";
  $("sMs").textContent=msEma?msEma.toFixed(1)+" ms":"–";
  $("sRate").textContent=rateEma?Math.round(rateEma)+" Myr/s":"–";
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
  if(view==="elev")L.innerHTML=`<span class="num">−8 km</span><span class="bar" style="background:${gradCss(ELEV_STOPS,-8,8)}"></span><span class="num">+8 km</span><span>Sea level is the blue–green break; dark lines are trenches.</span>`;
  else if(view==="rock")L.innerHTML=ROCKS.slice(1).map(r=>`<span><i class="swatch" style="background:${r.c}"></i>${r.name}</span>`).join("")+`<span>Darker means older.</span>`;
  else if(view==="l1")L.innerHTML=ROCKS.slice(1).map(r=>`<span><i class="swatch" style="background:${r.c}"></i>${r.name}</span>`).join("")+`<span>Smoothed L0 elevation and formations, cut by stand-in fine relief (noise, not a real L1 model). Updates every few seconds while running.</span>`;
  else if(view==="faults")L.innerHTML=`<span><i class="swatch" style="background:var(--conv)"></i>Thrust</span><span><i class="swatch" style="background:var(--div)"></i>Normal</span><span><i class="swatch" style="background:var(--trans)"></i>Strike-slip</span><span>Stronger colour means more total slip; faded means inactive.</span>`;
  else if(view==="age")L.innerHTML=`<span class="num">0</span><span class="bar" style="background:${gradCss(AGE_STOPS,0,200)}"></span><span class="num">200 Myr</span><span>Continents shown grey.</span>`;
  else if(view==="thick")L.innerHTML=`<span class="num">0</span><span class="bar" style="background:${gradCss(THICK_STOPS,0,80)}"></span><span class="num">80 km</span>`;
  else L.innerHTML=`<span><i class="swatch" style="background:var(--conv)"></i>Convergent</span><span><i class="swatch" style="background:var(--div)"></i>Divergent</span><span><i class="swatch" style="background:var(--trans)"></i>Transform</span><span>Light fill is continental crust.</span>`;
}
function frame(now){
  let need=false;
  if(running){const n=+$("speed").value,t0=performance.now();for(let k=0;k<n;k++)world.step();
    const ms=(performance.now()-t0)/n;msEma=msEma?msEma*0.9+ms*0.1:ms;
    const fdt=(now-lastFrame)/1000;if(fdt>0&&fdt<1){const rate=n*world.dt/fdt;rateEma=rateEma?rateEma*0.9+rate*0.1:rate;}
    colorsDirty=true;need=true;if(handoff)staleHandoff();updateUi();if(selCell>=0)showColumn();}
  if(proj==="globe"&&$("spin").checked&&!dragging){yaw+=0.0035;camDirty=true;need=true;}
  if(need)render();lastFrame=now;requestAnimationFrame(frame);
}
$("play").onclick=()=>{running=!running;$("play").textContent=running?"Pause":"Run";};
$("stepBtn").onclick=()=>{const t0=performance.now();world.step();msEma=performance.now()-t0;colorsDirty=true;if(handoff)staleHandoff();render();updateUi();if(selCell>=0)showColumn();};
$("speed").oninput=e=>$("speedVal").textContent=e.target.value;
$("nPlates").oninput=e=>$("nPlatesVal").textContent=e.target.value;
$("contT").oninput=e=>$("contTVal").textContent=e.target.value+"%";
$("contT").onchange=build;
$("reset").onclick=build;
$("randSeed").onclick=()=>{const w=["basalt","craton","ophiolite","gneiss","tethys","rodinia","zircon","karst","obsidian","pelagic","moho","flysch","horst","graben","olivine","andesite"];
  $("seed").value=w[Math.floor(Math.random()*w.length)]+"-"+Math.floor(Math.random()*1000);build();};
$("seed").addEventListener("keydown",e=>{if(e.key==="Enter")build();});
["mode","grid"].forEach(id=>$(id).onchange=build);
["arrows","hots","traces"].forEach(id=>$(id).onchange=()=>render());
$("nnr").onchange=()=>{world.opts.nnr=$("nnr").checked;if(world.opts.nnr){world.removeNetRotation();world.updateStats(true);updateUi();render();}};
document.querySelectorAll("#views button").forEach(b=>b.onclick=()=>{view=b.dataset.v;
  document.querySelectorAll("#views button").forEach(x=>x.setAttribute("aria-pressed",x===b?"true":"false"));colorsDirty=true;updateLegend();render();});
document.querySelectorAll("#proj button").forEach(b=>b.onclick=()=>{proj=b.dataset.p;
  document.querySelectorAll("#proj button").forEach(x=>x.setAttribute("aria-pressed",x===b?"true":"false"));
  wrap.className="mapwrap "+proj;$("hint").textContent=proj==="globe"?"Drag to turn the globe. Tap a spot to inspect it.":"Tap a spot to inspect it.";
  pixIdx=null;render();});
let dragging=false,downX=0,downY=0,lastX=0,lastY=0,moved=false;
wrap.addEventListener("pointerdown",e=>{downX=lastX=e.clientX;downY=lastY=e.clientY;moved=false;if(proj==="globe"){dragging=true;wrap.setPointerCapture(e.pointerId);}});
wrap.addEventListener("pointermove",e=>{if(!dragging)return;const dx=e.clientX-lastX,dy=e.clientY-lastY;lastX=e.clientX;lastY=e.clientY;
  if(Math.abs(e.clientX-downX)+Math.abs(e.clientY-downY)>5)moved=true;
  const k=3/wrap.clientWidth;yaw-=dx*k;pitch=Math.max(-1.45,Math.min(1.45,pitch+dy*k));camDirty=true;render();});
wrap.addEventListener("pointerup",e=>{dragging=false;if(!moved)inspect(e);});
wrap.addEventListener("pointercancel",()=>{dragging=false;});
function inspect(e){
  const r=canvas.getBoundingClientRect(),px=Math.floor((e.clientX-r.left)/r.width*pixW),py=Math.floor((e.clientY-r.top)/r.height*pixH);
  if(px<0||py<0||px>=pixW||py>=pixH)return;const i=pixIdx[py*pixW+px];if(i<0)return;
  selCell=i;showColumn();
  const w=world,e2=w.elev[i],box=$("inspect");box.style.display="block";
  box.textContent=`${cellLoc(i)}: ${e2>=0?"elevation":"depth"} ${Math.abs(e2).toFixed(1)} km. Column shown below the map.`;
  clearTimeout(box._t);box._t=setTimeout(()=>box.style.display="none",3000);
}
function cellLoc(i){const p=world.cube.pos,lat=Math.asin(p[3*i+2])*180/Math.PI,lon=Math.atan2(p[3*i+1],p[3*i])*180/Math.PI;
  return `${Math.abs(lat).toFixed(0)}°${lat>=0?"N":"S"} ${Math.abs(lon).toFixed(0)}°${lon>=0?"E":"W"}`;}
function fmtAgo(a){a=Math.max(0,a);return a>=1000?(a/1000).toFixed(2)+" Gyr":Math.round(a)+" Myr";}
function fmtKm(v){const a=Math.abs(v);return a>=100?Math.round(a).toLocaleString("en-US"):a.toFixed(a>=10?0:1);}
const COMPASS=["north","northeast","east","southeast","south","southwest","west","northwest"];
function showColumn(){
  const body=$("colBody"),empty=$("colEmpty");
  if(selCell<0||!world){body.hidden=true;empty.hidden=false;return;}
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
  body.innerHTML=html;body.hidden=false;empty.hidden=true;
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
  handoff=world.extractFaults();handoff.l1=world.buildL1();l1=handoff.l1;l1Time=performance.now();renderSummary();render();};
$("download").onclick=async()=>{const d=await downloadsP;if(!d||!handoff)return;
  const data=JSON.stringify(world.exportHandoff(handoff));
  const name=`l0-handoff-${($("seed").value||"genesis").replace(/[^a-z0-9-]/gi,"_")}-${Math.round(handoff.t)}myr.json`;
  try{await d.save({filename:name,data});$("dlNote").textContent="Saved.";}
  catch(err){$("dlNote").textContent=err&&err.code==="declined"?"Download cancelled.":"The download didn’t go through.";}};
downloadsP=(window.claude&&typeof window.claude.use==="function")?window.claude.use("downloads").catch(()=>null):Promise.resolve(null);
window.addEventListener("resize",()=>{pixIdx=null;render();});
build();updateLegend();requestAnimationFrame(frame);
</script>
</body>
</html>
