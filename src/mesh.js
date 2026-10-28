// ===== MESH-START =====
// The marker core for L0 (L0_MESH_SPEC.md, Stage 1): markers carry the crust; a Delaunay triangulation of their
// positions on the sphere is kept up to date as they move; plates are the connected sets of markers between interface
// markers. Stage 1 keeps the lattice core's kinematics and surface rules, so the two can be compared on one seed.
function MESHmodule(SIM){

// ---------- Delaunay triangulation on the unit sphere ----------
// Vertices are indices into X, Y, Z (owned here, grown together with the caller's marker arrays). Triangles are
// counter-clockwise seen from outside. Edge k of a triangle is the one opposite its vertex k; tn holds the triangle
// across it. On the sphere the Delaunay triangulation is the convex hull of the points, so
//   orient(a,b,c) = det[a,b,c] > 0  for a counter-clockwise triangle, and
//   d lies inside the circumcircle of (a,b,c)  iff  det[b-a, c-a, d-a] > 0  (d is beyond the plane of a, b, c).
// It is kept Delaunay under insertion, deletion and small motions by Lawson flips.
const IC_EPS=1e-18,ON_EPS=1e-13;
class SphereTri{
  constructor(vcap,tcap){this.X=new Float64Array(vcap);this.Y=new Float64Array(vcap);this.Z=new Float64Array(vcap);this.vt=new Int32Array(vcap).fill(-1);
    this.tv=new Int32Array(3*tcap);this.tn=new Int32Array(3*tcap);this.tAlive=new Uint8Array(tcap);this.freeT=[];this.nT=0;this.tcap=tcap;this.vcap=vcap;
    this.stamp=new Int32Array(tcap);this.stampN=0;this.walkSeed=1;this.stack=[];
    // a cache per triangle, filled by its user: [area, and three numbers per edge]; cleared whenever the triangle is rewritten
    this.tCache=new Float64Array(4*tcap);this.tOk=new Uint8Array(tcap);}
  growV(n){if(n<=this.vcap)return;const c=Math.max(n,Math.ceil(this.vcap*1.5));
    for(const k of ["X","Y","Z"]){const a=new Float64Array(c);a.set(this[k]);this[k]=a;}const vt=new Int32Array(c).fill(-1);vt.set(this.vt);this.vt=vt;this.vcap=c;}
  growT(){const c=Math.ceil(this.tcap*1.5)+16,tv=new Int32Array(3*c),tn=new Int32Array(3*c),al=new Uint8Array(c),st=new Int32Array(c);
    tv.set(this.tv);tn.set(this.tn);al.set(this.tAlive);st.set(this.stamp);const tc=new Float64Array(4*c),ok=new Uint8Array(c);tc.set(this.tCache);ok.set(this.tOk);
    Object.assign(this,{tv,tn,tAlive:al,stamp:st,tCache:tc,tOk:ok,tcap:c});}
  newT(){if(this.freeT.length)return this.freeT.pop();if(this.nT>=this.tcap)this.growT();return this.nT++;}
  freeTri(t){this.tAlive[t]=0;this.freeT.push(t);}
  setT(t,a,b,c,na,nb,nc){const o=3*t,tv=this.tv,tn=this.tn;tv[o]=a;tv[o+1]=b;tv[o+2]=c;tn[o]=na;tn[o+1]=nb;tn[o+2]=nc;this.tAlive[t]=1;this.tOk[t]=0;
    this.vt[a]=t;this.vt[b]=t;this.vt[c]=t;}
  orient(a,b,c){const X=this.X,Y=this.Y,Z=this.Z;
    return X[a]*(Y[b]*Z[c]-Z[b]*Y[c])-Y[a]*(X[b]*Z[c]-Z[b]*X[c])+Z[a]*(X[b]*Y[c]-Y[b]*X[c]);}
  orientP(a,b,px,py,pz){const X=this.X,Y=this.Y,Z=this.Z;
    return X[a]*(Y[b]*pz-Z[b]*py)-Y[a]*(X[b]*pz-Z[b]*px)+Z[a]*(X[b]*py-Y[b]*px);}
  incircle(a,b,c,d){const X=this.X,Y=this.Y,Z=this.Z,ax=X[a],ay=Y[a],az=Z[a];
    const bx=X[b]-ax,by=Y[b]-ay,bz=Z[b]-az,cx=X[c]-ax,cy=Y[c]-ay,cz=Z[c]-az,dx=X[d]-ax,dy=Y[d]-ay,dz=Z[d]-az;
    return bx*(cy*dz-cz*dy)-by*(cx*dz-cz*dx)+bz*(cx*dy-cy*dx);}
  // index of the edge of triangle u that faces triangle t
  edgeTo(u,t){const o=3*u,tn=this.tn;return tn[o]===t?0:tn[o+1]===t?1:tn[o+2]===t?2:-1;}
  relink(u,from,to){if(u<0)return;const o=3*u,tn=this.tn;if(tn[o]===from)tn[o]=to;else if(tn[o+1]===from)tn[o+1]=to;else if(tn[o+2]===from)tn[o+2]=to;}
  // the triangle containing point p, by a visibility walk from hint (random edge order, so it cannot cycle); a point
  // within ON_EPS of an edge counts as inside both triangles
  locate(px,py,pz,hint){const tv=this.tv,tn=this.tn;let t=hint>=0&&this.tAlive[hint]?hint:this.anyT();
    for(let it=0;it<1e6;it++){let moved=false;let x=this.walkSeed;x^=x<<13;x^=x>>>17;x^=x<<5;this.walkSeed=x|0||1;const r=((x>>>0)%3);
      for(let q=0;q<3;q++){const k=(q+r)%3,o=3*t,a=tv[o+(k+1)%3],b=tv[o+(k+2)%3];
        if(this.orientP(a,b,px,py,pz)<-ON_EPS){t=tn[o+k];moved=true;break;}}
      if(!moved)return t;}
    return -1;}
  anyT(){for(let t=0;t<this.nT;t++)if(this.tAlive[t])return t;return -1;}
  // flip the edge k of t (opposite its vertex a); returns false if the quadrilateral is not convex
  flip(t,k){const tv=this.tv,tn=this.tn,o=3*t,u=tn[o+k];if(u<0)return false;const j=this.edgeTo(u,t),p=3*u;
    const a=tv[o+k],b=tv[o+(k+1)%3],c=tv[o+(k+2)%3],d=tv[p+j];
    if(this.orient(a,b,d)<=0||this.orient(a,d,c)<=0)return false;
    const nTb=tn[o+(k+1)%3],nTc=tn[o+(k+2)%3],nUc=tn[p+(j+1)%3],nUb=tn[p+(j+2)%3];
    // u was (d,c,b): opposite c is edge (b,d) -> nUc; opposite b is edge (d,c) -> nUb
    this.setT(t,a,b,d,nUc,u,nTc);this.setT(u,a,d,c,nUb,nTb,t);
    this.relink(nUc,u,t);this.relink(nTb,t,u);return true;}
  // Lawson flips until every edge of the queued triangles (and of the triangles flips create) is locally Delaunay
  legalize(queue){const tv=this.tv,tn=this.tn,st=this.stack;st.length=0;this.stampN++;
    for(const t of queue)if(this.tAlive[t]){st.push(t);}
    let flips=0;
    while(st.length){const t=st.pop();if(!this.tAlive[t])continue;
      for(let k=0;k<3;k++){const o=3*t,u=tn[o+k];if(u<0)continue;const j=this.edgeTo(u,t),d=tv[3*u+j];
        if(this.incircle(tv[o],tv[o+1],tv[o+2],d)>IC_EPS){if(this.flip(t,k)){flips++;st.push(t,u);if(this.touched)this.touched.push(t,u);break;}}}}
    return flips;}
  // insert vertex v (coordinates already set); returns false if it coincides with an existing vertex. A point on an
  // edge splits the edge (two triangles become four), so no zero-area triangle is made.
  insert(v,hint){const X=this.X,Y=this.Y,Z=this.Z,px=X[v],py=Y[v],pz=Z[v],t=this.locate(px,py,pz,hint);if(t<0)return false;
    const tv=this.tv,tn=this.tn,o=3*t,a=tv[o],b=tv[o+1],c=tv[o+2];
    for(const q of [a,b,c]){const dx=X[q]-px,dy=Y[q]-py,dz=Z[q]-pz;if(dx*dx+dy*dy+dz*dz<1e-16)return false;}
    for(let k=0;k<3;k++){const e1=tv[o+(k+1)%3],e2=tv[o+(k+2)%3];if(Math.abs(this.orientP(e1,e2,px,py,pz))<ON_EPS){this.splitEdge(t,k,v);return true;}}
    const na=tn[o],nb=tn[o+1],nc=tn[o+2],t1=this.newT(),t2=this.newT();
    this.setT(t,v,b,c,na,t1,t2);this.setT(t1,v,c,a,nb,t2,t);this.setT(t2,v,a,b,nc,t,t1);
    this.relink(nb,t,t1);this.relink(nc,t,t2);this.vt[a]=t1;this.vt[b]=t;this.vt[c]=t;this.vt[v]=t;
    this.legalize([t,t1,t2]);return true;}
  // p on edge k of t (the edge (b,c) opposite a); u across it holds (d,c,b)
  splitEdge(t,k,p){const tv=this.tv,tn=this.tn,o=3*t,u=tn[o+k],j=this.edgeTo(u,t),q=3*u;
    const a=tv[o+k],b=tv[o+(k+1)%3],c=tv[o+(k+2)%3],d=tv[q+j];
    const nTb=tn[o+(k+1)%3],nTc=tn[o+(k+2)%3],nUc=tn[q+(j+1)%3],nUb=tn[q+(j+2)%3];
    const T1=t,T3=u,T2=this.newT(),T4=this.newT();
    this.setT(T1,p,a,b,nTc,T3,T2);this.setT(T2,p,c,a,nTb,T1,T4);this.setT(T3,p,b,d,nUc,T4,T1);this.setT(T4,p,d,c,nUb,T2,T3);
    this.relink(nTb,t,T2);this.relink(nUb,u,T4); // nTc keeps t (= T1), nUc keeps u (= T3)
    this.vt[a]=T1;this.vt[b]=T1;this.vt[c]=T2;this.vt[d]=T3;this.vt[p]=T1;
    this.legalize([T1,T2,T3,T4]);}
  // the ring of v: triangles and opposite edges in counter-clockwise order around v
  ring(v,outTris,outVerts){outTris.length=0;outVerts.length=0;const tv=this.tv,tn=this.tn,t0=this.vt[v];if(t0<0)return;let t=t0;
    for(let it=0;it<1000;it++){const o=3*t,i=tv[o]===v?0:tv[o+1]===v?1:2;outTris.push(t);outVerts.push(tv[o+(i+1)%3]);
      t=tn[o+(i+1)%3];if(t===t0)return;} // across edge (v, c) opposite b: the next triangle counter-clockwise
    throw new Error("ring walk did not close at vertex "+v);}
  // remove vertex v and fill the hole with Delaunay ears
  remove(v){const tris=[],P=[];this.ring(v,tris,P);if(!tris.length)return;
    const tv=this.tv,tn=this.tn,k=tris.length,outT=new Array(k),outE=new Array(k);
    for(let q=0;q<k;q++){const t=tris[q],o=3*t,i=tv[o]===v?0:tv[o+1]===v?1:2,u=tn[o+i];outT[q]=u;outE[q]=u>=0?this.edgeTo(u,t):-1;}
    for(const t of tris)this.freeTri(t);this.vt[v]=-1;
    // polygon P[q] -> P[q+1] has outside triangle outT[q] (its edge outE[q]) — the ring edge is (b,c) of triangle q
    const pv=P.slice(),pT=outT.slice(),pE=outE.slice(),made=[];
    const link=(T,e,u,ue)=>{this.tn[3*T+e]=u;if(u>=0)this.tn[3*u+ue]=T;};
    while(pv.length>3){const n=pv.length;let best=-1;
      for(let q=0;q<n&&best<0;q++){const a=pv[q],b=pv[(q+1)%n],c=pv[(q+2)%n];if(this.orient(a,b,c)<=0)continue;
        let ok=true;for(let r=0;r<n;r++){const d=pv[r];if(d===a||d===b||d===c)continue;if(this.incircle(a,b,c,d)>IC_EPS){ok=false;break;}}
        if(ok)best=q;}
      if(best<0)for(let q=0;q<n;q++){const a=pv[q],b=pv[(q+1)%n],c=pv[(q+2)%n];if(this.orient(a,b,c)>0){best=q;break;}}
      if(best<0){best=0;this.broken=true;} // a folded hole: the caller rebuilds (MeshWorld checks this flag)
      const q=best,q1=(q+1)%n,a=pv[q],b=pv[q1],c=pv[(q+2)%n],T=this.newT();
      this.setT(T,a,b,c,-1,-1,-1);link(T,0,pT[q1],pE[q1]);link(T,2,pT[q],pE[q]);made.push(T);
      // the new polygon edge (a,c) faces T across T's edge 1 (opposite b)
      pT[q]=T;pE[q]=1;pv.splice(q1,1);pT.splice(q1,1);pE.splice(q1,1);}
    const T=this.newT(),a=pv[0],b=pv[1],c=pv[2];this.setT(T,a,b,c,-1,-1,-1);
    // edge opposite a is (b,c) = polygon edge 1; opposite b is (c,a) = polygon edge 2; opposite c is (a,b) = polygon edge 0
    link(T,0,pT[1],pE[1]);link(T,1,pT[2],pE[2]);link(T,2,pT[0],pE[0]);made.push(T);
    this.legalize(made);}
  // build from scratch: a tetrahedron from four well-spread vertices, then insertion in the given order
  build(verts){this.tAlive.fill(0);this.freeT=[];this.nT=0;for(const v of verts)this.vt[v]=-1;
    const X=this.X,Y=this.Y,Z=this.Z,dirs=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],pick=[];
    for(const d of dirs){let best=-1,bd=-2;for(const v of verts){const s=(X[v]*d[0]+Y[v]*d[1]+Z[v]*d[2]);if(s>bd&&!pick.includes(v)){bd=s;best=v;}}pick.push(best);}
    // the four faces, each turned counter-clockwise seen from outside (the origin lies inside the tetrahedron)
    const [a,b,c,d]=pick,F=[[a,b,c],[a,b,d],[a,c,d],[b,c,d]].map(f=>this.orient(f[0],f[1],f[2])>0?f:[f[0],f[2],f[1]]);
    const T=F.map(()=>this.newT());for(let i=0;i<4;i++){const [p,q,r]=F[i];this.setT(T[i],p,q,r,-1,-1,-1);}
    this.rebuildAdjacency(T);
    let hint=T[0];for(const v of verts){if(pick.includes(v))continue;this.insert(v,hint);if(this.vt[v]>=0)hint=this.vt[v];}}
  // adjacency of a set of triangles by matching edges (used only for the start tetrahedron)
  rebuildAdjacency(ts){const m=new Map(),tv=this.tv;for(const t of ts)for(let k=0;k<3;k++){const a=tv[3*t+(k+1)%3],b=tv[3*t+(k+2)%3];m.set(a+","+b,[t,k]);}
    for(const t of ts)for(let k=0;k<3;k++){const a=tv[3*t+(k+1)%3],b=tv[3*t+(k+2)%3],r=m.get(b+","+a);this.tn[3*t+k]=r?r[0]:-1;}}
  // the spherical area of triangle t (Eriksson's formula)
  area(t){const tv=this.tv,o=3*t,a=tv[o],b=tv[o+1],c=tv[o+2],X=this.X,Y=this.Y,Z=this.Z;
    const det=this.orient(a,b,c),ab=X[a]*X[b]+Y[a]*Y[b]+Z[a]*Z[b],bc=X[b]*X[c]+Y[b]*Y[c]+Z[b]*Z[c],ca=X[c]*X[a]+Y[c]*Y[a]+Z[c]*Z[a];
    return 2*Math.atan2(det,1+ab+bc+ca);}
  // full check: orientation, adjacency symmetry, vertex links, Euler's formula and (optionally) the Delaunay property
  check(nVerts,delaunay=true){const tv=this.tv,tn=this.tn,errs=[];let F=0,minO=Infinity,badIC=0;
    for(let t=0;t<this.nT;t++){if(!this.tAlive[t])continue;F++;const o=3*t,ori=this.orient(tv[o],tv[o+1],tv[o+2]);if(ori<minO)minO=ori;if(ori<=0)errs.push("inverted "+t);
      for(let k=0;k<3;k++){const u=tn[o+k];if(u<0||!this.tAlive[u]){errs.push("open edge "+t);continue;}const j=this.edgeTo(u,t);if(j<0)errs.push("asymmetric "+t+"/"+u);
        else{const a=tv[o+(k+1)%3],b=tv[o+(k+2)%3],c=tv[3*u+(j+1)%3],d=tv[3*u+(j+2)%3];if(a!==d||b!==c)errs.push("edge mismatch "+t+"/"+u);
          if(delaunay&&this.incircle(tv[o],tv[o+1],tv[o+2],tv[3*u+j])>1e-12)badIC++;}}}
    if(F!==2*nVerts-4)errs.push(`Euler: ${F} triangles for ${nVerts} vertices`);
    let badVt=0;for(let v=0;v<this.vcap;v++){const t=this.vt[v];if(t<0)continue;if(!this.tAlive[t]||(tv[3*t]!==v&&tv[3*t+1]!==v&&tv[3*t+2]!==v))badVt++;}
    if(badVt)errs.push(badVt+" vertices with a stale incident triangle");
    if(badIC)errs.push(badIC+" non-Delaunay edges");
    return{ok:errs.length===0,errs:errs.slice(0,10),F,minOrient:minO};}
}

// ---------- the marker world (Stage 1) ----------
const {World,C,ROCK,R_KM,FS,elevOf,kindOf,randUnit,norm3,cross}=SIM;
const INTERIOR=1,IFACE=2;                         // marker kinds (0 = free slot)
const TRANSFORM=0,RIDGE=1,CONVERGENT=2,JUNCTION=3; // interface kinds
// marker arrays and their element types; each is grown together
const MFIELDS={mId:Uint32Array,mKind:Uint8Array,mPlate:Int16Array,mA:Float64Array,mT:Float32Array,mE:Float32Array,mB:Float32Array,mS:Float32Array,
  mK:Uint8Array,mCol:Int32Array,mBorn:Float32Array,arcAdd:Float64Array,collAdd:Float64Array,mBn:Float32Array,
  iN:Uint8Array,iPlate:Int16Array,iSX:Float64Array,iSY:Float64Array,iSZ:Float64Array,iKind:Uint8Array,iWin:Int8Array,
  X0:Float64Array,Y0:Float64Array,Z0:Float64Array,dT:Float64Array,stamp:Int32Array,
  sA:Float64Array,sTh:Float64Array,sOc:Float64Array,sCc:Float64Array,sFed:Float32Array,sFc:Float32Array,mU:Float64Array,iNrm:Float64Array,iLen:Float64Array};
const MWIDTH={mBn:3,iPlate:3,iSX:3,iSY:3,iSZ:3,iNrm:3};
function hash2(i,salt){let h=Math.imul(i,0x27d4eb2d)^salt;h=Math.imul(h^(h>>>16),0x7feb352d);h=Math.imul(h^(h>>>15),0x846ca68b);h^=h>>>16;return[(h&0xffff)/65536-0.5,(h>>>16)/65536-0.5];}
function rotM(W,ang){const w=Math.hypot(W[0],W[1],W[2]);if(w<=0||ang===0)return[1,0,0,0,1,0,0,0,1];
  const ux=W[0]/w,uy=W[1]/w,uz=W[2]/w,c=Math.cos(ang),s=Math.sin(ang),t=1-c;
  return[c+ux*ux*t,ux*uy*t-uz*s,ux*uz*t+uy*s, uy*ux*t+uz*s,c+uy*uy*t,uy*uz*t-ux*s, uz*ux*t-uy*s,uz*uy*t+ux*s,c+uz*uz*t];}

// ---------- Stage 2 physics (L0_MESH_SPEC.md §4): columns, slabs, and the constants of the force balance ----------
// SI units throughout this block (m, kg, s, N); times in Myr only where noted. Each constant is labelled as in the spec:
// [measured] a literature value or mid-range, [derived] follows from others, [rule] a modelling choice, [calibrated] fitted.
const PH={g:9.8,rhoM:3300,rhoW:1030,rhoOc:2900,rhoCc:2800,    // [measured] mantle, sea water, oceanic and continental crust
  alpha:3e-5,dTm:1300,kappa:1e-6,hMax:125e3,                   // [measured] expansivity, mantle - surface temperature, diffusivity; plate-model cap on thickness
  eRidge:-2600,rootH:150e3,                                    // [rule] ridge crest depth (as the elevation rule); continental root thickness
  drhoEcl:250,lEcl:80e3,lMax:700e3,                            // [measured] eclogite excess density, depth reached before it forms, slab length cap (660 km)
  lPlug:150e3,                                                 // [rule] longest buoyant plug of continental crust at the top of a slab
  etaUM:3e20,rBend:200e3,cGeom:6,                              // [measured] upper-mantle viscosity, bending radius; [calibrated, provisionally] slab
                                                               // shape factor: gave RMS 4-5 cm/yr at n = 48 before welding and slab attachment were
                                                               // made independent of resolution. At n = 96 the mean speed is set by how few plates
                                                               // survive welding, not by this (L0_MESH_SPEC.md §11.2)
  Cb:3e14,kKeel:5,                                             // [measured range, mid] basal drag eta/h (Pa s/m); [conjecture] factor under continents
  cSuc:0.2,                                                    // [calibrated in the spec's range 0.1-0.3, not fitted] trench suction
  cRidge:3e20,cTrench:1e21,cColl:1e22,                         // [conjecture: orders of magnitude, not fitted] interface coupling at low slip rates (Pa s)
  sRidge:3e11,sTrench:2e12,sColl:6e12,                         // [measured orders] its frictional limit (N/m): about 20 MPa over a 15 km fault at ridges and
                                                               // transforms; 15 MPa over 100-150 km of megathrust; the force Tibet transmits (5-7e12)
  weldCellKm2:4*Math.PI*R_KM*R_KM/(6*96*96),                   // [rule] the weld counter's unit of collided area: a cell of the n = 96 grid
  tauDetach:10,                                                // [rule, Myr] a slab no longer fed detaches with this e-folding time (stands in for §4.6 break-off)
  MYR:3.15576e13};
// thermal boundary layer thickness at an age in Myr: half-space cooling, capped (plate model)
function hth(age){return Math.min(PH.hMax,2.32*Math.sqrt(PH.kappa*(age>0?age:0)*PH.MYR));}
const HALF_ARD=0.5*PH.alpha*PH.rhoM*PH.dTm;                    // mean excess density of the boundary layer, about 64 kg/m3
// column buoyancy per area relative to the asthenosphere (kg/m2; times g for a pressure), §4.1; continents' thermal
// term is taken as cancelled by their depleted roots
function buoy(K,Tkm,age){return K?(PH.rhoM-PH.rhoCc)*Tkm*1e3:(PH.rhoM-PH.rhoOc)*Tkm*1e3-HALF_ARD*hth(age);}
// Gravitational potential energy of a column (N/m), §4.2: the vertical integral of lithostatic pressure from the surface
// (sea level where there is water) to a depth below every root, with the column floated isostatically. Layers: water,
// crust, a lithospheric mantle whose excess density falls linearly to zero at its base, then asthenosphere. The
// asthenosphere's free level is set so a ridge sits at eRidge; the continental root's excess mass so that 33 km of
// crust sits at sea level, as the elevation rule has it [derived from that rule]. Differences between columns are forces
// per length: ridge push at 100 Myr comes out 3.9e12 N/m (Turcotte & Schubert give about 3.9e12). Returned relative to
// a ridge column.
const GPE_D=300e3,EA=(PH.eRidge*(PH.rhoM-PH.rhoW)-(PH.rhoM-PH.rhoOc)*7e3)/PH.rhoM,ROOT=(PH.rhoM-PH.rhoCc)*C.CONT_BASE*1e3+PH.rhoM*EA;
function gpeRaw(K,Tkm,age){const g=PH.g,rm=PH.rhoM,T=Tkm*1e3;let rc,H,lith;
  if(K){rc=PH.rhoCc;H=PH.rootH;lith=ROOT;}else{rc=PH.rhoOc;H=hth(age);lith=HALF_ARD*H;}
  const mu=(rm-rc)*T-lith;let e=EA+mu/rm;if(e<0)e=(rm*EA+mu)/(rm-PH.rhoW);
  // layers top down; for a layer of thickness h whose density goes linearly from a to b, with pressure P on top:
  // integral of pressure = P h + g (a h^2/2 + (b - a) h^2/6), and the pressure below it = P + g (a + b) h / 2
  let P=0,I=0,z=e<0?0:-e;
  if(e<0){const h=-e,a=PH.rhoW;I+=g*a*h*h/2;P+=g*a*h;z+=h;}
  I+=P*T+g*rc*T*T/2;P+=g*rc*T;z+=T;
  if(H>0){const a=rm+2*lith/H;I+=P*H+g*(a*H*H/2+(rm-a)*H*H/6);P+=g*(a+rm)/2*H;z+=H;}
  const h=GPE_D-z;I+=P*h+g*rm*h*h/2;return I;}
const GPE0=gpeRaw(0,7,0);
function gpe(K,Tkm,age){return gpeRaw(K,Tkm,age)-GPE0;}
// symmetric positive-definite solve by Cholesky (n small: three unknowns per plate)
function cholSolve(A,b,n){for(let j=0;j<n;j++){let d=A[j*n+j];for(let k=0;k<j;k++)d-=A[j*n+k]*A[j*n+k];d=Math.sqrt(Math.max(d,1e-300));A[j*n+j]=d;
    for(let i=j+1;i<n;i++){let s=A[i*n+j];for(let k=0;k<j;k++)s-=A[i*n+k]*A[j*n+k];A[i*n+j]=s/d;}}
  const x=Float64Array.from(b);for(let i=0;i<n;i++){let s=x[i];for(let k=0;k<i;k++)s-=A[i*n+k]*x[k];x[i]=s/A[i*n+i];}
  for(let i=n-1;i>=0;i--){let s=x[i];for(let k=i+1;k<n;k++)s-=A[k*n+i]*x[k];x[i]=s/A[i*n+i];}return x;}

// The lattice core builds the starting state (so both cores start identically); toMarkers() then puts one marker
// in each fixed cell, jittered within it, and drops the plates' lattices. From then on the markers are the crust.
class MeshWorld extends World{
  constructor(o){super(o);this.core="mesh";this.forces=o.kin!=="walk";this.toMarkers();
    // Stage 2: velocities come from the force balance, so the start's random rotations are dropped; with no slabs yet,
    // the first motions come from gravitational potential energy alone (no spin-up)
    if(this.forces){for(const P of this.plates)P.W.fill(0);this.classifyInterfaces();this.solveVelocities();this.updateStats(true);}}
  allocMarkers(cap){this.mcap=cap;for(const k in MFIELDS){const w=MWIDTH[k]||1;this[k]=new MFIELDS[k](cap*w);}this.freeM=[];this.nSlots=0;this.nextId=1;}
  growMarkers(){const cap=Math.ceil(this.mcap*1.5);for(const k in MFIELDS){const w=MWIDTH[k]||1,a=new MFIELDS[k](cap*w);a.set(this[k]);this[k]=a;}this.mcap=cap;this.tri.growV(cap);}
  newMarker(kind){let m;if(this.freeM.length)m=this.freeM.pop();else{if(this.nSlots>=this.mcap)this.growMarkers();m=this.nSlots++;}
    this.mKind[m]=kind;this.mId[m]=this.nextId++;this.arcAdd[m]=0;this.collAdd[m]=0;this.mBn[3*m]=this.mBn[3*m+1]=this.mBn[3*m+2]=0;this.iN[m]=0;this.mBorn[m]=this.t;
    this.sA[m]=this.sTh[m]=this.sOc[m]=this.sCc[m]=0;this.sFed[m]=this.sFc[m]=-1e9;return m;}
  freeMarker(m){this.tri.remove(m);this.mKind[m]=0;this.freeM.push(m);}
  setPos(m,x,y,z){const l=Math.hypot(x,y,z),T=this.tri;T.X[m]=x/l;T.Y[m]=y/l;T.Z[m]=z/l;}
  toMarkers(){
    const {N,cube,plates,owner,wT,wB,wS,wK,wCol}=this,pos=cube.pos;
    this.h0=Math.sqrt(4*Math.PI/N);this.cellA=4*Math.PI*R_KM*R_KM/N;
    const cap=Math.ceil(N*1.3);this.tri=new SphereTri(cap,2*cap+16);this.allocMarkers(cap);
    const order=[];
    for(let i=0;i<N;i++){const m=this.newMarker(INTERIOR),x=pos[3*i],y=pos[3*i+1],z=pos[3*i+2];
      // jitter within the cell (a fixed hash): breaks the grid's co-circular quadruples and its lattice directions
      const e1=norm3(Math.abs(z)<0.9?cross([0,0,1],[x,y,z]):cross([1,0,0],[x,y,z])),e2=cross([x,y,z],e1),[a,b]=hash2(i,0x51ab),j=0.35*this.h0;
      this.setPos(m,x+j*(a*e1[0]+b*e2[0]),y+j*(a*e1[1]+b*e2[1]),z+j*(a*e1[2]+b*e2[2]));
      this.mPlate[m]=owner[i];this.mA[m]=cube.area[i]*R_KM*R_KM;this.mT[m]=wT[i];this.mB[m]=wB[i];this.mS[m]=wS[i];this.mK[m]=wK[i];this.mCol[m]=wCol[i];
      this.mBorn[m]=0;order.push(m);}
    for(const P of plates){P.T=P.B=P.S=P.K=P.col=P.clr=P.seen=null;}
    this.tri.build(order);
    this.ledger={};this.mStats={rebuilds:0,substeps:0,consumed:0,created:0,ifaceIns:0,ifaceDel:0};
    this.maintainInterfaces();this.buildAdjacency();this.assignOceanAreas();this.computeMarkerElev();this.computeView();this.updateStats(true);}
  addPlate(W){if(!this.mKind)return super.addPlate(W);const id=this.plates.length;
    const P={id,Q:new Float64Array([1,0,0,0,1,0,0,0,1]),W:W.slice(),w0:Math.max(this.wmin*3,Math.min(this.wmax*0.75,Math.hypot(W[0],W[1],W[2]))),area:0,contArea:0,maxArea:0,c:[1,0,0],alive:true,hue:(id*137.508+20)%360};
    this.plates.push(P);return P;}
  led(k,v){this.ledger[k]=(this.ledger[k]||0)+v;}
  alive(m){return this.mKind[m]!==0;}
  // ---------- neighbours ----------
  // CSR adjacency from the triangles (each directed edge a->b once, with half the cotangent of the angle opposite it),
  // then the material neighbours of each interior marker: its interior neighbours, plus, through each interface
  // marker beside it, the interior markers of the other plates around that interface marker (a "virtual" edge
  // across the boundary, weight 1/4 each way), so crust processes cross plate boundaries as they do on the lattice.
  buildAdjacency(){const T=this.tri,tv=T.tv,ns=this.nSlots,X=T.X,Y=T.Y,Z=T.Z,kind=this.mKind,pl=this.mPlate;
    const off=this.off&&this.off.length>=ns+1?this.off:new Int32Array(ns+1+1024);off.fill(0,0,ns+1);
    for(let t=0;t<T.nT;t++){if(!T.tAlive[t])continue;const o=3*t;off[tv[o]+1]++;off[tv[o+1]+1]++;off[tv[o+2]+1]++;}
    for(let i=0;i<ns;i++)off[i+1]+=off[i];const ne=off[ns];
    if(!this.adj||this.adj.length<ne){this.adj=new Int32Array(ne+4096);this.cw=new Float32Array(ne+4096);}
    if(!this._fill||this._fill.length<ns)this._fill=new Int32Array(ns+1024);
    const adj=this.adj,cw=this.cw,fill=this._fill;fill.set(off.subarray(0,ns));
    // Each triangle's area and the half-cotangents of its angles are cached on the triangulation: a triangle whose
    // three markers ride one rigid plate keeps its shape, so only rewritten triangles and those straddling plates or
    // touching an interface marker (invalidated after each motion) are recomputed.
    const tc=T.tCache,tok=T.tOk;
    for(let t=0;t<T.nT;t++){if(!T.tAlive[t])continue;const o=3*t,q=4*t;
      if(!tok[t]){tc[q]=T.area(t);
        for(let k=0;k<3;k++){const a=tv[o+k],b=tv[o+(k+1)%3],c=tv[o+(k+2)%3];
          // half the cotangent of the angle at c, opposite the directed edge a -> b
          const ax=X[a]-X[c],ay=Y[a]-Y[c],az=Z[a]-Z[c],bx=X[b]-X[c],by=Y[b]-Y[c],bz=Z[b]-Z[c],d=ax*bx+ay*by+az*bz;
          const cx=ay*bz-az*by,cy=az*bx-ax*bz,cz=ax*by-ay*bx;tc[q+1+k]=0.5*d/(Math.sqrt(cx*cx+cy*cy+cz*cz)||1e-30);}
        tok[t]=1;}
      const a=tv[o],b=tv[o+1],c=tv[o+2];adj[fill[a]]=b;cw[fill[a]++]=tc[q+1];adj[fill[b]]=c;cw[fill[b]++]=tc[q+2];adj[fill[c]]=a;cw[fill[c]++]=tc[q+3];}
    this.off=off;
    // material neighbours
    const st=this.stamp,mo=this.mo&&this.mo.length>=ns+1?this.mo:new Int32Array(ns+1+1024);let ma=this.ma&&this.ma.length>=2*ne?this.ma:new Int32Array(2*ne+4096),mw=this.mw&&this.mw.length>=ma.length?this.mw:new Float32Array(ma.length);
    let n=0,stampN=this._stampN||0;
    for(let m=0;m<ns;m++){mo[m]=n;if(kind[m]!==INTERIOR)continue;stampN++;st[m]=stampN;
      for(let e=off[m];e<off[m+1];e++){const q=adj[e];
        if(kind[q]===INTERIOR){if(st[q]!==stampN){st[q]=stampN;if(n>=ma.length){const a=new Int32Array(2*ma.length);a.set(ma);ma=a;const w=new Float32Array(2*mw.length);w.set(mw);mw=w;}ma[n]=q;mw[n++]=cw[e];}}
        else if(kind[q]===IFACE){for(let f=off[q];f<off[q+1];f++){const r=adj[f];if(kind[r]!==INTERIOR||pl[r]===pl[m]||st[r]===stampN)continue;st[r]=stampN;
          if(n>=ma.length){const a=new Int32Array(2*ma.length);a.set(ma);ma=a;const w=new Float32Array(2*mw.length);w.set(mw);mw=w;}ma[n]=r;mw[n++]=0.25;}}}}
    mo[ns]=n;this._stampN=stampN;this.mo=mo;this.ma=ma;this.mw=mw;}
  dualArea(v){const T=this.tri,tris=[],vs=[];T.ring(v,tris,vs);let a=0;for(const t of tris)a+=T.area(t)/3;return a;}
  ringOf(v){const out=[];const o=this.off;for(let e=o[v];e<o[v+1];e++)out.push(this.adj[e]);return out;}
  dist2(a,b){const X=this.tri.X,Y=this.tri.Y,Z=this.tri.Z,dx=X[a]-X[b],dy=Y[a]-Y[b],dz=Z[a]-Z[b];return dx*dx+dy*dy+dz*dz;}
  // the interior marker of `plate` nearest to point p, searched two rings out from the triangle holding p
  nearestOf(px,py,pz,plate,hint){const T=this.tri,t=T.locate(px,py,pz,hint>=0?T.vt[hint]:-1);if(t<0)return -1;
    const X=T.X,Y=T.Y,Z=T.Z,kind=this.mKind,pl=this.mPlate,seen=new Set(),front=[T.tv[3*t],T.tv[3*t+1],T.tv[3*t+2]];let best=-1,bd=Infinity;
    for(let ring=0;ring<3&&front.length;ring++){const next=[];
      for(const v of front){if(seen.has(v))continue;seen.add(v);
        if(kind[v]===INTERIOR&&pl[v]===plate){const d=(X[v]-px)**2+(Y[v]-py)**2+(Z[v]-pz)**2;if(d<bd){bd=d;best=v;}}
        const tris=[],vs=[];T.ring(v,tris,vs);for(const q of vs)if(!seen.has(q))next.push(q);}
      if(best>=0&&ring>=1)break;front.length=0;front.push(...next);}
    return best;}
  // ---------- interfaces ----------
  // Interface markers are boundary nodes with no material and no area: the crust's area is carried by the interior
  // markers alone, which tile the sphere (the audit's area closure). In the view, the nearest interior marker covers them.
  // Every edge between interior markers of two plates gets an interface marker at its midpoint, with one side per
  // plate. Interface markers whose neighbourhood no longer holds two plates are removed, as are those crowding
  // another interface marker; sides follow the plates around them (a third side at a triple junction).
  maintainInterfaces(){const T=this.tri,kind=this.mKind,pl=this.mPlate,h=this.h0,plates=this.plates,tris=[],vs=[],cp=new Int32Array(16),cc=new Int32Array(16),lim=(0.4*h)**2;
    // remove or re-side interface markers (once)
    for(let m=0;m<this.nSlots;m++){if(kind[m]!==IFACE)continue;T.ring(m,tris,vs);let nc=0,crowd=false;
      for(const v of vs){if(kind[v]===INTERIOR&&plates[pl[v]].alive){const p=pl[v];let k=0;while(k<nc&&cp[k]!==p)k++;if(k===nc){if(nc<16){cp[nc]=p;cc[nc++]=0;}else continue;}cc[k]++;}
        else if(kind[v]===IFACE&&v<m&&this.dist2(v,m)<lim)crowd=true;}
      if(nc<2||crowd){this.freeMarker(m);this.mStats.ifaceDel++;continue;}
      // the three best-represented plates become the sides (kept in plate order); a side already there keeps its position
      const want=[];for(let k=0;k<nc;k++)want.push(k);want.sort((a,b)=>cc[b]-cc[a]||cp[a]-cp[b]);const ps=want.slice(0,3).map(k=>cp[k]).sort((a,b)=>a-b);
      const n=this.iN[m];let same=n===ps.length;for(let s2=0;s2<n&&same;s2++)if(this.iPlate[3*m+s2]!==ps[s2])same=false;if(same)continue;
      const old=[];for(let s2=0;s2<n;s2++){const q=3*m+s2;old.push([this.iPlate[q],this.iSX[q],this.iSY[q],this.iSZ[q]]);}
      this.iN[m]=ps.length;for(let s2=0;s2<ps.length;s2++){const q=3*m+s2,o=old.find(e=>e[0]===ps[s2]);this.iPlate[q]=ps[s2];
        this.iSX[q]=o?o[1]:T.X[m];this.iSY[q]=o?o[2]:T.Y[m];this.iSZ[q]=o?o[3]:T.Z[m];}}
    // insert at edges joining two plates' interior markers (insertion can expose new ones nearby, hence the passes)
    for(let pass=0;pass<4;pass++){const tv=T.tv,todo=[],ta=T.tAlive,nT=T.nT;
      const edge=(a,b)=>{if(a<b&&kind[a]===INTERIOR&&kind[b]===INTERIOR&&pl[a]!==pl[b])todo.push(a,b);};
      for(let t=0;t<nT;t++){if(!ta[t])continue;const o=3*t,a=tv[o],b=tv[o+1],c=tv[o+2],pa=pl[a];
        if(kind[a]===INTERIOR&&kind[b]===INTERIOR&&kind[c]===INTERIOR&&pl[b]===pa&&pl[c]===pa)continue; // inside one plate: most triangles
        edge(b,c);edge(c,a);edge(a,b);}
      if(!todo.length)break;
      for(let q=0;q<todo.length;q+=2){const a=todo[q],b=todo[q+1];if(!kind[a]||!kind[b])continue;
        const m=this.newMarker(IFACE);this.setPos(m,T.X[a]+T.X[b],T.Y[a]+T.Y[b],T.Z[a]+T.Z[b]);
        if(!T.insert(m,T.vt[a])){this.mKind[m]=0;this.freeM.push(m);continue;}
        const ps=[pl[a],pl[b]].sort((x,y)=>x-y);this.iN[m]=2;this.iKind[m]=TRANSFORM;this.iWin[m]=-1;
        for(let s2=0;s2<2;s2++){this.iPlate[3*m+s2]=ps[s2];this.iSX[3*m+s2]=T.X[m];this.iSY[3*m+s2]=T.Y[m];this.iSZ[3*m+s2]=T.Z[m];}
        this.mStats.ifaceIns++;}}}
  // the frontier marker of a side: the nearest interior marker of that plate around the interface marker
  frontier(m,p){const o=this.off,adj=this.adj;let best=-1,bd=Infinity;
    for(let e=o[m];e<o[m+1];e++){const v=adj[e];if(this.mKind[v]!==INTERIOR||this.mPlate[v]!==p)continue;const d=this.dist2(v,m);if(d<bd){bd=d;best=v;}}return best;}
  // the lattice core's overlap rule: continental beats oceanic, then thicker continent, then younger ocean
  // (with forces, L0_MESH_SPEC.md §4.6: between two oceanic columns the more buoyant one overrides, b of §4.1)
  wins(a,b){const K=this.mK;if(K[a]!==K[b])return K[a]>K[b];if(K[a])return this.mT[a]>=this.mT[b];
    if(this.forces)return buoy(0,this.mT[a],this.t-this.mB[a])>=buoy(0,this.mT[b],this.t-this.mB[b]);return this.mB[a]>=this.mB[b];}
  // kinds of interface at the start of a step, from the plates' relative velocity across the boundary
  classifyInterfaces(){const T=this.tri,X=T.X,Y=T.Y,Z=T.Z,plates=this.plates,o=this.off,adj=this.adj,pl=this.mPlate,kind=this.mKind,thr=2/R_KM;
    for(let m=0;m<this.nSlots;m++){if(kind[m]!==IFACE)continue;
      this.classifyOne(m);
      // sides start each step at the marker and ride their plates through it
      for(let s=0;s<this.iN[m];s++){const q=3*m+s;this.iSX[q]=X[m];this.iSY[q]=Y[m];this.iSZ[q]=Z[m];}}
    // One polarity per stretch of convergent boundary: interface markers of one plate pair that touch form a segment,
    // and the segment takes its markers' majority verdict (neighbours that disagreed would ride different plates
    // and run into each other).
    const ns=this.nSlots,seg=new Int32Array(ns).fill(-1),find=a=>{while(seg[a]!==a){seg[a]=seg[seg[a]];a=seg[a];}return a;};
    const conv=this._conv=[];for(let m=0;m<ns;m++)if(kind[m]===IFACE&&this.iKind[m]===CONVERGENT){seg[m]=m;conv.push(m);} // the list the substeps scan
    for(let m=0;m<ns;m++){if(seg[m]<0)continue;for(let e=o[m];e<o[m+1];e++){const v=adj[e];if(seg[v]<0||this.iPlate[3*v]!==this.iPlate[3*m]||this.iPlate[3*v+1]!==this.iPlate[3*m+1])continue;
      const a=find(m),b=find(v);if(a!==b)seg[Math.max(a,b)]=Math.min(a,b);}}
    const vote=new Map();for(let m=0;m<ns;m++){if(seg[m]<0)continue;const r=find(m);vote.set(r,(vote.get(r)||0)+(this.iWin[m]===0?1:-1));}
    for(let m=0;m<ns;m++){if(seg[m]<0)continue;const v=vote.get(find(m));this.iWin[m]=v>0?0:v<0?1:this.iWin[find(m)];}}
  classifyOne(m){const T=this.tri,X=T.X,Y=T.Y,Z=T.Z,plates=this.plates,o=this.off,adj=this.adj,pl=this.mPlate,kind=this.mKind,thr=2/R_KM;
    {if(this.iN[m]!==2){this.iKind[m]=JUNCTION;this.iWin[m]=-1;return;}
      const p=this.iPlate[3*m],q=this.iPlate[3*m+1];let cp=[0,0,0],cq=[0,0,0],np=0,nq=0;
      for(let e=o[m];e<o[m+1];e++){const v=adj[e];if(kind[v]!==INTERIOR)continue;if(pl[v]===p){cp[0]+=X[v];cp[1]+=Y[v];cp[2]+=Z[v];np++;}else if(pl[v]===q){cq[0]+=X[v];cq[1]+=Y[v];cq[2]+=Z[v];nq++;}}
      if(!np||!nq){this.iKind[m]=TRANSFORM;this.iWin[m]=-1;return;}
      const x=[X[m],Y[m],Z[m]];let n=[cp[0]/np-cq[0]/nq,cp[1]/np-cq[1]/nq,cp[2]/np-cq[2]/nq];const d=n[0]*x[0]+n[1]*x[1]+n[2]*x[2];n=norm3([n[0]-d*x[0],n[1]-d*x[1],n[2]-d*x[2]]);
      this.iNrm[3*m]=n[0];this.iNrm[3*m+1]=n[1];this.iNrm[3*m+2]=n[2]; // from plate q toward plate p
      const Wp=plates[p].W,Wq=plates[q].W,v=cross([Wp[0]-Wq[0],Wp[1]-Wq[1],Wp[2]-Wq[2]],x),r=v[0]*n[0]+v[1]*n[1]+v[2]*n[2];
      if(r>thr){this.iKind[m]=RIDGE;this.iWin[m]=-1;}
      else if(r<-thr){this.iKind[m]=CONVERGENT;const fp=this.frontier(m,p),fq=this.frontier(m,q);this.iWin[m]=fp<0?1:fq<0?0:(this.wins(fp,fq)?0:1);}
      else{this.iKind[m]=TRANSFORM;this.iWin[m]=-1;}}}
  // ---------- motion ----------
  // Rigid rotation of each plate's markers; interface sides ride their plates, the marker itself sits on the
  // winner's side at a convergent boundary and midway elsewhere. The step is cut into substeps so that boundary
  // markers move at most a quarter spacing relative to each other; a substep that turns a triangle over is
  // halved. Flips keep the triangulation Delaunay; loser markers reaching a convergent interface are consumed.
  placeAt(f){const T=this.tri,X=T.X,Y=T.Y,Z=T.Z,np=this.plates.length;
    // every plate's rotation for fraction f of the step, flat (9 numbers a plate)
    const R=this._Rf&&this._Rf.length>=9*np?this._Rf:(this._Rf=new Float64Array(9*np+90));
    for(let p=0;p<np;p++){const P=this.plates[p];if(!P.alive)continue;R.set(rotM(P.W,f*Math.hypot(P.W[0],P.W[1],P.W[2])*this.dt),9*p);}
    const kind=this.mKind,pl=this.mPlate,X0=this.X0,Y0=this.Y0,Z0=this.Z0,S0=this._S0,iP=this.iPlate,iN=this.iN,ns=this.nSlots,iKind=this.iKind,iWin=this.iWin,iSX=this.iSX,iSY=this.iSY,iSZ=this.iSZ;
    const alive=this._palive&&this._palive.length>=np?this._palive:(this._palive=new Uint8Array(np+64));for(let p=0;p<np;p++)alive[p]=this.plates[p].alive?1:0;
    for(let m=0;m<ns;m++){const k=kind[m];if(!k)continue;
      if(k===INTERIOR){const o=9*pl[m],x=X0[m],y=Y0[m],z=Z0[m];X[m]=R[o]*x+R[o+1]*y+R[o+2]*z;Y[m]=R[o+3]*x+R[o+4]*y+R[o+5]*z;Z[m]=R[o+6]*x+R[o+7]*y+R[o+8]*z;continue;}
      const n=iN[m];let sx=0,sy=0,sz=0;const w=iKind[m]===CONVERGENT?iWin[m]:-1;
      for(let s=0;s<n;s++){const p=iP[3*m+s];if(!alive[p])continue;const o=9*p,q=3*m+s,x=S0[3*q],y=S0[3*q+1],z=S0[3*q+2];
        const nx=R[o]*x+R[o+1]*y+R[o+2]*z,ny=R[o+3]*x+R[o+4]*y+R[o+5]*z,nz=R[o+6]*x+R[o+7]*y+R[o+8]*z;iSX[q]=nx;iSY[q]=ny;iSZ[q]=nz;
        if(w<0||w===s){sx+=nx;sy+=ny;sz+=nz;}}
      const l=Math.hypot(sx,sy,sz)||1;X[m]=sx/l;Y[m]=sy/l;Z[m]=sz/l;}}
  rebuildMesh(){const T=this.tri,vs=[];for(let m=0;m<this.nSlots;m++)if(this.mKind[m])vs.push(m);T.build(vs);T.broken=false;this.mStats.rebuilds++;
    for(let m=0;m<this.nSlots;m++)if(this.mKind[m]&&T.vt[m]<0){this.mKind[m]=0;this.freeM.push(m);this.mStats.lostInRebuild=(this.mStats.lostInRebuild||0)+1;}}
  mixedTris(){const T=this.tri,tv=T.tv,kind=this.mKind,pl=this.mPlate,out=[];
    for(let t=0;t<T.nT;t++){if(!T.tAlive[t])continue;const o=3*t,a=tv[o],b=tv[o+1],c=tv[o+2];
      if(kind[a]===INTERIOR&&kind[b]===INTERIOR&&kind[c]===INTERIOR&&pl[a]===pl[b]&&pl[a]===pl[c])continue;out.push(t);}
    return out;}
  move(){const T=this.tri,ns=this.nSlots,kind=this.mKind,plates=this.plates,dt=this.dt,h=this.h0;
    this.X0.set(T.X.subarray(0,ns));this.Y0.set(T.Y.subarray(0,ns));this.Z0.set(T.Z.subarray(0,ns));
    const S0=this._S0=new Float64Array(9*ns);for(let q=0;q<3*ns;q++){S0[3*q]=this.iSX[q];S0[3*q+1]=this.iSY[q];S0[3*q+2]=this.iSZ[q];}
    // the largest relative motion across an interface this step, in radians
    let maxRel=0;
    for(let m=0;m<ns;m++){if(kind[m]!==IFACE)continue;const x=[T.X[m],T.Y[m],T.Z[m]],n=this.iN[m];
      for(let a=0;a<n;a++)for(let b=a+1;b<n;b++){const Wa=plates[this.iPlate[3*m+a]].W,Wb=plates[this.iPlate[3*m+b]].W;
        const v=cross([Wa[0]-Wb[0],Wa[1]-Wb[1],Wa[2]-Wb[2]],x);maxRel=Math.max(maxRel,Math.hypot(v[0],v[1],v[2])*dt);}}
    let nsub=Math.max(1,Math.ceil(maxRel/(0.25*h))),fDone=0,df=1/nsub,depth=0;
    while(fDone<1-1e-12){const f=Math.min(1,fDone+df);const mixed=this.mixedTris();this.placeAt(f);
      // every triangle turned over by this substep, and the repair for each (found at the trial positions, applied at
      // the last good ones, all at once)
      let bad=-1;const tv=T.tv,fixes=[];for(const t of mixed){if(!T.tAlive[t])continue;const o=3*t;if(T.orient(tv[o],tv[o+1],tv[o+2])<=0){if(bad<0)bad=t;const fx=this.culprit(t);if(fx)fixes.push(fx);}}
      if(bad>=0){this.placeAt(fDone);
        if(fixes.length){for(const fx of fixes)fx();this.mStats.culprits=(this.mStats.culprits||0)+fixes.length;continue;}
        if(df>1/512){df/=2;depth++;continue;}
        // last resort: rebuild the triangulation at the new positions
        this.placeAt(f);const vs=[];for(let m=0;m<this.nSlots;m++)if(kind[m])vs.push(m);T.build(vs);this.mStats.rebuilds++;}
      // the triangles that may join two plates now: those that did, and those the flips rewrote (every new one)
      let list=null;if(bad<0){T.touched=[];T.legalize(mixed);list=mixed.concat(T.touched);T.touched=null;}
      this.mStats.substeps++;fDone=f;this.consumeAtTrenches(list);}
    this._S0=null;this._R=null;
    if(T.broken)this.rebuildMesh();
    for(const t of this.mixedTris())T.tOk[t]=0; // their shapes changed
    // the plates' orientations (fault strikes are kept in plate frames)
    for(const P of plates){if(!P.alive)continue;const R=rotM(P.W,Math.hypot(P.W[0],P.W[1],P.W[2])*dt),Q=P.Q,M=new Float64Array(9);
      for(let r=0;r<3;r++)for(let k=0;k<3;k++)M[r*3+k]=R[r*3]*Q[k]+R[r*3+1]*Q[3+k]+R[r*3+2]*Q[6+k];
      if(this.steps%64===0){for(let r=0;r<3;r++){for(let q=0;q<r;q++){const d=M[r*3]*M[q*3]+M[r*3+1]*M[q*3+1]+M[r*3+2]*M[q*3+2];for(let k=0;k<3;k++)M[r*3+k]-=d*M[q*3+k];}
        const l=Math.hypot(M[r*3],M[r*3+1],M[r*3+2]);for(let k=0;k<3;k++)M[r*3+k]/=l;}}
      P.Q=M;}}
  // An inverted triangle (at the trial positions) usually has a culprit that should not survive the substep: two
  // interface markers crowding each other (one goes), or a marker of the losing plate at a convergent interface (it is
  // consumed). Returns the repair, to be applied at the last good positions, or null.
  culprit(t){const tv=this.tri.tv,o=3*t,vs=[tv[o],tv[o+1],tv[o+2]],kind=this.mKind,h=this.h0;
    for(let a=0;a<3;a++)for(let b=a+1;b<3;b++){const p=vs[a],q=vs[b];
      if(kind[p]===IFACE&&kind[q]===IFACE&&this.dist2(p,q)<(0.35*h)**2){const v=Math.max(p,q);return()=>{if(kind[v]===IFACE){this.freeMarker(v);this.mStats.ifaceDel++;}};}}
    for(const I of vs){if(kind[I]!==IFACE||this.iKind[I]!==CONVERGENT)continue;const W=this.iPlate[3*I+this.iWin[I]];
      for(const v of vs)if(kind[v]===INTERIOR&&this.mPlate[v]!==W){let side=false;for(let s=0;s<this.iN[I];s++)if(this.iPlate[3*I+s]===this.mPlate[v])side=true;
        if(side)return()=>{if(kind[v]===INTERIOR)this.consume(v,W,I);};}}
    for(let a=0;a<3;a++)for(let b=a+1;b<3;b++){const p=vs[a],q=vs[b];
      if(kind[p]===INTERIOR&&kind[q]===INTERIOR&&this.mPlate[p]!==this.mPlate[q]){const pw=this.wins(p,q),v=pw?q:p,W=this.mPlate[pw?p:q],n=pw?p:q;return()=>{if(kind[v]===INTERIOR)this.consume(v,W,n);};}}
    // otherwise an interface marker is in the way (a ridge's moving at the mean velocity, or one at a boundary that
    // opens at one end and closes at the other): it goes, and upkeep puts interface markers back where plates touch
    let I=-1;for(const v of vs)if(kind[v]===IFACE&&(I<0||this.iKind[v]!==CONVERGENT))I=v;
    if(I>=0)return()=>{if(kind[I]===IFACE){this.freeMarker(I);this.mStats.ifaceDel++;}};
    return null;}
  // loser markers within half a spacing of a convergent interface marker are consumed; so is any interior marker
  // that comes within a third of a spacing of another plate's interior marker (triple junctions, stray overlaps)
  consumeAtTrenches(list){const T=this.tri,kind=this.mKind,pl=this.mPlate,h=this.h0,lim=(0.3*h)**2,lim2=(0.3*h)**2,victims=[];
    const tris=[],vs=[];
    for(const m of this._conv){if(kind[m]!==IFACE||this.iKind[m]!==CONVERGENT)continue;const w=this.iWin[m],W=this.iPlate[3*m+w];
      T.ring(m,tris,vs);for(const v of vs)if(kind[v]===INTERIOR&&pl[v]!==W&&this.dist2(v,m)<lim)victims.push(v,W,m);}
    const tv=T.tv;
    for(const t of list||this.mixedTris()){if(!T.tAlive[t])continue;for(let k=0;k<3;k++){const a=tv[3*t+(k+1)%3],b=tv[3*t+(k+2)%3];
      if(a<b&&kind[a]===INTERIOR&&kind[b]===INTERIOR&&pl[a]!==pl[b]&&this.dist2(a,b)<lim2){const aw=this.wins(a,b);victims.push(aw?b:a,aw?pl[a]:pl[b],aw?a:b);}}}
    for(let q=0;q<victims.length;q+=3){const v=victims[q];if(kind[v]!==INTERIOR)continue;this.consume(v,victims[q+1],victims[q+2]);}}
  // consumption of one loser marker by winner plate Wp, near marker `near`: the lattice core's overlap side effects
  consume(m,Wp,near){const T=this.tri,X=T.X,Y=T.Y,Z=T.Z,x=X[m],y=Y[m],z=Z[m],L=this.plates[this.mPlate[m]],Wn=this.plates[Wp],h=this.h0;
    let r=cross([L.W[0]-Wn.W[0],L.W[1]-Wn.W[1],L.W[2]-Wn.W[2]],[x,y,z]);let rl=Math.hypot(r[0],r[1],r[2]);
    if(rl<1e-12){r=[X[near]-x,Y[near]-y,Z[near]-z];rl=Math.hypot(r[0],r[1],r[2])||1;}r=[r[0]/rl,r[1]/rl,r[2]/rl];
    const A=this.mA[m],at=k=>this.nearestOf(x+r[0]*h*k,y+r[1]*h*k,z+r[2]*h*k,Wp,near),lk=this.cube.lookup;
    if(!this.mK[m]){
      const a=(C.ARC+0.4*this.mS[m])*A;for(let k=1;k<4;k++){const tm=at(k);if(tm>=0)this.arcAdd[tm]+=a*C.ARCW[k];}
      this.trench[lk(x-r[0]*h,y-r[1]*h,z-r[2]*h)]=1;
      const age=this.t-this.mB[m],c=lk(x+r[0]*h*C.SLAB_INBOARD,y+r[1]*h*C.SLAB_INBOARD,z+r[2]*h*C.SLAB_INBOARD);
      this.slab[this.slabQ*this.N+c]+=A*Math.min(1,Math.sqrt(Math.max(0,age)/80));
      if(this.forces)this.feedSlab(at(0),A,this.mT[m],hth(age),0);
    } else {
      const bn=[this.mBn[3*m],this.mBn[3*m+1],this.mBn[3*m+2]],nl=Math.hypot(bn[0],bn[1],bn[2]);
      const f=nl>0.5?Math.abs(r[0]*bn[0]+r[1]*bn[1]+r[2]*bn[2])/nl:1,V=this.mT[m]*A;
      for(let k=0;k<4;k++){const w=f*C.COLLW[k]+(k===0?1-f:0),tm=at(k);if(tm>=0)this.collAdd[tm]+=V*w;else this.led("lostCollisionNoTarget",V*w);}
      this.led("contCollisionOut",-V);if(this.forces)this.feedSlab(at(0),A,this.mT[m],0,1);
      // the weld counter counts collided area in cells of the page's default grid (n = 96, where its thresholds were
      // set), so welding does not depend on resolution
      const key=Math.min(Wp,L.id)*4096+Math.max(Wp,L.id);this.pairs.set(key,(this.pairs.get(key)||0)+A/PH.weldCellKm2);this.pairLoc.set(key,[x,y,z]);}
    this.freeMarker(m);this.mStats.consumed++;}
  // Gaps: where an interface marker that is not at a convergent boundary has no interior marker of one of its plates
  // within GAP x spacing, a new oceanic marker of that plate is born between them (at a ridge, as the plates open; at
  // a continental rift, between the continents). Its area is set from the geometry, like all oceanic markers'.
  fillGaps(){const T=this.tri,X=T.X,Y=T.Y,Z=T.Z,kind=this.mKind,pl=this.mPlate,h=this.h0,lim=(1.3*h)**2,born=[],tris=[],vs=[],tris2=[],vs2=[];
    for(let m=0;m<this.nSlots;m++){if(kind[m]!==IFACE||this.iKind[m]===CONVERGENT)continue;const n=this.iN[m];
      // nearest interior marker of each side's plate, two rings out
      const best=[-1,-1,-1],bd=[Infinity,Infinity,Infinity];T.ring(m,tris,vs);
      const look=v=>{if(kind[v]!==INTERIOR)return;for(let s=0;s<n;s++)if(pl[v]===this.iPlate[3*m+s]){const d=this.dist2(v,m);if(d<bd[s]){bd[s]=d;best[s]=v;}}};
      for(const v of vs)look(v);
      let far=false;for(let s=0;s<n;s++)if(best[s]<0||bd[s]>lim)far=true;
      if(far)for(const v of vs){T.ring(v,tris2,vs2);for(const q of vs2)look(q);} // two rings out only when one ring did not settle it
      for(let s=0;s<n;s++){const f=best[s];if(f<0||bd[s]<=lim)continue;
        const d=Math.sqrt(bd[s]),k=Math.min(0.6*h,d/2)/d;born.push([this.iPlate[3*m+s],X[m]+(X[f]-X[m])*k,Y[m]+(Y[f]-Y[m])*k,Z[m]+(Z[f]-Z[m])*k,m]);}}
    for(const [p,x,y,z,near] of born){if(kind[near]!==IFACE)continue;
      const m=this.newMarker(INTERIOR);this.setPos(m,x,y,z);
      if(!T.insert(m,T.vt[near])){this.mKind[m]=0;this.freeM.push(m);continue;}
      this.mPlate[m]=p;this.mA[m]=0;this.mT[m]=C.OCEAN_T;this.mB[m]=this.t;this.mS[m]=0;this.mK[m]=0;this.mCol[m]=this.pool.make(ROCK.OCEAN,this.t,C.OCEAN_T,p);
      this.mStats.created++;}}
  // Oceanic crust is made and destroyed by definition, so its area is simply the space it occupies: each oceanic
  // marker's share of the triangles around it (a third of each), plus an equal share of the area around each interface
  // marker beside it. Continental markers keep their carried area (their volume is conserved, through the ledger).
  assignOceanAreas(){const T=this.tri,tv=T.tv,ns=this.nSlots,kind=this.mKind,off=this.off,adj=this.adj,R2=R_KM*R_KM;
    const dual=this._dual&&this._dual.length>=ns?this._dual:(this._dual=new Float64Array(ns+1024));dual.fill(0,0,ns);
    for(let t=0;t<T.nT;t++){if(!T.tAlive[t])continue;const a=T.tCache[4*t]/3,o=3*t;dual[tv[o]]+=a;dual[tv[o+1]]+=a;dual[tv[o+2]]+=a;} // areas cached by buildAdjacency
    for(let m=0;m<ns;m++){if(kind[m]!==IFACE)continue;let c=0;for(let e=off[m];e<off[m+1];e++)if(kind[adj[e]]===INTERIOR)c++;if(!c)continue;
      const sh=dual[m]/c;for(let e=off[m];e<off[m+1];e++)if(kind[adj[e]]===INTERIOR)dual[adj[e]]+=sh;}
    for(let m=0;m<ns;m++)if(kind[m]===INTERIOR&&!this.mK[m])this.mA[m]=dual[m]*R2;}
  // ---------- forces and the torque balance (Stage 2, L0_MESH_SPEC.md §4, §5.1) ----------
  // Slab records hang beneath the overriding plate: each consumed marker adds its area, thermal deficit and crust to the
  // record of the nearest interior marker of the plate that consumed it (that marker rides the overriding plate, as the
  // trench does). A record is the slab beneath about one marker spacing of trench.
  // Continental lithosphere's mantle root is taken as neutrally buoyant (§4.1); its crust is recorded in sCc. Two clocks:
  // sFed, when ocean last went down (or the record was last found behind an active oceanic trench); sFc, when
  // continental crust last went down.
  feedSlab(f,Akm2,Tkm,h,cont){if(f<0)return;const A=Akm2*1e6;this.sA[f]+=A;
    if(cont){this.sCc[f]+=A*Tkm*1e3;this.sFc[f]=this.t;}else{this.sTh[f]+=A*h;this.sOc[f]+=A*Tkm*1e3;this.sFed[f]=this.t;}}
  // the slab's pull (N), §4.3: the thermal deficit, plus oceanic crust that is eclogite below lEcl and buoyant above it,
  // minus the buoyancy of continental crust dragged down with an arriving margin. Consumption moves continental crust's
  // volume into the overriding plate's orogen (the ledger), so the slab is not given that crust to keep: it counts only
  // the plug a margin forms at the top of the slab before it is scraped off and exhumed, at most lPlug long [rule;
  // ultra-high-pressure rocks show continental crust reaching 100-200 km, measured]. Length is area over one spacing.
  slabForce(f){const a=this.sA[f];if(a<=0)return 0;const w=this.h0*R_KM*1e3,L=a/w,fd=L>PH.lEcl?1-PH.lEcl/L:0,plug=Math.min(this.sCc[f],C.CONT_BASE*1e3*PH.lPlug*w);
    return PH.g*(HALF_ARD*this.sTh[f]+this.sOc[f]*(fd*PH.drhoEcl-(1-fd)*(PH.rhoM-PH.rhoOc))-(PH.rhoM-PH.rhoCc)*plug);}
  // Slabs longer than lMax pass their deepest part on (scaled, keeping their mean make-up). Break-off [rule standing in
  // for §4.6]: once no ocean has gone down behind a record (nor has it been found behind an active oceanic trench) for
  // a step, its oceanic part detaches with e-folding time tauDetach; a buoyant plug of continental crust left behind
  // is exhumed on the same time scale once no more continental crust goes down. So a collision runs: the margin's
  // plug builds against the pull, the oceanic slab breaks off, the plug is left pushing back, then it too is gone.
  slabUpkeep(){const w=this.h0*R_KM*1e3,cap=PH.lMax*w,dec=Math.exp(-this.dt/PH.tauDetach),t=this.t,dt=this.dt,kind=this.mKind;
    for(let m=0;m<this.nSlots;m++){if(kind[m]!==INTERIOR||this.sA[m]<=0)continue;let f=1,fo=1,fc=1;
      if(this.sA[m]>cap)f=cap/this.sA[m];
      if(t-this.sFed[m]>1.5*dt)fo=dec;if(t-this.sFc[m]>1.5*dt)fc=dec;
      if(f===1&&fo===1&&fc===1)continue;
      // area goes with its parts: the oceanic share (by crust volume over 7 km) with fo, the rest with fc
      const aO=Math.min(this.sA[m],this.sOc[m]/(C.OCEAN_T*1e3)),aC=this.sA[m]-aO;
      this.sA[m]=f*(aO*fo+aC*fc);this.sTh[m]*=f*fo;this.sOc[m]*=f*fo;this.sCc[m]*=f*fc;
      if(this.sA[m]<1e-3*w*w)this.sA[m]=this.sTh[m]=this.sOc[m]=this.sCc[m]=0;}}
  // the length of boundary each interface marker stands for: half the distance between its two nearest interface
  // neighbours (which straightens the zigzag of edge midpoints: a great circle comes out within a few per cent)
  ifaceLengths(){const off=this.off,adj=this.adj,kind=this.mKind,h=this.h0,Rm=R_KM*1e3;
    for(let m=0;m<this.nSlots;m++){if(kind[m]!==IFACE)continue;let d1=Infinity,d2=Infinity,v1=-1,v2=-1;
      for(let e=off[m];e<off[m+1];e++){const v=adj[e];if(kind[v]!==IFACE)continue;const d=this.dist2(v,m);if(d<d1){d2=d1;v2=v1;d1=d;v1=v;}else if(d<d2){d2=d;v2=v;}}
      let w=v2>=0?Math.sqrt(this.dist2(v1,v2))/2:v1>=0?Math.sqrt(d1):h;w=Math.min(2*h,Math.max(0.2*h,w));this.iLen[m]=w*Rm;}}
  // One angular velocity per plate from the balance of torques (Stokes flow: no inertia), §5.1. Every force is local
  // and linear in the velocities except suction, which uses the previous step's velocities:
  //   basal drag (per marker area, stronger under continents), slab resistance (per trench interface marker, on the
  //   subducting plate's velocity into the trench), interface coupling (on the relative velocity of the two sides);
  //   driving: slab pull (at the trench, toward the overriding plate), the gradient of GPE (a body force over every
  //   triangle), and trench suction on the overriding plate where the slab sinks faster than the plate feeds it.
  // The system is symmetric positive definite, 3 unknowns a plate; the frame is the mantle's (no net rotation removed).
  solveVelocities(){const plates=this.plates,np=plates.length,T=this.tri,X=T.X,Y=T.Y,Z=T.Z,kind=this.mKind,pl=this.mPlate,ns=this.nSlots,t=this.t;
    const R=R_KM*1e3,R2=R*R,MY=PH.MYR,off=this.off,adj=this.adj,U=this.mU;
    const idx=new Int32Array(np).fill(-1),cnt=new Int32Array(np);for(let m=0;m<ns;m++)if(kind[m]===INTERIOR)cnt[pl[m]]++;
    let n=0;for(let p=0;p<np;p++)if(plates[p].alive&&cnt[p]>0)idx[p]=n++;const N3=3*n;let K=new Float64Array(N3*N3),b=new Float64Array(N3);
    if(!n)return;
    // diagnostics (this.budget): each term's own matrix or torque vector, to say afterwards what balanced what
    const bud=this.budget?{}:null,term=nm=>{if(!bud)return;if(!bud[nm])bud[nm]=nm[0]==="K"?new Float64Array(N3*N3):new Float64Array(N3);if(nm[0]==="K")K=bud[nm];else b=bud[nm];};
    const addProj=(i,j,c,x,y,z)=>{const o=3*i*N3+3*j;K[o]+=c*(1-x*x);K[o+1]-=c*x*y;K[o+2]-=c*x*z;K[o+N3]-=c*y*x;K[o+N3+1]+=c*(1-y*y);K[o+N3+2]-=c*y*z;
      K[o+2*N3]-=c*z*x;K[o+2*N3+1]-=c*z*y;K[o+2*N3+2]+=c*(1-z*z);};
    const addOuter=(i,c,a)=>{const o=3*i*N3+3*i;for(let r=0;r<3;r++)for(let q=0;q<3;q++)K[o+r*N3+q]+=c*a[r]*a[q];};
    const addTq=(i,f,x,y,z)=>{b[3*i]+=f*x;b[3*i+1]+=f*y;b[3*i+2]+=f*z;};
    // basal drag (per plate: sum of c (I - r r^T) over its markers, from the moments of c), and each column's GPE
    term("Kbasal");const mom=new Float64Array(7*n),cb=PH.Cb*1e6*R2,mK=this.mK,mT=this.mT,mB=this.mB,mA=this.mA;
    for(let m=0;m<ns;m++){if(kind[m]!==INTERIOR)continue;const K_=mK[m];U[m]=gpe(K_,mT[m],t-mB[m]);const i=idx[pl[m]];if(i<0)continue;
      const c=cb*(K_?PH.kKeel:1)*mA[m],x=X[m],y=Y[m],z=Z[m],o=7*i;mom[o]+=c;mom[o+1]+=c*x*x;mom[o+2]+=c*y*y;mom[o+3]+=c*z*z;mom[o+4]+=c*x*y;mom[o+5]+=c*x*z;mom[o+6]+=c*y*z;}
    for(let i=0;i<n;i++){const o=7*i,d=3*i*N3+3*i,c=mom[o];
      K[d]+=c-mom[o+1];K[d+N3+1]+=c-mom[o+2];K[d+2*N3+2]+=c-mom[o+3];K[d+1]-=mom[o+4];K[d+N3]-=mom[o+4];K[d+2]-=mom[o+5];K[d+2*N3]-=mom[o+5];K[d+N3+2]-=mom[o+6];K[d+2*N3+1]-=mom[o+6];}
    for(let m=0;m<ns;m++){if(kind[m]!==IFACE)continue;let s=0,c=0;for(let e=off[m];e<off[m+1];e++){const v=adj[e];if(kind[v]===INTERIOR){s+=U[v];c++;}}U[m]=c?s/c:0;}
    // GPE: on each triangle the force is -grad U times its area, = -(R/2) sum_k U_k (nhat x opposite edge_k); its torque
    // goes a third to each corner's plate (an interface corner shares its third among its sides)
    term("bgpe");const tv=T.tv,G=[0,0,0];let gpeF=0;
    for(let q=0;q<T.nT;q++){if(!T.tAlive[q])continue;const o=3*q,a=tv[o],bb=tv[o+1],c=tv[o+2];
      const abx=X[bb]-X[a],aby=Y[bb]-Y[a],abz=Z[bb]-Z[a],acx=X[c]-X[a],acy=Y[c]-Y[a],acz=Z[c]-Z[a];
      let nx=aby*acz-abz*acy,ny=abz*acx-abx*acz,nz=abx*acy-aby*acx;const nl=Math.sqrt(nx*nx+ny*ny+nz*nz)||1;nx/=nl;ny/=nl;nz/=nl;
      // e_a = c-b, e_b = a-c, e_c = b-a; sum U_k e_k
      const ua=U[a],ub=U[bb],uc=U[c];
      const ex=ua*(X[c]-X[bb])+ub*(X[a]-X[c])+uc*(X[bb]-X[a]),ey=ua*(Y[c]-Y[bb])+ub*(Y[a]-Y[c])+uc*(Y[bb]-Y[a]),ez=ua*(Z[c]-Z[bb])+ub*(Z[a]-Z[c])+uc*(Z[bb]-Z[a]);
      const fx=-0.5*R*(ny*ez-nz*ey),fy=-0.5*R*(nz*ex-nx*ez),fz=-0.5*R*(nx*ey-ny*ex);gpeF+=Math.sqrt(fx*fx+fy*fy+fz*fz);
      let cx=X[a]+X[bb]+X[c],cy=Y[a]+Y[bb]+Y[c],cz=Z[a]+Z[bb]+Z[c];const cl=Math.sqrt(cx*cx+cy*cy+cz*cz)||1;cx/=cl;cy/=cl;cz/=cl;
      G[0]=R*(cy*fz-cz*fy)/3;G[1]=R*(cz*fx-cx*fz)/3;G[2]=R*(cx*fy-cy*fx)/3;
      if(kind[a]===INTERIOR&&kind[bb]===INTERIOR&&kind[c]===INTERIOR&&pl[a]===pl[bb]&&pl[a]===pl[c]){const i=idx[pl[a]];if(i>=0)addTq(i,3,G[0],G[1],G[2]);continue;}
      for(let k3=0;k3<3;k3++){const v=tv[o+k3];if(kind[v]===INTERIOR){const i=idx[pl[v]];if(i>=0)addTq(i,1,G[0],G[1],G[2]);}
        else{const k=this.iN[v];let live=0;for(let s2=0;s2<k;s2++)if(idx[this.iPlate[3*v+s2]]>=0)live++;
          for(let s2=0;s2<k;s2++){const i=idx[this.iPlate[3*v+s2]];if(i>=0)addTq(i,1/live,G[0],G[1],G[2]);}}}}
    // interfaces
    this.ifaceLengths();
    const trench=[],tq=new Map(),cpl=[]; // cpl: coupling entries (two plates, coefficient, limit, length, position)
    for(const P of plates){P.perim=0;P.subLen=0;}
    for(let m=0;m<ns;m++){if(kind[m]!==IFACE)continue;const k=this.iN[m],w=this.iLen[m],x=X[m],y=Y[m],z=Z[m],conv=this.iKind[m]===CONVERGENT&&k===2;
      let coll=false,fs=null;
      if(conv){fs=[this.frontier(m,this.iPlate[3*m]),this.frontier(m,this.iPlate[3*m+1])];coll=fs[0]>=0&&fs[1]>=0&&this.mK[fs[0]]&&this.mK[fs[1]];}
      const cI=conv?(coll?PH.cColl:PH.cTrench):PH.cRidge,sI=conv?(coll?PH.sColl:PH.sTrench):PH.sRidge;
      for(let s1=0;s1<k;s1++){const P1=plates[this.iPlate[3*m+s1]];P1.perim+=w/1e3;
        for(let s2=s1+1;s2<k;s2++){const i=idx[this.iPlate[3*m+s1]],j=idx[this.iPlate[3*m+s2]];if(i<0||j<0)continue;cpl.push(i,j,cI,sI,w,x,y,z);}}
      if(!conv)continue;
      const win=this.iWin[m],Lp=this.iPlate[3*m+1-win],Wp=this.iPlate[3*m+win],sg=win===0?1:-1;
      const nx=sg*this.iNrm[3*m],ny=sg*this.iNrm[3*m+1],nz=sg*this.iNrm[3*m+2];          // from the subducting side toward the overriding side
      const a=[y*nz-z*ny,z*nx-x*nz,x*ny-y*nx],fL=fs[1-win],fW=fs[win];
      const hb=fL>=0?(this.mK[fL]?PH.hMax:hth(t-this.mB[fL])):hth(0),Cs=PH.cGeom*PH.etaUM*(1+hb/PH.rBend)**3;
      term("Kresist");if(idx[Lp]>=0)addOuter(idx[Lp],Cs*w*R2,a);
      if(fL>=0&&!this.mK[fL])plates[Lp].subLen+=w/1e3;
      tq.set(m,trench.length);trench.push(m,Lp,Wp,fW,Cs,coll?1:0,nx,ny,nz);}
    // Each slab record pulls at the nearest trench interface marker of its own (overriding) plate within two rings;
    // a record behind a trench where ocean goes down is still attached (behind a collision it is left to break off)
    const Fm=new Float64Array(trench.length/9),seen=this._seenR&&this._seenR.length>=ns?this._seenR:(this._seenR=new Int32Array(ns+1024));let mark=(this._seenN||0);
    for(let r=0;r<ns;r++){if(kind[r]!==INTERIOR||!(this.sA[r]>0))continue;const P=pl[r];let best=-1,bd=Infinity;mark++;seen[r]=mark;
      const look=v=>{const q=tq.get(v);if(q===undefined||trench[q+2]!==P)return;const d=this.dist2(v,r);if(d<bd){bd=d;best=q;}};
      for(let e=off[r];e<off[r+1];e++){const v=adj[e];seen[v]=mark;if(kind[v]===IFACE)look(v);}
      if(best<0)for(let e=off[r];e<off[r+1];e++){const v=adj[e];for(let f2=off[v];f2<off[v+1];f2++){const u=adj[f2];if(seen[u]===mark)continue;seen[u]=mark;if(kind[u]===IFACE)look(u);}}
      if(best<0)continue;Fm[best/9]+=this.slabForce(r);if(!trench[best+5])this.sFed[r]=t;}
    this._seenN=mark;
    // slab pull on the subducting plate, suction on the overriding plate (from last step's velocity into the trench)
    let pullF=0,sucF=0;
    for(let q=0;q<trench.length;q+=9){const m=trench[q],Lp=trench[q+1],Wp=trench[q+2],Cs=trench[q+4],F=Fm[q/9];
      if(F===0)continue;const x=X[m],y=Y[m],z=Z[m],nx=trench[q+6],ny=trench[q+7],nz=trench[q+8];
      const ax=y*nz-z*ny,ay=z*nx-x*nz,az=x*ny-y*nx;pullF+=F;
      term("bpull");if(idx[Lp]>=0)addTq(idx[Lp],R*F,ax,ay,az);term("bsuction");
      const WL=plates[Lp].W,vn=R*(WL[0]*ax+WL[1]*ay+WL[2]*az)/MY,ex=F-Cs*this.iLen[m]*vn;
      if(ex>0&&idx[Wp]>=0){const Fs=PH.cSuc*ex;sucF+=Fs;addTq(idx[Wp],-R*Fs,ax,ay,az);}}
    if(bud){K=new Float64Array(N3*N3);b=new Float64Array(N3);for(const nm in bud){const src=bud[nm];if(nm[0]==="K")for(let i=0;i<N3*N3;i++)K[i]+=src[i];else for(let i=0;i<N3;i++)b[i]+=src[i];}}
    // Interface coupling with a frictional limit (§4.8): traction c dv up to S, i.e. a secant coefficient
    // min(c, S/|dv|), found by Picard iteration from the previous step's velocities (as L1 does for its power law)
    const Kb=K,nc=cpl.length/8,ce=new Float64Array(nc),wv=p=>plates[p].W;
    const slip=(oi,oj,q)=>{const x=cpl[q+5],y=cpl[q+6],z=cpl[q+7],dx=oi[0]-oj[0],dy=oi[1]-oj[1],dz=oi[2]-oj[2],d=dx*x+dy*y+dz*z;return R*Math.sqrt(Math.max(0,dx*dx+dy*dy+dz*dz-d*d));};
    const pid=new Int32Array(n);for(let p=0;p<np;p++)if(idx[p]>=0)pid[idx[p]]=p;
    for(let c=0;c<nc;c++){const q=8*c,dv=slip(wv(pid[cpl[q]]).map(v=>v/MY),wv(pid[cpl[q+1]]).map(v=>v/MY),q);ce[c]=Math.min(cpl[q+2],dv>0?cpl[q+3]/dv:Infinity);}
    const addCoupling=(Kt,scale)=>{K=Kt;for(let c=0;c<nc;c++){const q=8*c,i=cpl[q],j=cpl[q+1],v=ce[c]*cpl[q+4]*R2*scale,x=cpl[q+5],y=cpl[q+6],z=cpl[q+7];
      addProj(i,i,v,x,y,z);addProj(j,j,v,x,y,z);addProj(i,j,-v,x,y,z);addProj(j,i,-v,x,y,z);}};
    let om=null,K0=null,iters=0;
    for(let it=0;it<8;it++){iters++;const Kt=Kb.slice();addCoupling(Kt,1);
      let dmax=0;for(let i=0;i<N3;i++)dmax=Math.max(dmax,Kt[i*N3+i]);for(let i=0;i<N3;i++)Kt[i*N3+i]+=1e-12*dmax;
      K0=this.checkSolve?Kt.slice():null;om=cholSolve(Kt,b,N3);
      let change=0;for(let c=0;c<nc;c++){const q=8*c,dv=slip(om.subarray(3*cpl[q],3*cpl[q]+3),om.subarray(3*cpl[q+1],3*cpl[q+1]+3),q),v=Math.min(cpl[q+2],dv>0?cpl[q+3]/dv:Infinity);
        change=Math.max(change,Math.abs(v-ce[c])/ce[c]);ce[c]=v;}
      if(change<0.02)break;}
    this.couplingIters=iters;
    if(bud){bud.Kcouple=new Float64Array(N3*N3);addCoupling(bud.Kcouple,1);}
    // invariant I7 (tests): the torques balance, |K omega - b| / |b|
    if(K0){let r2=0,b2=0;for(let i=0;i<N3;i++){let s=-b[i];for(let j=0;j<N3;j++)s+=K0[i*N3+j]*om[j];r2+=s*s;b2+=b[i]*b[i];}this.solveResidual=b2>0?Math.sqrt(r2/b2):0;}
    for(let p=0;p<np;p++){const i=idx[p];if(i<0)continue;const W=plates[p].W;W[0]=om[3*i]*MY;W[1]=om[3*i+1]*MY;W[2]=om[3*i+2]*MY;}
    // diagnostics: per converging plate pair, boundary length (km), mean convergence rate (cm/yr, length-weighted) and
    // the length where both sides are continental (a collision); and the forces in play (N)
    const pairs=new Map();
    for(let q=0;q<trench.length;q+=9){const m=trench[q],Lp=trench[q+1],Wp=trench[q+2],x=X[m],y=Y[m],z=Z[m],nx=trench[q+6],ny=trench[q+7],nz=trench[q+8];
      const a=plates[Lp].W,c=plates[Wp].W,dx=a[0]-c[0],dy=a[1]-c[1],dz=a[2]-c[2],v=((dy*z-dz*y)*nx+(dz*x-dx*z)*ny+(dx*y-dy*x)*nz)*R_KM/10;
      const key=Math.min(Lp,Wp)*4096+Math.max(Lp,Wp),e=pairs.get(key)||[0,0,0],l=this.iLen[m]/1e3;e[0]+=l;e[1]+=v*l;if(trench[q+5])e[2]+=l;pairs.set(key,e);}
    for(const e of pairs.values())e[1]/=e[0];
    this.fstat={pairs,pullN:pullF,suctionN:sucF,gpeN:gpeF};
    // the budget: per plate, each term's torque along the plate's own rotation axis (N m; driving > 0, resisting < 0)
    if(bud){const out={};for(const nm in bud){const v=new Float64Array(n),src=bud[nm];
        for(let i=0;i<n;i++){const wl=Math.hypot(om[3*i],om[3*i+1],om[3*i+2])||1;let s2=0;
          for(let r=0;r<3;r++){let tq;if(nm[0]==="K"){tq=0;for(let j=0;j<N3;j++)tq-=src[(3*i+r)*N3+j]*om[j];}else tq=src[3*i+r];s2+=tq*om[3*i+r]/wl;}v[i]=s2;}
        out[nm.slice(1)]=v;}
      out.plate=[];for(let p=0;p<np;p++)if(idx[p]>=0)out.plate[idx[p]]=p;this.fstat.budget=out;}}
  walkVelocities(){if(!this.forces)super.walkVelocities();}
  // ---------- one step ----------
  step(){const dt=this.dt,kind=this.mKind;
    {const b=Math.floor(this.t/C.SLAB_DT),q=((b%C.SLAB_BINS)+C.SLAB_BINS)%C.SLAB_BINS;if(this.slabT[q]!==b*C.SLAB_DT){this.slabT[q]=b*C.SLAB_DT;this.slab.fill(0,q*this.N,(q+1)*this.N);}this.slabQ=q;}
    for(let i=0;i<this.N;i++)this.trench[i]*=0.85;
    this.classifyInterfaces();if(this.forces){this.slabUpkeep();this.solveVelocities();}this.move();this.fillGaps();this.maintainInterfaces();if(this.tri.broken)this.rebuildMesh();this.buildAdjacency();this.assignOceanAreas();
    this.surfaceProcesses();this.recordFaults();
    this.t+=dt;this.steps++;
    this.computeMarkerElev();if(this.deferView)this.viewStale=true;else this.computeView();
    this.updateStats(this.steps%4===0);this.hotspotEvents();this.plateEvents();
    if(this.steps%32===0)this.gc();}
  // elevOf from the lattice core, inlined: isostasy for continents; for ocean, depth with the square root of age,
  // plus thickness above normal crust and sediment
  computeMarkerElev(){const t=this.t,ns=this.nSlots,kind=this.mKind,K=this.mK,T=this.mT,B=this.mB,S=this.mS,E=this.mE,ISO=C.ISO,CB=C.CONT_BASE,OT=C.OCEAN_T;
    for(let m=0;m<ns;m++){if(kind[m]!==INTERIOR)continue;
      if(K[m]){E[m]=ISO*(T[m]-CB);continue;}const a=t-B[m],e=-2.6-0.33*Math.sqrt(a>0?a:0)+0.5*(T[m]-OT)+0.3*S[m];E[m]=e<-7?-7:e;}}
  // ---------- surface processes: the lattice core's rules, per marker ----------
  surfaceProcesses(){const dt=this.dt,t=this.t,pool=this.pool,ns=this.nSlots,kind=this.mKind,K=this.mK,T=this.mT,A=this.mA,S=this.mS,E=this.mE,pl=this.mPlate,mo=this.mo,ma=this.ma,mw=this.mw,r=this.rand;
    // arcs and collisions gathered during motion
    const got=new Uint8Array(ns);
    for(let m=0;m<ns;m++){if(kind[m]!==INTERIOR)continue;
      if(this.arcAdd[m]>0){got[m]=1;const a=this.arcAdd[m]/A[m],before=K[m]?A[m]*T[m]:0;T[m]+=a;if(K[m]?a>0.25:T[m]>C.OCEAN_T+4)pool.push(this.mCol[m],ROCK.ARC,t,a,pl[m]);if(K[m])this.led("contArc",A[m]*T[m]-before);this.arcAdd[m]=0;}
      if(this.collAdd[m]>0){const c=this.collAdd[m]/A[m];const wasK=K[m],before=wasK?A[m]*T[m]:0;T[m]+=c;if(!K[m]&&T[m]>C.ARC_CONV)K[m]=1;
        if(c>0.4)pool.push(this.mCol[m],ROCK.OROGEN,t,c*0.4,pl[m]);
        if(K[m])this.led(wasK?"contCollisionIn":"contAccretion",A[m]*T[m]-before);this.collAdd[m]=0;}}
    // hotspots: the marker over each, and its neighbours (weights scaled so the total matches the lattice's four)
    for(const h of this.hotspots){const p=this.cube.pos,c=this.nearestAny(p[3*h.i],p[3*h.i+1],p[3*h.i+2]);if(c<0)continue;
      const add=(m,dz)=>{const before=K[m]?A[m]*T[m]:0;T[m]+=dz;pool.push(this.mCol[m],ROCK.HOTSPOT,t,dz,pl[m]);if(K[m])this.led("contHotspot",A[m]*T[m]-before);};
      if(K[c])add(c,0.1*dt);else{add(c,C.HOT*dt);const deg=mo[c+1]-mo[c]||1;for(let e=mo[c];e<mo[c+1];e++){const j=ma[e];if(!K[j])add(j,0.3*C.HOT*dt*4/deg);}}}
    this.computeMarkerElev();
    const fb=Math.pow((this.contFrac||this.contTarget)/this.contTarget,6),sedConv=C.SED_CONV*fb,collapseP=0.05*dt/fb,kd=C.KD*dt,hk=Math.sqrt(this.cellA);
    const dV=this.dT;dV.fill(0,0,ns);
    for(let m=0;m<ns;m++){if(kind[m]!==INTERIOR)continue;
      if(K[m]){
        let drop=0;
        for(let e=mo[m];e<mo[m+1];e++){const j=ma[e];if(K[j]){const fl=kd*hk*hk*mw[e]*(T[j]-T[m]);dV[m]+=fl;dV[j]-=fl;}if(E[j]<E[m])drop+=E[m]-E[j];}
        if(E[m]>0&&drop>0){const loss=C.KE*E[m]*dt,V=loss*A[m];dV[m]-=V;pool.erode(this.mCol[m],loss);
          for(let e=mo[m];e<mo[m+1];e++){const j=ma[e];if(E[j]>=E[m])continue;const sh=V*(E[m]-E[j])/drop;
            if(K[j]){dV[j]+=sh;pool.push(this.mCol[j],ROCK.CLASTIC,t,sh/A[j],pl[j]);}else{S[j]+=0.5*sh/A[j];this.led("contErodedToOcean",-sh);}}}
      } else if(T[m]>C.OCEAN_T){const d=(T[m]-C.OCEAN_T)*C.VOLC_DECAY*dt;dV[m]-=d*A[m];pool.erode(this.mCol[m],d);}}
    for(let m=0;m<ns;m++){if(kind[m]!==INTERIOR)continue;let v=T[m]+dV[m]/A[m];
      if(K[m]){const raw=v;if(v>C.CAP)v=C.CAP;if(v<10)v=10;if(v!==raw)this.led("contClamp",(v-raw)*A[m]);}
      else{
        if(S[m]>sedConv){const add=Math.max(0.5,24-v);K[m]=1;v+=add;pool.push(this.mCol[m],ROCK.SHELF,t,add,pl[m]);S[m]=0;this.led("contFromShelf",v*A[m]);}
        else if(v>C.ARC_CONV*fb&&got[m]){K[m]=1;this.led("contFromArc",v*A[m]);}
        if(!K[m]){if(v<5)v=5;if(v>30)v=30;}}
      T[m]=v;}
    for(let m=0;m<ns;m++){if(kind[m]!==INTERIOR||!K[m]||T[m]<C.COLLAPSE||r()>collapseP)continue;const deg=mo[m+1]-mo[m];if(!deg)continue;
      const c=ma[mo[m]+Math.floor(r()*deg)];
      if(!K[c]&&pl[c]===pl[m]){const base=pool.R[this.mCol[m]*SIM.U];K[c]=1;T[c]=28;S[c]=0;this.mCol[c]=pool.make(base,t,26,pl[m]);pool.push(this.mCol[c],ROCK.CLASTIC,t,2,pl[m]);
        T[m]-=8;this.led("contCollapse",28*A[c]-8*A[m]);}}
    if(fb>1){const pz=0.1*dt*(fb-1);
      for(let m=0;m<ns;m++){if(kind[m]!==INTERIOR||!K[m]||E[m]>0.3||r()>pz)continue;let coast=false;for(let e=mo[m];e<mo[m+1];e++)if(!K[ma[e]])coast=true;
        if(coast){this.led("contToOcean",-T[m]*A[m]);K[m]=0;T[m]=C.OCEAN_T+1;this.mB[m]=t-60;S[m]=2;this.mCol[m]=pool.make(ROCK.OCEAN,t-60,T[m],pl[m]);}}}
    for(let m=0;m<ns;m++)if(kind[m]===INTERIOR)pool.reconcile(this.mCol[m],T[m]);
  }
  nearestAny(px,py,pz){const T=this.tri,t=T.locate(px,py,pz,this._hintAny||-1);if(t<0)return -1;this._hintAny=t;let best=-1,bd=Infinity;
    for(let k=0;k<3;k++){const v=T.tv[3*t+k];if(this.mKind[v]!==INTERIOR)continue;const d=(T.X[v]-px)**2+(T.Y[v]-py)**2+(T.Z[v]-pz)**2;if(d<bd){bd=d;best=v;}}
    if(best<0){const v=T.tv[3*t];best=this.nearestOf(px,py,pz,this.iN[v]?this.iPlate[3*v]:0,v);}return best;}
  // boundary kinematics into each marker's fault record (the lattice core's recordFaults, over material neighbours)
  recordFaults(){const T=this.tri,X=T.X,Y=T.Y,Z=T.Z,plates=this.plates,pool=this.pool,fpool=this.fpool,dt=this.dt,t=this.t,mo=this.mo,ma=this.ma,pl=this.mPlate;
    for(let m=0;m<this.nSlots;m++){if(this.mKind[m]!==INTERIOR)continue;const P=plates[pl[m]],x=X[m],y=Y[m],z=Z[m];
      let cnt=0,ac=0,al=0,nx=0,ny=0,nz=0;
      for(let e=mo[m];e<mo[m+1];e++){const j=ma[e],q=pl[j];if(q===pl[m])continue;const Qp=plates[q];
        let ex=X[j]-x,ey=Y[j]-y,ez=Z[j]-z;const d=ex*x+ey*y+ez*z;ex-=d*x;ey-=d*y;ez-=d*z;const el=Math.hypot(ex,ey,ez)||1;ex/=el;ey/=el;ez/=el;
        const dx=Qp.W[0]-P.W[0],dy=Qp.W[1]-P.W[1],dz=Qp.W[2]-P.W[2];const vx=dy*z-dz*y,vy=dz*x-dx*z,vz=dx*y-dy*x;
        const tx=y*ez-z*ey,ty=z*ex-x*ez,tz=x*ey-y*ex;ac+=-(vx*ex+vy*ey+vz*ez);al+=vx*tx+vy*ty+vz*tz;nx+=ex;ny+=ey;nz+=ez;cnt++;}
      const bn=this.mBn;if(!cnt){bn[3*m]=bn[3*m+1]=bn[3*m+2]=0;continue;}
      const nl=Math.hypot(nx,ny,nz);if(nl<1e-6){bn[3*m]=bn[3*m+1]=bn[3*m+2]=0;continue;}
      nx/=nl;ny/=nl;nz/=nl;bn[3*m]=nx;bn[3*m+1]=ny;bn[3*m+2]=nz;ac=ac/cnt*R_KM;al=al/cnt*R_KM;
      const c=this.mCol[m];let f=pool.flt[c];if(f<0){f=fpool.alloc();pool.flt[c]=f;}
      const D=fpool.D,o=f*FS,e=o+8,ki=kindOf(ac,al);
      if(isNaN(D[e])){D[e]=t;D[e+1]=t;D[e+2]=0;D[e+3]=0;D[e+4]=0;D[e+5]=0;D[e+6]=0;}
      else{let shift=t-D[e+1]>C.FAULT_GAP;
        if(!shift){if(kindOf(D[e+2],D[e+3])!==ki&&Math.hypot(D[e+2],D[e+3])>5){D[o]+=dt;if(D[o]>C.FAULT_SWITCH)shift=true;}else D[o]=0;}
        if(shift){for(let k=0;k<7;k++)D[o+1+k]=D[e+k];D[e]=t;D[e+1]=t;D[e+2]=0;D[e+3]=0;D[e+4]=0;D[e+5]=0;D[e+6]=0;D[o]=0;}}
      const sx=y*nz-z*ny,sy=z*nx-x*nz,sz=x*ny-y*nx,Q=P.Q;
      let bx=Q[0]*sx+Q[3]*sy+Q[6]*sz,by=Q[1]*sx+Q[4]*sy+Q[7]*sz,bz=Q[2]*sx+Q[5]*sy+Q[8]*sz;
      if(bx*D[e+4]+by*D[e+5]+bz*D[e+6]<0){bx=-bx;by=-by;bz=-bz;}
      const mm=Math.hypot(ac,al)*dt;D[e+1]=t;D[e+2]+=ac*dt;D[e+3]+=al*dt;D[e+4]+=bx*mm;D[e+5]+=by*mm;D[e+6]+=bz*mm;}}
  // ---------- the fixed-grid view the page reads ----------
  // Each fixed cell shows the nearest interior marker of the triangle holding its centre (or of its corners' rings):
  // plate, crust, rock column, elevation. A view only: nothing is written back.
  // The marker a fixed cell shows: greedy descent on the Delaunay graph from any start reaches the vertex nearest the
  // cell's centre; the cell shows that vertex if it is an interior marker, else the nearest interior marker within two
  // rings of it. Only the graph is used (no triangle walk), so a helper worker holding a copy can do it too.
  viewMarker(px,py,pz,v){const X=this.tri.X,Y=this.tri.Y,Z=this.tri.Z,kind=this.mKind,off=this.off,adj=this.adj;
    if(v<0||!kind[v]){v=this._anyV>=0&&kind[this._anyV]?this._anyV:-1;if(v<0){for(let m=0;m<kind.length;m++)if(kind[m]){v=m;break;}this._anyV=v;}}
    let d=(X[v]-px)**2+(Y[v]-py)**2+(Z[v]-pz)**2;
    for(let moved=true;moved;){moved=false;let nv=v;const e1=off[v+1];for(let e=off[v];e<e1;e++){const q=adj[e],dq=(X[q]-px)**2+(Y[q]-py)**2+(Z[q]-pz)**2;if(dq<d){d=dq;nv=q;moved=true;}}v=nv;}
    this._near=v;if(kind[v]===INTERIOR)return v;
    let best=-1,bd=Infinity;
    for(let e=off[v];e<off[v+1];e++){const q=adj[e];if(kind[q]===INTERIOR){const dq=(X[q]-px)**2+(Y[q]-py)**2+(Z[q]-pz)**2;if(dq<bd){bd=dq;best=q;}}}
    if(best<0)for(let e=off[v];e<off[v+1];e++){const r=adj[e];for(let f=off[r];f<off[r+1];f++){const q=adj[f];if(kind[q]!==INTERIOR)continue;const dq=(X[q]-px)**2+(Y[q]-py)**2+(Z[q]-pz)**2;if(dq<bd){bd=dq;best=q;}}}
    return best;}
  markerForCell(i){const p=this.cube.pos;return this.viewMarker(p[3*i],p[3*i+1],p[3*i+2],-1);}
  computeView(){const pos=this.cube.pos,N=this.N,kind=this.mKind,X=this.tri.X,Y=this.tri.Y,Z=this.tri.Z,off=this.off,adj=this.adj;
    if(!this._vcache||this._vcache.length!==N){this._vcache=new Int32Array(N).fill(-1);this._vcid=new Uint32Array(N);}const vc=this._vcache,vid=this._vcid,mid=this.mId;let prev=-1;
    for(let i=0;i<N;i++){const px=pos[3*i],py=pos[3*i+1],pz=pos[3*i+2];let v=vc[i],best;
      // start from last step's vertex if it is still the same marker (slots are reused far away), else from the
      // previous cell's answer (cells are visited in spatial order); then viewMarker's rule, inlined for speed
      if(v<0||!kind[v]||mid[v]!==vid[i])v=prev>=0&&kind[prev]?prev:-1;
      if(v<0){best=this.viewMarker(px,py,pz,-1);v=this._near;}
      else{let d=(X[v]-px)**2+(Y[v]-py)**2+(Z[v]-pz)**2;
        for(let moved=true;moved;){moved=false;let nv=v;const e1=off[v+1];for(let e=off[v];e<e1;e++){const q=adj[e],dq=(X[q]-px)**2+(Y[q]-py)**2+(Z[q]-pz)**2;if(dq<d){d=dq;nv=q;moved=true;}}v=nv;}
        best=kind[v]===INTERIOR?v:this.viewMarker(px,py,pz,v);}
      vc[i]=v;vid[i]=mid[v];prev=v;
      this.owner[i]=this.mPlate[best];this.bodyIdx[i]=best;this.wT[i]=this.mT[best];this.wB[i]=this.mB[best];this.wS[i]=this.mS[best];this.wK[i]=this.mK[best];
      this.wCol[i]=this.mCol[best];this.elev[i]=this.mE[best];}
    this.viewStale=false;this.surface();}
  ensureView(){if(this.viewStale)this.computeView();}
  // the island event reads the markers (the view may be built elsewhere): the marker the hotspot's cell shows
  hotspotEvents(){if(!this.mKind)return super.hotspotEvents();
    for(const h of this.hotspots){const m=this.markerForCell(h.i),up=this.mE[m]>0;
      if(up&&!h.wasUp&&!this.mK[m]&&this.t-h.last>15){this.log("island","A volcanic island emerges over a hotspot",h.v);h.last=this.t;}h.wasUp=up;}}
  // ---------- statistics, plates, events ----------
  updateStats(centroids){if(!this.mKind)return super.updateStats(centroids);
    const {plates,pool,t}=this,X=this.tri.X,Y=this.tri.Y,Z=this.tri.Z,ca=this.cellA;
    for(const P of plates){P.area=0;P.contArea=0;P._c=[0,0,0];P.nM=0;P._v2=0;}
    let land=0,cont=0,tot=0,oa=0,on=0,units=0,nc=0;
    for(let m=0;m<this.nSlots;m++){if(this.mKind[m]!==INTERIOR)continue;const P=plates[this.mPlate[m]],a=this.mA[m];P.area+=a/ca;P.nM++;tot+=a;
      if(this.mK[m]){P.contArea+=a/ca;cont+=a;units+=pool.cnt[this.mCol[m]];nc++;}else{oa+=(t-this.mB[m])*a;on+=a;}if(this.mE[m]>0)land+=a;
      if(centroids){P._c[0]+=X[m]*a;P._c[1]+=Y[m]*a;P._c[2]+=Z[m]*a;const W=P.W,d=W[0]*X[m]+W[1]*Y[m]+W[2]*Z[m];P._v2+=a*(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]-d*d);}}
    for(const P of plates){if(!P.alive)continue;if(P.area>P.maxArea)P.maxArea=P.area;if(centroids&&P.area>0)P.c=norm3(P._c);if(P.nM===0)P.area=0;}
    this.landFrac=land/tot;this.contFrac=cont/tot;this.oceanAge=on?oa/on:0;this.meanUnits=nc?units/nc:0;this.markerArea=tot;
    // speeds: the root mean square over the crust's area of |omega x r| (cm/yr), per plate and overall
    if(centroids){this.netRot=Math.hypot(...this.netRotation())*180/Math.PI;let ss=0;
      for(const P of plates)if(P.alive){ss+=P._v2;P.vrms=P.area>0?Math.sqrt(P._v2/(P.area*ca))*R_KM/10:0;}this.rmsSpeed=Math.sqrt(ss/(tot||1))*R_KM/10;}}
  // net rotation over the markers' carried areas (the lattice core integrates over fixed cells)
  netRotation(){if(!this.mKind)return super.netRotation();const X=this.tri.X,Y=this.tri.Y,Z=this.tri.Z,M=this.cube.Minv,s=1/(R_KM*R_KM);let sx=0,sy=0,sz=0;
    for(let m=0;m<this.nSlots;m++){if(this.mKind[m]!==INTERIOR)continue;const W=this.plates[this.mPlate[m]].W,x=X[m],y=Y[m],z=Z[m],d=W[0]*x+W[1]*y+W[2]*z,a=this.mA[m]*s;
      sx+=a*(W[0]-d*x);sy+=a*(W[1]-d*y);sz+=a*(W[2]-d*z);}
    return[M[0]*sx+M[1]*sy+M[2]*sz,M[3]*sx+M[4]*sy+M[5]*sz,M[6]*sx+M[7]*sy+M[8]*sz];}
  // the lattice core's split, on markers: a noisy great circle through the plate's centroid
  split(P){const {rand:r,noise:nz}=this,X=this.tri.X,Y=this.tri.Y,Z=this.tri.Z;
    const c=P.c;let n=randUnit(r);const d=n[0]*c[0]+n[1]*c[1]+n[2]*c[2];n=norm3([n[0]-d*c[0],n[1]-d*c[1],n[2]-d*c[2]]);
    const a=norm3(cross(c,n)),Qp=this.addPlate(P.W);Qp.Q.set(P.Q);Qp.parent=P.id;const moved=[];let contA=0,contB=0,rest=0;
    for(let m=0;m<this.nSlots;m++){if(this.mKind[m]!==INTERIOR||this.mPlate[m]!==P.id)continue;const x=X[m],y=Y[m],z=Z[m];
      const s=x*n[0]+y*n[1]+z*n[2]+(nz(x+40,y,z,2.2,4)-0.5)*0.5;if(s>0){moved.push(m);if(this.mK[m])contB++;}else{rest++;if(this.mK[m])contA++;}}
    if(moved.length<60||rest<60){Qp.alive=false;Qp.area=0;return;}
    for(const m of moved)this.mPlate[m]=Qp.id;
    const s=(1.5+2.5*r())*10/R_KM;for(let k=0;k<3;k++){Qp.W[k]=P.W[k]+a[k]*s/2;P.W[k]-=a[k]*s/2;}
    this.log("rift",`${contA>40&&contB>40?"Continental rift":"Rift"}: plate ${P.id} splits and plate ${Qp.id} drifts away`,c);
    this.maintainInterfaces();this.buildAdjacency();this.updateStats(true);}
  // the lattice core's weld: the smaller plate's markers join the larger; interface markers between them dissolve
  merge(a,b){const purge=()=>{for(const k of [...this.pairs.keys()]){const x=Math.floor(k/4096),y=k%4096;if(x===a.id||y===a.id||x===b.id||y===b.id)this.pairs.delete(k);}};
    if(!a.alive||!b.alive){purge();return;}
    const [Wn,L]=a.area>=b.area?[a,b]:[b,a];
    for(let m=0;m<this.nSlots;m++)if(this.mKind[m]===INTERIOR&&this.mPlate[m]===L.id)this.mPlate[m]=Wn.id;
    const tot=Wn.area+L.area;for(let k=0;k<3;k++)Wn.W[k]=(Wn.W[k]*Wn.area+L.W[k]*L.area)/tot;Wn.area=tot;Wn.contArea+=L.contArea;
    this.log("suture",`Plates ${Wn.id} and ${L.id} weld into one along a suture`,Wn.c);
    L.mergedInto=Wn.id;this.kill(L);purge();this.maintainInterfaces();this.buildAdjacency();}
  kill(P){P.alive=false;P.area=0;}
  gc(){const {pool,fpool}=this;pool.mark.fill(0);fpool.mark.fill(0);
    for(let m=0;m<this.nSlots;m++)if(this.mKind[m]===INTERIOR&&this.mCol[m]>=0)pool.mark[this.mCol[m]]=1;
    for(let c=0;c<pool.cap;c++){if(!pool.used[c])continue;if(!pool.mark[c]){const f=pool.flt[c];pool.release(c);if(f>=0){fpool.release(f);pool.flt[c]=-1;}}else if(pool.flt[c]>=0)fpool.mark[pool.flt[c]]=1;}
    for(let f=0;f<fpool.cap;f++)if(fpool.used[f]&&!fpool.mark[f])fpool.release(f);}
  restore(snap){super.restore(snap);Object.setPrototypeOf(this.tri,SphereTri.prototype);}
  // marker-level checks used by the tests: area closure, mesh validity, identity
  audit(){const T=this.tri,kind=this.mKind;let nv=0,A=0,ids=new Set(),dupIds=0,cross=0,ifaceA=0;
    for(let m=0;m<this.nSlots;m++){if(!kind[m])continue;nv++;if(ids.has(this.mId[m]))dupIds++;ids.add(this.mId[m]);if(kind[m]===INTERIOR)A+=this.mA[m];}
    const dual=new Float64Array(this.nSlots);const tv=T.tv;
    for(let t=0;t<T.nT;t++){if(!T.tAlive[t])continue;const a=T.area(t)/3,o=3*t;dual[tv[o]]+=a;dual[tv[o+1]]+=a;dual[tv[o+2]]+=a;
      for(let k=0;k<3;k++){const p=tv[o+(k+1)%3],q=tv[o+(k+2)%3];if(kind[p]===INTERIOR&&kind[q]===INTERIOR&&this.mPlate[p]!==this.mPlate[q])cross++;}}
    for(let m=0;m<this.nSlots;m++)if(kind[m]===IFACE)ifaceA+=dual[m]*R_KM*R_KM;
    const sphere=4*Math.PI*R_KM*R_KM,chk=T.check(nv,false);
    return{markers:nv,closure:(A-sphere)/sphere,interiorArea:A/sphere,ifaceArea:ifaceA/sphere,dupIds,crossEdges:cross/2,mesh:chk};}
}

// The worker side of the marker core. The page builds a Web Worker from SIMmodule's and MESHmodule's sources plus this
// function's, and talks to it by messages: it owns a MeshWorld, steps it (yielding between steps so a pause or a
// question gets through), and after each step posts the fixed-grid view (fresh copies, transferred, not shared) with
// the plates, statistics and new events; on request it posts the inspector's column, the handoff, or the export.
// With a view helper (init's `helper`), the step worker does not build the view: it hands a copy of the markers to a
// second worker and starts the next step at once; the view comes back and goes to the page with that step's summary.
// While the helper is busy, steps are not queued: when it is free, it gets the newest state.
function workerMain(){
  let w=null,running=false,until=0,handoff=null,markers=false,lastEv=0,ms=0,helper=null,busy=false,snap=null,evBuf=[];
  const ch=new MessageChannel();ch.port1.onmessage=()=>tick();const later=()=>ch.port2.postMessage(0);
  const VIEW=["owner","wK","wT","wB","wS","wR","wRA","elev","bodyIdx"];
  const summary=()=>({t:w.t,steps:w.steps,landFrac:w.landFrac,contFrac:w.contFrac,oceanAge:w.oceanAge,meanUnits:w.meanUnits,netRot:w.netRot,rmsSpeed:w.rmsSpeed,
    msPerStep:ms,plates:w.plates.map(P=>({id:P.id,alive:P.alive,hue:P.hue,W:P.W.slice(),c:P.c.slice(),area:P.area,contArea:P.contArea,Q:Array.from(P.Q),parent:P.parent,mergedInto:P.mergedInto})),
    nEvents:w.nEvents,trench:w.trench.slice()});
  const takeEvents=()=>{for(const e of w.events)if(e.n>lastEv)evBuf.push(e);lastEv=w.nEvents;const out=evBuf;evBuf=[];return out;};
  const markerCopy=()=>{const n=w.nSlots;return{X:w.tri.X.slice(0,n),Y:w.tri.Y.slice(0,n),Z:w.tri.Z.slice(0,n),K:w.mKind.slice(0,n),I:w.iKind.slice(0,n)};};
  const buffers=o=>{const b=[];for(const k in o)if(o[k]&&o[k].buffer)b.push(o[k].buffer);return b;};
  // without a helper: the view built here
  function stateHere(type){w.ensureView();const out={type,running,...summary(),events:takeEvents(),view:{}};for(const k of VIEW)out.view[k]=w[k].slice();
    const tr=buffers(out.view);tr.push(out.trench.buffer);if(markers){out.markers=markerCopy();tr.push(...buffers(out.markers));}self.postMessage(out,tr);}
  // with a helper: a copy of what the view needs goes out; the summary waits for it
  function sendSnap(){const n=w.nSlots,pool=w.pool,M=markerCopy();
    const m={steps:w.steps,n:w.n,N:w.N,X:M.X,Y:M.Y,Z:M.Z,kind:M.K,I:M.I,mId:w.mId.slice(0,n),plate:w.mPlate.slice(0,n),T:w.mT.slice(0,n),B:w.mB.slice(0,n),S:w.mS.slice(0,n),
      K:w.mK.slice(0,n),E:w.mE.slice(0,n),col:w.mCol.slice(0,n),off:w.off.slice(0,n+1),adj:w.adj.slice(0,w.off[n]),cnt:pool.cnt.slice(),pT:pool.T.slice(),pR:pool.R.slice(),pB:pool.B.slice(),markers};
    snap={steps:w.steps,summary:summary()};busy=true;helper.postMessage(m,buffers(m));}
  function onView(v){busy=false;const s=snap,last=!running&&w.steps===s.steps;
    const out={type:last?"stopped":"state",running:!last,...s.summary,events:takeEvents(),view:v.view},tr=buffers(out.view);tr.push(out.trench.buffer);
    if(v.markers){out.markers=v.markers;tr.push(...buffers(v.markers));}self.postMessage(out,tr);
    if(w.steps!==s.steps)sendSnap();}
  function report(type){if(!helper)return stateHere(type);if(!busy)sendSnap();} // a busy helper sends the newest state when it is done
  function tick(){if(!running)return;const t0=performance.now();w.step();ms=performance.now()-t0;if(w.t>=until-1e-9)running=false;report(running?"state":"stopped");if(running)later();}
  function startHelper(){try{const src=`${SIMmodule}\nconst SIM=SIMmodule();\n${MESHmodule}\nconst MESH=MESHmodule(SIM);\n(${MESH.viewHelperMain})();`;
      const h=new Worker(URL.createObjectURL(new Blob([src],{type:"text/javascript"})));
      h.onmessage=e=>onView(e.data);h.onerror=e=>{if(e.preventDefault)e.preventDefault();helper=null;busy=false;w.deferView=false;w.viewStale=true;report(running?"state":"stopped");};
      helper=h;w.deferView=true;w.viewStale=true;}catch(e){helper=null;}}
  self.onmessage=e=>{const m=e.data;try{
    if(m.type==="init"){if(helper){helper.terminate();helper=null;}w=new MESH.MeshWorld(m.opts);lastEv=0;evBuf=[];handoff=null;running=false;busy=false;
      if(m.helper&&typeof Worker!=="undefined")startHelper();
      self.postMessage({type:"ready",helper:!!helper,hotspots:w.hotspots.map(h=>({v:h.v,i:h.i})),dt:w.dt,cellKm:w.cellKm,contTarget:w.contTarget});report("stopped");}
    else if(m.type==="run"){until=m.until;if(!running&&w.t<until-1e-9){running=true;later();}}
    else if(m.type==="pause"){if(running){running=false;report("stopped");}}
    else if(m.type==="step"){running=false;const t0=performance.now();w.step();ms=performance.now()-t0;report("stopped");}
    else if(m.type==="markers"){markers=!!m.on;if(!helper||!busy)report(running?"state":"stopped");}
    else if(m.type==="column"){w.ensureView();self.postMessage({type:"column",i:m.i,steps:w.steps,...w.columnInfo(m.i)});}
    else if(m.type==="extract"){w.ensureView();handoff=w.extractFaults();handoff.l1=w.buildL1();self.postMessage({type:"extract",steps:w.steps,handoff});}
    else if(m.type==="buildL1"){w.ensureView();self.postMessage({type:"buildL1",steps:w.steps,l1:w.buildL1()});}
    else if(m.type==="export"){w.ensureView();self.postMessage({type:"export",json:handoff?JSON.stringify(w.exportHandoff(handoff)):null});}
    else if(m.type==="hash"){w.ensureView();self.postMessage({type:"hash",steps:w.steps,t:w.t,audit:w.audit()});}
  }catch(err){self.postMessage({type:"error",message:String(err&&err.stack||err)});}};
}
// The view helper: given a copy of the markers, their adjacency and the rock columns, it builds the fixed-grid view
// with MeshWorld's own computeView (on an object holding just those arrays) and sends it back.
function viewHelperMain(){
  const o=Object.create(MESH.MeshWorld.prototype);
  self.onmessage=e=>{const m=e.data;
    if(!o.cube||o.N!==m.N){o.cube=SIM.getCube(m.n);o.N=m.N;o._vcache=null;}
    Object.assign(o,{tri:{X:m.X,Y:m.Y,Z:m.Z},mKind:m.kind,mId:m.mId,off:m.off,adj:m.adj,mPlate:m.plate,mT:m.T,mB:m.B,mS:m.S,mK:m.K,mCol:m.col,mE:m.E,
      pool:{cnt:m.cnt,T:m.pT,R:m.pR,B:m.pB}});
    const N=m.N;Object.assign(o,{owner:new Int16Array(N),bodyIdx:new Int32Array(N),wT:new Float32Array(N),wB:new Float32Array(N),wS:new Float32Array(N),wK:new Uint8Array(N),
      wCol:new Int32Array(N),elev:new Float32Array(N),wR:new Uint8Array(N),wRA:new Float32Array(N)});
    o.computeView();
    const view={owner:o.owner,wK:o.wK,wT:o.wT,wB:o.wB,wS:o.wS,wR:o.wR,wRA:o.wRA,elev:o.elev,bodyIdx:o.bodyIdx},tr=[];for(const k in view)tr.push(view[k].buffer);
    const out={steps:m.steps,view};if(m.markers){out.markers={X:m.X,Y:m.Y,Z:m.Z,K:m.kind,I:m.I};for(const k of ["X","Y","Z","kind","I"])tr.push(m[k].buffer);}
    self.postMessage(out,tr);};
}

return {SphereTri,MeshWorld,workerMain,viewHelperMain,PH,gpe,hth,buoy};
}
const MESH=MESHmodule(SIM);
// ===== MESH-END =====
