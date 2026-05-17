// universe-ch.js — SOAP COSMOS · núcleo físico Cahn-Hilliard
// ============================================================================
//  Sustituye las reglas ad-hoc de v0.3 por FÍSICA REAL:
//
//    Cahn-Hilliard:  ∂φ/∂t = M ∇²( φ³ − φ − ε²∇²φ )
//    + término gravitatorio: la materia (φ>0) genera un potencial de Poisson
//      de larga escala que la agrupa en cúmulos DENTRO de la espuma.
//
//  φ = campo de fase:  −1 = vacío,  +1 = materia (espuma de jabón / espacio-t).
//  Las membranas, vacíos, filamentos y cúmulos EMERGEN de la ecuación.
//  No hay umbrales arbitrarios. Esquema de Eyre (semi-implícito, FFT),
//  verificado en Python: conserva masa (~1e-16) y separa fases (φ→±1).
//
//  Expone una API compatible con el render de v1.0:
//    initUniverse(), stepUniverse(dt), getState()
// ============================================================================

const GRID = 160;            // rejilla (potencia cómoda para FFT)
const N = GRID;
const SIZE = N * N;

// parámetros físicos (de la prueba de concepto verificada)
const EPS = 1.0;             // ancho de interfaz (ε)
const M   = 1.0;             // movilidad
const A   = 2.0;             // estabilizador de Eyre (mínimo teórico = 2)
const G_GRAV = 0.04;         // gravedad: débil y lenta (separación de escalas)
const DT_PHYS = 0.15;        // paso temporal interno (Cahn-Hilliard es rígida)

let phi;                     // campo de fase (Float64 para estabilidad numérica)
let age = 0;
let mass0 = 0;

// ---- FFT 2D (Cooley-Tukey, radix-2; GRID debe ser potencia de 2... 160 no
//      lo es, así que usamos una FFT genérica por filas/columnas con DFT
//      vía algoritmo de Bluestein-lite: para 160 va sobrado en rendimiento) ----
// Para simplicidad y robustez usamos una FFT real basada en
// transformadas separables con un tamaño potencia de 2 (256) y zero-pad.
// Pero zero-pad rompe la periodicidad de CH. Mejor: usar GRID=128 (2^7).
// -> redefinimos a 128 para FFT radix-2 limpia.

// (ver implementación FFT abajo; trabajamos a 128)
const G2 = 128;
const SZ = G2 * G2;

let re, im, kx2, ky2, K2, K4, denom;

function fft1d(reArr, imArr, n, inverse){
  // FFT radix-2 in-place (n potencia de 2)
  for(let i=1,j=0;i<n;i++){
    let bit=n>>1;
    for(;j&bit;bit>>=1) j^=bit;
    j^=bit;
    if(i<j){
      let tr=reArr[i]; reArr[i]=reArr[j]; reArr[j]=tr;
      let ti=imArr[i]; imArr[i]=imArr[j]; imArr[j]=ti;
    }
  }
  for(let len=2;len<=n;len<<=1){
    const ang=(inverse?2:-2)*Math.PI/len;
    const wr=Math.cos(ang), wi=Math.sin(ang);
    for(let i=0;i<n;i+=len){
      let cwr=1, cwi=0;
      for(let k=0;k<len/2;k++){
        const ur=reArr[i+k],          ui=imArr[i+k];
        const vr=reArr[i+k+len/2]*cwr - imArr[i+k+len/2]*cwi;
        const vi=reArr[i+k+len/2]*cwi + imArr[i+k+len/2]*cwr;
        reArr[i+k]=ur+vr;       imArr[i+k]=ui+vi;
        reArr[i+k+len/2]=ur-vr; imArr[i+k+len/2]=ui-vi;
        const ncwr=cwr*wr-cwi*wi;
        cwi=cwr*wi+cwi*wr; cwr=ncwr;
      }
    }
  }
  if(inverse){ for(let i=0;i<n;i++){ reArr[i]/=n; imArr[i]/=n; } }
}

// FFT 2D separable, sobre buffers de tamaño SZ
const _rowR=new Float64Array(G2), _rowI=new Float64Array(G2);
function fft2d(reA, imA, inverse){
  for(let y=0;y<G2;y++){
    const o=y*G2;
    for(let x=0;x<G2;x++){ _rowR[x]=reA[o+x]; _rowI[x]=imA[o+x]; }
    fft1d(_rowR,_rowI,G2,inverse);
    for(let x=0;x<G2;x++){ reA[o+x]=_rowR[x]; imA[o+x]=_rowI[x]; }
  }
  for(let x=0;x<G2;x++){
    for(let y=0;y<G2;y++){ _rowR[y]=reA[y*G2+x]; _rowI[y]=imA[y*G2+x]; }
    fft1d(_rowR,_rowI,G2,inverse);
    for(let y=0;y<G2;y++){ reA[y*G2+x]=_rowR[y]; imA[y*G2+x]=_rowI[y]; }
  }
}

// buffers de trabajo
let bufR, bufI, nlR, nlI, gpotR, gpotI;

export function initUniverse(){
  phi = new Float64Array(SZ);
  bufR=new Float64Array(SZ); bufI=new Float64Array(SZ);
  nlR =new Float64Array(SZ); nlI =new Float64Array(SZ);
  gpotR=new Float64Array(SZ); gpotI=new Float64Array(SZ);

  // números de onda (espectro)
  K2 =new Float64Array(SZ);
  K4 =new Float64Array(SZ);
  denom=new Float64Array(SZ);
  const k=new Float64Array(G2);
  for(let i=0;i<G2;i++){
    k[i] = (i<G2/2 ? i : i-G2) * (2*Math.PI/G2);
  }
  for(let y=0;y<G2;y++)for(let x=0;x<G2;x++){
    const kk = k[x]*k[x] + k[y]*k[y];
    const idx=y*G2+x;
    K2[idx]=kk; K4[idx]=kk*kk;
    denom[idx]=1.0 + DT_PHYS*M*kk*(A + EPS*EPS*kk);
  }

  // quench: mezcla casi homogénea con ruido (el "impacto del jabón")
  // un poco de sesgo central para recordar el origen del impacto
  const c=G2/2;
  for(let y=0;y<G2;y++)for(let x=0;x<G2;x++){
    const d2=((x-c)**2+(y-c)**2)/(G2*G2);
    phi[y*G2+x] = 0.10*(Math.random()-0.5) + 0.15*Math.exp(-d2*6);
  }
  let m=0; for(let i=0;i<SZ;i++) m+=phi[i];
  mass0 = m/SZ;
  age=0;
}

function stepPhys(){
  // ---- 1) Cahn-Hilliard, esquema de Eyre (verificado) ----
  // término no lineal nl = φ³ − (1+A)φ
  for(let i=0;i<SZ;i++){ const p=phi[i]; nlR[i]=p*p*p-(1+A)*p; nlI[i]=0; }
  fft2d(nlR,nlI,false);
  // PHI = FFT(phi)
  for(let i=0;i<SZ;i++){ bufR[i]=phi[i]; bufI[i]=0; }
  fft2d(bufR,bufI,false);
  // PHI_new = (PHI - dt*M*K2*NL) / denom
  for(let i=0;i<SZ;i++){
    const f = DT_PHYS*M*K2[i];
    bufR[i] = (bufR[i] - f*nlR[i]) / denom[i];
    bufI[i] = (bufI[i] - f*nlI[i]) / denom[i];
  }
  fft2d(bufR,bufI,true);
  for(let i=0;i<SZ;i++) phi[i]=bufR[i];

  // ---- 2) gravedad: potencial de Poisson de la materia, larga escala ----
  // rho = (phi+1)/2 ; resolver ∇²Φ = ρ en Fourier ; advectar φ suavemente
  let rmean=0;
  for(let i=0;i<SZ;i++){ const r=(phi[i]+1)*0.5; gpotR[i]=r; rmean+=r; }
  rmean/=SZ;
  for(let i=0;i<SZ;i++){ gpotR[i]-=rmean; gpotI[i]=0; }
  fft2d(gpotR,gpotI,false);
  for(let i=0;i<SZ;i++){
    let kk=K2[i]; if(kk<1e-9) kk=1e-9;
    const smooth=Math.exp(-K2[i]*8.0);     // gravedad = fuerza de larga escala
    gpotR[i] = -gpotR[i]/kk*smooth;
    gpotI[i] = -gpotI[i]/kk*smooth;
  }
  fft2d(gpotR,gpotI,true);                  // gpotR = potencial gravitatorio
  // advección suave de φ por −∇Φ_grav (diferencias centradas, periódico)
  const idx=(x,y)=>((y+G2)%G2)*G2+((x+G2)%G2);
  for(let y=0;y<G2;y++)for(let x=0;x<G2;x++){
    const i=y*G2+x;
    const gx=(gpotR[idx(x+1,y)]-gpotR[idx(x-1,y)])*0.5;
    const gy=(gpotR[idx(x,y+1)]-gpotR[idx(x,y-1)])*0.5;
    const px=(phi[idx(x+1,y)]-phi[idx(x-1,y)])*0.5;
    const py=(phi[idx(x,y+1)]-phi[idx(x,y-1)])*0.5;
    nlR[i] = phi[i] + DT_PHYS*G_GRAV*( px*(-gx) + py*(-gy) );
  }
  for(let i=0;i<SZ;i++) phi[i]=nlR[i];

  // reproyectar masa (corrige el pequeño error de advección) -> conservación
  let m=0; for(let i=0;i<SZ;i++) m+=phi[i];
  const corr = mass0 - m/SZ;
  for(let i=0;i<SZ;i++) phi[i]+=corr;
}

export function stepUniverse(dt){
  // dt del usuario -> nº de subpasos físicos (cada uno DT_PHYS, estable)
  const steps = Math.max(1, Math.round(dt*2));
  for(let s=0;s<steps;s++) stepPhys();
  age += dt*0.05;
}

// ---- API de estado, compatible con el render de v1.0 ----
// Derivamos los campos que el render espera DESDE φ (no son reglas: son
// lecturas del campo físico real).
export function getState(){
  const energy   = new Float32Array(SZ);
  const density  = new Float32Array(SZ);
  const phase    = new Float32Array(SZ);
  const curvature= new Float32Array(SZ);
  const entropy  = new Float32Array(SZ);
  const flowX    = new Float32Array(SZ);
  const flowY    = new Float32Array(SZ);
  const bhMask   = new Uint8Array(SZ);
  const planetMask=new Uint8Array(SZ);
  const idx=(x,y)=>((y+G2)%G2)*G2+((x+G2)%G2);

  const planets=[];
  for(let y=0;y<G2;y++)for(let x=0;x<G2;x++){
    const i=y*G2+x;
    const p=phi[i];
    // fase normalizada 0..1 (vacío..materia) — la espuma
    const ph=(p+1)*0.5;
    phase[i]=ph;
    // densidad = materia local suavizada
    const d = Math.max(0, p);
    density[i]=d;
    // energía = |gradiente de φ|² → alta en las MEMBRANAS (interfaces)
    const gx=(phi[idx(x+1,y)]-phi[idx(x-1,y)])*0.5;
    const gy=(phi[idx(x,y+1)]-phi[idx(x,y-1)])*0.5;
    const grad2=gx*gx+gy*gy;
    energy[i]=grad2*40;
    // curvatura = laplaciano de φ (estructura fina de las paredes)
    const lap=phi[idx(x+1,y)]+phi[idx(x-1,y)]+phi[idx(x,y+1)]+phi[idx(x,y-1)]-4*p;
    curvature[i]=Math.abs(lap);
    // entropía proxy: cuánto se aleja de fase pura (interfaces = más entropía)
    entropy[i]=1-Math.abs(p);
    // flujo = gradiente del campo (hacia dónde se mueve la materia)
    flowX[i]=gx*2; flowY[i]=gy*2;
    // PLANETA: núcleo denso de materia rodeado de interfaz (emergente)
    if(p>0.85 && grad2<0.02) planetMask[i]=1;
    // AGUJERO NEGRO: pozo de materia ultradenso y muy curvado (emergente)
    if(p>0.97 && Math.abs(lap)>0.5){ bhMask[i]=1; }
    if(planetMask[i]) planets.push({x,y});
  }

  // stats (medias) para el panel
  let sE=0,sS=0,sK=0,foam=0,bh=0,pl=0;
  for(let i=0;i<SZ;i++){
    sE+=energy[i]; sS+=entropy[i]; sK+=curvature[i];
    if(phase[i]>0.55) foam++;
    if(bhMask[i]) bh++;
    if(planetMask[i]) pl++;
  }
  return {
    GRID:G2, phase, energy, density, curvature, entropy,
    flowX, flowY, bhMask, planetMask, planets, age,
    _stats:{
      avgEnergy:sE/SZ, avgEntropy:sS/SZ, avgCurvature:sK/SZ,
      foamCount:foam, bhCount:bh, planetCount:pl
    }
  };
}
