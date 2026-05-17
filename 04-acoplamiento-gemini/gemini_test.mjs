import { initUniverse, stepUniverse, getState } from './universe-ch.js';
import fs from 'fs';
initUniverse();
for(let i=0;i<7;i++) stepUniverse(25);   // Cahn-Hilliard desarrollado

const st0=getState(); const G=st0.GRID;
const N=35000;
// particulas de gas libre (las que Gemini hace orbitar)
const part=[];
for(let i=0;i<N;i++) part.push({x:Math.random()*G, y:Math.random()*G, vx:0, vy:0});

// === MOTOR DE GEMINI, calibrado a nuestros numeros reales ===
function applyCosmicGravity(particles, state){
  const { GRID, density } = state;
  const G_GRAV   = 0.06;    // calibrado (Gemini puso 1.5, escala distinta)
  const softening= 0.8;
  const orbitF   = 1.2;
  const THRESH   = 0.8;     // CALIBRADO: nuestra densidad es 0-1, no 0-3 (Gemini puso 1.8)
  const R = 4;

  for(const p of particles){
    const gx=Math.floor(p.x), gy=Math.floor(p.y);
    if(gx<=R||gx>=GRID-R||gy<=R||gy>=GRID-R) continue;
    let fX=0,fY=0,tX=0,tY=0;
    for(let ny=-R;ny<=R;ny++)for(let nx=-R;nx<=R;nx++){
      const tx=gx+nx, ty=gy+ny;
      const d=density[ty*GRID+tx];
      if(d>THRESH){
        const dx=tx-p.x, dy=ty-p.y;
        const distSq=dx*dx+dy*dy+softening;
        const dist=Math.sqrt(distSq);
        const f=(G_GRAV*d)/distSq;
        fX+=(dx/dist)*f;  fY+=(dy/dist)*f;
        tX+=(-dy/dist)*f*orbitF;  tY+=(dx/dist)*f*orbitF;  // fuerza orbital
      }
    }
    p.vx=(p.vx+fX+tX)*0.98;       // friccion (condensa en planetas)
    p.vy=(p.vy+fY+tY)*0.98;
    p.x+=p.vx; p.y+=p.vy;
    // si sale del mapa, reaparece como gas nuevo
    if(p.x<2||p.x>=G-2||p.y<2||p.y>=G-2){ p.x=Math.random()*G; p.y=Math.random()*G; p.vx=p.vy=0; }
  }
}

const frames=[];
for(let f=0;f<4;f++){
  stepUniverse(8);
  const st=getState();
  for(let step=0;step<50;step++) applyCosmicGravity(part, st);
  const half=G/2, snap=[];
  for(let i=0;i<N;i+=2){
    const p=part[i];
    const gx=Math.floor(p.x), gy=Math.floor(p.y);
    const d=st.density[gy*G+gx]||0;
    const sp=Math.sqrt(p.vx*p.vx+p.vy*p.vy);
    let r,g,b;
    if(d>0.8){ r=1.0;g=0.80;b=0.38; }            // en grumo: planeta/sol ambar
    else if(sp>0.5){ r=0.9;g=0.6+sp*0.1;b=0.4; } // moviendose rapido: estela calida
    else { r=0.2;g=0.5;b=0.85; }                  // gas libre: azul
    snap.push([p.x-half, p.y-half, r,g,b, sp]);
  }
  frames.push(snap);
  // medir: ¿se estan agrupando? varianza de posiciones (baja = colapso a grumos)
  let mx=0,my=0; for(const p of part){mx+=p.x;my+=p.y;} mx/=N;my/=N;
  let v=0; for(const p of part) v+=(p.x-mx)**2+(p.y-my)**2; v/=N;
  console.log(`f${f}: dispersion=${Math.sqrt(v).toFixed(1)} (baja=agrupado) | planetas_fisica=${st._stats.planetCount}`);
}
fs.writeFileSync('/tmp/gemini.json',JSON.stringify(frames));
console.log('motor Gemini calibrado OK');
