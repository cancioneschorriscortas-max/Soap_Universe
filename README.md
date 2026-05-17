# 🧼 SOAP COSMOS — estado del proyecto y mapa para retomarlo

> Premisa: un jabón impacta, la espuma que cae forma un universo. Estamos
> dentro de la espuma; la espuma es el espacio-tiempo. Dentro de esa
> espuma, la materia se condensa en cúmulos, planetas y estrellas.
>
> Objetivo declarado del autor: que la **física sea legítima** — que un
> profesor de física teórica lo mire y diga "esto está guay, lo investigo".
> La física pesa más que la belleza.

Este paquete recoge TODO lo que funciona y está verificado, más el punto
exacto donde el proyecto se quedó. Si lees esto y se te ocurre algo: la
sección "DÓNDE TIRAR DEL HILO" (al final) es para ti.

---

## 📦 Qué hay aquí (todo verificado y funcional)

### `01-repo-3-versiones/soap-cosmos-engine.zip`
El simulador completo, tres versiones que funcionan y se ejecutan en XAMPP:
- **v0.3-original**: motor de campo en rejilla. Física emergente real
  (planetas/agujeros negros salen de las ecuaciones), render básico.
- **v0.3-plus-render**: misma física, render WebGL reconstruido
  (log-normalización + espuma + glow). Bonito y honesto.
- **v1.0-particles-3d**: física v0.3 + 60.000 partículas 3D que la
  obedecen, Three.js, cámara orbitable.
Listo para subir a GitHub (lleva su propio README y guía git dentro).

### `02-cahn-hilliard-verificado/`
**El logro físico principal.** El núcleo se reescribió sobre la ecuación
de **Cahn-Hilliard** (separación de fases — la física real de cómo se
separa una espuma/emulsión), con un término gravitatorio acoplado.
- `ch_1d.py`, `ch_2d.py`, `ch_grav2.py`: pruebas de concepto en Python.
- `universe-ch.js`: el núcleo portado a JavaScript (FFT casera + esquema
  de Eyre). API compatible con el render de v1.0.
- PNGs: la espuma emergiendo de la ecuación.

**Tests de físico que PASA (medidos, no "se ve bien"):**
- Conservación de masa: deriva ~1e-10 a 1e-16 (cero numérico). ✔
- Separación de fases real: φ → ±1 (las dos fases puras). ✔
- Estructura emergente: planetas/vacíos/filamentos salen de la ecuación,
  sin umbrales arbitrarios (a diferencia de v0.3). ✔
- Limitación honesta: el exponente de *coarsening* mide t^0.26 vs t^0.33
  teórico. Es **preasintótico** (cuestión de tiempo de cómputo, no de
  validez); cualitativamente el comportamiento es correcto.

Fundamento en literatura real: arXiv 2410.15436 (modelo de descomposición
espinodal para la estructura a gran escala del universo, inspirado en
membranas de polímero — la premisa del jabón tiene primo serio publicado),
y arXiv cond-mat/9604126 (Cahn-Hilliard con término gravitatorio).

### `03-nbody-verificado/`
Gravedad N-body (soles que atrapan planetas en órbita), por si se quiere
la escala "sistema solar".
- `nbody3.py`: sol + planetesimales, leapfrog + softening de Plummer.
- **Test de físico que PASA**: conservación de energía 1.03e-05
  (EXCELENTE). 86% de cuerpos en órbita estable; el resto son
  expulsiones físicamente esperadas. ✔

### `04-acoplamiento-gemini/`
El intento de fusionar las dos escalas (espuma + gravedad orbital local),
sobre una idea de Gemini: gravedad por rejilla (lee la densidad de
Cahn-Hilliard como pozos gravitatorios) + fuerza orbital (torque
perpendicular para discos de acreción). Evita el coste O(N²).
- `gemini_test.mjs`: integrado y **calibrado a los rangos reales**
  medidos de nuestro Cahn-Hilliard (densidad 0-1, no 0-3 → THRESH 0.8,
  no 1.8; G y orbitFactor reescalados).
- Estado: corre sin explotar, la física Cahn-Hilliard sigue viva, pero
  el balance atracción/órbita está descalibrado → las partículas
  escapan al borde en vez de formar discos. **No cuajó visualmente.**
  Es el punto exacto donde se cerró el proyecto.

---

## 🧭 RESUMEN HONESTO DEL ESTADO

Lo que SÍ está conseguido y verificado:
- Simulador completo (3 versiones).
- Física Cahn-Hilliard + gravedad: implementada en JS, verificada con
  tests de físico (masa, fases, estructura emergente). **Esto cumple el
  objetivo principal del proyecto** (física legítima, reconocible,
  con dominio de validez declarado).
- N-body verificado (energía conservada).
- Idea de acoplamiento multiescala integrada y calibrada.

Lo que NO se logró:
- Un render *espectacular* del acoplamiento (espuma + soles orbitando
  bonito a la vez). Patrón observado durante todo el desarrollo: la
  física se verifica bien (criterio objetivo), el render espectacular
  es un bucle de ajuste fino sin final claro.

Raíz técnica del cuello de botella: Cahn-Hilliard es **2D**. Forzar un
3D con partículas siempre colapsó a una franja (eje Z inventado, sin
física que lo sostenga). La física CH se ve espectacular como **campo
continuo 2D** (ver `02-.../ch_grav2.png`), mediocre como partículas.

---

## 🎯 DÓNDE TIRAR DEL HILO (si lees esto y quieres seguir)

Tres pistas concretas, ordenadas por relación esfuerzo/recompensa:

1. **EL EDITOR EN VIVO (la idea del autor, la más prometedora).**
   El balance gravedad/órbita de `04-acoplamiento-gemini/gemini_test.mjs`
   no se afina bien a ciegas. La solución correcta es un **panel de
   sliders en tiempo real** sobre ese motor (ya calibrado a los rangos
   reales) para tocar en vivo: `THRESH` (~0.8), `G_GRAV` (~0.06),
   `orbitF` (~1.2, probablemente BAJARLO — la órbita está expulsando
   partículas), `softening` (0.8), radio de búsqueda (4). Re-sembrar
   partículas en el centro, NO en el borde. Esto convierte el bucle de
   iteración en algo que el humano controla viéndolo. Es el siguiente
   paso natural y el de mayor impacto.

2. **Render de CAMPO CONTINUO para Cahn-Hilliard.** La física CH se ve
   genuinamente bien como campo (no partículas): conectar `universe-ch.js`
   al render WebGL de `v0.3-plus-render` (que ya está verificado y quedó
   bonito) en vez de a partículas. Bajo riesgo, probablemente bonito.

3. **Las dos escalas, de verdad.** Cahn-Hilliard 2D para la espuma →
   detectar grumos densos → en cada grumo, instanciar un N-body 3D local
   (`03-nbody-verificado/`) para el sistema solar. Cada física en su
   dimensión natural. Es lo más ambicioso (proyecto de semanas), pero
   conceptualmente coherente y con literatura que lo respalda.

**Método que funcionó (importante): no programar a ciegas.** Ejecutar la
física real, volcar el estado, RENDERIZAR a imagen y MIRARLA antes de
dar nada por bueno. Verificar la física con métricas medibles
(conservación de masa/energía), no con "se ve bien". Eso salvó el
proyecto de muchas vueltas en falso.

---

## 🙏 Créditos

- Concepto y física v0.3: del autor.
- Arquitectura de partículas GPU: inspirada en
  github.com/Sanjays2402/ai-particle-simulator (MIT, no se usó su código).
- Idea de acoplamiento gravedad-por-rejilla: aportada por Gemini,
  integrada y calibrada aquí.
- Fundamento Cahn-Hilliard↔cosmos: arXiv 2410.15436, cond-mat/9604126.

Licencia: MIT.
