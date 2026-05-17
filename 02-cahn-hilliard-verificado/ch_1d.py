import numpy as np
# Diagnostico 1D: ¿se separa phi en -1/+1? Aqui es facil de ver.
# Esquema correcto: term lineal -phi tambien implicito.
# d phi/dt = M d2/dx2 ( phi^3 - phi - eps^2 d2/dx2 phi )
# En Fourier: d PHI/dt = -M k^2 ( FFT(phi^3) - PHI + eps^2 k^2 PHI )
# Eyre: parte contractiva (phi^3 + a*phi) explicita, resto implicito.

N=512; L=100.0; dx=L/N; eps=1.0; M=1.0; dt=0.1; STEPS=20000
x=np.arange(N)*dx
k=np.fft.fftfreq(N, d=dx)*2*np.pi
k2=k*k
rng=np.random.default_rng(1)
phi=0.05*(rng.random(N)-0.5)
mass0=phi.mean()

# Eyre convex splitting: f(phi)=phi^3-phi.  f = fc' - fe'
#  fc = phi^4/4 + (a/2)phi^2  (convexa, implicita)
#  fe = (1+a)/2 phi^2         (concava, explicita)
# => mu = phi^3 + a phi   (impl. lineal en a)  - (1+a) phi (expl.)
a=2.0
denom = 1.0 + dt*M*k2*( a + eps**2*k2 )
for s in range(STEPS):
    nl = phi**3 - (1.0+a)*phi          # parte explicita
    NL = np.fft.fft(nl)
    PHI= np.fft.fft(phi)
    PHI_new = (PHI - dt*M*k2*NL) / denom
    phi = np.real(np.fft.ifft(PHI_new))
    if s%4000==0 or s==STEPS-1:
        print(f"t={ (s+1)*dt:7.0f}  phi[min={phi.min():+.3f} max={phi.max():+.3f}]  "
              f"mean={phi.mean():+.2e} (dmass={abs(phi.mean()-mass0):.1e})")

# cuantos puntos estan cerca de +1 o -1 (fases puras)?
pure = np.mean((phi>0.8)|(phi<-0.8))
print(f"\nFraccion del dominio en fase pura (|phi|>0.8): {pure*100:.1f}%")
print("Si es alto (>70%), la separacion de fases FUNCIONA.")
np.save('/tmp/ch1d.npy', phi)
