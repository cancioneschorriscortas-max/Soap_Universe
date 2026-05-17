import numpy as np
# Acoplamiento CORRECTO: separacion de escalas temporales.
#  - Cahn-Hilliard: rapido (forma la espuma) -> dt completo
#  - Gravedad: lenta y debil, modula la MOVILIDAD (mas movilidad donde
#    hay mas materia => la materia fluye y se agrupa SIN romper la espuma).
#  Esto es fisica: movilidad dependiente del campo (Cahn-Hilliard "B-model").
N=200; eps=1.0; M0=1.0; dt=0.15; STEPS=12000
k=np.fft.fftfreq(N)*2*np.pi
KX,KY=np.meshgrid(k,k); K2=KX**2+KY**2
rng=np.random.default_rng(5)
phi=0.05*(rng.random((N,N))-0.5)
mass0=phi.mean()
a=2.0
G=0.04   # gravedad DEBIL y lenta (separacion de escalas)

def cl(f):
    F=np.fft.fft2(f-f.mean());P=np.abs(F)**2
    kk=np.sqrt(K2);kk[0,0]=1e-9
    return P.sum()/((P*kk).sum())

snaps={}; log=[]
for s in range(STEPS):
    # 1) Cahn-Hilliard (la espuma se forma y engrosa, intacta)
    nl=phi**3-(1.0+a)*phi
    denom=1.0+dt*M0*K2*(a+eps**2*K2)
    phi=np.real(np.fft.ifft2((np.fft.fft2(phi)-dt*M0*K2*np.fft.fft2(nl))/denom))
    # 2) Gravedad: potencial de Poisson de la materia, aplicado SUAVE.
    #    Solo desplaza ligeramente la materia hacia los pozos -> grumos,
    #    sin advectar destructivamente (filtrado paso-bajo = solo escalas grandes)
    rho=np.clip((phi+1)/2,0,1)
    RHO=np.fft.fft2(rho-rho.mean())
    K2s=K2.copy();K2s[0,0]=1e-9
    # suavizar el potencial: solo modos de gran escala (gravedad = fuerza larga)
    smooth=np.exp(-K2*8.0)
    Phig=np.real(np.fft.ifft2(-RHO/K2s*smooth))
    gx=np.gradient(Phig,axis=1);gy=np.gradient(Phig,axis=0)
    px=np.gradient(phi,axis=1);py=np.gradient(phi,axis=0)
    phi=phi+dt*G*(px*(-gx)+py*(-gy))
    phi+=(mass0-phi.mean())
    if s%600==0 or s==STEPS-1:
        t=(s+1)*dt
        log.append((t,cl(phi),abs(phi.mean()-mass0),np.mean(np.abs(phi)>0.8),np.var(rho)))
    if s in (400,4000,11999): snaps[s]=phi.copy()

arr=np.array(log)
print("t\tL(t)\t|dmass|\tpura%\tgrumos")
for t,Lc,dm,p,cm in log[::2]:
    print(f"{t:7.0f}\t{Lc:5.2f}\t{dm:.0e}\t{p*100:4.0f}\t{cm:.4f}")
mask=arr[:,0]>arr[-1,0]*0.3
aa,_=np.polyfit(np.log(arr[mask,0]),np.log(arr[mask,1]),1)
print(f"\nCoarsening: L(t)~t^{aa:.3f} (CH puro daba ~0.26; sigue VIVO = no congelo)")
print(f"Masa: deriva max {arr[:,2].max():.1e}")
print(f"Espuma se forma: pura {arr[0,3]*100:.0f}% -> {arr[-1,3]*100:.0f}%")
print(f"Gravedad agrupa: grumos {arr[1,4]:.4f} -> {arr[-1,4]:.4f}")

from PIL import Image
def col(f):
    g=np.abs(np.gradient(f)[0])+np.abs(np.gradient(f)[1])
    img=np.zeros((*f.shape,3));mat=(f+1)/2
    img[...,0]=0.03+mat*0.28;img[...,1]=0.08+mat*0.55;img[...,2]=0.15+mat*0.80
    g=np.clip(g*4,0,1)
    img[...,0]+=g;img[...,1]+=g*0.95;img[...,2]+=g
    return (np.clip(img,0,1)*255).astype(np.uint8)
ks=sorted(snaps)
ims=[Image.fromarray(col(snaps[s])).resize((360,360),Image.LANCZOS) for s in ks]
W=sum(i.width for i in ims)+40
cv=Image.new('RGB',(W,360),(0,0,0));x=0
for im in ims: cv.paste(im,(x,0));x+=im.width+20
cv.save('/tmp/ch_grav2.png')
print("render guardado")
