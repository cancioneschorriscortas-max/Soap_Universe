import numpy as np
# Cahn-Hilliard 2D con Eyre convex-splitting (esquema verificado en 1D)
N=200; dx=1.0; eps=1.0; M=1.0; dt=0.15; STEPS=14000
k=np.fft.fftfreq(N,d=dx)*2*np.pi
KX,KY=np.meshgrid(k,k); K2=KX**2+KY**2
rng=np.random.default_rng(3)
phi=0.05*(rng.random((N,N))-0.5)        # quench 50/50
mass0=phi.mean()
a=2.0
denom=1.0+dt*M*K2*(a+eps**2*K2)

def cl(f):
    F=np.fft.fft2(f-f.mean()); P=np.abs(F)**2
    kk=np.sqrt(K2); kk[0,0]=1e-9
    return P.sum()/((P*kk).sum())

snaps={}
log=[]
for s in range(STEPS):
    nl=phi**3-(1.0+a)*phi
    PHI_new=(np.fft.fft2(phi)-dt*M*K2*np.fft.fft2(nl))/denom
    phi=np.real(np.fft.ifft2(PHI_new))
    if s%500==0 or s==STEPS-1:
        t=(s+1)*dt
        pure=np.mean(np.abs(phi)>0.8)
        log.append((t,cl(phi),abs(phi.mean()-mass0),pure))
    if s in (200, 2000, 13999):
        snaps[s]=phi.copy()

arr=np.array(log)
print("t\tL(t)\t|dmass|\tfase_pura%")
for t,Lc,dm,p in log[::3]:
    print(f"{t:7.0f}\t{Lc:6.2f}\t{dm:.1e}\t{p*100:5.1f}")
mask=arr[:,0]>arr[-1,0]*0.3
aa,bb=np.polyfit(np.log(arr[mask,0]),np.log(arr[mask,1]),1)
print(f"\n*** COARSENING medido: L(t) ~ t^{aa:.3f}   (teoria Cahn-Hilliard: t^0.333) ***")
print(f"*** Conservacion de masa: deriva max = {arr[:,2].max():.1e}  (cero numerico) ***")
print(f"*** Fase pura final: {arr[-1,3]*100:.1f}%  (la espuma se ha formado) ***")

# RENDER de los 3 momentos
from PIL import Image
def col(f):
    # phi=-1 vacio (negro), phi=+1 materia (claro), interfaz = membrana brillante
    g=np.abs(np.gradient(f)[0])+np.abs(np.gradient(f)[1])
    img=np.zeros((*f.shape,3))
    mat=(f+1)/2
    img[...,0]=0.05+mat*0.25
    img[...,1]=0.10+mat*0.55
    img[...,2]=0.18+mat*0.75
    g=np.clip(g*4,0,1)
    img[...,0]+=g*0.9; img[...,1]+=g*0.95; img[...,2]+=g*1.0
    return (np.clip(img,0,1)*255).astype(np.uint8)
keys=sorted(snaps)
ims=[Image.fromarray(col(snaps[s])).resize((360,360),Image.NEAREST) for s in keys]
W=sum(i.width for i in ims)+40
canvas=Image.new('RGB',(W,360),(0,0,0))
x=0
for im in ims:
    canvas.paste(im,(x,0)); x+=im.width+20
canvas.save('/tmp/ch_evolution.png')
print("\nrender de la evolucion guardado (3 momentos: caos -> espuma -> coarsening)")
