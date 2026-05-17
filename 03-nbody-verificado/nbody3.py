import numpy as np
rng=np.random.default_rng(42)
N=300; G=1.0; SOFT=0.25; DT=0.001; STEPS=15000  # mas softening, dt menor

mass=np.empty(N); mass[0]=200.0; mass[1:]=rng.uniform(0.1,0.8,N-1)
pos=np.zeros((N,3)); vel=np.zeros((N,3))
for i in range(1,N):
    r=rng.uniform(3.0,12.0); th=rng.uniform(0,2*np.pi)
    pos[i]=[r*np.cos(th),r*np.sin(th),rng.uniform(-0.3,0.3)]
    vo=np.sqrt(G*mass[0]/np.sqrt(r*r+SOFT*SOFT))  # vel circular con softening
    vel[i]=[-np.sin(th)*vo,np.cos(th)*vo,0]

def accel(p):
    d=p[None,:,:]-p[:,None,:]
    r2=(d**2).sum(2)+SOFT**2
    np.fill_diagonal(r2,1e18)
    return G*(d*(mass[None,:]/r2**1.5)[:,:,None]).sum(1)
def energy(p,v):
    ke=0.5*(mass[:,None]*v**2).sum()
    d=p[None,:,:]-p[:,None,:]; r=np.sqrt((d**2).sum(2)+SOFT**2)
    iu=np.triu_indices(N,1)
    return ke-G*((mass[:,None]*mass[None,:])[iu]/r[iu]).sum()

E0=energy(pos,vel); a=accel(pos); log=[]; snaps={}
for s in range(STEPS):
    vel+=0.5*DT*a; pos+=DT*vel; a=accel(pos); vel+=0.5*DT*a
    if s%2000==0 or s==STEPS-1:
        E=energy(pos,vel); dr=abs((E-E0)/E0)
        rr=np.sqrt((pos[1:]**2).sum(1))
        log.append((s*DT,dr,rr.mean(),np.mean(rr<30)))
    if s in (0,7000,14999): snaps[s]=pos.copy()
print("t\tderiva_E\tr_medio\t%ligados")
for t,d,rm,b in log: print(f"{t:6.1f}\t{d:.2e}\t{rm:5.2f}\t{b*100:.0f}%")
f=log[-1]
print(f"\nConservacion energia: {f[1]:.2e} -> "+("EXCELENTE" if f[1]<0.01 else "ACEPTABLE" if f[1]<0.05 else "FALLA"))
print(f"Sistema estable: {f[3]*100:.0f}% orbitando")
# render de las 3 instantaneas (vista cenital del sistema)
from PIL import Image, ImageFilter
ims=[]
for k in sorted(snaps):
    P=snaps[k]; S=480; img=np.zeros((S,S,3),np.float32)
    sx=(S/2+P[:,0]*16).astype(int); sy=(S/2+P[:,1]*16).astype(int)
    for i in range(N):
        x,y=sx[i],sy[i]
        if 0<=x<S and 0<=y<S:
            if i==0: 
                for ddx in range(-3,4):
                    for ddy in range(-3,4):
                        if 0<=x+ddx<S and 0<=y+ddy<S: img[y+ddy,x+ddx]+=[1.0,0.85,0.4]
            else: img[y,x]+=[0.5,0.7,1.0]
    pim=Image.fromarray((np.clip(img,0,1)*255).astype(np.uint8))
    gl=pim.filter(ImageFilter.GaussianBlur(2))
    ims.append(Image.fromarray(np.clip(np.array(pim)+np.array(gl)*0.6,0,255).astype(np.uint8)))
cv=Image.new('RGB',(480*3+40,480),(0,0,0)); x=0
for im in ims: cv.paste(im,(x,0)); x+=500
cv.save('/tmp/nbody_seq.png')
print("render: inicio -> medio -> final")
