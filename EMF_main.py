#! /usr/bin/python3
import sys
import math
import numpy as np
import scipy.stats
import scipy.special as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import warnings
warnings.filterwarnings("ignore")
import modulo as mtd


##############################
##      EXPERIMENTO 1       ##
##############################

##  FIGURA 7    ##

#obtener datos del .txt
#M7=np.loadtxt('f7x6a.txt')
M7=np.loadtxt('f7x6b.txt')
#M7=np.loadtxt('f7x6c.txt')


def func7a1r(b,x):
    r=20
    a=r
    return 1/r*((1+((x+a/2)/r)**2)**-1.5+(1+((x-a/2)/r)**2)**-1.5), b/r*((1+((x+a/2)/r)**2)**-1.5+(1+((x-a/2)/r)**2)**-1.5)

def func7ar2(b,x):
    r=20
    a=r/2
    return 1/r*((1+((x+a/2)/r)**2)**-1.5+(1+((x-a/2)/r)**2)**-1.5), b/r*((1+((x+a/2)/r)**2)**-1.5+(1+((x-a/2)/r)**2)**-1.5)

def func7a2r(b,x):
    r=20
    a=2*r
    return 1/r*((1+((x+a/2)/r)**2)**-1.5+(1+((x-a/2)/r)**2)**-1.5), b/r*((1+((x+a/2)/r)**2)**-1.5+(1+((x-a/2)/r)**2)**-1.5)


a1r=mtd.nolingen(np.column_stack((M7[:,0],M7[:,1])),1,func7a1r,1e-3)
ar2=mtd.nolingen(np.column_stack((M7[:,2],M7[:,3])),1,func7ar2,1e-3)
a2r=mtd.nolingen(np.column_stack((M7[:,4],M7[:,5])),1,func7a2r,1e-3)


plt.figure(1)
plt.plot(M7[:,0],M7[:,1],'o-',label=r'$\alpha=R$',markersize=5,lw=.3)
plt.plot(M7[:,2],M7[:,3],'o-',label=r'$\alpha=R/2$',markersize=5,lw=.3)
plt.plot(M7[:,4],M7[:,5],'o-',label=r'$\alpha=2R$',markersize=5,lw=.3)
plt.plot(M7[:,0],func7a1r(float(a1r[0]),M7[:,0])[1],color='C0',label='$%1.2f\,A$, $r^2=%1.2f$'%(a1r[0]/2/np.pi/154*100,a1r[1][0]))
plt.plot(M7[:,2],func7ar2(float(ar2[0]),M7[:,2])[1],color='C1',label='$%1.2f\,A$, $r^2=%1.2f$'%(ar2[0]/2/np.pi/154*100,ar2[1][0]))
plt.plot(M7[:,4],func7a2r(float(a2r[0]),M7[:,4])[1],color='C2',label='$%1.2f\,A$, $r^2=%1.2f$'%(a2r[0]/2/np.pi/154*100,a2r[1][0]))
plt.xlabel('$z$ [cm]')
plt.ylabel('$B_z$ [mT]')
plt.grid(ls=':',color='grey',alpha=.5)
plt.legend()


##  FIGURA 8    ##
#M8=np.loadtxt('f8x5a.txt')
M8=np.loadtxt('f8x5b.txt')
#M8=np.loadtxt('f8x5c.txt')

def func8a1r(b,x):
    r=20
    a=r
    return 1/r*((1+((x+a/2)/r)**2)**-1.5+(1+((x-a/2)/r)**2)**-1.5), b/r*((1+((x+a/2)/r)**2)**-1.5+(1+((x-a/2)/r)**2)**-1.5)

f8a1r=mtd.nolingen(np.column_stack((M8[:,0],M8[:,1])),1,func8a1r,1e-3)

label8=['$r=0\,cm$','$r=10\,cm$','$r=14\,cm$','$r=16\,cm$']

plt.figure(2)
for i in range(4):
    plt.plot(M8[:,0],M8[:,i+1],'o-',label=label8[i],markersize=5,lw=.3)
plt.plot(M8[:,0],func8a1r(float(f8a1r[0]),M8[:,0])[1],color='C0',label='$%1.2f\,A$, $r^2=%1.2f$'%(f8a1r[0]/2/np.pi/154*100,f8a1r[1][0]))
plt.xlabel('$z$ [cm]')
plt.ylabel('$B_z$ [mT]')
plt.grid(ls=':',color='grey',alpha=.5)
plt.legend()


##  FIGURA 9    ##
#M9=np.loadtxt('f9x5a.txt')
M9=np.loadtxt('f9x5b.txt')
#M9=np.loadtxt('f9x5c.txt')

label9=['$r=0\,cm$','$r=10\,cm$','$r=14\,cm$','$r=16\,cm$']

plt.figure(3)
for i in range(4):
    plt.plot(M9[:,0],M9[:,i+1],'o-',label=label9[i],markersize=5,lw=.3)
plt.xlabel('$z$ [cm]')
plt.ylabel('$B_r$ [mT]')
plt.grid(ls=':',color='grey',alpha=.5)
plt.legend()


##  FIGURA 10    ##
#M10=np.loadtxt('f10x3a.txt')
M10=np.loadtxt('f10x3b.txt')
#M10=np.loadtxt('f10x3c.txt')

def fk(M):
    a=20
    r=M
    th=np.pi/2
    return ((4*a*r*np.sin(th))/(a**2+r**2+2*a*r*np.sin(th)))**.5

#kk=fk(M10[1,1])
#print(kk)
#print(((2-kk**2)*sp.ellipk(kk)-2*sp.ellipe(kk))/kk**2)



plt.figure(4)
plt.plot(M10[:,0],M10[:,1],'o-',label='$Coil\,1,\;z=+R/2$',markersize=5,lw=.3)
plt.plot(M10[:,0],M10[:,2],'o-',label='$Coil\,2,\;z=-R/2$',markersize=5,lw=.3)
plt.plot(M10[:,0],M10[:,1]+M10[:,2],'--')
plt.xlabel('$r$ [cm]')
plt.ylabel('$B_r$ [mT]')
plt.grid(ls=':',color='grey',alpha=.5)
plt.legend()



##############################
##      EXPERIMENTO 2       ##
##############################

##  FIGURA 2    ##
M2=np.loadtxt('f2x2.txt')

def lmsq(M):
    qK=sum(M[:,0]*M[:,1])/sum(M[:,0]**2)
    ym=sum(M[:,1])/len(M)
    st=sum((M[:,1]-ym)**2)
    sr=sum((M[:,1]-qK*M[:,0])**2)
    return qK,((st-sr)/st)

a2=lmsq(M2)

plt.figure(5)
plt.plot(M2[:,0],M2[:,1],'o',markersize=5)
plt.plot(M2[:,0],a2[0]*M2[:,0],'--',color='C0',lw=.8)
plt.ylabel('$^h\,B_H$ [mT]')
plt.xlabel('$I_H$ [A]')
plt.grid(ls=':',color='grey',alpha=.5)
plt.xlim(0,3)
plt.ylim(0,2)
plt.text(2,.5,'$K=%1.3f$ $[\dfrac{mT}{A}]$\n$r^2=%1.4f$'%(a2[0],a2[1]))


##  FIGURA 4    ##

M4aux=np.array([[0,.1,.12,.15,.2,.25,.3,.35,.4,.67],[0,70,74,76,78,80,82,83,84,86]]).T

M4=np.copy(M4aux)
for i in range(len(M4)):
    M4[i,0]=np.tan(M4aux[i,1]*np.pi/180)
    M4[i,1]=M4aux[i,0]*a2[0]*1000

T4=np.array([[28,32],[20,56],[18,38]])

a4=lmsq(M4)

for i in range(len(T4)):
    print('Ensayo %g'%(i+1))
    vbe=np.tan(sum(T4[i])/2*np.pi/180)*a4[0]
    be=(vbe**2+a4[0]**2)**.5
    print('h^B_E= %1.2f [\u03BCT]\tv^B_E= %1.2f [\u03BCT]\tB_E= %1.2f [\u03BCT]'%(a4[0],vbe,be))

plt.figure(6)
plt.plot(M4[:,0],M4[:,1],'o',color='C2',markersize=5)
plt.plot(M4[:,0],a4[0]*M4[:,0],'--',color='C2',lw=.8)
plt.xlabel(r'$\tan(\alpha)$')
plt.ylabel(r'$I_H \cdot K$ [$\mu T$]')
plt.grid(ls=':',color='grey',alpha=.5)
plt.xlim(0,16)
plt.ylim(0,500)
plt.text(10,200,'$^h\,B_E=%1.2f$ [$\mu$T] \n $r^2=%1.4f$'%(a4[0],a4[1]))


## MAPA CALOR B_Z(z;r) ##

z=np.arange(-20,21)
#r=np.array([0,10,14,16])
r=np.array([-16,-14,10,0,10,14,16])
s=np.loadtxt('bz_zr2.txt')

fig,ax=plt.subplots()
im=ax.imshow(s.T,cmap='hot_r')
ax.set_yticks(np.arange(len(r)))
ax.set_xticks(np.arange(len(z)))
ax.set_yticklabels(r)
ax.set_xticklabels(z)
ax.set_ylabel('$r\;[cm]$')
ax.set_xlabel('$z\;[cm]$')
cbar=ax.figure.colorbar(im, ax=ax, orientation='horizontal')
cbar.ax.set_title('$B_z\;[mT]$')#,location='bottom')

plt.show()
