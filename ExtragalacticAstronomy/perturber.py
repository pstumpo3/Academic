#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from scipy import integrate
import math
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product,combinations
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


#####DEFINISICO I PARAMETRI DELLA SFERA
R = 40 #pc
M_sole = 1.99*10**30 #kg
Mtot = 1 #M_sole kg
a = 1
Npart = 30000

P = np.random.random(Npart)
P1 = np.random.random(Npart)
P2 = np.random.random(Npart)
######CAMPIONO LE POSIZIONI
##campiono theta
theta= np.arccos(2*P2-1)
plt.hist(theta,35,density=True,label='campionamento')
plt.xlabel('$\\theta$ $_{radianti}$')
plt.ylabel('p($\\theta$)')
thetareale=np.linspace(0,np.pi,10000)
i1=(integrate.simps((np.sin(thetareale)),thetareale))
plt.plot(thetareale,(np.sin(thetareale))/i1, label='p($\\theta$)$_{teorico})$')

plt.legend()
plt.title('$\\theta$$_{r}$')
#plt.savefig('plummeriniziale_theta.png')
plt.show()

#campiono phi
phi=2*np.pi*P
normalizzazione=np.empty(len(phi))
for i in range (len(phi)):
    normalizzazione[i]=1/(2*np.pi)
plt.hist(phi,30,density=True, label='campionamento')
plt.plot(phi,normalizzazione,label='p($\phi$)$_{teorico})$')
plt.xlabel('$\phi$ $_{radianti}$')
plt.ylabel('p($\phi$)')

plt.legend()
plt.title('$\phi$$_{r}$')
#plt.savefig('plummeriniziale_phi.png')
plt.show()

##campiono R
#definisco prima la funzione di densità
a = 1 #*3 *10**16
rreale = np.linspace(0,R,30000)
densità = 3*a**2 / (4*np.pi) * Mtot /((a**2 + rreale**2)**(5/2))
densità2= 3/(4*np.pi*a**3) * Mtot / ((1 + (rreale**2/a**2)/3)**(5/2)) 
def densità3(variabile,massa,parametro):
    return 3/(4*np.pi*parametro**3) * massa / ((1 + (variabile**2/parametro**2)/3)**(5/2)) * variabile**2 

maximusrho = max (3/(4*np.pi*a**3) * Mtot / ((1 + (rreale**2/a**2)/3)**(5/2)) *rreale**2 )
print (maximusrho)
integrale=(integrate.simps((densità2*rreale**2), rreale))
normalizzazione10 = integrale
#plt.plot(rreale,densità2*rreale**2)
plt.show()

distanza = []
i = 0
while i < (Npart):
    k = float(np.random.random(1)*R)
    j = np.random.random(1)*1.1 * maximusrho
    #print (k)
    #print (j)
    #print (densità3(k,Mtot,a))
    if j <= (densità3(k,Mtot,a)):
        distanza.append(k)
        i = i + 1 
        #print ('\n',distanza,i,'\n')

plt.hist(distanza,150,density=True,label='campionamento')
plt.xlabel('r $_{parsec}$')
plt.ylabel('densità di corpi per r$^2$')
plt.plot(rreale,((densità2) * rreale**2)/normalizzazione10,label='densità$_{teorica}$')

plt.xlim(0,R)
plt.legend()
plt.title('Distribuzione dei raggi')
#plt.savefig('plummeriniziale_r.png')
plt.show()

distanza.sort()
#integrale5=(integrate.simps((densità2*x**2)/normalizzazione, x))

#cambio in coordinate cartesiane
xc= distanza * np.sin(theta)*np.cos(phi)
yc= distanza * np.sin(theta)*np.sin(phi)
zc= distanza * np.cos(theta)


#grafico la mia situazione iniziale

fig = plt.figure()
ax = fig.gca(projection='3d')
zdata=zc#/(3 *10**16)
ydata=yc#/(3 *10**16)
xdata=xc#/(3 *10**16)
ax.set_xlim3d(-R,R)
ax.set_ylim3d(-R,R)
ax.set_zlim3d(-R,R)
ax.set_zlabel('z [pc]')
ax.set_xlabel('x [pc]')
ax.set_ylabel('y [pc]')
#ax.set_title('Condizione iniziale')
ax.scatter3D(xdata,ydata,zdata,c='red',s = 0.1)
plt.savefig('plummerinizialep.png')

#salvo il mio raggio lagrangiano iniziale
distanza.sort()
print (distanza[9000])


# In[3]:


#####CAMPIONO LE VELOCITà
######CAMPIONO gli angoli della VELOCITA'
PV = np.random.random(Npart)
PV1 = np.random.random(Npart)
PV2 = np.random.random(Npart)

get_ipython().run_line_magic('matplotlib', 'inline')
#campiono theta
theta1= np.arccos(2*PV2-1)
plt.hist(theta1,30,density=True,label = 'campionamento')
plt.xlabel('$\\theta$$_{radianti}$')
plt.ylabel('p($\\theta$)')
thetareale1=np.linspace(0,np.pi,10000)
plt.plot(thetareale1,np.sin(thetareale1)/2, label = 'p($\\theta$)$_{teorico})$')

plt.legend()
plt.title('$\\theta$$_{v}$')
#plt.savefig('plummerinizialev_theta')
plt.show()

#campiono phi
phi1=2*np.pi*PV1
normalizzazione6=np.empty(len(phi1))
for i in range (len(phi1)):
    normalizzazione6[i]=1/(2*np.pi)
plt.hist(phi1,30,density=True,label = 'campionamento')
plt.plot(phi1,normalizzazione6,label= 'p($\phi$)$_{teorico}$')
plt.xlabel('$\phi$ $_{radianti}$')
plt.ylabel('p($\phi$) $')
plt.legend()
plt.title('$\phi$$_{v}$')
#plt.savefig('plummerinizialev_phi.png')
plt.show()


# In[4]:


####CAMPIONO LA VELOCITà
distanza = np.array(distanza)

np.random.seed(2)
v_fuga = np.zeros([Npart])
psi = np.zeros([Npart])

for i in range (Npart):
    psi[i] = np.sqrt((Mtot)/distanza[i]**2+(a)**2)
    v_fuga[i] = np.sqrt(2*(psi[i]))


P4 = np.random.uniform(0,1,150000)
P5 = np.random.uniform(0,1,150000) 
q = np.zeros([150000])
X5 = 0.1*P5
prova = np.zeros([150000])

g = (P4**2)*(1-P4**2)**(3.5)
for i in range (150000):
    if X5[i] < g[i]:
        q[i] = P4[i]
        prova[i] = P5[i]
        i = i+1
    else:
        i = i+1

q = q[np.nonzero(q)]
prova = prova[np.nonzero(prova)]
q = q[:30000]
prova = prova[:30000]
ves = np.sqrt(2)*q*(1+distanza**2)**(-0.25)

v_x= ves * np.sin(theta1)*np.cos(phi1)
v_y= ves * np.sin(theta1)*np.sin(phi1)
v_z= ves * np.cos(theta1)

#grafico la mia situazione iniziale delle velocità

fig = plt.figure()
ax = fig.gca(projection='3d')
zdata=v_x
ydata=v_y
xdata=v_z
ax.set_xlim3d(-1,1)
ax.set_ylim3d(-1,1)
ax.set_zlim3d(-1,1)
ax.set_zlabel('z')
ax.set_xlabel('x')
ax.set_ylabel('y')
#ax.set_title('Condizione iniziale')
ax.scatter3D(xdata,ydata,zdata,c='red',s = 1)
#plt.savefig('plummeriniziale.png')


# In[5]:


###SCRIVO IL FILE TXT
m = Mtot / Npart
masse = np.zeros(Npart)
i = 0
for i in range (Npart):
    masse[i] = m

np.savetxt('plummerperturbermasseprova.txt',np.column_stack([masse]))
np.savetxt('plummerperturbercoordinateprova.txt',np.column_stack([xc,yc,zc]))
np.savetxt('plummerperturbervelocitàprova.txt',np.column_stack([v_x,v_y,v_z]))


# In[6]:


fname='plummercondizioni.txt' 
fp = open(fname,'w') 
fp.write('30000') 
fp.write('\n') 
fp.write('3') 
fp.write('\n') 
fp.write('0') 
fp.write('\n') 
fp.close()

filenames = ['plummercondizioni.txt','plummerperturbermasseprova.txt', 'plummerperturbercoordinateprova.txt','plummerperturbervelocitàprova.txt']
with open('plummer_inizio.txt', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)


# In[8]:


####ANALIZZO I MIEI RISULTATI
import pandas as pd
from scipy import integrate
import math
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product,combinations

Npart = 30000

fp=open('plummer_fine.txt','r')
Nlines=len(fp.readlines())
Nsnapshot=(Nlines)//(3*Npart+3)
print (Nlines,Nsnapshot)
t = np.zeros(Nsnapshot)
posp = np.zeros([Nsnapshot,Npart,3])
vfp = np.zeros([Nsnapshot,Npart,3])
Npart = 30000
mf = np.zeros([Nsnapshot,Npart])

i=0

for i in range (Nsnapshot):
    initial_line=i*(3*Npart+3)   #è il numero di linee che passano tra uno snapshot e l'altro
   
    #mf[i]=np.genfromtxt('perturberplummer_fine.txt',skip_header=initial_line+3,skip_footer=Nlines-initial_line-4)
    t[i]=np.genfromtxt('plummer_fine.txt',skip_header=initial_line+2, skip_footer=Nlines-initial_line-3)
    posp[i]=np.genfromtxt('plummer_fine.txt',skip_header=initial_line+30003, skip_footer=Nlines-initial_line-60003)
    vfp[i]=np.genfromtxt('plummer_fine.txt',skip_header=initial_line+60003, skip_footer=Nlines-initial_line-90003)
    #print(t[i],mf[i],pos[i],vf[i])


# In[11]:


z_2=np.zeros(Npart)
x_2=np.zeros(Npart)
y_2=np.zeros(Npart)
r2 = np.zeros([Nsnapshot,Npart])

for u in range (Nsnapshot):
    for i in range (Npart):
        z_2[i]=posp[u,i,2]
        x_2[i]=posp[u,i,0]
        y_2[i]=posp[u,i,1]
        #print (z_2)
    for w in range (Npart):
        r2[u,w]=((x_2[w]**2+y_2[w]**2+z_2[w]**2)**(1/2))
    

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlim3d(-R,R)
    ax.set_ylim3d(-R,R)
    ax.set_zlim3d(-R,R)
    ax.set_zlabel('z')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    #ax.set_title('5 t$_{dyn}$')

    ax.scatter3D(x_2,y_2,z_2,color='red', s = 0.7)
    #plt.savefig('2tempo_collasso.png')
    plt.show()  


#plt.plot(t,r2)


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlim3d(-R,R)
ax.set_ylim3d(-R,R)
ax.set_zlim3d(-R,R)
ax.set_zlabel('z [m]')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title('Situazione finale 360 My')

ax.scatter3D(x_2,y_2,z_2,color='red', s = 0.1)
plt.savefig('situazionefinaleplummer.png')
plt.show()  


# In[13]:


######CALCOLO IL CENTRO DI MASSA PERCHè IL MIO SISTEMA SI STA SPOSTANDO PROBABILMENTE,E VERIFICO
Mtot = 1 #* 2 *10**30 #kg
x_cm = np.zeros(Nsnapshot)
y_cm = np.zeros(Nsnapshot)
z_cm = np.zeros(Nsnapshot)
r_cm = np.zeros(Nsnapshot)

xnew = np.zeros([Nsnapshot,Npart])
ynew = np.zeros([Nsnapshot,Npart])
znew = np.zeros([Nsnapshot,Npart])
r2new = np.zeros([Nsnapshot,Npart])

for u in range (Nsnapshot):
    x_cm[u] = (sum(Mtot/Npart * posp[u,:,0])/(Mtot))
    y_cm[u] = (sum(Mtot/Npart * posp[u,:,1])/Mtot)
    z_cm[u] = (sum(Mtot/Npart * posp[u,:,2])/Mtot)
    r_cm[u] = np.sqrt(x_cm[u]**2+y_cm[u]**2+z_cm[u]**2)
    
        
   
   
plt.plot(t,x_cm, c = 'r',label='x_cm')
plt.plot(t,y_cm, c = 'b',label= 'y_cm')
plt.plot(t,z_cm, c = 'g',label = 'z_cm')
plt.plot(t,r_cm, c = 'black', ls='--', label = 'posizione del centro di massa')
plt.xlabel('t [s]')
plt.ylabel('')
plt.legend()
#plt.savefig('movcentrodimassa_h.png')


# In[15]:


###riscalo tutto nel centro di massa
x_2new=np.zeros([Npart])
y_2new=np.zeros([Npart])
z_2new=np.zeros([Npart])
print (R)

for u in range (Nsnapshot):
    for i in range (Npart):
        z_2new[i]=posp[u,i,2]-z_cm[u]
        x_2new[i]=posp[u,i,0]-x_cm[u]
        y_2new[i]=posp[u,i,1]-y_cm[u]
    for w in range (Npart):
        r2new[u,w]=((x_2new[w]**2+y_2new[w]**2+z_2new[w]**2)**(1/2))
        
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlim3d(-R,R)
    ax.set_ylim3d(-R,R)
    ax.set_zlim3d(-R,R)
    ax.set_zlabel('z [m]')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title('Situazione finale')

    ax.scatter3D(x_2new,y_2new,z_2new,color='red',s=1)
    plt.savefig('hernquist_finale.png')
    plt.show()  


# In[23]:


Mtot = 1 #* 2*10**30
rreale = np.linspace(0,R,30000)
integrale=(integrate.simps((densità2*rreale**2), rreale))

plt.hist(r2new[6,:],1000,density=True,label='campionamento finale')

print (R)
plt.plot(rreale,(((densità2)*rreale**2)/normalizzazione10),label='densità$_{teorica}$')
plt.xlim(0,R)
plt.xlabel('r [m]')
plt.ylabel('densità per r$^{2}$')
plt.title('Campionamento finale')
plt.legend()
plt.savefig('evoluzionefinaleperturber.png')


# In[24]:


r2new.sort()
tMy = t#/(1*10**6*365*24*60*60)
r2newpc = r2new#/( 3 *10**16)
plt.plot(tMy,r2newpc[:,393],label='r10')
plt.plot(tMy,r2newpc[:,676],label='r30')
plt.plot(tMy,r2newpc[:,663],label='r60')
plt.plot(tMy,r2newpc[:,2780],label='r90')
plt.legend()
plt.xlabel('t [My]')
plt.ylabel('r [pc]')
plt.show()


# In[25]:


Mr = np.zeros([Nsnapshot,Npart])
r2new.sort()
for u in range (Nsnapshot):
    for i in range (Npart):
            Mr[u,i] = Mtot * r2new[u,i]*3 * (r2new[u,i]**2+a**2)**(-(3/2))
            if 0.08<Mr[u,i]/Mtot<0.12:
                print ('r10')
                print (Mr[u,i], u, i)
            if 0.29<Mr[u,i]/Mtot<0.31:
                print ('r30')
                print (Mr[u,i], u, i)
            if 0.59<Mr[u,i]/Mtot<0.61:
                print ('r60')
                print (Mr[u,i], u, i)
            if 0.89<Mr[u,i]/Mtot<0.91:
                print ('r90')
                print (Mr[u,i], u, i)
                
                
        


# In[26]:


####INSERISCO IL PERTURBER
G = 1
a = 1
Npart = 30000
m_perturber=0.03
print (m_perturber)
Mr = np.zeros(Npart)
R_peri_i = (((0.5**(-2/3))-1))**(-0.5)*a
R_apo_i = R_peri_i/3
v_iniziale_q = (G/4)*((1/(np.sqrt(R_apo_i**2+a**2))-(1/(np.sqrt(R_peri_i**2+a**2)))))
print(v_iniziale_q)
v_i = np.sqrt(v_iniziale_q)
print (v_i)
print (R_apo_i,R_peri_i)

###CORDINATE E VECLOCITà PERTURBER
xp = R_peri_i
yp = 0
zp = 0

vxp = 0
vyp = v_i
vzp = 0


# In[29]:


####SALVO LE CONDIZIONI INIZIALI

fname='plummercondizioni1.txt' 
fp = open(fname,'w') 
fp.write('30001') 
fp.write('\n') 
fp.write('3') 
fp.write('\n') 
fp.write('0') 
fp.write('\n') 
fp.close()

np.savetxt('m_perturber.txt',np.column_stack([m_perturber]))
np.savetxt('coordinate_perturber.txt',np.column_stack([xp,yp,zp]))
np.savetxt('perturber_velocità.txt',np.column_stack([vxp,vyp,vzp]))


filenames = ['plummercondizioni1.txt','m_perturber.txt','plummerperturbermasseprova.txt','coordinate_perturber.txt', 'plummerperturbercoordinateprova.txt','perturber_velocità.txt', 'plummerperturbervelocitàprova.txt']
with open('perturber_inizio_03_160.txt', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)


# In[34]:


####ANALIZZO I MIEI RISULTATI
import pandas as pd
from scipy import integrate
import math
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product,combinations

Npart = 30001

fp=open('perturber_fine_03_360.txt','r')
Nlines=len(fp.readlines())
Nsnapshot=(Nlines)//(3*Npart+3)
print (Nlines,Nsnapshot)
t = np.zeros(Nsnapshot)
pos03 = np.zeros([Nsnapshot,Npart,3])
vf = np.zeros([Nsnapshot,Npart,3])
Npart = 30001
mf = np.zeros([Nsnapshot,Npart])

fp=open('perturber_fine_01_361.txt','r')
Nlines=len(fp.readlines())
Nsnapshot=(Nlines)//(3*Npart+3)
print (Nlines,Nsnapshot)
t = np.zeros(Nsnapshot)
pos01 = np.zeros([Nsnapshot,Npart,3])
vf = np.zeros([Nsnapshot,Npart,3])
Npart = 30001
mf = np.zeros([Nsnapshot,Npart])


i=0

for i in range (Nsnapshot):
    initial_line=i*(3*Npart+3)   #è il numero di linee che passano tra uno snapshot e l'altro
    
    #mf[i]=np.genfromtxt('perturberplummer_fine.txt',skip_header=initial_line+3,skip_footer=Nlines-initial_line-4)
    t[i]=np.genfromtxt('perturber_fine_03_360.txt',skip_header=initial_line+2, skip_footer=Nlines-initial_line-3)
    pos03[i]=np.genfromtxt('perturber_fine_03_360.txt',skip_header=initial_line+30004, skip_footer=Nlines-initial_line-60005)
    #vf[i]=np.genfromtxt('perturber_fine62.txt',skip_header=initial_line+60005, skip_footer=Nlines-initial_line-90006)
    #print(t[i],mf[i],pos[i],vf[i])
    print (i)
    initial_line=i*(3*Npart+3)   #è il numero di linee che passano tra uno snapshot e l'altro
    
    #mf[i]=np.genfromtxt('perturberplummer_fine.txt',skip_header=initial_line+3,skip_footer=Nlines-initial_line-4)
    t[i]=np.genfromtxt('perturber_fine_01_361.txt',skip_header=initial_line+2, skip_footer=Nlines-initial_line-3)
    pos01[i]=np.genfromtxt('perturber_fine_01_361.txt',skip_header=initial_line+30004, skip_footer=Nlines-initial_line-60005)
    #vf[i]=np.genfromtxt('perturber_fine62.txt',skip_header=initial_line+60005, skip_footer=Nlines-initial_line-90006)
    #print(t[i],mf[i],pos[i],vf[i])
    print (i)


# In[21]:


r_bh = np.zeros(Nsnapshot)
x_bh = np.zeros(Nsnapshot)
y_bh = np.zeros(Nsnapshot)
z_bh = np.zeros(Nsnapshot)
r_bh03 = np.zeros(Nsnapshot)
x_bh03 = np.zeros(Nsnapshot)
y_bh03 = np.zeros(Nsnapshot)
z_bh03 = np.zeros(Nsnapshot)
for u in range(Nsnapshot):
    x_bh[u] = pos01[u,0,0]
    y_bh[u] = pos01[u,0,1]
    z_bh[u] = pos01[u,0,2]
    r_bh[u] = np.sqrt(x_bh[u]**2+y_bh[u]**2+z_bh[u]**2)
    #print (u)
    
    x_bh03[u] = pos03[u,0,0]
    y_bh03[u] = pos03[u,0,1]
    z_bh03[u] = pos03[u,0,2]
    r_bh03[u] = np.sqrt(x_bh03[u]**2+y_bh03[u]**2+z_bh03[u]**2)
    #print (pos[u,0,0],pos[u,0,1],pos[u,0,2])

plt.plot(t*conversione_tempo/(10**6*365*24*60*60),r_bh,label='0.01 Mtot')
plt.axhline(2*0.05,color = 'red',ls='--',label='risoluzione')
plt.plot(t*conversione_tempo/(10**6*365*24*60*60),r_bh03,label='0.03 Mtot')
plt.axvline(298,color= 'black',ls = '-.',label='$t_{df[01]}$')
plt.axvline(173,color= 'grey',ls = '-.',label='$t_{df[03]}$')
plt.xlabel('t [My]')
plt.ylabel('r [pc]')
plt.legend()
plt.savefig('raggio_bh.png')
plt.show()

plt.plot(x_bh,y_bh,label='0.01 Mtot')
plt.plot(x_bh,y_bh03, label='0.03 Mtot')
plt.xlabel('x [pc]')
plt.ylabel('y [pc]')
plt.legend()
#plt.savefig('percorso_bh.png')
plt.show()


# In[4]:


#####SCRIVO I PARAMETRI DI CONVERSIONE
conversione_tempo = 4.7*10**14
conversione_velocità = 0.065

#tdyn_pu = tdyn*conversione_tempo
My = 10**6*365*24*60*60

print ( conversione_velocità)
#print (tdyn_pu/My)


# In[133]:


####FACCIO IL ROTATORE
#CAMBIO VELOCITà
m = 1/30000
L0 = xc*v_y*0.067*m*1.99*10**30 + yc*v_x*0.067*m*1.99*10**30
plt.hist(L0,1000)
plt.xlabel('$\ell$ [pc$^2$ kg / My]')
plt.savefig('Lz_hist.png')
plt.show()
ves1 = np.copy(ves)
for i in range (Npart-1):
    if (xc[i]*v_y[i]*m+yc[i]*v_x[i]*m)<=0:
        ves1[i] = -ves[i]
    else:
        ves[i] = ves[i]
        
v_x1= ves1 * np.sin(theta1)*np.cos(phi1)
v_y1= ves1 * np.sin(theta1)*np.sin(phi1)
v_z1= ves1 * np.cos(theta1)

print (v_x1,v_y1,v_z1)

L1 = xc*v_y1*0.067*m*1.99*10**30 + yc*v_x1*0.067*m*1.99*10**30
plt.hist(L1,1000)
plt.xlabel('$\ell$ [pc$^2$ kg / My]')
plt.savefig('Lz_positivi.png')
plt.show()

np.savetxt('plummerperturbervelocitàmodificate.txt',np.column_stack([v_x1,v_y1,v_z1]))
####SALVO LE CONDIZIONI INIZIALI

fname='plummercondizioni1.txt' 
fp = open(fname,'w') 
fp.write('30001') 
fp.write('\n') 
fp.write('3') 
fp.write('\n') 
fp.write('0') 
fp.write('\n') 
fp.close()

m_perturber = 0.03*Mtot
np.savetxt('m_perturber.txt',np.column_stack([m_perturber]))
np.savetxt('coordinate_perturber.txt',np.column_stack([xp,yp,zp]))
np.savetxt('perturber_velocità.txt',np.column_stack([vxp,vyp,vzp]))


filenames = ['plummercondizioni1.txt','m_perturber.txt','plummermassep.txt','coordinate_perturber.txt', 'plummercoordinatep.txt','perturber_velocità.txt', 'plummerperturbervelocitàmodificate.txt']
with open('rotatore_inizio_01_60.txt', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)


# In[1]:


####ANALIZZO I MIEI RISULTATI
import pandas as pd
from scipy import integrate
import math
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product,combinations

Npart = 30001

fp=open('rotatore_fine_03_199.txt','r')
Nlines=len(fp.readlines())
Nsnapshot=(Nlines)//(3*Npart+3)
print (Nlines,Nsnapshot)
t = np.zeros(Nsnapshot)
posr = np.zeros([Nsnapshot,Npart,3])
vf = np.zeros([Nsnapshot,Npart,3])
Npart = 30001
mf = np.zeros([Nsnapshot,Npart])
fp=open('rotatore_fine_01_199.txt','r')
Nlines=len(fp.readlines())
Nsnapshot=(Nlines)//(3*Npart+3)
print (Nlines,Nsnapshot)
t = np.zeros(Nsnapshot)
pos01r = np.zeros([Nsnapshot,Npart,3])
vf = np.zeros([Nsnapshot,Npart,3])
Npart = 30001
mf = np.zeros([Nsnapshot,Npart])

i=0

for i in range (Nsnapshot):
    initial_line=i*(3*Npart+3)   #è il numero di linee che passano tra uno snapshot e l'altro
    
    #mf[i]=np.genfromtxt('perturberplummer_fine.txt',skip_header=initial_line+3,skip_footer=Nlines-initial_line-4)
    t[i]=np.genfromtxt('rotatore_fine_03_199.txt',skip_header=initial_line+2, skip_footer=Nlines-initial_line-3)
    posr[i]=np.genfromtxt('rotatore_fine_03_199.txt',skip_header=initial_line+30004, skip_footer=Nlines-initial_line-60005)
    #vf[i]=np.genfromtxt('perturber_fine62.txt',skip_header=initial_line+60005, skip_footer=Nlines-initial_line-90006)
    #print(t[i],mf[i],pos[i],vf[i])
    print (i)
    initial_line=i*(3*Npart+3)   #è il numero di linee che passano tra uno snapshot e l'altro
    
    #mf[i]=np.genfromtxt('perturberplummer_fine.txt',skip_header=initial_line+3,skip_footer=Nlines-initial_line-4)
    t[i]=np.genfromtxt('rotatore_fine_01_199.txt',skip_header=initial_line+2, skip_footer=Nlines-initial_line-3)
    pos01r[i]=np.genfromtxt('rotatore_fine_01_199.txt',skip_header=initial_line+30004, skip_footer=Nlines-initial_line-60005)
    #vf[i]=np.genfromtxt('perturber_fine62.txt',skip_header=initial_line+60005, skip_footer=Nlines-initial_line-90006)
    #print(t[i],mf[i],pos[i],vf[i])
    print (i)


# In[6]:


r_bhr = np.zeros(Nsnapshot)
x_bhr = np.zeros(Nsnapshot)
y_bhr = np.zeros(Nsnapshot)
z_bhr = np.zeros(Nsnapshot)
r_bh03r = np.zeros(Nsnapshot)
x_bh03r = np.zeros(Nsnapshot)
y_bh03r = np.zeros(Nsnapshot)
z_bh03r = np.zeros(Nsnapshot)
for u in range(Nsnapshot):
    x_bhr[u] = pos01r[u,0,0]
    y_bhr[u] = pos01r[u,0,1]
    z_bhr[u] = pos01r[u,0,2]
    r_bhr[u] = np.sqrt(x_bhr[u]**2+y_bhr[u]**2)
    #print (u)
    
    x_bh03r[u] = posr[u,0,0]
    y_bh03r[u] = posr[u,0,1]
    z_bh03r[u] = posr[u,0,2]
    r_bh03r[u] = np.sqrt(x_bh03r[u]**2+y_bh03r[u]**2)


plt.plot(t*conversione_tempo/(10**6*365*24*60*60),r_bhr,label='0.01 Mtot')
plt.axhline(2.5*0.05,color = 'red',ls='--',label='risoluzione')
plt.axvline(538, color = 'grey',ls = '-.',label = '$t_{df[03]}$')
plt.plot(t*conversione_tempo/(10**6*365*24*60*60),r_bh03r,label='0.03 Mtot')
plt.xlabel('t [My]')
plt.ylabel('r [pc]')
#plt.savefig('raggio_bhr.png')
plt.legend()
plt.show()

plt.plot(x_bhr,y_bhr,label='0.01 Mtot')
plt.plot(x_bhr,y_bh03r, label='0.03 Mtot')
plt.xlabel('x [pc]')
plt.ylabel('y [pc]')
#plt.savefig('percorso_bhr.png')
plt.show()


# In[59]:


print(24*conversione_tempo/(10**6*365*24*60*60))


# In[ ]:




