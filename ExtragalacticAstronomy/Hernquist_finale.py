#!/usr/bin/env python
# coding: utf-8

# In[70]:


import pandas as pd
from scipy import integrate
import math
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product,combinations


#####PER CAMPIONARE, IN QUESTO CASO GENERO UN VETTORE CON I RAGGI, A CUI ASSOCIO LA DENSITà E IL POTENZIALE

M_sole = 2*10**30
M_sfera = 1 #* M_sole
M_bh = 0.5 * M_sfera
Npart = 10000
rh = 0.1 #* 3*10**(16)
Rmax = 4 #* 3*10**(16)

#GENERO IL VETTORE CON RAGGI CHE VANNO DA UNO PIù PICCOLO DEL RAGGIO SCALA FINO A PIù DEL RAGGIO MASSIMO, CON UN BINNAGGIO MAGGIORE DI NPART
R = np.linspace(0.000001,5,50000)
R = list(R)
R.sort()
R.reverse()
R = np.array(R)


#ora voglio la densità e il potenziale per questi raggi
rho =((M_sfera /(2*np.pi) * rh/R * 1/(R+rh)**3)*R**2)


rho10 = np.zeros(len(R))
for i in range (len(R)):
    rho10[i] = math.log10(rho[i])

plt.plot(R,rho)
plt.show()

#plt.plot(R,(M_sfera*rh/(2*np.pi) * 1/(R*(1+R/rh)**3)*R**2))
#plt.show()

P1 = np.random.random(Npart)
P2 = np.random.random(Npart)
#campiono le particelle per raggio
i = 0
rhistor = np.linspace(0,4,10000)
rhisto=[]
while i < Npart:
    x1 = float(np.random.random(1)*4)
    y1 = np.random.random(1)*0.3
    if y1 <= ((M_sfera/(2*np.pi)) * (rh/x1) * (1/(x1+rh)**3) *x1**2):
        rhisto.append(x1)
        i = i + 1

R1 = np.linspace(0,Rmax,10000)

plt.hist(rhisto,100,density=True,label = 'Campionamento')
plt.plot(R1,((M_sfera*rh/(2*np.pi) * 1/(R1*(1+R1/rh)**3))*R1**2)/0.0000764874, label='Andamento teorico')
plt.xlabel('r [m]')
plt.ylabel('$\\rho$(r) * r$^{2}$')
plt.legend()
plt.savefig('hernquistiniziale_raggi.png')
plt.show()


##campiono theta
theta= np.arccos(2*P2-1)
plt.hist(theta,30,density=True,label='campionamento')
plt.xlabel('$\\theta$ $_{radianti}$')
plt.ylabel('p($\\theta$)')
thetareale=np.linspace(0,np.pi,10000)
i1=(integrate.simps((np.sin(thetareale)),thetareale))
plt.plot(thetareale,(np.sin(thetareale))/i1, label='p($\\theta$)$_{teorico})$')

plt.legend()
plt.title('$\\theta$$_{r}$')
plt.savefig('hernquistiniziale_theta.png')
plt.show()

#campiono phi
phi=2*np.pi*P1
normalizzazione=np.empty(len(phi))
for i in range (len(phi)):
    normalizzazione[i]=1/(2*np.pi)
plt.hist(phi,30,density=True, label='campionamento')
plt.plot(phi,normalizzazione,label='p($\phi$)$_{teorico})$')
plt.xlabel('$\phi$ $_{radianti}$')
plt.ylabel('p($\phi$)')

plt.legend()
plt.title('$\phi$$_{r}$')
plt.savefig('hernquistiniziale_phi.png')
plt.show()


# In[ ]:





# In[71]:


#SITUAZIONE INIZIALE
#cambio in coordinate cartesiane
x= rhisto * np.sin(theta)*np.cos(phi)
y= rhisto * np.sin(theta)*np.sin(phi)
z= rhisto * np.cos(theta)

#grafico la mia situazione iniziale

fig = plt.figure()
ax = fig.gca(projection='3d')
zdata=z
ydata=y
xdata=x
ax.set_xlim3d(-4.500,4.500)
ax.set_ylim3d(-4.500,4.500)
ax.set_zlim3d(-4.500,4.500)
ax.set_zlabel('z [m]')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title('Condizione iniziale')
ax.scatter3D(xdata,ydata,zdata,s=1, c = 'r')
ax.scatter3D(0,0,0, s=1, c = 'b')
plt.savefig('hernquistiniziale.png')


# In[32]:


#GENERO ALLO STESSO MODO IL POTENZIALE
pot=np.zeros(len(R))
G= 1 


rhisto = np.array(rhisto)
psi = G*M_sfera/(rhisto+rh) + G*M_bh/rhisto
psi10 = np.zeros(len(rhisto))
for i in range (len(rhisto)):
    psi10[i]= math.log10(psi[i])


plt.plot(rhisto,psi10, '.')
plt.show()


####è TUTTO UN GROSSO PUNTO DI DOMANDA
#genero il vettore potenziale
pot = -G*M_sfera/(R+rh) - G*M_bh/R

pot10 = np.zeros(len(R))
for i in range (len(R)):
    pot10[i] = math.log10(-pot[i])
    
plt.plot(R,pot10)
plt.show()



# In[33]:


#voglio trovare il vettore g(epsilon)


deltarho = np.zeros(len(R))
for i in range (len(R)):
    if i == 0:
        deltarho[i]= rho[i+1]-rho[i]
    elif i != 0 and i != (len(R)-1) :
        deltarho[i]= (rho[i+1]-rho[i-1])/2
    elif i == (len(R)-1):
        deltarho[i]= (rho[i]-rho[i-1])




# In[34]:


####PARALLELIZZO PER VELOCIZZARE

from queue import Queue
from threading import Thread


g = np.zeros(len(R))


def task(j,g,deltarho,pot):
    g[j]=np.sum(deltarho[0:(j)]/(np.sqrt(2*(-pot[j]+pot[0:(j)]))))
    
    
def execute (q):
    while True:
        try:
            j,g,deltarho,pot= q.get()
            task(j,g,deltarho,pot)
            q.task_done()
        except Queue.Empty:
            continue
            
q = Queue(maxsize=0)
for i in range (100):
    t = Thread(target = execute, args = (q,))
    #t.setDeamon(True)
    t.start()
    
for j in range (len(R)):
    q.put((j,g,deltarho,pot))
q.join()


# In[122]:





# In[21]:



eps = np.zeros(len(R))
for i in range (len(R)):
    eps[i]=-pot[i]



f = np.zeros(len(R))
f10 = np.zeros(len(R))
for i in range (len(R)):
    if i == 0:
        f[i]= (g[i+1]-g[i])/(eps[i+1]-eps[i])
    elif i != 0 and i != (len(R)-1) :
        f[i]= (g[i+1]-g[i-1])/(2*(eps[i+1]-eps[i-1]))
    elif i == (len(R)-1):
        f[i]= (g[i]-g[i-1])/(eps[i]-eps[i-1])


# In[18]:


print (R[19999],eps[19999], (np.sqrt(2*eps[19999])))


# In[23]:


plt.loglog()
plt.plot(eps,f)
plt.xlim(0.1,3500)
plt.show()

print (eps[0])


# In[56]:


#CALCOLO LA CUMULATIVA F(EPS). DEVO FARE L'INTEGRALE DA 0 A RAD(2PSI) DI f(EPS) IN dEPS
F = np.zeros(len(R))
a = np.zeros(len(R))
b = np.zeros(len(R))
c = np.zeros(len(R))
cnum = np.zeros(len(R))
cden = c = np.zeros(len(R))
F10 = np.zeros(len(R))
for i in range (len(R)):
    F[i]=(np.sqrt(2)/(32*np.pi**3))*(-((17-29*M_bh+14*M_bh**2-3*eps[i]+3*M_bh*eps[i]-2*eps[i]**2)*np.sqrt(eps[i]))/((2-M_bh+eps[i]))**3+2/np.sqrt(M_bh)*np.arctan(np.sqrt(M_bh*eps[i]))+((1-M_bh)*(1-4*M_bh*-2*M_bh**2-16*eps[i]-17*M_bh*eps[i]+4*eps[i]**2))/(2-M_bh+eps[i])**(7/2)*np.arctanh(np.sqrt((2-M_bh+eps[i])*eps[i])/(1+eps[i])))
    F10[i] = math.log10(F[i])
   


# In[119]:


plt.plot(eps10,F10)


# In[21]:


pot_iesimo = np.zeros(len(rhisto))

pot_iesimo10 = np.zeros(len(rhisto))
F_iesimo = np.zeros(len(rhisto))
rhisto.sort()
rhisto = list(rhisto)
rhisto.reverse()
rhisto = np.array(rhisto)

for i in range (len(rhisto)):
    pot_iesimo[i]= -G*M_sfera/(rhisto[i]+rh) - G*M_bh/rhisto[i]
    pot_iesimo10[i] = math.log10(-pot_iesimo[i])


    
        
plt.plot(rhisto,pot_iesimo10)
#plt.plot(R,pot10)


# In[38]:


epsnuovo = 1/(R+1) + (M_bh/R)
Fsimulato = np.zeros(len(R))

for i in range (len(R)-1):
    if i == 0:
        Fsimulato[i] += 4*np.pi*f[i]*np.sqrt(-2*pot[i])*(-pot[i]-epsnuovo[i])
    else:
        Fsimulato[i] += 4*np.pi*f[i]*np.sqrt(2*(-pot[i+1]+pot[i]))*(-pot[i]-epsnuovo[i])
    


# In[45]:


print (Fsimulato,epsnuovo)
Fsimulato.sort()
plt.loglog()
plt.plot(eps,Fsimulato)
plt.xlim(0.2,1000000)


# In[46]:


fv_rid = np.zeros([Npart])
v = np.zeros([Npart])
for i in range (Npart):
    for k in range(i):
        if -pot[i] >= epsnuovo[k]:
            fv_rid[i] = np.random.uniform(0,1)*Fsimulato[k-1]
    


# In[47]:





for i in range (Npart):
    if fv_rid[i] < Fsimulato[i]:
        v[i] = np.sqrt(2*np.abs(pot[i+1]-pot[i]))

    print (v[i],i)


# In[127]:


v_circolare = np.zeros([Npart])
G = 1
for i in range (Npart):
    v_circolare[i] = np.sqrt((G*M_sfera*rhisto[i]) / ((rhisto[i]+rh)**2))
    print (v_circolare[i])


# In[31]:




######CAMPIONO LE VELOCITA'
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
plt.savefig('hernquistinizialev_theta')
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
plt.savefig('hernquistinizialev_phi.png')
plt.show()


plt.plot(rhisto,v)
plt.show()
'''
v_x= v * np.sin(theta1)*np.cos(phi1)
v_y= v * np.sin(theta1)*np.sin(phi1)
v_z= v * np.cos(theta1)
'''
v_x1= v * np.sin(theta1)*np.cos(phi1)
v_y1= v * np.sin(theta1)*np.sin(phi1)
v_z1= v * np.cos(theta1)

plt.plot(rhisto,v_x1)
#plt.plot(rhisto,v_y1)
#plt.plot(rhisto,v_z1)
plt.show()
plt.hist(v,100,density=True)


# In[28]:


'''
fig = plt.figure()
ax = fig.gca(projection='3d')

ax.set_xlim3d(-0.005,0.005)
ax.set_ylim3d(-0.005,0.005)
ax.set_zlim3d(-0.005,0.005)
ax.set_zlabel('z')
ax.set_xlabel('x')
ax.set_ylabel('y')
#ax.set_title('Condizione iniziale')
ax.scatter3D(v_x,v_y,v_z,s=0.1, c = 'r')

#plt.savefig('plummeriniziale.png')
'''
fig = plt.figure()
ax = fig.gca(projection='3d')

ax.set_xlim3d(-0.005,0.005)
ax.set_ylim3d(-0.005,0.005)
ax.set_zlim3d(-0.005,0.005)
ax.set_zlabel('z')
ax.set_xlabel('x')
ax.set_ylabel('y')
#ax.set_title('Condizione iniziale')
ax.scatter3D(v_x1,v_y1,v_z1, c = 'r')

plt.savefig('hernquistinizialev.png')


# In[27]:


#####GENERO IL FILE TXT
m = M_sfera / Npart
masse = np.zeros(Npart)
i = 0
for i in range (Npart):
    masse[i] = m
#print (rhisto,x,y,z,v_x1,v_y1,v_z1)

np.savetxt('hernquistmasseF.txt',np.column_stack([masse]))
np.savetxt('hernquistcoordinateF.txt',np.column_stack([x,y,z]))
#np.savetxt('hernquistvelocità3.txt',np.column_stack([v_x1,v_y1,v_z1]))


# In[54]:


####ANALIZZO I MIEI RISULTATI
import pandas as pd
from scipy import integrate
import math
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product,combinations

Npart = 10000

fp=open('hernquist_fine30.txt','r')
Nlines=len(fp.readlines())
Nsnapshot=(Nlines)//(3*Npart+3)
print (Nlines,Nsnapshot)
t = np.zeros(Nsnapshot)
pos = np.zeros([Nsnapshot,Npart,3])
vf = np.zeros([Nsnapshot,Npart,3])
Npart = 10000
mf = np.zeros([Nsnapshot,Npart])

i=0

for i in range (Nsnapshot):
    initial_line=i*(3*Npart+3)   #è il numero di linee che passano tra uno snapshot e l'altro
    mf[i]=np.genfromtxt('hernquist_fine30.txt',skip_header=initial_line+3,skip_footer=Nlines-initial_line-4)
    t[i]=np.genfromtxt('hernquist_fine30.txt',skip_header=initial_line+2, skip_footer=Nlines-initial_line-3)
    pos[i]=np.genfromtxt('hernquist_fine30.txt',skip_header=initial_line+10003, skip_footer=Nlines-initial_line-20003)
    vf[i]=np.genfromtxt('hernquist_fine30.txt',skip_header=initial_line+20003, skip_footer=Nlines-initial_line-30003)
    #print(t[i],mf[i],pos[i],vf[i])


# In[ ]:


z_2=np.zeros(Npart)
x_2=np.zeros(Npart)
y_2=np.zeros(Npart)
r2 = np.zeros([Nsnapshot,Npart])

for u in range (Nsnapshot):
    for i in range (Npart):
        z_2[i]=pos[u,i,2]
        x_2[i]=pos[u,i,0]
        y_2[i]=pos[u,i,1]
    for w in range (Npart):
        r2[u,w]=((x_2[w]**2+y_2[w]**2+z_2[w]**2)**(1/2))

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlim3d(-4,4)
    ax.set_ylim3d(-4,4)
    ax.set_zlim3d(-4,4)
    ax.set_zlabel('z')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    #ax.set_title('5 t$_{dyn}$')

    ax.scatter3D(x_2,y_2,z_2,color='red', s = 0.1)
    #plt.savefig('2tempo_collasso.png')
    plt.show()  


plt.plot(t,r2)


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlim3d(-10,10)
ax.set_ylim3d(-10,10)
ax.set_zlim3d(-10,10)
ax.set_zlabel('z')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Situazione finale')

ax.scatter3D(x_2,y_2,z_2,color='red', s = 0.1)
#plt.savefig('situazione_finale.png')
plt.show()  


# In[21]:


#calcolo softening
eps = (m / (M_sfera*rh/(2*np.pi) * 1/(rh*(1+rh/rh)**3)))**(1/3)
print (eps)


# In[ ]:


#calcolo tdyn
tdyn = np.sqrt((3*np.pi)/(16*G*(M_sfera*rh/(2*np.pi) * 1/(rh*(1+rh/rh)**3)))
print ()


# In[73]:


#####CALCOLO IL CENTRO DI MASSA PERCHè IL MIO SISTEMA SI STA SPOSTANDO PROBABILMENTE,E VERIFICO
M_sfera = 1
M_bh = 0.5
Mtot=M_sfera+M_bh
x_cm = np.zeros(Nsnapshot)
y_cm = np.zeros(Nsnapshot)
z_cm = np.zeros(Nsnapshot)
r_cm = np.zeros(Nsnapshot)

xnew = np.zeros([Nsnapshot,Npart])
ynew = np.zeros([Nsnapshot,Npart])
znew = np.zeros([Nsnapshot,Npart])
r2new = np.zeros([Nsnapshot,Npart])

for u in range (Nsnapshot):
    x_cm[u] = (sum(Mtot/Npart * pos[u,:,0]))/(Mtot)
    y_cm[u] = (sum(Mtot/Npart * pos[u,:,1])/Mtot)
    z_cm[u] = (sum(Mtot/Npart * pos[u,:,2])/Mtot)
    print (x_cm[u],y_cm[u],z_cm[u])
    
        
   
   
plt.plot(t,x_cm, c = 'r',label='x_cm')
plt.plot(t,y_cm, c = 'b',label= 'y_cm')
plt.plot(t,z_cm, c = 'g',label = 'z_cm')
plt.xlabel('t [s]')
plt.ylabel('m')
plt.legend()
plt.savefig('movcentrodimassa_h.png')


# In[85]:


x_2new=np.zeros([Npart])
y_2new=np.zeros([Npart])
z_2new=np.zeros([Npart])


for u in range (Nsnapshot):
    for i in range (Npart):
        z_2new[i]=pos[u,i,2]-z_cm[u]
        x_2new[i]=pos[u,i,0]-x_cm[u]
        y_2new[i]=pos[u,i,1]-y_cm[u]
    for w in range (Npart):
        r2new[u,w]=((x_2new[w]**2+y_2new[w]**2+z_2new[w]**2)**(1/2))
        
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlim3d(-4,4)
    ax.set_ylim3d(-4,4)
    ax.set_zlim3d(-4,4)
    ax.set_zlabel('z [m]')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title('Situazione finale')

    ax.scatter3D(x_2new,y_2new,z_2new,color='red',s=1)
    plt.savefig('hernquist_finale.png')
    plt.show()  
print (r2new,t)


# In[75]:


print (r2new)
plt.plot(t,r2new)
plt.xlabel('t [s]')
plt.ylabel('r [m]')
plt.savefig('raggineltempo.png')
plt.show()


# In[27]:


###calcolo softening
m = M_sfera/Npart
eps = (m/(M_sfera*rh/(2*np.pi) * 1/(rh*(1+rh/rh)**3) ))**(1/3)
print (eps)


# In[7]:


###calcolo tdyn
G=6.67*10**(-11)
M_sfera = 1.5 * 2 * 10**(30)
tdyn = (np.sqrt((3*np.pi)/(16*G*(M_sfera*rh/(2*np.pi) * 1/(4*(1+4/rh)**3)))))
print (tdyn)


# In[3]:


M = np.zeros([Nsnapshot,Npart])
r_90 = np.zeros(Nsnapshot)
r_60 = np.zeros(Nsnapshot)
r_30 = np.zeros(Nsnapshot)
r_15 = np.zeros(Nsnapshot)
r_bh = np.zeros(Nsnapshot)
r2new.sort()
rh = 0.1

for u in range (Nsnapshot):
    for i in range (Npart):
        M[u,i] = (M_sfera * (r2new[u,i]/(r2new[u,i]+rh))**2)
        if 0.14 < M[u,i]/M_sfera < 0.16:
            r_15[u] = r2new[u,i]
        elif 0.29< M[u,i]/M_sfera < 0.31:
            r_30[u] = r2new[u,i]
        elif 0.59 <M[u,i]/M_sfera <0.61:
            r_60[u] = r2new[u,i]
        elif 0.89 <M[u,i]/M_sfera < 0.91:
            r_90[u] = r2new[u,i]
        r_bh[u] = r2new[u,0]
print (Mtot)
        
plt.plot(t,r_15, label = 'r$_{lagrangiano15}$')
plt.plot(t,r_30, label = 'r$_{lagrangiano30}$')
plt.plot(t,r_60, label = 'r$_{lagrangiano60}$')
plt.plot(t,r_90, label = 'r$_{lagrangiano90}$')
plt.plot(t,r_bh, label = 'r$_{bh}$')
plt.xlabel('t [s]')
plt.ylabel('r [m]')
plt.legend()
#plt.savefig('raggiobh.png')


# In[69]:


plt.hist(r2new[5,:],30,density = True,label = 'campionamento finale')
R1 = np.linspace(0,4,10000)
rh = 0.1
plt.plot(R1,((M_sfera*rh/(2*np.pi) * 1/(R1*(1+R1/rh)**3))*R1**2)/0.0000764874,label = 'andamento teorico')
plt.xlabel('r [m]')
plt.ylabel('$\\rho$(r) * r$^{2}$')
plt.legend()
plt.xlim(-0.01,4)
plt.savefig('rfinale.png')


# In[4]:


v_circolare2 = np.zeros([Npart])
G = 1
for i in range (Npart):
    v_circolare2[i] = np.sqrt((G*M_sfera*rhisto[i]) / ((rhisto[i]+rh)**2))
    print (v_circolare2[i])


# In[2]:


import pandas as pd
from scipy import integrate
import math
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product,combinations


# In[ ]:




