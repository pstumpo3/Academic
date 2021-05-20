#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product,combinations
get_ipython().run_line_magic('matplotlib', 'inline')

#estraggo coordinate casuali per i miei punti
Npart=100
#distribuzione di phi
P=np.random.random(Npart)
P2=np.random.random(Npart)
P3=np.random.random(Npart)

#distribuzione di phi

phi=2*np.pi*P
normalizzazione=np.empty(len(phi))
for i in range (len(phi)):
    normalizzazione[i]=1/(2*np.pi)
plt.hist(phi,10,density=True)
plt.plot(phi,normalizzazione)
plt.xlabel('$\phi$')
plt.ylabel('p($\phi$)')
#plt.savefig('distribuzione_phi.png')
plt.show()

#distibuzione theta
theta= np.arccos(2*P2-1)
plt.hist(theta,10,density=False)
#plt.plot(theta,np.sin(theta),'ro')
plt.xlabel('$\\theta$')
plt.ylabel('#corpi')
#plt.savefig('distribuzione_theta.png')
plt.show()
   
#distribuzione raggio
R_max=4
m=0.03
r=R_max*P2**(1/3)
plt.hist(r,10,density=False)
plt.xlabel('r')
plt.ylabel('#corpi')
#plt.savefig('distribuzione_r.png')
plt.show()
#verifico che la distibuzione nel raggio sia conforme teoria (sfera omogonea): 
#la densità deve rimanere costante per ogni guscio sfeirco
r.sort()
R=0
part_shell=[]
dr=[]
raggio=[]
rho_shell=[]
step=0.05
rho_tot=0

while R<R_max:
    npart = 0
    for i in range (len(r)):
        if R<r[i]<R+step:
            npart=npart + 1
        part_shell.append(npart)
        diffr=(R+step)**3-R**3
        raggio.append(R)
        dr.append(diffr)

    R = R + step

dr=np.array(dr)
part_shell=np.array(part_shell)
rho_shell= 3*part_shell*m/(4*np.pi)*(1/dr)
rho_0=np.empty(len(rho_shell))
for i in range(len(raggio)):
    rho_0[i]=3*m*Npart/(4*np.pi)*(1/R_max**3)
rho_tot=0
for i in range(len(rho_shell)):
    rho_tot=rho_tot+rho_shell[i]
rho_medio= rho_tot/(len(rho_shell))

plt.plot(raggio,rho_shell, label='$\\rho$(r)')
plt.plot(raggio,rho_0,c='red', label='$\\rho_0$')
plt.axhline(rho_medio,c='black', label='$\\rho$$^-$ calcolato con le shell')
plt.xlabel('r')
plt.ylabel('$\\rho$')
plt.legend()
#plt.savefig('densità_shell.png')
plt.show()

#il codice ha bisogno delle coordinate cartesiane masse e velocità, e un file txt
x= r * np.sin(theta)*np.cos(phi)
y= r * np.sin(theta)*np.sin(phi)
z= r * np.cos(theta)

v=np.zeros(Npart)
m_part=np.empty(Npart)
for i in range (0,Npart):
    m_part[i]=m

#np.savetxt('condizioni_inziali.txt',np.column_stack([m_part,x,y,z,v,v,v]))
    
#grafico la mia situazione iniziale


fig = plt.figure()
ax = fig.gca(projection='3d')
zdata=z
ydata=y
xdata=x
ax.set_xlim3d(-6,6)
ax.set_ylim3d(-6,6)
ax.set_zlim3d(-6,6)
ax.set_zlabel('z')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Condizione iniziale')

#ax.scatter3D(xdata,ydata,zdata,c='red')
plt.savefig('condizione.png')

plt.show()






# In[141]:


#grafico la mia situazione iniziale


fig = plt.figure()
ax = fig.gca(projection='3d')
zdata=z
ydata=y
xdata=x
ax.set_xlim3d(-6,6)
ax.set_ylim3d(-6,6)
ax.set_zlim3d(-6,6)
ax.set_zlabel('z')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Condizione iniziale')

#ax.scatter3D(xdata,ydata,zdata,c='red')
plt.savefig('condizione.png')

plt.show()


# In[2]:


#######EVOLVO TEMPORALMENTE
#devo leggere i miei dati e sistemarl in array

fp=open('fine.out.txt','r')
Nlines=len(fp.readlines())
Nsnapshot=(Nlines)//(Npart+2)
pos=np.empty([Nsnapshot,Npart,7])

t=np.empty(Nsnapshot)

for i in range (Nsnapshot):
    initial_line=i*(Npart+2)   #è il numero di linee che passano tra uno snapshot e l'altro
    t[i]=np.genfromtxt('fine.out.txt',skip_header=initial_line+1, skip_footer=Nlines-initial_line-2)
    pos[i]=np.genfromtxt('fine.out.txt',skip_header=initial_line+2, skip_footer=Nlines-initial_line-102)
    print (i)


# In[157]:


##vedo come varia il raggio dei miei corpi
raggio = np.empty(Nsnapshot)

T_collasso=np.sqrt((3*np.pi)/(32*(m*Npart)/(4/3*np.pi*R_max**3)))
for i in range (Npart):
    for n in range (Nsnapshot):
        raggio[n]= np.sqrt((pos[n,i,1])**2 + (pos[n,i,2])**2 + (pos[n,i,3])**2)
    plt.plot(t,raggio)
plt.axvline(T_collasso,color='black',linestyle='--')
plt.xlabel('t$_{I.U}$')
plt.ylabel('r$_{I.U.}$')
#lt.savefig('raggio(t).png')


# In[143]:


#voglio vedere anche qual è la situazione alla fine del collasso, e al doppio del tempo di collasso della mia sfera omogenea

z_collasso=np.empty(Npart)
x_collasso=np.empty(Npart)
y_collasso=np.empty(Npart)

for i in range (Npart):
    z_collasso[i]=pos[Nsnapshot//2,i,3]
    x_collasso[i]=pos[Nsnapshot//2,i,1]
    y_collasso[i]=pos[Nsnapshot//2,i,2]
      
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlim3d(-6,6)
ax.set_ylim3d(-6,6)
ax.set_zlim3d(-6,6)
ax.set_zlabel('z')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('t$_{collasso}$')

ax.scatter3D(x_collasso,y_collasso,z_collasso,color='red')
plt.savefig('tempo_collasso.png')
plt.show()

z_2collasso=np.empty(Npart)
x_2collasso=np.empty(Npart)
y_2collasso=np.empty(Npart)

for i in range (Npart):
    z_2collasso[i]=pos[Nsnapshot-1,i,3]
    x_2collasso[i]=pos[Nsnapshot-1,i,1]
    y_2collasso[i]=pos[Nsnapshot-1,i,2]
      
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlim3d(-6,6)
ax.set_ylim3d(-6,6)
ax.set_zlim3d(-6,6)
ax.set_zlabel('z')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('2 t$_{collasso}$')

ax.scatter3D(x_2collasso,y_2collasso,z_2collasso,color='red')
plt.savefig('2tempo_collasso.png')
plt.show()


# In[136]:


T_collasso=np.sqrt((3*np.pi)/(32*(m*Npart)/(4/3*np.pi*R_max**3)))
print (T_collasso)


# In[164]:


import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product,combinations
get_ipython().run_line_magic('matplotlib', 'inline')

#########rifaccio l'esercizio per più corpi
#estraggo coordinate casuali per i miei punti
Npart=10000
#distribuzione di phi
P=np.random.random(Npart)
P2=np.random.random(Npart)
P3=np.random.random(Npart)

#distribuzione di phi

phi=2*np.pi*P
normalizzazione=np.empty(len(phi))
for i in range (len(phi)):
    normalizzazione[i]=1/(2*np.pi)
plt.hist(phi,10,density=True)
plt.plot(phi,normalizzazione)
plt.xlabel('$\phi$')
plt.ylabel('p($\phi$)')
plt.savefig('distribuzione_phi2.png')
plt.show()

#distibuzione theta
theta= np.arccos(2*P2-1)
plt.hist(theta,10,density=False)
#plt.plot(theta,np.sin(theta),'ro')
plt.xlabel('$\\theta$')
plt.ylabel('#corpi')
plt.savefig('distribuzione_theta2.png')
plt.show()
   
#distribuzione raggio
R_max=4
m=3/Npart
r=R_max*P2**(1/3)
plt.hist(r,10,density=False)
plt.xlabel('r')
plt.ylabel('#corpi')
plt.savefig('distribuzione_r2.png')
plt.show()
#verifico che la distibuzione nel raggio sia conforme teoria (sfera omogonea): 
#la densità deve rimanere costante per ogni guscio sfeirco
r.sort()
R=0
part_shell=[]
dr=[]
raggio=[]
rho_shell=[]
step=0.05
rho_tot=0

while R<R_max:
    npart = 0
    for i in range (len(r)):
        if R<r[i]<R+step:
            npart=npart + 1
        part_shell.append(npart)
        diffr=(R+step)**3-R**3
        raggio.append(R)
        dr.append(diffr)

    R = R + step

dr=np.array(dr)
part_shell=np.array(part_shell)
rho_shell= 3*part_shell*m/(4*np.pi)*(1/dr)
rho_0=np.empty(len(rho_shell))
for i in range(len(raggio)):
    rho_0[i]=3*m*Npart/(4*np.pi)*(1/R_max**3)
rho_tot=0
for i in range(len(rho_shell)):
    rho_tot=rho_tot+rho_shell[i]
rho_medio= rho_tot/(len(rho_shell))

plt.plot(raggio,rho_shell, label='$\\rho$(r)')
plt.plot(raggio,rho_0,c='red', label='$\\rho_0$')
plt.axhline(rho_medio,c='black', label='$\\rho$$^-$ calcolato con le shell')
plt.xlabel('r')
plt.ylabel('$\\rho$')
plt.legend()
plt.savefig('densità_shell2.png')
plt.show()

#il codice ha bisogno delle coordinate cartesiane masse e velocità, e un file txt
x= r * np.sin(theta)*np.cos(phi)
y= r * np.sin(theta)*np.sin(phi)
z= r * np.cos(theta)

v=np.zeros(Npart)
m_part=np.empty(Npart)
for i in range (0,Npart):
    m_part[i]=m

np.savetxt('condizioni_inziali2.txt',np.column_stack([m_part,x,y,z,v,v,v]))
    
#grafico la mia situazione iniziale


fig = plt.figure()
ax = fig.gca(projection='3d')
zdata=z
ydata=y
xdata=x
ax.set_xlim3d(-6,6)
ax.set_ylim3d(-6,6)
ax.set_zlim3d(-6,6)
ax.set_zlabel('z')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Condizione iniziale')

ax.scatter3D(xdata,ydata,zdata,c='red')
plt.savefig('condizione2.png')

plt.show()


# In[ ]:




