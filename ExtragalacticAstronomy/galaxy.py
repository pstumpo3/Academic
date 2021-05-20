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


######APRO IL FILE DI TESTO DELLA GALASSIA

Npart = 20002
fp=open('SingleGalaxy.in.txt','r')
Nlines=len(fp.readlines())
Nsnapshot=(Nlines)//(3*Npart+3)
print (Nlines,Nsnapshot)
t = np.zeros(Nsnapshot)
pos = np.zeros([Nsnapshot,Npart,3])
vf = np.zeros([Nsnapshot,Npart,3])

mf = np.zeros([Nsnapshot,Npart])

i=0

for i in range (Nsnapshot):
    initial_line=i*(3*Npart+3)   #è il numero di linee che passano tra uno snapshot e l'altro
   
    mf[i]=np.genfromtxt('SingleGalaxy.in.txt',skip_header=initial_line+3,skip_footer=Nlines-initial_line-20005)
    t[i]=np.genfromtxt('SingleGalaxy.in.txt',skip_header=initial_line+2, skip_footer=Nlines-initial_line-3)
    pos[i]=np.genfromtxt('SingleGalaxy.in.txt',skip_header=initial_line+20005, skip_footer=Nlines-initial_line-40007)
    vf[i]=np.genfromtxt('SingleGalaxy.in.txt',skip_header=initial_line+40007, skip_footer=Nlines-initial_line-60009)
    #print(t[i],mf[i],pos[i],vf[i])


# In[3]:


#grafico la mia situazione iniziale


fig = plt.figure()
ax = fig.gca(projection='3d')
zdata=pos[0,:,2]
ydata=pos[0,:,1]
xdata=pos[0,:,0]
ax.set_xlim3d(-4,4)
ax.set_ylim3d(-4,4)
ax.set_zlim3d(-4,4)
ax.set_zlabel('z')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Condizione iniziale')
ax.scatter3D(xdata,ydata,zdata,c='red',s = 0.1)
plt.savefig('galassia.png')

#salvo il mio raggio lagrangiano iniziale



# In[4]:


####MODIFICO VELOCITà E POSIZIONI DELLE PARTICELLE e DUPLICO LA GALASSIA

#Galaxy 1
pos1= np.copy(pos)
v1 = np.copy(vf)

pos1[0,:,0] = pos[0,:,0] + 20 
v1[0,:,1] = v1[0,:,1] + 0.3
    
#Galaxy 2
pos2= np.copy(pos)
v2 = np.copy(vf)

pos2[0,:,0] = pos[0,:,0] - 20 
v2[0,:,1] = v1[0,:,1] - 0.3
    


# In[5]:


#grafico la mia situazione iniziale


fig = plt.figure()
ax = fig.gca(projection='3d')
zdata1=pos1[0,:,2]
ydata1=pos1[0,:,1]
xdata1=pos1[0,:,0]

zdata2=pos2[0,:,2]
ydata2=pos2[0,:,1]
xdata2=pos2[0,:,0]
ax.set_xlim3d(-24,24)
ax.set_ylim3d(-24,24)
ax.set_zlim3d(-24,24)
ax.set_zlabel('z [pc]')
ax.set_xlabel('x [pc]')
ax.set_ylabel('y [pc]')
#ax.set_title('Condizione iniziale')
ax.scatter3D(xdata2,ydata2,zdata2,c='blue',s = 0.01)
#plt.savefig('plummerinizialep.png')

#salvo il mio raggio lagrangiano iniziale

#ax.set_title('Condizione iniziale')
ax.scatter3D(xdata1,ydata1,zdata1,c='red',s = 0.01)
#plt.savefig('plummerinizialep.png')

#salvo il mio raggio lagrangiano iniziale



# In[6]:


#####SALVO LE MASSE,POSIZIONI E VELOCITà DELLE MIE DIVERSE COMPONENTI: BH, MATTER, DARK MATTER.
#GALAXY1
#BH
m_bh1 = mf[0,0]
x_bh1 = pos1[0,0,0]
y_bh1 = pos1[0,0,1]
z_bh1 = pos1[0,0,2]

#MATTER
m_m1 = mf[0,1:10002]
x_m1 = pos1[0,1:10002,0]
y_m1 = pos1[0,1:10002,1]
z_m1 = pos1[0,1:10002,2]


#DARK MATTER
m_dm1 = mf[0,10002:20002]
x_dm1 = pos1[0,10002:20002,0]
y_dm1 = pos1[0,10002:20002,1]
z_dm1 = pos1[0,10002:20002,2]

#GALAXY2
#BH
m_bh2 = mf[0,0]
x_bh2 = pos2[0,0,0]
y_bh2 = pos2[0,0,1]
z_bh2 = pos2[0,0,2]

#MATTER
m_m2 = mf[0,1:10002]
x_m2 = pos2[0,1:10002,0]
y_m2 = pos2[0,1:10002,1]
z_m2 = pos2[0,1:10002,2]

#DARK MATTER
m_dM2 = mf[0,10002:20002]
x_dm2 = pos2[0,10002:20002,0]
y_dm2 = pos2[0,10002:20002,1]
z_dm2 = pos2[0,10002:20002,2]

print(x_m2.shape)


# In[7]:


#grafico la mia situazione iniziale
get_ipython().run_line_magic('matplotlib', 'inline')

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.set_xlim3d(-30,30)
ax.set_ylim3d(-30,30)
ax.set_zlim3d(-30,30)
ax.set_zlabel('z [pc]')
ax.set_xlabel('x [pc]')
ax.set_ylabel('y [pc]')
ax.set_title('Condizione iniziale')

ax.scatter3D(x_m1,y_m1,z_m1,c='blue',s = 1,label='M$_{1}$')
ax.scatter3D(x_dm1,y_dm1,z_dm1,c= 'lightblue',s=0.1, label = 'DM$_{1}$')
#ax.scatter3D(x_bh1,y_bh1,z_bh1,c='black',s = 10,label='M$_{1}$')

ax.scatter3D(x_m2,y_m2,z_m2,c='red',s = 1, label='M$_{2}$')
ax.scatter3D(x_dm2,y_dm2,z_dm2,c='coral',s = 0.1, label='DM$_{2}$')
#ax.scatter3D(x_bh2,y_bh2,z_bh2,c='black',s = 10,label='M$_{1}$')
plt.legend()
plt.savefig('initial_spiriling3.png')

#salvo il mio raggio lagrangiano iniziale



# In[ ]:


####MI CALCOLO LA DENSITà INIZIALE.
r_m1 = np.sqrt(x_m1**2+y_m1**2+z_m1**2)
r_dm1 = np.sqrt(x_dm1**2+y_dm1**2+z_dm1**2)
r_m1.sort()
r_dm1.sort()

r_m2 = np.sqrt(x_m2**2+y_m2**2+z_m2**2)
r_dm2 = np.sqrt(x_dm2**2+y_dm2**2+z_dm2**2)
r_m2.sort()
r_dm2.sort()

R = 0
R_max=200
r_step=0.1
a1 = (2*np.pi*r)
V2 = (4/3 * np.pi * r**3)

while R < R_max:
    npart = 0 
    if R < r_m1[i] < R+r_step:
        npart = npart+1
        


# In[21]:


####GENERO IL FILE TXT

Npart_f= 40004
np.savetxt('m1.txt',np.column_stack([mf[0,:]]))
np.savetxt('m2.txt',np.column_stack([mf[0,:]]))
np.savetxt('galaxy1pos.txt',np.column_stack([pos1[0,:,0],pos1[0,:,1],pos1[0,:,2]]))
np.savetxt('galaxy2pos.txt',np.column_stack([pos2[0,:,0],pos2[0,:,1],pos2[0,:,2]]))
np.savetxt('galaxy1v.txt',np.column_stack([v1[0,:,0],v1[0,:,1],v1[0,:,2]]))
np.savetxt('galaxy2v.txt',np.column_stack([v2[0,:,0],v2[0,:,1],v2[0,:,2]]))
    

fname='galaxycondition.txt' 
fp = open(fname,'w') 
fp.write('40004') 
fp.write('\n') 
fp.write('3') 
fp.write('\n') 
fp.write('0') 
fp.write('\n') 
fp.close()

filenames = ['galaxycondition.txt','m1.txt','m2.txt', 'galaxy1pos.txt','galaxy2pos.txt','galaxy1v.txt','galaxy2v.txt']
with open('DoubleGalaxy.in.txt', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)


# In[107]:


######ANALIZZO LA SIMULAZIONE

Npart = 40004
fp=open('DoubleGalaxy.out.txt','r')
Nlines=len(fp.readlines())
Nsnapshot=(Nlines)//(3*Npart+3)
print (Nlines,Nsnapshot)
t = np.zeros(Nsnapshot)
pos_f = np.zeros([Nsnapshot,Npart,3])
v_f = np.zeros([Nsnapshot,Npart,3])


i=0

for i in range (Nsnapshot):
    initial_line=i*(3*Npart+3)   #è il numero di linee che passano tra uno snapshot e l'altro
   
   # mf[i]=np.genfromtxt('SingleGalaxy.in.txt',skip_header=initial_line+3,skip_footer=Nlines-initial_line-20005)
    t[i]=np.genfromtxt('DoubleGalaxy.out.txt',skip_header=initial_line+2, skip_footer=Nlines-initial_line-3)
    pos_f[i]=np.genfromtxt('DoubleGalaxy.out.txt',skip_header=initial_line+40007, skip_footer=Nlines-initial_line-80011)
    #v_f[i]=np.genfromtxt('DoubleGalaxy.out.txt',skip_header=initial_line+80011, skip_footer=Nlines-initial_line-60009)
    #print(t[i],mf[i],pos[i],vf[i])
    print (i)


# In[156]:


###SALVO LE COORDINATE DELLE 2 GALASSIE

#galaxy 1

x_bh1_f = pos_f[:,0,0]
y_bh1_f = pos_f[:,0,1]
z_bh1_f = pos_f[:,0,2]

x_m1_f = pos_f[:,1:10002,0]
y_m1_f = pos_f[:,1:10002,1]
z_m1_f = pos_f[:,1:10002,2]

x_dm1_f = pos_f[:,10002:20002,0]
y_dm1_f = pos_f[:,10002:20002,1]
z_dm1_f = pos_f[:,10002:20002,2]



#galaxy 2
x_bh2_f = pos_f[:,20002,0]
y_bh2_f = pos_f[:,20002,1]
z_bh2_f = pos_f[:,20002,2]

x_m2_f = pos_f[:,20003:30004,0]
y_m2_f = pos_f[:,20003:30004,1]
z_m2_f = pos_f[:,20003:30004,2]

x_dm2_f =  pos_f[:,30004:40004,0]
y_dm2_f =  pos_f[:,30004:40004,1]
z_dm2_f = pos_f[:,30004:40004,2]



# In[164]:


###GRAFICO L'EVOLUZIONE 
get_ipython().run_line_magic('matplotlib', 'inline')
for u in range (Nsnapshot):
    pictures = np.zeros([Nsnapshot])
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.set_xlim3d(-30,30)
    ax.set_ylim3d(-30,30)
    ax.set_zlim3d(-30,30)
    ax.set_zlabel('z [pc]')
    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    
    ax.scatter3D(x_bh1_f[u],y_bh1_f[u],z_bh1_f[u],c='black',s = 10,label='M$_{bh}$')
    ax.scatter3D(x_m1_f[u],y_m1_f[u],z_m1_f[u],c='blue',s = 1,label='M$_{1}$')
    ax.scatter3D(x_dm1_f[u],y_dm1_f[u],z_dm1_f[u],c= 'lightblue',s=0.1, label = 'DM$_{1}$')
    #ax.scatter3D(x_bh1,y_bh1,z_bh1,c='black',s = 10,label='M$_{1}$')
    
    ax.scatter3D(x_bh2_f[u],y_bh2_f[u],z_bh2_f[u],c='black',s = 10,label='M$_{bh}$')
    ax.scatter3D(x_m2_f[u],y_m2_f[u],z_m2_f[u],c='red',s = 1, label='M$_{2}$')
    ax.scatter3D(x_dm2_f[u],y_dm2_f[u],z_dm2_f[u],c='coral',s = 0.1, label='DM$_{2}$')
    #ax.scatter3D(x_bh2,y_bh2,z_bh2,c='black',s = 10,label='M$_{1}$')
    #plt.legend()
    plt.show()
    plt.savefig('evolution_spiriling[u].png')


# In[4]:


#ANALIZZO LO ZOOM NEL TEMPO
Npart_z = 40004
fp=open('provaout.txt','r')
Nlines_z=len(fp.readlines())
Nsnapshot_z=(Nlines_z)//(3*Npart_z+3)
print (Nlines_z,Nsnapshot_z)
tz = np.zeros(Nsnapshot_z)
pos_fz = np.zeros([Nsnapshot_z,Npart_z,3])
v_fz = np.zeros([Nsnapshot_z,Npart_z,3])
mf_z = np.zeros([Nsnapshot_z,Npart_z])

i=0

for i in range (Nsnapshot_z):
    initial_line=i*(3*Npart_z+3)   #è il numero di linee che passano tra uno snapshot e l'altro
   
    mf_z[i]=np.genfromtxt('provaout.txt',skip_header=initial_line+3,skip_footer=Nlines_z-initial_line-40007)
    tz[i]=np.genfromtxt('provaout.txt',skip_header=initial_line+2, skip_footer=Nlines_z-initial_line-3)
    pos_fz[i]=np.genfromtxt('provaout.txt',skip_header=initial_line+40007, skip_footer=Nlines_z-initial_line-80011)
    #v_fz[i]=np.genfromtxt('DoubleGalaxy.out.txt',skip_header=initial_line+80011, skip_footer=Nlines-initial_line-60009)
    #print(t[i],mf[i],pos[i],vf[i])
    print (i)


# In[6]:


from tempfile import TemporaryFile
m_t = TemporaryFile()
pos_t = TemporaryFile()
t_t = TemporaryFile()

np.save(m_t, mf_z)
np.save(pos_t,pos_fz)
np.save(t_t,tz)


# In[7]:


###SALVO LE COORDINATE DELLE 2 GALASSIE

#galaxy 1

x_bh1_fz = pos_fz[:,0,0]
y_bh1_fz = pos_fz[:,0,1]
z_bh1_fz = pos_fz[:,0,2]

x_m1_fz = pos_fz[:,1:10002,0]
y_m1_fz = pos_fz[:,1:10002,1]
z_m1_fz = pos_fz[:,1:10002,2]

x_dm1_fz = pos_fz[:,10002:20002,0]
y_dm1_fz = pos_fz[:,10002:20002,1]
z_dm1_fz = pos_fz[:,10002:20002,2]



#galaxy 2
x_bh2_fz = pos_fz[:,20002,0]
y_bh2_fz = pos_fz[:,20002,1]
z_bh2_fz = pos_fz[:,20002,2]

x_m2_fz = pos_fz[:,20003:30004,0]
y_m2_fz = pos_fz[:,20003:30004,1]
z_m2_fz = pos_fz[:,20003:30004,2]

x_dm2_fz =  pos_fz[:,30004:40004,0]
y_dm2_fz =  pos_fz[:,30004:40004,1]
z_dm2_fz =  pos_fz[:,30004:40004,2]

print('ciao')


# In[11]:


###GRAFICO L'EVOLUZIONE 
get_ipython().run_line_magic('matplotlib', 'inline')

for u in range (Nsnapshot_z):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.set_xlim3d(-35,35)
    ax.set_ylim3d(-35,35)
    ax.set_zlim3d(-35,35)
    ax.set_zlabel('z')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    ax.scatter3D(x_bh1_fz[u],y_bh1_fz[u],z_bh1_fz[u],c='black',s = 10,label='M$_{bh}$')
    ax.scatter3D(x_m1_fz[u],y_m1_fz[u],z_m1_fz[u],c='blue',s = 0.5,label='M$_{1}$')
    ax.scatter3D(x_dm1_fz[u],y_dm1_fz[u],z_dm1_fz[u],c= 'lightblue',s=0.01, label = 'DM$_{1}$')
    #ax.scatter3D(x_bh1,y_bh1,z_bh1,c='black',s = 10,label='M$_{1}$')
    
    ax.scatter3D(x_bh2_fz[u],y_bh2_fz[u],z_bh2_fz[u],c='black',s = 10,label='M$_{bh}$')
    ax.scatter3D(x_m2_fz[u],y_m2_fz[u],z_m2_fz[u],c='red',s = 0.5, label='M$_{2}$')
    ax.scatter3D(x_dm2_fz[u],y_dm2_fz[u],z_dm2_fz[u],c='coral',s = 0.01, label='DM$_{2}$')
    #ax.scatter3D(x_bh2,y_bh2,z_bh2,c='black',s = 10,label='M$_{1}$')
    #plt.legend()
    ax.set_title('t = %d' %(5*u))
    plt.savefig('evoluzione{u}.png'.format(u=u))
    plt.show()
    


# In[8]:


#####DEVO CALCOLARMI L'EVOLUZIONE DELLA DENSITà
#DEFINISCO I VETTORI DISTANZA DAL CENTRO, CONTENENTE L'INFO DELLA MASSA.
r = np.zeros([Nsnapshot_z,Npart_z])

for u in range (Nsnapshot_z):
    for i in range (Npart_z):
        r[u,i] = np.sqrt(pos_fz[u,i,0]**2+pos_fz[u,i,1]**2+pos_fz[u,i,2]**2)
        
r.sort()
print (r)


# In[ ]:





# In[11]:


#DENSITà EVOLUZIONE
r.sort()
r_step = 1
R = 0
cont_matter = 0
cont_dmatter = 0
cont_bh = 0
rho_shell = []
m1_shell = []
dm1_shell = []
bh1_shell = []
dr = []
raggio = []

while R<50:
    cont_matter = 0
    cont_dmatter = 0
    cont_bh = 0
    for i in range (Npart_z):
        if R < r[0,i] < R + r_step:
            print ('particella=',i)
            if mf_z[0,i] > 0.2:
                cont_bh = cont_bh + 1
                print ('Black hole= ', cont_bh, i)
            if 0.0005 < mf_z[0,i] < 0.2:
                cont_dmatter = cont_dmatter + 1
                print ('materia oscura=',i,cont_dmatter)
            if mf_z[0,i] < 0.0005:
                cont_matter = cont_matter + 1
                print ('materia=',cont_matter,i)
        
        
        
        m1_shell.append(cont_matter)
        dm1_shell.append(cont_dmatter)
        bh1_shell.append(cont_bh)
        diffr = (R+r_step)**3-R**3
        dr.append(diffr)
        raggio.append(dr)
     
    print ('R=',R)
    R = R + r_step

    
m_shell = np.array (m_shell)
dm_shell = np.array (dm_shell)
bh_shell = np.array (bh_shell)
dr = np.array(dr)



# In[11]:


###Massa totale
Mtot = 0
for i in range (Npart_z):
    Mtot = Mtot + mf_z[0,i]
print (Mtot)


# In[16]:


R = 0
r_step = 0.5


# In[45]:


#RAGGI LAGRANGIANI nel piano
r_m1 = np.sqrt(x_m1_fz**2+y_m1_fz**2+z_m1_fz**2)
r_dm1 = np.sqrt(x_dm1_fz**2+y_dm1_fz**2+z_dm1_fz**2)
r_bh1 = np.sqrt(x_bh1_fz**2+y_bh1_fz**2+z_bh1_fz**2)

r_m2 = np.sqrt(x_m2_fz**2+y_m2_fz**2+z_m2_fz**2)
r_dm2 = np.sqrt(x_dm2_fz**2+y_dm2_fz**2+z_dm2_fz**2)
r_bh2 = np.sqrt(x_bh2_fz**2+y_bh2_fz**2+z_bh2_fz**2)

r_m1.sort()
r_m2.sort()


r_mtot = np.zeros([Nsnapshot_z,20002])

for i in range (Nsnapshot_z):
    r_mtot[i] = np.union1d(r_m1[i,:],r_m2[i,:])
    


# In[52]:


#RAGGI LAGRANGIANI nel sistema del cdm delle galassie

print (x_bh1_fz.shape)

r1_m1= np.zeros([Nsnapshot_z,10001])
r2_m2= np.zeros([Nsnapshot_z,10001])


for u in range (Nsnapshot_z):
    for i in range (10001):
        r1_m1[u,i] = np.sqrt((x_m1_fz[u,i]-x_bh1_fz[u])**2+(y_m1_fz[u,i]-y_bh1_fz[u])**2+(z_m1_fz[u,i]-z_bh1_fz[u])**2)
        r2_m2[u,i] = np.sqrt((x_m2_fz[u,i]-x_bh2_fz[u])**2+(y_m2_fz[u,i]-y_bh2_fz[u])**2+(z_m2_fz[u,i]-z_bh2_fz[u])**2)

        
r1_m1.sort()
r2_m2.sort()






# In[44]:


plt.plot(tz,r1_m1[:,9000], label='R$_{90}$ galaxy1')
plt.plot(tz,r2_m2[:,9000],label='R$_{90}$ galaxy2')
plt.legend()
plt.savefig('esplosione90.png')
plt.show()


# In[51]:


plt.plot(tz,r_mtot[:,18000], label='R$_{90}$')

plt.legend()
plt.savefig('esplosione90_00.png')
plt.show()


# In[60]:


plt.plot(tz,r1_m1[:,1000],label='R$_{10}$')
plt.plot(tz,r1_m1[:,3000],label='R$_{30}$')
plt.plot(tz,r1_m1[:,5000],label='R$_{50}$')
plt.plot(tz,r1_m1[:,6000],label='R$_{60}$')
plt.plot(tz,r1_m1[:,7000],label='R$_{70}$')
plt.xlabel('t')
plt.ylabel('r')
plt.legend()
plt.savefig('raggi_lagr_galaxy1.png')

plt.show()

plt.plot(tz,r2_m2[:,1000],label='R$_{10}$')
plt.plot(tz,r2_m2[:,3000],label='R$_{30}$')
plt.plot(tz,r2_m2[:,5000],label='R$_{50}$')
plt.plot(tz,r2_m2[:,6000],label='R$_{60}$')
plt.plot(tz,r2_m2[:,7000],label='R$_{70}$')
plt.xlabel('t')
plt.ylabel('r')
plt.legend()
plt.savefig('raggi_lagr_galaxy2.png')
plt.show()


# In[61]:



plt.plot(tz,r_mtot[:,2000],label='R$_{10}$')
plt.plot(tz,r_mtot[:,6000],label='R$_{30}$')
plt.plot(tz,r_mtot[:,10000],label='R$_{50}$')
plt.plot(tz,r_mtot[:,12000],label='R$_{60}$')
plt.plot(tz,r_mtot[:,14000],label='R$_{70}$')
plt.xlabel('t')
plt.ylabel('r')
plt.legend()
plt.savefig('raggi_lagr_00.png')
plt.show()


# In[58]:


####FACCIO ANCHE PER LA MATERIA
r_dmtot = np.zeros([Nsnapshot_z,20000])

for i in range (Nsnapshot_z):
    r_dmtot[i] = np.union1d(r_dm1[i,:],r_dm2[i,:])
    

    


# In[62]:


plt.plot(tz,r_dmtot[:,2000])
plt.plot(tz,r_dmtot[:,6000])
plt.plot(tz,r_dmtot[:,10000])
plt.plot(tz,r_dmtot[:,12000])
plt.plot(tz,r_dmtot[:,14000])
plt.xlabel('t')
plt.ylabel('r')
plt.savefig('lagr_dm.png')
plt.show()


# In[73]:


for u in range (Nsnapshot_z):
    plt.plot(x_m1_fz[u,:],y_m1_fz[u,:],'.',c='blue')
    plt.plot(x_m2_fz[u,:],y_m2_fz[u,:],'.',c='red')
    plt.xlim(-50,50)
    plt.ylim(-50,50)
    plt.show()


# In[121]:


plt.hist2d(x_m1_fz[0,:],y_m1_fz[0,:],bins=(6000,6000),density=True)
plt.xlim(16,24)
plt.ylim(-4,4)


# In[ ]:





# In[149]:


from scipy.stats import gaussian_kde


# In[159]:


xy = np.vstack([x_m1_fz[0,:],y_m1_fz[0,:]])
z = gaussian_kde(xy)(xy)
xy2 = np.vstack([x_m2_fz[0,:],y_m2_fz[0,:]])
z2 = gaussian_kde(xy2)(xy2)

idx = z.argsort()
x_m1_fz[0,:], y_m1_fz[0,:], z = x_m1_fz[0,idx], y_m1_fz[0,idx], z[idx]
idx2 = z2.argsort()
x_m2_fz[0,:], y_m2_fz[0,:], z = x_m2_fz[0,idx2], y_m2_fz[0,idx2], z2[idx2]


fig, ax = plt.subplots()
ax.scatter(x_m1_fz[0,:], y_m1_fz[0,:], c=z, s=25, edgecolor='')
ax.scatter(x_m2_fz[0,:], y_m2_fz[0,:], c=z2, s=25, edgecolor='')
plt.xlim(-40,40)
plt.ylim(-40,40)
plt.savefig('densità_iniziale.png')
plt.show()


# In[160]:


xy = np.vstack([x_m1_fz[100,:],y_m1_fz[100,:]])
z = gaussian_kde(xy)(xy)
xy2 = np.vstack([x_m2_fz[100,:],y_m2_fz[100,:]])
z2 = gaussian_kde(xy2)(xy2)

idx = z.argsort()
x_m1_fz[100,:], y_m1_fz[100,:], z = x_m1_fz[100,idx], y_m1_fz[100,idx], z[idx]
idx2 = z2.argsort()
x_m2_fz[100,:], y_m2_fz[100,:], z = x_m2_fz[100,idx2], y_m2_fz[100,idx2], z2[idx2]


fig, ax = plt.subplots()
ax.scatter(x_m1_fz[100,:], y_m1_fz[100,:], c=z, s=25, edgecolor='')
ax.scatter(x_m2_fz[100,:], y_m2_fz[100,:], c=z2, s=25, edgecolor='')
plt.xlim(-40,40)
plt.ylim(-40,40)
plt.savefig('densità_finale.png')
plt.show()


# In[1]:


from sklearn.neighbors import KDTree

vect = np.concatenate((r1_m1, r2_m2), axis = 1)
print (vect.shape)
density = np.zeros([Nsnapshot_z,20002])
for i in range (Nsnapshot_z):
    tree = KDTree(vect[i,:],leaf_size = 2)
    dist,ind = tree.query(vect[i,:,:],k=5)
    density[i] = (6*mf_z[1])/(np.pi*dist[:,4]**2)


# In[ ]:




