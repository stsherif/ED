import matplotlib.pyplot as plt
import numpy as np


def lat_bonds(Nx,Ny):
    N = Nx*Ny
    bond=[]
    bonds = []
    for i in range(1,N):
        if i%Ny==1 and i<N-Ny:
            b=[i,i+1]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny]
            bonds.append(b)
            bond.append(b)
            b=[i,i+Ny+1]
            bonds.append(b)
            bond.append(b)
            b=[i,i+Ny-1]
            #bonds.append(b)
            bond.append(b)
            b=[i,i+2*Ny-1]
            #bonds.append(b)
            bond.append(b)
   
        elif (i%Ny==0 and i<N):
            b=[i,i+Ny]
            bond.append(b)
            bonds.append(b)

        elif i>N-Ny:
            b=[i,i+1]
            bond.append(b)
            bonds.append(b)
            if i==N-Ny+1:
                b=[i,i+Ny-1]
                bond.append(b)

        elif i%2==1:
            b=[i,i+1]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny-1]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny+1]
            bond.append(b)
            bonds.append(b)

        elif i%2==0:
            b=[i,i+1]
            bond.append(b)
            bonds.append(b)
            b=[i,i+Ny]
            bond.append(b)
            bonds.append(b)
    #I added this part to remove simillar elements from bond (only happens for Ny = 2)
    unique_bond=sorted(set(tuple(b) for b in bond))
    return unique_bond, bonds     #bond for cylindrical bc and bonds for open

"""
bond,_ = lat_bonds()

for i  in range(len(bond)):
    print(bond[i])
"""

#finding coordinates
def coordinates(Nx,Ny):
    N = Nx*Ny
    a1 = np.array([1,0])
    a2 = np.array([0.5,np.sqrt(3)/2.0])

    Rxy = {}

    for j in range(Ny,0,-1):
        for i in range(1,Nx+1):
            label = j+(i-1)*Ny
            if j == 6:
                Rxy[label] = (i-1)*a1
            elif j==5:
                Rxy[label] = (i-1)*a1 + a2
            elif j==4:
                Rxy[label] = (i-2)*a1 + 2*a2
            elif j==3:
                Rxy[label] = (i-2)*a1 + 3*a2
            elif j==2:
                Rxy[label] = (i-3)*a1 + 4*a2
            else:
                Rxy[label] = (i-3)*a1 + 5*a2
    return Rxy

'''
Nx = 2
Ny = 2
N = Nx*Ny
bond,bonds = lat_bonds(Nx,Ny)
Rxy = coordinates(Nx,Ny)
for k,bnd in enumerate(bonds):
    s1=bnd[0]
    s2=bnd[1]
    
    plt.plot([Rxy[s1][0],Rxy[s2][0]],[Rxy[s1][1],Rxy[s2][1]],lw=1.3,c='k')
    plt.text(Rxy[s1][0]+0.05,Rxy[s1][1],str(s1),fontsize=10)#,weight="bold",color='green')
#plt.text(Rxy[N][0]+0.05,Rxy[36][1],str(N),fontsize=10)

#-----------------------------------------------temp------------just_for_a_pic-----------------------------

for k,bnd in enumerate(bonds):
    s1=bnd[0]
    s2=bnd[1]
    plt.scatter(Rxy[s1][0],Rxy[s1][1],s=80,c='k')
ll=22    
#plt.scatter(Rxy[ll][0],Rxy[ll][1],s=80,c='white')

plt.axis("off")
#plt.savefig('XC_6x6_1h.pdf',dpi=600,bbox_inches='tight')
print(bond)
print(Rxy)
#plt.show()
'''
