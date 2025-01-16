import numpy as np
import math
import matplotlib.pyplot as pylab
import sys
from ED_necessary_functions import produce_states, fix_bonds, swap_bits, CdagC, is_hermitian
from bonds_n_coords import lat_bonds
from bonds_n_coords_diagonal import diagonal_lat_bonds

if len(sys.argv) < 7:
    print('Provide inputs: Nx Ny Nup Ndn string geometery string U')
    sys.exit(1)

Nx = int(sys.argv[1])
Ny = int(sys.argv[2])
Nsite = int(Nx*Ny)
Nup = int(sys.argv[3])
Ndn = int(sys.argv[4])
geometery = sys.argv[5]
U = sys.argv[6]
Np = Nup + Ndn
t=1

if geometery == 'diagonal':
    bonds, open_bonds = diagonal_lat_bonds(Nx,Ny)
else:    
    bonds, open_bonds = lat_bonds(Nx,Ny)

print(bonds)
bonds = fix_bonds(bonds)


states = produce_states(Nsite, Nup, Ndn, U)

#---------------------------------Make Hamiltoinan------------------------------
Hamiltonian = np.zeros((len(states),len(states)))
for i in range(len(states)):
    for j in range(len(bonds)):
        new_state , si = swap_bits(states[i],Nsite,bonds[j][0],bonds[j][1],U)
        if new_state != states[i]:
            for k in range(len(states)):
                if new_state == states[k]:
                     Hamiltonian[i][k] = -t * si
                     Hamiltonian[k][i] = -t * si

print('Is Hamiltonian hermitian ' ,is_hermitian(Hamiltonian))
eigenvalues, vectors = np.linalg.eigh(Hamiltonian)
print("Eigenvalues:", eigenvalues)
print("Eigenvalue of GS:", eigenvalues[0])
GS=vectors[:,0]
print("Ground state ket",'\n', vectors[:,0])

#-------------------------------Measure CdagupCup ----------------------------------------------------
#print('------------------------------------------------------')
#print('i   j   CdagupCup')
site_i_fixed_back =0
CdagupCup_final = []
site_i_final=[]
site_j_final=[]
for site_i in range(0,2*Nsite,2):
    site_i_fixed_back +=1
    site_j_fixed_back =0
    for site_j in range(0,2*Nsite,2):
           site_j_fixed_back +=1
           site_i_final.append(site_i_fixed_back)
           site_j_final.append(site_j_fixed_back)
           CdagupCup = np.zeros((len(states),len(states)))
           for i in range(len(states)):
                      new_state , si = CdagC(states[i],Nsite,site_i,site_j,U)
                      for k in range(len(states)):
                               if new_state == states[k]:
                                     CdagupCup[i][k] = si
           obs_CdagupCup = np.dot(GS.T,np.dot(CdagupCup,GS))
           CdagupCup_final.append(obs_CdagupCup) 
           #print(site_i_fixed_back,' ',site_j_fixed_back,' ',round(obs_CdagupCup,6))
           #print(CdagupCup)

#-------------------------------------Measure CdagdnCdn----------------------------------
#print("---------------------------------")
#print('i   j   CdagdnCdn')
site_i_fixed_back =0
CdagdnCdn_final = []
for site_i in range(1,2*Nsite,2):
    site_i_fixed_back +=1
    site_j_fixed_back =0
    for site_j in range(1,2*Nsite,2):
           site_j_fixed_back +=1
           CdagdnCdn = np.zeros((len(states),len(states)))
           for i in range(len(states)):
                      new_state , si = CdagC(states[i],Nsite,site_i,site_j,U)
                      for k in range(len(states)):
                               if new_state == states[k]:
                                     CdagdnCdn[i][k] = si
           obs_CdagdnCdn = np.dot(GS.T,np.dot(CdagdnCdn,GS))
           CdagdnCdn_final.append(obs_CdagdnCdn)
           #print(CdagupCup)
           #print(site_i_fixed_back,' ',site_j_fixed_back,' ',round(obs_CdagdnCdn,6))
print("-------------------------------------------------")
print('i   j   CdagupCup   CdagdnCdn')
for k in range(len(CdagupCup_final)):
    print(site_i_final[k],' ',site_j_final[k],' ',round(CdagupCup_final[k],6),'  ',round(CdagdnCdn_final[k],6))

    
   
