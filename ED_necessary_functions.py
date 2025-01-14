import numpy as np
import math
import matplotlib.pyplot as pylab
import sys

def produce_states(nsites, nup, ndn, U):
    states=[]
    for state in range(2 ** (2*nsites)):
        binary_string = bin(state)[2:].zfill(2*nsites) #remove the 0b and fill in zeros
        even_sites = sum(binary_string[j].count('1') for j in range(0,2*nsites,2))
        odd_sites = sum(binary_string[j].count('1') for j in range(1,2*nsites,2))
        if even_sites == nup and odd_sites == ndn:
            if U == "zero":
                 states.append(state)
                 print(state,' ',binary_string,' ',binary_string.count('1'))
            elif U == "inf":
                no_of_douple = 0
                for i in range(0,2*nsites,2):
                       if (binary_string[i] + binary_string[i+1] == bin(3)[2:]):
                              no_of_douple += 1
                if no_of_douple == 0:
                    states.append(state)
                    print(state,' ',binary_string,' ',binary_string.count('1'))
    return states                


def fix_bonds(bonds):
    sites_doubles_for_up_and_dn=[]
    fixed_bonds=[]
    for i in range(0,2*max(max(bonds)),2):
        sites_doubles_for_up_and_dn.append([i,i+1])
    for j in range(len(bonds)):
        fixed_bonds.append([sites_doubles_for_up_and_dn[bonds[j-1][0]-1][0],sites_doubles_for_up_and_dn[bonds[j-1][1]-1][0]])
        fixed_bonds.append([sites_doubles_for_up_and_dn[bonds[j-1][0]-1][1],sites_doubles_for_up_and_dn[bonds[j-1][1]-1][1]])
    print("fixed_bonds = " ,fixed_bonds)
    return fixed_bonds

def swap_bits(num,N, n, m, U):
    n=2*N-1-n
    m=2*N-1-m
    sign = -1 
    if m>n:
        n, m = m, n
    mask = ((1 << (n-m-1))-1) << (m+1)
    bit_n = (num >> n) & 1
    bit_m = (num >> m) & 1
    if bit_n != bit_m:
        if bit_n == 0:
            if (n%2 ==0 and ((num >> n+1) & 1)==0 and U == 'inf') or (n%2 ==1 and ((num >> n-1) & 1)==0 and U == 'inf') or (U == 'zero'):
                num ^= (1<<n)
                num ^= (1<<m)
                sign = sign ** bin(num & mask).count('1')
        if bit_m == 0:
            if (m%2 ==0 and ((num >> m+1) & 1)==0 and U == 'inf') or (m%2 ==1 and ((num >> m-1) & 1)==0 and U == 'inf') or (U == 'zero'):
                num ^= (1<<n)
                num ^= (1<<m)
                sign = sign ** bin(num & mask).count('1')
    return num, sign

def is_hermitian(matrix):
        return np.allclose(matrix, matrix.conj().T)

def CdagC(num, N, n, m, U):
    n = 2*N-1-n
    m = 2*N-1-m
    sign = -1
    mask = 0
    bit_n = (num >> n) & 1
    bit_m = (num >> m) & 1
    if n==m and bit_n==1:
        return num, 1
    elif n!=m and bit_n == 0:
        if n>m:
            if n-m-1 >= 0:
                mask = ((1 << (n-m-1))-1) << (m+1)
                if (n%2 ==0 and ((num >> n+1) & 1)==0 and U == 'inf') or (n%2 ==1 and ((num >> n-1) & 1)==0 and U == 'inf') or (U == 'zero'):
                     num ^= (1<<n)
                     num ^= (1<<m)
                     sign = sign ** bin(num & mask).count('1')
                     return num, sign
        if n<m:
            if m-n-1 >= 0:
               mask = ((1 << (m-n-1))-1) << (n+1)
               if (m%2 ==0 and ((num >> m+1) & 1)==0 and U == 'inf') or (m%2 ==1 and ((num >> m-1) & 1)==0 and U == 'inf') or (U == 'zero'):
                     num ^= (1<<n)
                     num ^= (1<<m)
                     sign = sign ** bin(num & mask).count('1')
                     return num, sign
    return num, 0                 
