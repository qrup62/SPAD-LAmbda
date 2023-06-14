#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 16:05:59 2022

@author: rotemliran

Numerical exercise 5 physical statistics
"""
import numpy as np
import functionNumerical5 as fun
import timeit
## Initial condition
# creat a lattice N*N were each particale have 
# spin up (1) or down(-1) randomly
N=36


## Dynamic of the lattice 
# each step, run over all the particals and by raffle determine 
# to change or stay.

# The raffle is done by - 
# the particle energy by the field and the nighbours
 
## record of the system 
# U the energy of the system 

# M the magnetization of the system

# c_V  The heat capacity of the system

# k_B* T = 1/Beta = constant
# eta = j*beta 
# h = mu * B * beta

start1 = timeit.default_timer()
MagneticField = 0
eta = np.linspace(0.42,0.46,9)
eta = 2

M0_All = np.zeros([np.size(MagneticField) , np.size(eta) ])
U0_All = np.zeros([np.size(MagneticField) , np.size(eta) ])
U0_squareAll=  np.zeros([np.size(MagneticField) , np.size(eta) ])
M0SquareAll =  np.zeros([np.size(MagneticField) , np.size(eta) ])

for i in range(np.size(MagneticField)):
    for j in range(np.size(eta)):
        Lattice = fun.randomspinLattice(N) # initial lattice 
        M0_Square = 0
        M0_K = 0 
        M_Khalf = 0
        U0 = 0
        U0_square = 0 
        B  = MagneticField
        J = eta
        K=1*(10**5)
        Delta = 10**(-3)

        flip_N1 = 0
    ## forget initail condition
        while flip_N1<(K/2):  
            [Lattice , flip_N] = fun.promotation(Lattice ,J,B)
            flip_N1 =  flip_N1  + flip_N
    ## premote half K
        print('Flip number th first: ' , flip_N1)
        flip_N1 = 0
        while flip_N1<(K/2): 
            flip_N2 = 0
            for k5 in range(5):
                [Lattice , flip_N] = fun.promotation(Lattice ,J,B)
                flip_N2 = flip_N2 + flip_N
            M_number = sum(sum(Lattice))  
            M_Khalf = M_Khalf  +M_number  
            flip_N1 =  flip_N1  + flip_N2
             
        M_Khalf = M_Khalf/(flip_N1)
        
    ## premote  until the settsfay condition
        
        stillNot = True
        while stillNot:
            flip_N1 = 0 
            DIVIDER = 0 
            while flip_N1 < K:
                flip_N2 = 0
                
                
                for k5 in range(5): 
                    [Lattice , flip_N] = fun.promotation(Lattice ,J,B)
                    flip_N2 = flip_N2 + flip_N
                
                DIVIDER = DIVIDER+1
                M_number = sum(sum(Lattice))  
                M0_K = M0_K  +M_number  
                M0_Square = M0_Square + (M_number**2)
                energy = fun.TotalEnergy(Lattice,J,B)            
                U0 = U0 + energy
                U0_square =  U0_square + (energy**2)
                flip_N1 =  flip_N1  + flip_N2
            

            M0_K = M0_K/DIVIDER
            M0_Square = M0_Square/DIVIDER
            U0 = U0/DIVIDER
            U0_square = U0_square/DIVIDER
            condition = abs(M0_K - M_Khalf)/abs(M0_K)
            if condition < Delta or K>=(10**8):
                M0_All[i,j] = M0_K
                M0SquareAll[i,j] = M0_Square
                U0_All[i,j] = U0
                U0_squareAll[i,j] = U0_square


                print('M_0 = ' , M0_K , '    i = '  , i , '\n j  = ' , j )
                stillNot = False
            else:
                
                K = 2*K

                print('Condition : ', condition , '\n M_0 = ' , M0_K 
                      ,'\n M0_Khalf =  ',M_Khalf ,  ' \n Doubled K to: ' , K)
                M_Khalf = M0_K
                M0_K = 0 
                M0_Square = 0
                U0 = 0
                U0_square = 0 
    
            end1 = timeit.default_timer()
            print ('time: ', end1 - start1)
            start1 = timeit.default_timer()




np.savetxt("M0_AllJ.csv", M0_All, delimiter=",")
np.savetxt("M0SquareAllJ.csv", M0SquareAll, delimiter=",")
np.savetxt("U0_AllJ.csv", U0_All, delimiter=",")
np.savetxt("U0_squareAll.csv", U0_squareAll, delimiter=",")
np.savetxt("etaJ.csv", eta, delimiter=",")
