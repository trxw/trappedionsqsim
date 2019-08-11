#############################################################################################################
#      Simulate Single and Multi Qubit Gates Setting Experimental Parameters of a Multi-Chromatic Field
#      Author: Daniel Murphy
#              Georgia Institute of Technology
#############################################################################################################


import numpy as np
from numpy import pi
from qutip import *
import matplotlib.pyplot as plt
import operators.operators as operators
import experiment
import groundcoupling

##############Raman Dynamics Simulation Yb171+ E2 (2S1/2 ==> 2P1/2, 2P3/2)#####################
###############################################################################################

######Level Structure Info (ground in basis of |0,0>, |1,-1>, |1,0>, |1,1>)###################################################
mg = [0, -1, 0, 1]
fg = [0, 1, 1, 1]
fint = [0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
mint = [0, -1, 0, 1, -1, 0, 1, -2, -1, 0, 1, 2]
dipole = [3.49, 3.49, 3.49, 3.49, 4.66, 4.66, 4.66, 4.66, 4.66, 4.66, 4.66, 4.66]
ground = [0, 12640, 12640, 12640]
detuning = [-33000000, -33000000, -33000000, -33000000, 66000000, 66000000, 66000000, 66000000, 66000000, 66000000, 66000000, 66000000]
##############################################################################################################################


###############################################################################################################################
##################A simulation of single qubit Raman Dynamics in Bichromatic Field#############################################
###############################################################################################################################
###############################################################################################################################

def SQramanDynamicsBC(Expm):
    
    if (Expm.B_C):
        det = Expm.Dets
        EP = Expm.E[0][0]
        ES = Expm.E[1][0]
        tmax = Expm.tmax
        tstep = Expm.tstep
        dPhase = Expm.P
        
        gc = groundcoupling.GroundCoupling(detuning, fint, mint, fg, mg, dipole)
        
        H = []
        
        #########Calculate Second Order Stark Shift#########
        ss = []
        
        ss.append(gc.getCoupling(0, 0, EP, EP))
        
        for i in range(1, 4):
            
            ss.append(gc.getCoupling(i, i, ES, ES))
        
        ####################################################
        
        #########Calculate Couplings################
        O = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
        
        for i in range(1, 4):
            O[0][i] = gc.getCoupling(0, i, EP, ES)
            O[i][0] = np.conjugate(O[0][i])
        
        for i in range(1, 4):
            for j in range(1, 4):
                O[i][j] = gc.getCoupling(i, j, ES, ES)
        ############################################

        ##########Calculate Detunings From Degenerate Levels###############
        dets = [0, 0, 0]

        for i in range(1, 4):

            dets[i - 1] = det + (ground[2] - ground[i]) + (ss[2] - ss[i])

        ##################################################################

        for i in range(1,4):
            for j in range(1,4):
                if(i != j):
                    func = ('exp(1j*(%f*t))' % (det[i-1] - det[j-1]))
                    H.append([O[i][j]*basis(4, i)*basis(4, j).dag(), func])

        for i in range(1,4):
            func = ('exp(1j*(%f*t))' % (det[i-1]))
            H.append([O[i][0]*basis(4, i)*basis(4, 0).dag(), func])


        for i in range(1,4):
            func = ('exp(-1j*(%f*t))' % (det[i-1]))
            H.append([O[i][0]*basis(4, i)*basis(4, 0).dag(), func])


        psi0 = basis(4, 0)

        t = np.linspace(0, tmax, tstep)
        sol = mesolve(H, psi0, t, [], [], progress_bar = True)

        return sol

    else:

        print("Not valid Bichromatic experiment")



##########################################################################################################################################
##################A simulation of Phonon Coupled Raman Dynamics in a Bichromatic Field#################################################
###########################################################################################################################################
def MQramanDynamicsBC(Expm):
    
    
    
    modes = Expm.N_m
    trunc = Expm.D_f
    vibF = Expm.interactionFrequencies()
    Phase = Expm.P
    tstep = Expm.tstep
    tmax = Expm.tmax
    ES = Expm.F_s[1]
    EP = Expm.F_s[0][0]
    detBR = Expm.Dets
    eta = Expm.Dicke
    
    
    gc = groundcoupling.GroundCoupling(detuning, fint, mint, fg, mg, dipole)
    
    if (Expm.B_C):
        
        #######Calculate 2nd Order Stark Effect######
        ss = [0, 0, 0, 0]
        ss[0] = ss[0] + gc.getCoupling(0, 0, EP, EP)
        for c in range(0,2):
                for i in range(1,4):
                        ss[i] = ss[i] + gc.getCoupling(i, i, ES[c], ES[c])
        #############################################


        O = [[[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]],
             [[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]]]
        ####Store Effective Rabi Couplings################
        for c in range(0, 2):
            for i in range(1, 4):
                    O[c][0][i][0] += gc.getCoupling(i, 0, ES[c], EP)
                    O[c][0][0][i] += np.conjugate(O[c][0][i][0])

        for c1 in range(0, 2):
            for c2 in range(0, 2):
                for i in range(1, 4):
                    for j in range(1, 4):
                        O[c1][c2][j][i] = gc.getCoupling(j, i, ES[c1], ES[c2])

        print(O[0][0][2][0])
        ######################################################

        det = [[0, 0, 0], [0, 0, 0]]
        #Calculate Detuning from resonance for all States up to 2nd order stark shift#
        for c in range(0,2):
            for i in range(1, 4):
                det[c][i - 1] = detBR[c] + ((ground[2] - ground[i]) + (ss[2] - ss[i]))

        det = np.real(det)
        ##############################################################################



        H = []
        HM = []
        op = operators.Operators(2, 4, modes, trunc)

        #######Create Motional Hamiltonian Components##########
        for m in range(0, modes):
            HM.append([op.id(), eta[m]*1j*op.ad(m+1), eta[m]*1j*op.a(m+1)])

        ######################################################

        #####First Ion Hamiltonian#####
        for m in range(0, modes):
            for e in range(0, 3):
                for c in range(0, 2):
                    for i in range(0, 3):
                        func = ('exp(1j*(%f*t + %f))' % (det[c][i]+vibF[m][e], Phase[c]))
                        H.append([O[c][0][i + 1][0]*op.coupling([[-1, 0]], [[-1, i + 1]])*HM[m][e], func])


        for m in range(0, modes):
            for e in range(0, 3):
                for c in range(0, 2):
                    for i in range(0, 3):
                        func = ('exp(-1j*(%f*t + %f))' % (det[c][i]+vibF[m][e], Phase[c]))
                        H.append([O[c][0][0][i + 1]*op.coupling([[]-1, i + 1]], [[-1, 0]])*HM[m][e].dag(), func])

        for c1 in range(0, 2):
            for c0 in range(0, 2):
                for i in range(1, 4):
                    for j in range(1, 4):
                        if (i != j):
                            func = ('exp(1j*(%f*t + %f))' % (det[c1][i-1] - det[c0][j-1], Phase[c0] - Phase[c1]))
                            H.append([O[c][0][i][j]*op.coupling([[-1, j]], [[-1, i]])*op.id(), func])

        #######################################

        ####Second Ion Hamiltonian#####
        for m in range(0, modes):
            for e in range(0, 3):
                for c in range(0, 2):
                    for i in range(0, 3):
                        func = ('exp(1j*(%f*t + %f))' % (det[c][i]+vibF[m][e], Phase[c]))
                        H.append([O[c][0][i + 1][0]*op.coupling([[0, -1]], [[i + 1, -1]])*HM[m][e], func])

        for m in range(0, modes):
            for e in range(0, 3):
                for c in range(0, 2):
                    for i in range(0, 3):
                        func = ('exp(-1j*(%f*t + %f))' % (det[c][i]+vibF[m][e], Phase[c]))
                        H.append([O[c][0][0][i + 1]*op.coupling([[i + 1, -1]], [[0, -1]])*HM[m][e].dag(), func])

        for c1 in range(0, 2):
            for c0 in range(0, 2):
                for i in range(1, 4):
                    for j in range(1, 4):
                        if (i != j):
                            func = ('exp(1j*(%f*t + %f))' % (det[c1][i-1] - det[c0][j-1], Phase[c0] - Phase[c1]))
                            H.append([O[c][0][i][j]*op.coupling([[j, -1]], [[i, -1]])*op.id(), func])

        #########################################

        #######Simulate Dynamics########
        t = np.linspace(0, tmax, tstep)
        ###Start in Ground State of All Modes###
        psi0 = tensor(basis(4, 0), basis(4, 0), tensor([basis(trunc, 0) for i in range(0, modes)]))
        ########################################
        data = mesolve(H, psi0, t, [], [], progress_bar = True)

        return data

    else:

        print("Not Valid Bichromatic Experiment")

        return null


###########################################################################################################









####Define Experiment####
Expm = experiment.Experiment(B_C = True, N_m = 1, D_f = 20, F_f = [3], Dets = [-2.85, 2.85],
                             F_s = [[[15000, 0, 0]], [[0, 15000, 0], [0, 15000, 0]]],
                             P = [0, 0], Dicke = [.1], t = 200, tstep = 40000)
#########################
























    












