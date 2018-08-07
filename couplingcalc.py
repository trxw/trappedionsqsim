
import numpy as np
import clebschgordon

class CouplingCalc:
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        Object for computing Electric Field Coupling Between States of F1 and F2 given F2 and reduced dipole element by Electric Field E
        
        Author: Daniel Murphy
                Georgia Institute of Technology
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    

    
    def __init__(self, F1, m1, E = [0, 0, 0]):

        self.E = E
        self.F1 = F1
        self.m1 = m1
        self.cg = clebschgordon.ClebschGordonCalc()

    ####Transforms Carteisian Vector to Spherical Basis#####
    def toSpherical(self, v):

        vp = np.empty([3], dtype = complex)
        
        vp[0] = (v[0] + 1j*v[2])/(np.sqrt(2))
        vp[1] = v[2]
        vp[2] = -1*(v[0] - 1j*v[1])/(np.sqrt(2))
        
        return vp

    ####Get Couplings Between States of m1 and m2####
    def getCoupling(self, F2, m2, dipole):

        dipVec = []

        for m in range(-1, 2):
            
            dipVec.append(dipole*self.cg.CGCoeff(self.F1, self.m1, F2, m2, m)/(np.sqrt(2*self.F1 + 1)))


        return np.vdot(self.toSpherical(self.E), dipVec)




