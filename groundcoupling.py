
import numpy as np
import couplingcalc


class GroundCoupling:

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        Calculates Two-Photon Coupling Between Ground States Given Field
        and Given States
        
        
        Author: Daniel Murphy
                Georgia Institute of Technology
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    def __init__(self, detuning, fint, mint, fg, mg, dipole):

        self.detuning = detuning
        self.fint = fint
        self.mint = mint
        self.fg = fg
        self.mg = mg
        self.dipole = dipole
    
    def intermediateRabiArray(self, i, E):
        
        v = []
        
        m = self.mg[i]
        f = self.fg[i]
        
        cc = couplingcalc.CouplingCalc(f, m, E)
        
        for j in range(0, len(self.fint)):
            
            v.append(cc.getCoupling(self.fint[j], self.mint[j], self.dipole[j]))
    
        return v
            
            
    def intermediateSum(self, v1, v2):

        sum = 0
        
        for e in range(0, len(self.detuning)):
            
            sum += (-1/(4*self.detuning[e]))*v1[e]*np.conjugate(v2[e])

        return sum


    def getCoupling(self, i, j, E1, E2):

        v1 = self.intermediateRabiArray(i, E1)
        v2 = self.intermediateRabiArray(j, E2)

        return self.intermediateSum(v1, v2)







