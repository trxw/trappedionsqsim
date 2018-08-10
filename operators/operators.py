from __future__ import division, absolute_import, print_function, unicode_literals


from qutip import *
import numpy as np
from scipy import *

class Operators(object):
    """
    Create necessary operators for constructing a Hamiltonian of an ion chain with arbitrary Hilbert space dimensions 
    number_of_ions,
    dim_of_electronic_states_space,
    number_of_motional_modes,
    dim_of_each_Fock_space
    """
    def __init__(self, number_of_ions = 1,
                dim_of_electronic_states_space = 2,
                number_of_motional_modes = 1,
                dim_of_each_Fock_space = 20):
        self.N_e = number_of_ions
        self.D_F = dim_of_each_Fock_space
        self.N_F = number_of_motional_modes
        self.D_e = dim_of_electronic_states_space


    def a_func(self, i):
        """Define an annihiliation operator on the motional state of i-th ion.
        i = 1, 2, 3, ... number_of_ions
        """
        if i <= self.D_F:
            return tensor( [qeye(self.D_e) 
                            for i in range(1, self.N_e+1)] + 
                            [ destroy(self.D_F) if j == i else qeye(self.D_F) 
                            for j in range(1, self.N_F+1) ] )
        else: 
            print("Mode number should be between 1 and "+str(self.N_F))

    def ad_func(self, i):
        """Define a creation operator on the motional state of i-th ion.
        i = 1, 2, 3, ... number_of_ions
        """
        if i <= self.D_F:
            return tensor( [qeye(self.D_e) for i in range(1, self.N_e+1)] 
                           + [ create(self.D_F) if j == i else qeye(self.D_F) 
                            for j in range(1, self.N_F+1) ] )
        else: 
            print("Mode number should be between 1 and "+str(self.N_F))

    
    def sm_func(self, i):
        """Define an annihiliation operator on the electronic state of i-th ion.
        i = 1, 2, 3, ... number_of_ions
        """
        return tensor( [destroy(self.D_e) if j == i else qeye(self.D_e) 
                        for j in range(1, self.N_e+1)]
                        + [ qeye(self.D_F) for j in range(1, self.N_F+1) ] )
    
    def sp_func(self, i):
        """Define a creation operator on the electronic state of i-th ion.
        i = 1, 2, 3, ... number_of_ions
        """
    
        return tensor( [create(self.D_e) if j == i else qeye(self.D_e) 
                        for j in range(1, self.N_e+1)] 
                        + [ qeye(self.D_F) for j in range(1, self.N_F+1) ] )

    def ket(self, *state):
        """Define a pure quantum state, with first self.N enries the electronic states of ions,
        and the second self.N entries the motional states of ions.
        Note: Use projection() to redefine this method.
        """

        if len(state) == self.N_e+self.N_F:

            return tensor([basis(self.D_e, state[i-1]) 
                            for i in range(1, self.N_e+1)]
                            +[basis(self.D_F, state[i-1]) 
                            for i in range(self.N_e+1, self.N_e+self.N_F+1)])
        else:
            print("Number of arguments should be "+str(2*self.N))

    def thermal_dm(self, *args):
        average_phonon_num = args[0]
        state              = args[1:]
        if len(state) != self.N_e:
            print("Number of states must be equal to dimension of Hilbert space of electronic states.")

        elif average_phonon_num<self.D_F:
            return tensor( [basis(self.D_e, state[i-1])*basis(self.D_e, state[i-1]).dag() 
                            for i in range(1, self.N_e+1)] +
                           [thermal_dm(self.D_F, average_phonon_num) 
                            for i in range(1, self.N_F+1)] )
        else:
            print("Argument average_phonon_num must be less than Dimension of each mode's Fock space.")

            
    def coherent_dm(self, *args): 
        if len(args) != self.N_F+self.N_e:
            print("Number of argument must be number of motional modes plus number of ions.")
        else:
            alpha = args[0:self.N_F]
            state = args[self.N_F:]
            if len(state) != self.N_e:
                print("Number of states must be equal to number of ions.")

            else:
                return tensor( [basis(self.D_e, state[i-1])*basis(self.D_e, state[i-1]).dag() 
                                for i in range(1, self.N_e+1)] +
                               [coherent_dm(self.D_F, alpha[i-1]) 
                                for i in range(1, self.N_F+1)] )

            
    def eprojection(self, state):
        """Projection operator |state><state|
        where state is the state number 
        going from 0 to self.N_e-1.
        """
        if len(state) == self.N_e:

            return tensor(tensor([basis(self.D_e, state[i-1])*basis(self.D_e, state[i-1]).dag() 
                            for i in range(1, self.N_e+1)]), tensor([qeye(self.D_F) for i in range(self.N_F)]))
        else:
            print("Number of arguments should be ")

    def coupling(self, state1, state2):
        """Projection operator |state1><state2|
        where state is the state number 
        going from 0 to self.N_e-1.
        use state -1 for projection of identity onto correpsonding ion
        """
        if len(state1) == self.N_e and len(state2) == self.N_e:

            return tensor([basis(self.D_e, state1[i-1])*basis(self.D_e, state2[i-1]).dag()
                                  if ((state1[i - 1] >= 0) and (state2[i - 1] >=0))
                                  else qeye(self.D_e) for i in range(1, self.N_e + 1)] +
                                  [qeye(self.D_F) for i in range(1, self.N_F + 1)])
        else:
            print("Number of arguments should be ")

     
        

    def sigmaz_func(self, i):
        """Return sigmaz() operator acting on the electronic state of 
        i-th ion, where i = 1, 2, ... self.N
        Note: Use projection() to rewrite this.
        """
        if i>0 and i<=self.N_e:

            return tensor( [sigmaz() if j == i else qeye(2) 
                            for j in range(1, self.N_e+1) ] + 
                            [ qeye(self.D_F) for j in range(1, self.N_F+1)] )
        else:
            print("Error: Ion number must be between 1 and number of ions")

    def sigmax_func(self, i):
        """Return sigmax() operator acting on the electronic state of 
        i-th ion, where i = 1, 2, ... self.N
        Note: Use projection() to rewrite this.
        """
        if i>0 and i<=self.N_e:
            return tensor( [sigmax() if j == i else qeye(2) 
                            for j in range(1, self.N_e+1) ] + 
                            [ qeye(self.D_F) for j in range(1, self.N_F+1)] )
        else:
            print("Error: Ion number must be between 1 and number of ions")

    def sigmay_func(self, i):
        """Return sigmay() operator acting on the electronic state of 
        i-th ion, where i = 1, 2, ... self.N
        Note: Use projection() to rewrite this.
        """
        if i>0 and i<=self.N_e:
            return tensor( [sigmay() if j == i else qeye(2) 
                            for j in range(1, self.N_e+1) ] + 
                            [ qeye(self.D_F) for j in range(1, self.N_F+1)] )
        else:
            print("Error: Ion number must be between 1 and number of ions")


    def id(self):
        """Return identity operator acting on the entire Hilbert space 
        of the problem.
        """
        return tensor( [qeye(self.D_e) for j in range(1, self.N_e+1) ] + 
                   [qeye(self.D_F) for j in range(1, self.N_F+1)] )


