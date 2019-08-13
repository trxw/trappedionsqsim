from __future__ import division, absolute_import, print_function, unicode_literals

import qutip as qtp
import numpy as np
from scipy import *

#A library for automatical generation of quantum operators 
#Created by Omid Khosravani
#Duke University and Georgia Institute of Technology 
#2017-2019 
 

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
                number_of_motional_modes = 0,
                dim_of_each_Fock_space = 0):

        self.N_e = number_of_ions
        self.D_F = dim_of_each_Fock_space
        self.N_F = number_of_motional_modes
        self.D_e = dim_of_electronic_states_space
        #Verify the range of valid values:
        self._error(self.N_e, 1, float("inf"), "Number of ions")
        self._error(self.D_e, 2, float("inf"), "Dimension of Electronic Hilbert space of each ion")
        self._error(self.N_F, 0, self.N_e, "Number of Fock states")
        if self.N_F>0:
            self._error(self.D_F, 2, float("inf"), "Dimension of Fock space of the motional state of each ion")

    def a(self, mode_num):
        """
        params
        mode_num = 1, 2, 3, ... number_of_motional_modes

        return an annihiliation operator on the mode_num-th Fock space.
        """
        if self.N_F == 0:
            raise ValueError("Fock space is not defined")
        self._error(mode_num, 1, self.N_F, "Mode number")

        return qtp.tensor( [qtp.qeye(self.D_e)
            for i in range(self.N_e)] + 
            [ qtp.destroy(self.D_F) if j == mode_num else qtp.qeye(self.D_F) for j in range(1, self.N_F+1) ] )

    def ad(self, mode_num):
        """
        params
        mode_num = 1, 2, 3, ... number_of_motional_modes

        return an creation operator on the mode_num-th Fock space.
        """
        return self.a(mode_num).dag()

    def _error(self, i,min_lim, max_lim, st):

        if not isinstance(i,int) or i<min_lim or i>max_lim:
            raise ValueError(st+" "+"must be greater than or equal to "+str(min_lim)+ " and " + "less than or equal to "+str(max_lim))


    def sm(self, ion_num):
        """Return an annihiliation operator on the electronic state of i-th ion.
        i = 1, 2, 3, ... number_of_ions
        """

        self._error(ion_num, 1, self.N_e, "Ion number")

        
        return qtp.tensor( [qtp.destroy(self.D_e) if j == ion_num else qtp.qeye(self.D_e) 
                        for j in range(1, self.N_e+1)]
                        + [ qtp.qeye(self.D_F) for j in range(self.N_F) ] )

    
    def sp(self, ion_num):
        """Return a creation operator on the electronic state of i-th ion.
        i = 1, 2, 3, ... number_of_ions
        """
        self._error(ion_num, 1, self.N_e, "Ion number")

        
        return qtp.tensor( [qtp.create(self.D_e) if j == ion_num else qtp.qeye(self.D_e) 
                        for j in range(1, self.N_e+1)]
                        + [ qtp.qeye(self.D_F) for j in range(self.N_F) ] )


  

    def dm_pure(self, e_state_list, motional_state_list=[]):
        """Return a density matrix of a pure quantum state, with first self.N enries the electronic states of ions,
        and the second self.N entries the motional states of ions.
        """
        return qtp.tensor(self.state.ket(e_state_list, motional_state_list), self.state.ket(e_state_list, motional_state_list).dag())
    
    def ket(self, e_state_list, motional_state_list=[]):
        """Return a pure quantum state, with first self.N enries the electronic states of ions,
        and the second self.N entries the motional states of ions.
        """
        if len(e_state_list) != self.N_e:
            raise ValueError("Length of electronic_state_list must be equal to " + str(self.N_e) )
        if len(motional_state_list) != self.N_F:
            raise ValueError("Length of motional_state_list must be equal to " + str(self.N_F))

        self._error(len(motional_state_list), 0, self.N_F, "Length of motional state list")

        for i in e_state_list:
            self._error(i, 0, self.D_e-1, "Electronic state number")
        for j in motional_state_list:
            self._error(j, 0, self.D_F-1, "Motional state number")

        return qtp.tensor([qtp.basis(self.D_e, e_state_list[i]) 
                        for i in range(self.N_e)]  
                        +[qtp.basis(self.D_F, motional_state_list[j]) 
                        for j in range(self.N_F)])


    def eprojection(self, e_state_list, motional_state_list): 
        """Projection operator |state><state|
        where state is the states number list 
        go from 0 to self.D_e-1
        and 
        go from 0 to self.D_F-1
        """
        return qtp.tensor( self.ket(e_state_list, motional_state_list),  self.ket(e_state_list, motional_state_list).dag() )
    
    def coupling(self, state1, state2): #Not tested
        """Projection operator |state1><state2|
         where state1&2 each are the states number list 
        go from 0 to self.D_e-1
        and 
        go from 0 to self.D_F-1
        """
        if len(state1)==2 and len(state2) ==2:
            e_state_list1, motional_state_list1 = state1[0], state1[1]
            e_state_list2, motional_state_list2 = state2[0], state2[1]
        elif len(state1)==1 and len(state2) ==1:
            e_state_list1, motional_state_list1 = state1[0], []
            e_state_list2, motional_state_list2 = state2[0], []

        return qtp.tensor( self.ket(e_state_list1, motional_state_list1),  self.ket(e_state_list2, motional_state_list2).dag() )


    def sx(self, ion_num):
        """Return sigmax() operator acting on the electronic state of 
        i-th ion, where i = 1, 2, ... self.N
        """
        self._error(ion_num, 1, self.N_e, "Ion number")
        return qtp.tensor( [qtp.sigmax() if j == ion_num else qtp.qeye(2) for j in range(1, self.N_e+1) ] + 
                        [ qtp.qeye(self.D_F) for j in range(1, self.N_F+1)] )

    def sy(self, ion_num):
        """Return sigmay() operator acting on the electronic state of 
        i-th ion, where i = 1, 2, ... self.N
        """
        self._error(ion_num, 1, self.N_e, "Ion number")
        return qtp.tensor( [qtp.sigmay() if j == ion_num else qtp.qeye(2) for j in range(1, self.N_e+1) ] + 
                        [ qtp.qeye(self.D_F) for j in range(1, self.N_F+1)] )

    def sz(self, ion_num):
        """Return sigmaz() operator acting on the electronic state of 
        i-th ion, where i = 1, 2, ... self.N
        """
        self._error(ion_num, 1, self.N_e, "Ion number")
        return qtp.tensor( [qtp.sigmaz() if j == ion_num else qtp.qeye(2) for j in range(1, self.N_e+1) ] + 
                        [ qtp.qeye(self.D_F) for j in range(1, self.N_F+1)] )


    def id(self, ion_num=1):
        """Return identity operator acting on the entire Hilbert space 
        of the problem.
        """
        return qtp.tensor( [qtp.qeye(self.D_e) for j in range(self.N_e) ] + 
                   [qtp.qeye(self.D_F) for j in range(self.N_F)] )


class States(Operators):
    def __init__(self, number_of_ions = 1,
                dim_of_electronic_states_space = 2,
                number_of_motional_modes = 0,
                dim_of_each_Fock_space = 0):
        super(Operators, self).__init__(number_of_ions,
                dim_of_electronic_states_space,
                number_of_motional_modes,
                dim_of_each_Fock_space)



    def ket(self, e_state_list, motional_state_list=[]):
        """Return a pure quantum state, with first self.N enries the electronic states of ions,
        and the second self.N entries the motional states of ions.
        """
        if len(e_state_list) != self.N_e:
            raise ValueError("Length of electronic_state_list must be equal to " + str(self.N_e) )
        if len(motional_state_list) != self.N_F:
            raise ValueError("Length of motional_state_list must be equal to " + str(self.N_F))

        self._error(len(motional_state_list), 0, self.N_F, "Length of motional state list")

        for i in e_state_list:
            self._error(i, 0, self.D_e-1, "Electronic state number")
        for j in motional_state_list:
            self._error(j, 0, self.D_F-1, "Motional state number")

        return qtp.tensor([qtp.basis(self.D_e, e_state_list[i]) 
                        for i in range(self.N_e)]  
                        +[qtp.basis(self.D_F, motional_state_list[j]) 
                        for j in range(self.N_F)])

    '''

       
   
    def coherent_dm(self, *args): #Not tested
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

    '''