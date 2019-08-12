from __future__ import division, absolute_import, print_function, unicode_literals

from qutip import *
import numpy as np
from scipy import *
import types
from trappedionsqsim.utils.operators import Operators
#A library for simulation of arbitrarily long ion chains
#Created by Omid Khosravani, 
#Duke University and Georgia Institute of Technology 
#2017-2019 


class Simulation(Operators):
    def __init__(self, number_of_ions = 1,
                dim_of_electronic_states_space = 2,
                number_of_motional_modes = 0,
                dim_of_each_Fock_space = 0, 
                LD_param = .0,
                LD_order = .0):
        self.LD_param = LD_param
        self.LD_order = LD_order
        super(Simulation, self).__init__(number_of_ions,
                dim_of_electronic_states_space,
                number_of_motional_modes,
                dim_of_each_Fock_space)
        self.Hamiltonian = Hamiltonian(self.N_e,
                self.D_e,
                self.N_F,
                self.D_F)

    def add_hamiltonian(self, hamiltonian):
        self.Hamiltonian()


    def one_qubit_rotation(self, f: types.LambdaType, g: types.LambdaType, h: types.LambdaType=None):
        return (f(0))


class Hamiltonian(Operators):
    def __init__(self, number_of_ions = 1,
                dim_of_electronic_states_space = 2,
                number_of_motional_modes = 0,
                dim_of_each_Fock_space = 0): 
        super().__init__(number_of_ions,
                        dim_of_electronic_states_space,
                        number_of_motional_modes,
                        dim_of_each_Fock_space)

    def Hamiltonian(self, channel):
        pass
