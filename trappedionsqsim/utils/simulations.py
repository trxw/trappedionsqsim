from __future__ import division, absolute_import, print_function, unicode_literals

import qutip as qtp
import numpy as np
from scipy import *
import types
from .operators import Operators 
#A library for simulation of arbitrarily long ion chains
#Created by Omid Khosravani, okhosravani@gatech.edu
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

        #self.Hamiltonian = Hamiltonian(self.N_e,
        #        self.D_e,
        #        self.N_F,
        #        self.D_F)
        self.gate = Gate(self.N_e,
                self.D_e,
                self.N_F,
                self.D_F)
        self.Hamiltonian_list = []
        self.time_evolve_list = []

        self.reset()
    """
    def add_hamiltonian(self, hamiltonian, t_arr):
        '''
        params:
        hamiltonian is an instance of Hamiltonian class
        t_arr is an array of time points during which the given Hamiltonian will be applied

        '''
        self.Hamiltonian_list += [hamiltonian]
        if len(self.time_evolve_list) != 0 and t_arr[0] != self.time_evolve_list[-1][-1]:
            raise ValueError("The first time stamp t1 should be equal to the last time stamp (t2) of previous t_list in time_evolve_list")  

        self.time_evolve_list += [t_arr]
    """

    def set_curr_state(self, psi): #tested
        self.curr_state = psi
    def get_curr_state(self): #tested
        return self.curr_state

    def evolve_spline(self, Hamiltonian_list, t_int, N_steps=None, save_to_states_list=True, c_ops = [], observable_list=[]): #tested
        '''
        params
        Hamiltonian_list is a list of  [ [Hamiltonian1, coef_function1], [Hamiltonian2, coef_function2], ...] some of which is the total Hamiltonian applied in t_arr
        t_int is the time-interval. If it hsa length 2 it specifies the first and last time [t1,t2]
        otherwise it must be an array will all the time stamps
        N_steps is the number of time steps to divide the time_interval corresponding to this Hamiltonian 
        turn each function into a spline, simulate time evolution and return output, a qutip object 
        verbose will print out the average value of time_step_isValid 

        '''
        if not (isinstance(t_int, list) or isinstance(t_int, np.ndarray)) or len(t_int) <2:
            raise ValueError("Time interval must be a list or array of length at least 2.")
        elif len(t_int) <2:
            raise ValueError("Time interval must be a list or array of length at least 2.")

        elif len(t_int) == 2:
            if N_steps is None:
                raise ValueError("Number of time steps must be specified.")
            elif isinstance(N_steps, int):
                if N_steps <2:
                    print("N_steps must be at list 2.")
                else:
                    t_arr = np.array(t_int[0], t_int[1], N_steps)

        
        elif len(t_int) <100:
            raise ValueError("Set t_int to at least 100 time steps.")

        elif len(t_int)>=100:
            print("t_arr is set")
            t_arr = t_int

        if t_arr[0] != self.curr_t:
            raise ValueError("First time stamp does not match the current time stamp curr_t.")


        if len(Hamiltonian_list)<1:
            raise ValueError("Set at least one Hamiltonian in Hamiltonian_list: [ [Hamiltonian1, coef_function1],..].")

        for H in Hamiltonian_list: 
            if not self.time_step_isValid(H[0], self.curr_state, t_arr): 
                raise ValueError("Time steps are too large.")

        if self.curr_state is None:
            raise ValueError("Initial state not set.")

        #self.states_list += [self.curr_state]


        output = qtp.mcsolve(Hamiltonian_list, self.curr_state, t_arr, c_ops, observable_list)
        self.set_curr_state(output.states[-1])


        if  save_to_states_list:
            if len(self.states_list) != 0 :
                self.states_list = self.states_list[:-1] + output.states
            else:
                self.states_list = output.states

            if len( self.time_arr )  != 0:
                self.time_arr = np.concatenate( (self.time_arr[:-1], t_arr ) )
            else: 
                self.time_arr = t_arr

            self.output_list += [output]

        self.curr_t = t_arr[-1]


        return output

    
    def time_step_isValid(self, Hamiltonian, psi, t_arr):#t1,t2, steps_num): #tested
        '''Check to see if time step is small enough 
        '''
        print("Checking time step in time-evolution:")
        print("Transition matrix element after first time-step: ", qtp.expect(Hamiltonian, psi)*(t_arr[1]-t_arr[0]) /np.pi) 
        return 50*int(abs(qtp.expect(Hamiltonian, psi)*(t_arr[1]-t_arr[0]) /np.pi) ) < 1.


    def apply_gate(self, gate_string,  save_to_states_list=True): #tested
        '''
        params
        set qubit gate S operator from Pauli basis {X, Y, Z, I} for 1 qubit, 
        and S operator from {XX, XY, XZ, XI, YX, YY, YZ, YI, ...} for 2 qubits,
        and so on. 
        '''

        self.curr_state = self.gate.get_pauli(gate_string) * self.curr_state
        if  save_to_states_list:
            self.states_list += [self.curr_state]


    def reset(self):
        #self.set_curr_state(None)
        self.states_list = []
        self.output_list = []
        self.time_arr = np.array([])
        self.curr_t = 0. #Current time in simulation


class Gate(Operators): #tested
    def __init__(self, number_of_ions = 1,
                dim_of_electronic_states_space = 2,
                number_of_motional_modes = 0,
                dim_of_each_Fock_space = 0): 
        #super().__init__(number_of_ions,
        #                dim_of_electronic_states_space,
        #                number_of_motional_modes,
        #                dim_of_each_Fock_space)
        self.ops = Operators(number_of_ions,
                        dim_of_electronic_states_space,
                        number_of_motional_modes,
                        dim_of_each_Fock_space)
        self.gates_dic = {'X': self.ops.sx, 'Y':self.ops.sy, 'Z':self.ops.sz, 'I':self.ops.id
                            }

    def get_pauli(self, gate_string):
        '''
        params
        gate_string is a string of Pauli gates applied to each ion 1,2,3...
        All upper-case as in gates_dic
        return corresponding multiqubit Pauli gate
        '''
        for g in gate_string:
            if g not in self.gates_dic:
                raise ValueError("Wrong gate specified.")

        if len(gate_string) != self.ops.N_e:
            raise ValueError("Length of gate_string must be equal to the number of qubit") 
        #print( self.gates_dic[gate_string[0]] )
        pauli = self.ops.id()
        for i_num in range(len(gate_string)):
            pauli *= self.gates_dic[gate_string[i_num]](i_num+1) 
        return pauli


"""
class Hamiltonian(Operators): #Not tested
    def __init__(self, number_of_ions = 1,
                dim_of_electronic_states_space = 2,
                number_of_motional_modes = 0,
                dim_of_each_Fock_space = 0): 
        super().__init__(number_of_ions,
                        dim_of_electronic_states_space,
                        number_of_motional_modes,
                        dim_of_each_Fock_space)
        self.rotation_functions = []

    def set_pauli_rotation_funcs(self, *func_list):
        '''
        params
        set qubit rotations with functions fi as coefficients  
        for S operator from Pauli basis {X, Y, Z, I} for 1 qubit, 
            --> func_list list of 4 functions expected
        and S operator from {XX, XY, XZ, XI, YX, YY, YZ, YI, ...} for 2 qubits,
            --> func_list list of 16 functions expected
        and so on. 
        '''
        if len(func_list) != (self.D_e**2)**self.N_e:
            raise ValueError("All Paulis must be specified")
      
        for f in func_list:
            self.rotation_functions += [f]


    def set_Hamiltonian(self):
        pass



class Channel(Operators):
    def __init__(self, number_of_ions = 1,
                dim_of_electronic_states_space = 2,
                number_of_motional_modes = 0,
                dim_of_each_Fock_space = 0): 
        super().__init__(number_of_ions,
                        dim_of_electronic_states_space,
                        number_of_motional_modes,
                        dim_of_each_Fock_space)
"""