
class Experiment:

    """""""Define Experimental Paramaters for Raman Gates""""""""""
    B_C = Boolean Identifying whether Bichromatic Field
    N_m = Number of Modes
    Dicke = Ld Parameter of Modes [eta1, eta2, ...]
    D_f = Dimension of Fock Spaces Considered
    F_f = Frequency of oscillator spaces (MHz)
    Dets = Red Blue Detuning (MHz)
    F_s = Field (MHz) [Pump, Stokes]
    P = Relative Phase of Stokes Beams [P1, P2] if BC
    F_C = Boolean Identifying whether Frequency Comb
    F_t = Blue Red tone difference (positive)
    S = width of pulse
    T = Pulse Peri of comb
    P = Relative phase of Combs if comb
    tstep = Integration Step
    tmax = Simulation Time Length
    
    
    Author: Daniel Murphy
            Georgia Institute of Technology
    """""""""""""""""""""""""""""""""""""""""""""""""""""

    def __init__(self, B_C = False, N_m = 0, D_f = 0, F_f = [0], Dets = [0, 0],
                 F_s = [[[0, 0, 0]], [[0, 0, 0],[0, 0, 0]]],
                 P = 0, Dicke = [0], F_C = False, F_t = 0, S = 0, T = 0,
                 tstep = 0, t = 0):

                self.B_C = B_C
                self.N_m = N_m
                self.D_f = D_f
                self.Dets = Dets
                self.F_s = F_s
                self.F_f = F_f
                self.P = P
                self.Dicke = Dicke
                self.F_C = F_C
                self.F_t = F_t
                self.S = S
                self.T = T
                self.tstep = tstep
                self.tmax = t



    def interactionFrequencies(self):

        v = []

        for i in range(0, self.N_m):

            v.append([0, self.F_f[i], -1*self.F_f[i]])

        return v


    


    
    
