
import numpy as np


class ClebschGordonCalc:

    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    Computes the clebsch gordon coefficient <F1, m1| 1, m, F2, m2>
    Used for computing matrix elements of rank 1 spherical tensor operators from reduced
    matrix elements
    
    Author: Daniel Murphy
            Georgia Institute of Technology
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    def __init__(self):
        
        self.data = 0
    
    def CGCoeff(self, f1, m1, f2, m2, m):
        
        if ((f1 == 0) and (f2 == 0)):
            return 0
        elif (np.abs(f1-f2) > 1):
            return 0
        elif (m1 != (m2 + m)) :
            return 0
        elif ((f1 == f2+1) and (m == 1)):
            return np.sqrt(((f2 + m1)*(f2 + m1 + 1))/((2*f2 + 1)*(2*f2 + 2)))
        elif ((f1 == f2+1) and (m == 0)):
            return np.sqrt(((f2 - m1 + 1)*(f2 + m1 + 1))/((2*f2 + 1)*(f2 + 1)))
        elif ((f1 == f2+1) and (m == -1)):
            return np.sqrt(((f2 - m1)*(f2 - m1 + 1))/((2*f2 + 1)*(2*f2 + 2)))/(np.sqrt(2*f1 + 1))
        elif ((f1 == f2) and (m == 1)):
            return -1*np.sqrt(((f2 + m1)*(f2 - m1 + 1))/(2*f2*(f2 + 1)))
        elif ((f1 == f2) and (m == 0)):
            return (m1/(np.sqrt(f2*(f2+1))))/(np.sqrt(2*f1 + 1))
        elif ((f1 == f2) and (m == -1)):
            return np.sqrt(((f2 - m1)*(f2 + m1 + 1))/(2*f2*(f2 + 1)))/(np.sqrt(2*f1 + 1))
        elif ((f1 == f2-1) and (m == 1)):
            return np.sqrt(((f2 - m1)*(f2 - m1 + 1))/(2*f2*(2*f2 + 1)))/(np.sqrt(2*f1 + 1))
        elif ((f1 == f2-1) and (m == 0)):
            return -1*np.sqrt(((f2 - m1)*(f2 + m1))/(f2*(2*f2 + 1)))/(np.sqrt(2*f1 + 1))
        elif ((f1 == f2-1) and (m == -1)):
            return np.sqrt(((f2 + m1)*(f2 + m1 + 1))/(2*f2*(2*f2 + 1)))/(np.sqrt(2*f1 + 1))


        

