import pennylane as qml
import numpy as np

from adiaframe.utils import fam_code_add

class PFrame:
    def __init__(self, n:int):
        assert n>0 and isinstance(n, int), "n must be a positive integer."
        self.n = n
        
        # frame
        
        # (IZII, IXII) ...
        num = 1
        frame = []
        for i in range(n):
            frame.append([[num, 0], [0, num]])
            num = 2*num
        self.frame = frame
        self.gates = []
    def reset_gate_log(self):
        self.gates = []
    def to_pennylane_circuit(self, replace_rz:None):
        pass
    def CX(self, i, j):
        s_i, sd_i = self.frame[i]
        s_j, sd_j = self.frame[j]
        
        sd_i_new = fam_code_add()
        self.frame[i] = [s_i, sd_i_new]
    def H(self, i):
        s_i, sd_i = self.frame[i]
        self.frame[i] = [sd_i, s_i]
    def S(self, i):
        s_i, sd_i = self.frame[i]
        
        
            
                

def initial_frame(n:int):
    """Initial frame generating function

    Args:
        n (int): Number of wires
    """
    frame = []    
    for i in range(n):
        zstr = (i)*"I" + "Z" + (n-i-1)*"I"
        xstr = (i)*"I" + "X" + (n-i-1)*"I"
        frame.append(zstr, xstr)
    return frame
def apply_term_to_circuit(frame, pauli):
    
    #...
    
    return frame 