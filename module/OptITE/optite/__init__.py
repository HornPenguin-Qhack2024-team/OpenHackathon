import Pennylane as qml



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