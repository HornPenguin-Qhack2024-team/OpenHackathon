from adiaframe import *

# default bit functions
def krons(*oper_list): # Operator Kronecker delta
    if len(oper_list) == 1:
        oper_list = oper_list[0]
    return reduce(np.kron, oper_list)
def frobenius_inner(A, B): # Frobenius inner product.
    n, n2 = A.shape
    return np.trace((A.conj().T)@B)/(n)

def commute_reggio(pa:Tuple[int, int], pb:Tuple[int, int]):
    #Reggio et al, Fast Partitioning of Pauli Strings into Commuting Families for Optimal Expectation Value Measurements of Dense Operators, 2023-06-07
    nx_a, nz_a = pa
    nx_b, nz_b = pb
    
    a = bin(nx_a & nz_b).count("1")%2
    b = bin(nx_b & nz_a).count("1")%2
    return a==b
    
def commute_reggio_df(s):
    a = bin(s.iloc[0] & s.iloc[3]).count("1")%2
    b = bin(s.iloc[1] & s.iloc[2]).count("1")%2
    return a == b

def integer_order_map(int_list):
    sorted_unique = np.unique(np.array(int_list))
    return {num: idx for idx, num in enumerate(sorted_unique)}

def get_coef(x_str, z_str): 
    # i coefficient in construction of general pauli-element from XZ elements.
    # Use this function in python module
    # The below bitwise implementations are slower than the current function.
    # They are just for further C implementation.
    n = len(x_str)
    x_str = x_str.replace("X", "1")
    x_str = x_str.replace("I", "0")
    z_str = z_str.replace("Z", "1")
    z_str = z_str.replace("I", "0")
    
    x_int = int(x_str, 2)
    z_int = int(z_str, 2)
    return get_coef_bin(x_int, z_int, n)
def get_coef_bin(x_int:int, z_int:int, n:str):
    y_pos = x_int&z_int
    y_pos = format(x_int&z_int, f"0{n}b")
    z_pos = format((x_int|z_int) - x_int, f"0{n}b")
    x_pos = format((x_int|z_int) - z_int, f"0{n}b")

    g_str = []
    for x,y,z in zip(x_pos, y_pos, z_pos):
        if x==y and y==z:
            g_str.append("I")
        elif x== "1":
            g_str.append("X")
        elif y == "1":
            g_str.append("Y")
        else:
            g_str.append("Z")
    return 1j**y_pos.count("1"), "".join(g_str)

# Bit operators=========================================
# pauli coef calculation
def pauli_xz_product_coef(x_int, z_int):
    return 1j**(bit_count(x_int&z_int))
def bit_count(n:int):
    #Brian Kernighan’s Algorithm.
    num = 0
    while n:
        n &= n-1
        num+=1
    return num

# Calculate string from ints
int_pchar = ["I", "Z", "X", "Y"]
def pstr_from_xz(x_int, z_int):
    z_modi = insert_zeros_in_gaps(z_int)
    x_modi = insert_zeros_in_gaps(x_int)
    x_modi <<= 1

    p_int = x_modi + z_modi 
    # In binary representation: (00)(10)(10)(11) form
    # 00:I, 10: X 01: Z, 11: Y

    # Get length of str
    len_p = 0
    tem = p_int
    while tem:
        len_p +=1
        tem >>=1
    len_p += len_p&1
    pstr = len_p*['']
    i = 1
    while p_int >0:
        p = p_int & 3
        p_int >>= 2
        pstr[-i] = int_pchar[p]
        i+=1
    return "".join(pstr)
def insert_zeros_in_gaps(n):
    result = 0
    bit_position = 0

    while n > 0:
        # Isolate the rightmost bit
        rightmost_bit = n & 1
        # Shift the bit to its new position
        result |= rightmost_bit << (bit_position << 1)
        # Move to the next bit
        n >>= 1
        bit_position += 1

    return result