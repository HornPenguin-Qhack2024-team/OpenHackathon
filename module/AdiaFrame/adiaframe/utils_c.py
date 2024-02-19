import ctypes
from pathlib import Path
from platform import system
# Find a path of shared library considering os dependency. 

platform_os = system()
libname = "libutilsc"
extname = ".dll" if "Windows" in platform_os else ".so"
path = Path(__file__).parent/"c_modules"/(libname+extname)
lib = ctypes.CDLL(str(path.absolute()))

# Define the argument and return types for each C function
lib.bit_count.argtypes = [ctypes.c_uint]
lib.bit_count.restype = ctypes.c_uint

lib.xz_coef_pow.argtypes = [ctypes.c_uint, ctypes.c_uint]
lib.xz_coef_pow.restypes = ctypes.c_uint

lib.insert_zeros.argtypes = [ctypes.c_uint]
lib.insert_zeros.restype = ctypes.c_uint

lib.encode_xzcode.argtypes = [ctypes.c_uint, ctypes.c_uint]
lib.encode_xzcode.restype = ctypes.c_uint

lib.decode_pcode.argtypes = [ctypes.c_uint, ctypes.c_size_t, ctypes.POINTER(ctypes.c_uint), ctypes.POINTER(ctypes.c_uint)]
lib.decode_pcode.restype = None

# Python wrappers for the C functions
def bit_count(n):
    return lib.bit_count(n)
def pauli_xz_product_coef(x, z):
    return 1j**lib.xz_coef_pow(x,z)
def insert_zeros(n):
    return lib.insert_zeros(n)
def encode_xzcode(x, z):
    return lib.encode_xzcode(x, z)
def decode_pcode(p_int, p_len):
    x = ctypes.c_uint()
    z = ctypes.c_uint()
    lib.decode_pcode(p_int, p_len, ctypes.byref(x), ctypes.byref(z))
    return x.value, z.value