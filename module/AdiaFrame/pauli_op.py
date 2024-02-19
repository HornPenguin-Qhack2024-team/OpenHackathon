from collections import OrderedDict
from itertools import combinations

import pennylane as qml
from pennylane import qchem
from pennylane import numpy as np
import pandas as pd


def hamiltonian_processing(molecule_name, bondlength, directory="dataset_processed", basis = "STO-3G"):
    data = qml.data.load("qchem", molname=molecule_name, basis=basis, bondlength = bondlength)[0]
    hamiltonian = data.hamiltonian
    
    h_basis = hamiltonian.terms()[1]
    qubits = len(hamiltonian.wires)


    pauli_strs = []
    for p in h_basis:
        pauli_strs.append(Pauli2String(p, qubits))
    
    z_f  = []
    x_f = []

    p_strs = {}
    for ps in pauli_strs:
        p_strs[ps] = PauliDecompose(ps)
        type_i, z_str, x_str = PauliDecompose(ps)

        if type_i ==0:
            z_f.append(z_str)
        elif type_i == 1:
            x_f.append(x_str)
        else:
            z_f.append(z_str)
            x_f.append(x_str) 
    
    df = pd.DataFrame.from_dict(p_strs, orient="index", columns = ["type", "Z", "X"])
    df= df.sort_values("type")
    df = df.replace("", None)
    df["Znum"] = df["Z"].apply(lambda x: int(x.replace("Z", "1").replace("I", "0"),2) if x is not None else 0).astype(int)
    df["Xnum"] = df["X"].apply(lambda x: int(x.replace("X", "1").replace("I", "0"),2) if x is not None else 0).astype(int)
    df= df.sort_values(by=["type", "Znum", "Xnum"])
    df.reset_index(inplace=True, names="Pstring")

    zdf = df[["Znum"]].drop_duplicates()
    xdf = df[["Xnum"]].drop_duplicates()

    zdf = zdf.reset_index(drop=True).reset_index(names="zindex")
    xdf = xdf.reset_index(drop=True).reset_index(names="xindex")

    df = df.merge(zdf, left_on="Znum", right_on='Znum')
    df = df.merge(xdf, left_on="Xnum", right_on='Xnum')
    df= df.sort_values(by=["type", "Znum", "Xnum"]).reset_index()

    # Done Df
    df.to_csv(f"{directory}/{molecule_name}-{bondlength}_node.csv",encoding="UTF-8")

    edge_df = pd.DataFrame(combinations(df["Pstring"].values ,2), columns=['source', 'target'])

    edge_df = edge_df.merge(df[["Pstring", "Znum", "Xnum"]], how="left", left_on="source", right_on='Pstring').drop("Pstring", 1)
    edge_df = edge_df.merge(df[["Pstring", "Znum", "Xnum"]], how="left", left_on="target", right_on='Pstring').drop("Pstring", 1)
    edge_df["commute"] = edge_df[["Znum_x", "Xnum_x", "Znum_y", "Xnum_y"]].apply(_commute, axis=1)
    edge_df["weight"] = edge_df[["Znum_x", "Xnum_x", "Znum_y", "Xnum_y"]].apply(_weight, axis=1)/qubits

    edge_df.to_csv(f"{directory}/{molecule_name}-{bondlength}_edge.csv",encoding="UTF-8")

def Pauli2String(pauli, num_wires=None): 
    # Convert the given Pennylane PauliBasis to Pauli string
    types = "".join(pauli.name).replace("Pauli", "").replace("Identity", "I")
    
    acting_wires = list(pauli.wires)
    if num_wires is None:
        num_wires = pauli.num_wires
        
    pstring = num_wires*["I"]
    for i, p in enumerate(acting_wires):
        pstring[p] = types[i]
    return "".join(pstring)

def PauliDecompose(pauli_string): 
    # Decompose product string into x, z family string 
    # 0: Z family
    # 1: X family
    # 2: product string
    st_type = 2
    st1 = "" # z
    st2 = "" # x
    pro_st = ""
    if "Y" not in pauli_string:
        if "X" in pauli_string and "Z" in pauli_string:
            pro_st = pauli_string
        elif "X" in pauli_string:
            st_type = 1
            st2 = pauli_string
        else:
            st_type = 0
            st1 = pauli_string
    else:
        pro_st = pauli_string
    
    # Decompose pro_st
    if st_type ==2:
        z_st = []
        x_st = []
        for s in pro_st:
            if s =="Y":
                x_st.append("X")
                z_st.append("Z")
            elif s =="X":
                x_st.append("X")
                z_st.append("I")
            elif s =="Z":
                x_st.append("I")
                z_st.append("Z")
            elif s =="I":
                x_st.append("I")
                z_st.append("I")
        st1 = "".join(z_st)
        st2 = "".join(x_st)
    return (st_type, st1, st2)

def _commute(s):
    a = bin(s[0] & s[3]).count("1")%2
    b = bin(s[1] & s[2]).count("1")%2
    return a == b
def _weight(s):
    binary_z = bin(abs((s[0] & s[1])^ (s[2]&s[3])))
    binary_x = bin(abs((s[1]^s[3])))
    return (binary_z + binary_x).count("1")
