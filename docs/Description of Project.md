



## XZ family code representation of Pauli term.


In Reggio et al, a Pauli-term is represented as (x, z) int tuple
to determine the commutative of two given terms.

Adopting such convention, in the project we implemented 
Pauili-group algebraic structure with well-designed binary operation.
This is an imitation of $Z_3$ group strucutre with 2 bit binary
and expansion.

This code represents location of the term on Latin matrix.
For example, length 1 Pauli-string terms, "I, X, Y, Z" are in XZ code.

- "I" = (0, 0)
- "X" = (1, 0)
- "Y" = (1, 1)
- "Z" = (0, 1)

| I | Z|
|:--:|:--:|
|X| Y|

for n-length Pauli strings, 

```
"IXYYZ" = (IIZZZ)+(IXXXI)
0b(IXYYZ) = (0b0001111110) = (0b|00|10|11|11|01)
=(0b-0-0-1-1-1) + (0b0-1-1-1-0-)
=(0b00111) + (0b01110)
```
### Synthesis

The synthesis of two term simply constructed as XOR bit operation 
of two coordinate.

```
p1 = [n1x, n1z]
p2 = [n2x, n2z]
p3 = [n3x, n3z] = p1 + p2 = [n1z^n2z, n1x^n2x],
```

### Commutation

Reggio et al, method

```
p1 = [n1x, n1z]
p2 = [n2x, n2z]
a = bin(nx_1 & nz_2).count("1")%2
b = bin(nx_2 & nz_1).count("1")%2

a==b #A of: Do they commute?
```