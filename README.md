# OpenHackathon

2024 Qhack Open Hackathon Challenge

Team name: HornPenguin

Members: name(nickname)

* HYUNSEONG Kim(Hornpenguin, Raspbie)
* Guillaume(Guillaume)
* Pech(Pech)
* Sharma(Yamagi1729)


## Project Name

## Challenge

Target challenge: Bridging the gap, Seeing the future

- A matter of taste: How to find a phase transition of the given spin system provided by Pennylane.
- Bridging the gap: Find a way to estimate spectral gap of the given system.
- Preparing for battle: Search or come up an innovative technique to build the given initial state using amplitude embedding.
- Seeing the future: Fancy way to simulate evolution process of a quantum system.
- The sound of silence: Find a way to reduce the number of qubits in quantum algorithm.

## Proejct description

Imaginary time approach implementation on Pennylane framework (Guillaume proposed).


Additional comment: 
The imaginary time evolution is an non-unitary operator. 
However, we can simulate the non-unitary operator with larger unitary opeartor, 
since sub-matrix of unitary matrix is not unitary. 
Therefore, $\exp(-H \tau) \approx \exp(-i A t)$ and restrict the result. 
Most of the implementation based on VQE 

### Questions

1. How imaginary time approach work and why it is better than the standard one.
2. 

