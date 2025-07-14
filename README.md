Polytime simulator of the Clifford modulo phase group.

### How to build this
```
make
make interface
```

### Input format
Number of qubits in the state, followed by a sequence of gates. Qubits are numbered from 0 to n-1.
```
n # number of qubits
H targ # apply H to targ
X targ # apply X to targ
Z targ # apply Z to targ
S targ # apply S to targ
CNOT src targ # apply CNOT from src to targ
SWAP src targ # apply SWAP from src to targ
CZ src targ # apply CZ from src to targ
```

### Example interaction
The -v flag prints the state in a verbose (human-readable) format.

```
./build/interface -v
2
H 0
(EOF)
Phase polynomial matrix:
0 0 
0 0 
Affine pair:
A:
0 1 
B: 0 
```

```
./build/interface -v
2
H 0
(EOF)
(2,00,0000,01,0)
```
