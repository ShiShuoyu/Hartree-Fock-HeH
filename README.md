# Hartree-Fock-HeH
The second homework of advanced quantum mechanics assigned by Prof. Li (xzli@pku.edu.cn).

### Method
Calculate the ground state energy of helium hydride ion (HeH+) by **Hartree-Fock** method. The interatomic distance set to the precise value R = 1.4632 a.u.

### Basis set
STO-3G minimal basis set from [Basis Set Exchange](https://www.basissetexchange.org/)
``` NWChem Format
#----------------------------------------------------------------------
# Basis Set Exchange
# Version 0.11
# https://www.basissetexchange.org
#----------------------------------------------------------------------
#   Basis set: STO-3G
# Description: STO-3G Minimal Basis (3 functions/AO)
#        Role: orbital
#     Version: 1  (Data from Gaussian09)
#----------------------------------------------------------------------


BASIS "ao basis" SPHERICAL PRINT
#BASIS SET: (3s) -> [1s]
H    S
      0.3425250914E+01       0.1543289673E+00
      0.6239137298E+00       0.5353281423E+00
      0.1688554040E+00       0.4446345422E+00
#BASIS SET: (3s) -> [1s]
He    S
      0.6362421394E+01       0.1543289673E+00
      0.1158922999E+01       0.5353281423E+00
      0.3136497915E+00       0.4446345422E+00
END
```
which uses three gaussians to fit the 1S orbit of H and He.