

## Cas cible 1 - many_spehre.sh

### 16 MPI processes

| solver | preconditionner | converged | iterations | residu | time |
| :-------------|:---------------|:---------------:|:---------------:|:---------------:| ---------------: |
| HyprePCG | ANY | &#10004; | 1 | 2.96608e-11 | 2.28593e+09 |
| HyprePCG | HypreBoomerAMG | &#10004; | 2 | 4.33152e-14 | 3.58924e+09 |
| HyprePCG | HypreILU | &#10004; | 2 | 4.30336e-14 | 3.98486e+09 |
| HyprePCG | HypreEuclid | &#10004; | 1 | 2.94184e-11 | 1.58152e+11 |
| HyprePCG | HypreDiagScale | &#10004; | 2 | 4.33098e-14 | 2.18682e+09 |
| HypreFGMRES | ANY | &cross; | &cross; | &cross; | &cross; |
| HypreFGMRES | HypreBoomerAMG | &#10004; | 1 | 3.01325e-11 | 2.68551e+09 |
| HypreFGMRES | HypreILU | &#10004; | 1 | 3.22595e-11 | 4.58519e+09 |
| HypreFGMRES | HypreParaSails | &#10004; | 1 | 3.20752e-11 | 3.99365e+09 |
| HypreFGMRES | HypreDiagScale | &#10004; | 4 | 1.51959e-11 | 4.80448e+10 |
| HypreGMRES | ANY | &cross; | &cross; | &cross; | &cross; |
| HypreGMRES | HypreBoomerAMG | &#10004; | 1 | 3.01323e-11 | 2.65321e+09 |
| HypreGMRES | HypreILU | &#10004; | 1 | 3.22599e-11 | 4.55605e+09 |
| HypreGMRES | HypreParaSails | &#10004; | 1 | 3.20753e-11 | 3.84917e+09 |
| HypreGMRES | HypreDiagScale | &#10004; | 4 | 5.21728e-14 | 4.60127e+10 |
| MUMPSSolver | ANY | &#10004; | 6 | 4.32443e-14 | 7.23665e+10 |
| CGSolver | ANY | &#10004; | 1 | 1.15378e-12 | 1.61118e+09 |
| CGSolver | HypreParaSails | &cross; | &cross; | &cross; | &cross; |
| CGSolver | HypreDiagScale | &#10004; | 1 | 7.57415e-12 | 1.16162e+09 |
| BiCGSTABSolver | ANY | &#10004; | 1 | 1.42026e-12 | 4.61881e+09 |
| BiCGSTABSolver | HypreBoomerAMG | &#10004; | 1 | 7.8214e-13 | 2.60989e+09 |
| BiCGSTABSolver | HypreILU | &#10004; | 1 | 1.13276e-12 | 3.62895e+09 |
| BiCGSTABSolver | HypreEuclid | &#10004; | 1 | 2.00048e-12 | 1.58506e+11 |
| BiCGSTABSolver | HypreParaSails | &#10004; | 1 | 7.17439e-13 | 2.77627e+09 |
| BiCGSTABSolver | HypreDiagScale | &#10004; | 1 | 2.12864e-12 | 2.4907e+09 |
| MINRESSolver | ANY | &#10004; | 1 | 4.74096e-12 | 1.69667e+09 |
| MINRESSolver | HypreBoomerAMG | &#10004; | 1 | 4.38625e-12 | 1.88573e+09 |
| MINRESSolver | HypreILU | &#10004; | 1 | 5.51621e-12 | 2.1444e+09 |
| MINRESSolver | HypreEuclid | &#10004; | 1 | 2.81607e-12 | 1.56406e+11 |
| MINRESSolver | HypreParaSails | &cross; | &cross; | &cross; | &cross; |
| MINRESSolver | HypreDiagScale | &#10004; | 1 | 8.01382e-12 | 1.22562e+09 |

### 128 MPI processes

| solver | preconditionner | converged | iterations | residu | time |
| :-------------|:---------------|:---------------:|:---------------:|:---------------:| ---------------: |
| HyprePCG | ANY | &#10004; | 1 | 3.19848e-11 | 4.38711e+08 |
| HyprePCG | HypreBoomerAMG | &#10004; | 2 | 4.40356e-14 | 8.71224e+08 |
| HyprePCG | HypreILU | &#10004; | 2 | 4.42119e-14 | 3.46134e+08 |
| HyprePCG | HypreEuclid | &#10004; | 2 | 4.42259e-14 | 7.5854e+10 |
| HyprePCG | HypreDiagScale | &#10004; | 2 | 4.4217e-14 | 2.64795e+08 |
| HypreFGMRES | ANY | &cross; | &cross; | &cross; | &cross; |
| HypreFGMRES | HypreBoomerAMG | &#10004; | 1 | 2.84829e-11 | 7.88338e+08 |
| HypreFGMRES | HypreILU | &#10004; | 1 | 3.0718e-11 | 5.93981e+08 |
| HypreFGMRES | HypreParaSails | &#10004; | 1 | 3.20754e-11 | 9.52163e+08 |
| HypreFGMRES | HypreDiagScale | &#10004; | 3 | 2.83674e-11 | 4.86338e+09 |
| HypreGMRES | ANY | &cross; | &cross; | &cross; | &cross; |
| HypreGMRES | HypreBoomerAMG | &#10004; | 1 | 2.84829e-11 | 7.8814e+08 |
| HypreGMRES | HypreILU | &#10004; | 1 | 3.07182e-11 | 5.97464e+08 |
| HypreGMRES | HypreParaSails | &#10004; | 1 | 3.20755e-11 | 5.94519e+08 |
| HypreGMRES | HypreDiagScale | &#10004; | 3 | 9.45393e-13 | 4.83148e+09 |
| MUMPSSolver | ANY | &#10004; | 1 | 1.8629e-13 | 9.23528e+09 |
| CGSolver | ANY | &#10004; | 1 | 1.18624e-12 | 2.00643e+08 |
| CGSolver | HypreParaSails | &cross; | &cross; | &cross; | &cross; |
| CGSolver | HypreDiagScale | &#10004; | 1 | 7.67158e-12 | 1.45991e+08 |
| BiCGSTABSolver | ANY | &#10004; | 1 | 1.62537e-12 | 6.99936e+08 |
| BiCGSTABSolver | HypreBoomerAMG | &#10004; | 1 | 1.00329e-12 | 5.45661e+08 |
| BiCGSTABSolver | HypreILU | &#10004; | 1 | 6.33664e-13 | 2.65186e+08 |
| BiCGSTABSolver | HypreEuclid | &#10004; | 1 | 4.57805e-13 | 3.813e+10 |
| BiCGSTABSolver | HypreParaSails | &#10004; | 1 | 6.68232e-13 | 3.825e+08 |
| BiCGSTABSolver | HypreDiagScale | &#10004; | 1 | 5.62894e-12 | 2.79097e+08 |
| MINRESSolver | ANY | &#10004; | 1 | 4.48072e-12 | 2.02166e+08 |
| MINRESSolver | HypreBoomerAMG | &#10004; | 1 | 4.02522e-12 | 4.0452e+08 |
| MINRESSolver | HypreILU | &#10004; | 1 | 5.65164e-12 | 1.80559e+08 |
| MINRESSolver | HypreEuclid | &#10004; | 1 | 2.11107e-12 | 3.79546e+10 |
| MINRESSolver | HypreParaSails | &cross; | &cross; | &cross; | &cross; |
| MINRESSolver | HypreDiagScale | &#10004; | 1 | 7.77112e-12 | 1.46421e+08 |

### 256 processes

| solver | preconditionner | converged | iterations | residu | time |
| :-------------|:---------------|:---------------:|:---------------:|:---------------:| ---------------: |
| HyprePCG | ANY | &#10004; | 1 | 3.16367e-11 | 6.31277e+08 |
| HyprePCG | HypreBoomerAMG | &#10004; | 2 | 4.43327e-14 | 6.86705e+08 |
| HyprePCG | HypreILU | &#10004; | 2 | 4.43256e-14 | 2.27496e+08 |
| HyprePCG | HypreEuclid | &#10004; | 2 | 4.45541e-14 | 5.12471e+10 |
| HyprePCG | HypreDiagScale | &#10004; | 1 | 3.22487e-11 | 8.85343e+07 |
| HypreFGMRES | ANY | &cross; | &cross; | &cross; | &cross; |
| HypreFGMRES | HypreBoomerAMG | &#10004; | 1 | 3.02046e-11 | 7.39224e+08 |
| HypreFGMRES | HypreILU | &#10004; | 1 | 3.19017e-11 | 4.63415e+08 |
| HypreFGMRES | HypreParaSails | &#10004; | 1 | 3.20756e-11 | 9.31205e+08 |
| HypreFGMRES | HypreDiagScale | &#10004; | 4 | 1.08478e-12 | 7.496e+09 |
| HypreGMRES | ANY | &cross; | &cross; | &cross; | &cross; |
| HypreGMRES | HypreBoomerAMG | &#10004; | 1 | 3.02048e-11 | 7.8395e+08 |
| HypreGMRES | HypreILU | &#10004; | 1 | 3.19024e-11 | 4.65424e+08 |
| HypreGMRES | HypreParaSails | &#10004; | 1 | 3.20757e-11 | 5.42546e+08 |
| HypreGMRES | HypreDiagScale | &#10004; | 3 | 2.77391e-11 | 5.48894e+09 |
| MUMPSSolver | ANY | &#10004; | 1 | 1.79382e-13 | 1.20456e+10 |
| CGSolver | ANY | &#10004; | 1 | 1.16704e-12 | 1.39094e+08 |
| CGSolver | HypreParaSails | &cross; | &cross; | &cross; | &cross; |
| CGSolver | HypreDiagScale | &#10004; | 1 | 7.48684e-12 | 9.28112e+07 |
| BiCGSTABSolver | ANY | &#10004; | 1 | 1.24806e-12 | 4.68141e+08 |
| BiCGSTABSolver | HypreBoomerAMG | &#10004; | 1 | 8.17793e-13 | 5.1469e+08 |
| BiCGSTABSolver | HypreILU | &#10004; | 1 | 1.00223e-12 | 1.82264e+08 |
| BiCGSTABSolver | HypreEuclid | &#10004; | 1 | 3.41084e-13 | 2.57518e+10 |
| BiCGSTABSolver | HypreParaSails | &#10004; | 1 | 1.09135e-12 | 2.74539e+08 |
| BiCGSTABSolver | HypreDiagScale | &#10004; | 1 | 1.07804e-12 | 1.90592e+08 |
| MINRESSolver | ANY | &#10004; | 1 | 4.49087e-12 | 1.39058e+08 |
| MINRESSolver | HypreBoomerAMG | &#10004; | 1 | 4.43451e-12 | 3.6128e+08 |
| MINRESSolver | HypreILU | &#10004; | 1 | 5.45301e-12 | 1.16021e+08 |
| MINRESSolver | HypreEuclid | &#10004; | 1 | 2.29658e-12 | 2.55112e+10 |
| MINRESSolver | HypreParaSails | &cross; | &cross; | &cross; | &cross; |
| MINRESSolver | HypreDiagScale | &#10004; | 1 | 7.77763e-12 | 9.14023e+07 |


### 512 processes

| solver | preconditionner | converged | iterations | residu | time |
| :-------------|:---------------|:---------------:|:---------------:|:---------------:| ---------------: |
| HyprePCG | ANY | &#10004; | 1 | 3.15931e-11 | 1.46977e+09 | 
| HyprePCG | HypreBoomerAMG | &#10004; | 2 | 4.48997e-14 | 9.36578e+08 | 
| HyprePCG | HypreILU | &#10004; | 2 | 4.4992e-14 | 2.12171e+08 | 
| HyprePCG | HypreEuclid | &#10004; | 2 | 4.44e-14 | 2.62706e+10 | 
| HyprePCG | HypreDiagScale | &#10004; | 2 | 4.48628e-14 | 1.41588e+08 | 
| HypreFGMRES | ANY | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreBoomerAMG | &#10004; | 1 | 3.18201e-11 | 9.33016e+08 | 
| HypreFGMRES | HypreILU | &#10004; | 1 | 3.19821e-11 | 4.3614e+08 | 
| HypreFGMRES | HypreParaSails | &#10004; | 1 | 3.20748e-11 | 1.39837e+09 | 
| HypreFGMRES | HypreDiagScale | &#10004; | 4 | 7.73158e-12 | 1.05439e+10 | 
| HypreGMRES | ANY | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreBoomerAMG | &#10004; | 1 | 3.18199e-11 | 9.56455e+08 | 
| HypreGMRES | HypreILU | &#10004; | 1 | 3.19826e-11 | 4.64469e+08 | 
| HypreGMRES | HypreParaSails | &#10004; | 1 | 3.20748e-11 | 5.54849e+08 | 
| HypreGMRES | HypreDiagScale | &#10004; | 6 | 7.18327e-12 | 1.4865e+10 | 
| MUMPSSolver | ANY | &#10004; | 1 | 1.71603e-13 | 1.14057e+10 | 
| CGSolver | ANY | &#10004; | 1 | 1.18068e-12 | 1.08672e+08 | 
| CGSolver | HypreParaSails | &cross; | &cross; | &cross; | &cross; | 
| CGSolver | HypreDiagScale | &#10004; | 1 | 7.3736e-12 | 7.17475e+07 | 
| BiCGSTABSolver | ANY | &#10004; | 1 | 1.36678e-12 | 4.05911e+08 | 
| BiCGSTABSolver | HypreBoomerAMG | &#10004; | 1 | 9.99699e-13 | 6.23385e+08 | 
| BiCGSTABSolver | HypreILU | &#10004; | 1 | 8.11712e-13 | 1.23463e+08 | 
| BiCGSTABSolver | HypreEuclid | &#10004; | 1 | 9.67875e-13 | 1.32186e+10 | 
| BiCGSTABSolver | HypreParaSails | &#10004; | 1 | 5.34142e-13 | 2.06533e+08 | 
| BiCGSTABSolver | HypreDiagScale | &#10004; | 1 | 9.2699e-13 | 1.36466e+08 | 
| MINRESSolver | ANY | &#10004; | 1 | 4.52263e-12 | 1.10697e+08 | 
| MINRESSolver | HypreBoomerAMG | &#10004; | 1 | 4.43764e-12 | 4.4563e+08 | 
| MINRESSolver | HypreILU | &#10004; | 1 | 5.61177e-12 | 8.29284e+07 | 
| MINRESSolver | HypreEuclid | &#10004; | 1 | 3.26885e-12 | 1.30963e+10 | 
| MINRESSolver | HypreParaSails | &cross; | &cross; | &cross; | &cross; | 
| MINRESSolver | HypreDiagScale | &#10004; | 1 | 7.76551e-12 | 6.73967e+07 | 
