# Feedback on Topaze

## Simulation 1 : many_spheres

Refinement level : 1 


| solver | preconditionner | number of mpi | converged | iterations | residu | time |
| :-------------|:---------------|:---------------:|:---------------:|:---------------:|:---------------:|---------------:|
| HypreFGMRES | HypreDiagScale | 1 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 2 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 32 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 64 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 128 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 256 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 512 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 1024 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 2048 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 4096 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 8192 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 16384 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreBoomerAMG | 1 | &#10004; | 1 | 1.62331e-11 | 5.07337e+11 | 
| HypreFGMRES | HypreBoomerAMG | 2 | &#10004; | 1 | 1.54386e-11 | 4.06457e+11 | 
| HypreFGMRES | HypreBoomerAMG | 32 | &#10004; | 2 | 4.13734e-14 | 5.49606e+11 | 
| HypreFGMRES | HypreBoomerAMG | 64 | &#10004; | 1 | 1.62221e-11 | 6.3186e+10 | 
| HypreFGMRES | HypreBoomerAMG | 128 | &#10004; | 1 | 1.62396e-11 | 2.57268e+10 | 
| HypreFGMRES | HypreBoomerAMG | 256 | &#10004; | 1 | 1.61277e-11 | 1.54386e+10 | 
| HypreFGMRES | HypreBoomerAMG | 512 | &#10004; | 1 | 1.61658e-11 | 7.0965e+09 | 
| HypreFGMRES | HypreBoomerAMG | 1024 | &#10004; | 1 | 1.54603e-11 | 8.75541e+09 | 
| HypreFGMRES | HypreBoomerAMG | 2048 | &#10004; | 1 | 1.62507e-11 | 1.06409e+10 | 
| HypreFGMRES | HypreBoomerAMG | 4096 | &#10004; | 1 | 1.62481e-11 | 1.84465e+10 | 
| HypreFGMRES | HypreBoomerAMG | 8192 | &#10004; | 2 | 4.16581e-14 | 4.06142e+10 | 
| HypreFGMRES | HypreBoomerAMG | 16384 | &#10004; | 2 | 5.43858e-13 | 4.59165e+10 | 
| HypreFGMRES | HypreILU | 1 | &#10004; | 1 | 1.57027e-11 | 2.75508e+11 | 
| HypreFGMRES | HypreILU | 2 | &#10004; | 1 | 1.57911e-11 | 2.84016e+11 | 
| HypreFGMRES | HypreILU | 32 | &#10004; | 2 | 4.13694e-14 | 3.58887e+11 | 
| HypreFGMRES | HypreILU | 64 | &#10004; | 2 | 8.93405e-13 | 1.87877e+11 | 
| HypreFGMRES | HypreILU | 128 | &#10004; | 1 | 1.59811e-11 | 4.2703e+10 | 
| HypreFGMRES | HypreILU | 256 | &#10004; | 1 | 1.55932e-11 | 1.14633e+10 | 
| HypreFGMRES | HypreILU | 512 | &#10004; | 2 | 4.17141e-14 | 1.50265e+10 | 
| HypreFGMRES | HypreILU | 1024 | &#10004; | 3 | 5.7957e-12 | 1.68581e+10 | 
| HypreFGMRES | HypreILU | 2048 | &#10004; | 5 | 1.18959e-12 | 2.85679e+10 | 
| HypreFGMRES | HypreILU | 4096 | &#10004; | 7 | 4.51387e-12 | 4.663e+10 | 
| HypreFGMRES | HypreILU | 8192 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreILU | 16384 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | ANY | 1 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | ANY | 2 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | ANY | 32 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | ANY | 64 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | ANY | 128 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | ANY | 256 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | ANY | 512 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | ANY | 1024 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | ANY | 2048 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | ANY | 4096 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | ANY | 8192 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | ANY | 16384 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 1 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 2 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 32 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 64 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 128 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 256 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 512 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 1024 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 2048 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 4096 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 8192 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 16384 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreBoomerAMG | 1 | &#10004; | 1 | 1.62336e-11 | 5.10664e+11 | 
| HypreGMRES | HypreBoomerAMG | 2 | &#10004; | 1 | 1.5439e-11 | 4.16752e+11 | 
| HypreGMRES | HypreBoomerAMG | 32 | &#10004; | 2 | 5.69757e-14 | 5.5562e+11 | 
| HypreGMRES | HypreBoomerAMG | 64 | &#10004; | 1 | 1.62226e-11 | 6.35948e+10 | 
| HypreGMRES | HypreBoomerAMG | 128 | &#10004; | 1 | 1.62399e-11 | 2.57708e+10 | 
| HypreGMRES | HypreBoomerAMG | 256 | &#10004; | 1 | 1.61281e-11 | 1.50037e+10 | 
| HypreGMRES | HypreBoomerAMG | 512 | &#10004; | 1 | 1.61662e-11 | 6.70454e+09 | 
| HypreGMRES | HypreBoomerAMG | 1024 | &#10004; | 1 | 1.54606e-11 | 8.48375e+09 | 
| HypreGMRES | HypreBoomerAMG | 2048 | &#10004; | 1 | 1.62509e-11 | 1.04408e+10 | 
| HypreGMRES | HypreBoomerAMG | 4096 | &#10004; | 1 | 1.62484e-11 | 1.73871e+10 | 
| HypreGMRES | HypreBoomerAMG | 8192 | &#10004; | 2 | 4.16572e-14 | 4.05755e+10 | 
| HypreGMRES | HypreBoomerAMG | 16384 | &#10004; | 2 | 4.15217e-14 | 4.67668e+10 | 
| HypreGMRES | HypreILU | 1 | &#10004; | 1 | 1.57034e-11 | 2.7907e+11 | 
| HypreGMRES | HypreILU | 2 | &#10004; | 1 | 1.5792e-11 | 2.89124e+11 | 
| HypreGMRES | HypreILU | 32 | &#10004; | 2 | 4.14632e-14 | 3.78568e+11 | 
| HypreGMRES | HypreILU | 64 | &#10004; | 2 | 1.56779e-13 | 1.86129e+11 | 
| HypreGMRES | HypreILU | 128 | &#10004; | 1 | 1.59816e-11 | 4.30084e+10 | 
| HypreGMRES | HypreILU | 256 | &#10004; | 1 | 1.55937e-11 | 1.14949e+10 | 
| HypreGMRES | HypreILU | 512 | &#10004; | 2 | 4.16661e-14 | 1.44532e+10 | 
| HypreGMRES | HypreILU | 1024 | &#10004; | 3 | 7.86323e-12 | 1.66078e+10 | 
| HypreGMRES | HypreILU | 2048 | &#10004; | 4 | 2.48858e-13 | 2.23566e+10 | 
| HypreGMRES | HypreILU | 4096 | &#10004; | 5 | 6.21654e-12 | 3.25909e+10 | 
| HypreGMRES | HypreILU | 8192 | &#10004; | 10 | 6.99954e-12 | 7.16895e+10 | 
| HypreGMRES | HypreILU | 16384 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | ANY | 1 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | ANY | 2 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | ANY | 32 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | ANY | 64 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | ANY | 128 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | ANY | 256 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | ANY | 512 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | ANY | 1024 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | ANY | 2048 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | ANY | 4096 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | ANY | 8192 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | ANY | 16384 | &cross; | &cross; | &cross; | &cross; | 
| HyprePCG | HypreDiagScale | 1 | &#10004; | 2 | 4.10891e-14 | 2.72545e+11 | 
| HyprePCG | HypreDiagScale | 2 | &#10004; | 2 | 4.10426e-14 | 1.59387e+11 | 
| HyprePCG | HypreDiagScale | 32 | &#10004; | 2 | 4.13774e-14 | 2.49085e+10 | 
| HyprePCG | HypreDiagScale | 64 | &#10004; | 2 | 4.20666e-14 | 1.34006e+10 | 
| HyprePCG | HypreDiagScale | 128 | &#10004; | 2 | 4.16713e-14 | 5.04386e+09 | 
| HyprePCG | HypreDiagScale | 256 | &#10004; | 2 | 4.15446e-14 | 1.65512e+09 | 
| HyprePCG | HypreDiagScale | 512 | &#10004; | 2 | 4.17033e-14 | 9.15308e+08 | 
| HyprePCG | HypreDiagScale | 1024 | &#10004; | 2 | 4.23039e-14 | 6.53304e+08 | 
| HyprePCG | HypreDiagScale | 2048 | &#10004; | 2 | 4.1713e-14 | 4.85061e+08 | 
| HyprePCG | HypreDiagScale | 4096 | &#10004; | 2 | 4.24909e-14 | 4.96495e+08 | 
| HyprePCG | HypreDiagScale | 8192 | &#10004; | 2 | 4.16776e-14 | 5.02405e+08 | 
| HyprePCG | HypreDiagScale | 16384 | &#10004; | 2 | 4.14963e-14 | 6.14634e+08 | 
| HyprePCG | HypreBoomerAMG | 1 | &#10004; | 2 | 4.11339e-14 | 3.72123e+11 | 
| HyprePCG | HypreBoomerAMG | 2 | &#10004; | 2 | 4.09397e-14 | 2.14618e+11 | 
| HyprePCG | HypreBoomerAMG | 32 | &#10004; | 2 | 4.1299e-14 | 3.53857e+10 | 
| HyprePCG | HypreBoomerAMG | 64 | &#10004; | 2 | 4.20937e-14 | 1.87617e+10 | 
| HyprePCG | HypreBoomerAMG | 128 | &#10004; | 2 | 4.16361e-14 | 7.99635e+09 | 
| HyprePCG | HypreBoomerAMG | 256 | &#10004; | 2 | 4.15697e-14 | 3.20422e+09 | 
| HyprePCG | HypreBoomerAMG | 512 | &#10004; | 2 | 4.15905e-14 | 2.65435e+09 | 
| HyprePCG | HypreBoomerAMG | 1024 | &#10004; | 2 | 4.24515e-14 | 2.43651e+09 | 
| HyprePCG | HypreBoomerAMG | 2048 | &#10004; | 2 | 4.16921e-14 | 2.40574e+09 | 
| HyprePCG | HypreBoomerAMG | 4096 | &#10004; | 2 | 4.23583e-14 | 2.50261e+09 | 
| HyprePCG | HypreBoomerAMG | 8192 | &#10004; | 2 | 4.1831e-14 | 3.22667e+09 | 
| HyprePCG | HypreBoomerAMG | 16384 | &#10004; | 2 | 4.17649e-14 | 2.83406e+09 | 
| HyprePCG | HypreILU | 1 | &#10004; | 2 | 4.11507e-14 | 2.39584e+11 | 
| HyprePCG | HypreILU | 2 | &#10004; | 2 | 4.1051e-14 | 1.63333e+11 | 
| HyprePCG | HypreILU | 32 | &#10004; | 2 | 4.13539e-14 | 2.65402e+10 | 
| HyprePCG | HypreILU | 64 | &#10004; | 2 | 4.21919e-14 | 1.58807e+10 | 
| HyprePCG | HypreILU | 128 | &#10004; | 2 | 4.16929e-14 | 7.6931e+09 | 
| HyprePCG | HypreILU | 256 | &#10004; | 2 | 4.16206e-14 | 3.19577e+09 | 
| HyprePCG | HypreILU | 512 | &#10004; | 2 | 4.16579e-14 | 1.03516e+09 | 
| HyprePCG | HypreILU | 1024 | &#10004; | 2 | 4.23183e-14 | 6.40727e+08 | 
| HyprePCG | HypreILU | 2048 | &#10004; | 2 | 4.17542e-14 | 4.92351e+08 | 
| HyprePCG | HypreILU | 4096 | &#10004; | 2 | 4.23467e-14 | 4.07229e+08 | 
| HyprePCG | HypreILU | 8192 | &#10004; | 2 | 4.17402e-14 | 4.31898e+08 | 
| HyprePCG | HypreILU | 16384 | &#10004; | 2 | 4.15727e-14 | 5.7153e+08 | 
| HyprePCG | ANY | 1 | &#10004; | 1 | 1.56261e-11 | 2.48668e+11 | 
| HyprePCG | ANY | 2 | &#10004; | 1 | 1.41618e-11 | 1.36043e+11 | 
| HyprePCG | ANY | 32 | &#10004; | 1 | 1.54485e-11 | 2.44395e+10 | 
| HyprePCG | ANY | 64 | &#10004; | 1 | 1.46064e-11 | 1.31294e+10 | 
| HyprePCG | ANY | 128 | &#10004; | 2 | 4.16581e-14 | 9.94432e+09 | 
| HyprePCG | ANY | 256 | &#10004; | 1 | 1.43737e-11 | 1.55256e+09 | 
| HyprePCG | ANY | 512 | &#10004; | 1 | 1.46915e-11 | 8.9548e+08 | 
| HyprePCG | ANY | 1024 | &#10004; | 1 | 1.60938e-11 | 6.1382e+08 | 
| HyprePCG | ANY | 2048 | &#10004; | 1 | 1.46525e-11 | 5.01002e+08 | 
| HyprePCG | ANY | 4096 | &#10004; | 1 | 1.49761e-11 | 5.039e+08 | 
| HyprePCG | ANY | 8192 | &#10004; | 2 | 4.17086e-14 | 1.00263e+09 | 
| HyprePCG | ANY | 16384 | &#10004; | 1 | 1.54911e-11 | 7.51577e+08 | 
| CGSolver | HypreDiagScale | 1 | &#10004; | 1 | 7.15066e-12 | 1.41992e+11 | 
| CGSolver | HypreDiagScale | 2 | &#10004; | 1 | 5.85546e-12 | 8.17612e+10 | 
| CGSolver | HypreDiagScale | 32 | &#10004; | 1 | 5.43121e-12 | 1.29406e+10 | 
| CGSolver | HypreDiagScale | 64 | &#10004; | 1 | 5.48768e-12 | 6.74793e+09 | 
| CGSolver | HypreDiagScale | 128 | &#10004; | 1 | 7.42059e-12 | 2.57615e+09 | 
| CGSolver | HypreDiagScale | 256 | &#10004; | 1 | 6.06739e-12 | 8.56503e+08 | 
| CGSolver | HypreDiagScale | 512 | &#10004; | 1 | 5.63854e-12 | 4.67084e+08 | 
| CGSolver | HypreDiagScale | 1024 | &#10004; | 1 | 5.87078e-12 | 3.35604e+08 | 
| CGSolver | HypreDiagScale | 2048 | &#10004; | 1 | 6.34586e-12 | 2.83345e+08 | 
| CGSolver | HypreDiagScale | 4096 | &#10004; | 1 | 6.68831e-12 | 2.50161e+08 | 
| CGSolver | HypreDiagScale | 8192 | &#10004; | 1 | 5.44574e-12 | 2.50983e+08 | 
| CGSolver | HypreDiagScale | 16384 | &#10004; | 1 | 5.70419e-12 | 4.09048e+08 | 
| CGSolver | ANY | 1 | &#10004; | 1 | 1.88249e-12 | 2.66211e+11 | 
| CGSolver | ANY | 2 | &#10004; | 1 | 1.84882e-12 | 1.39247e+11 | 
| CGSolver | ANY | 32 | &#10004; | 1 | 1.87214e-12 | 2.3948e+10 | 
| CGSolver | ANY | 64 | &#10004; | 1 | 1.86216e-12 | 1.36878e+10 | 
| CGSolver | ANY | 128 | &#10004; | 1 | 1.84577e-12 | 5.06449e+09 | 
| CGSolver | ANY | 256 | &#10004; | 1 | 1.84771e-12 | 1.59398e+09 | 
| CGSolver | ANY | 512 | &#10004; | 1 | 1.8225e-12 | 9.18757e+08 | 
| CGSolver | ANY | 1024 | &#10004; | 1 | 1.8938e-12 | 6.66258e+08 | 
| CGSolver | ANY | 2048 | &#10004; | 1 | 1.83262e-12 | 5.75772e+08 | 
| CGSolver | ANY | 4096 | &#10004; | 1 | 1.88339e-12 | 5.62388e+08 | 
| CGSolver | ANY | 8192 | &#10004; | 1 | 1.86897e-12 | 5.431e+08 | 
| CGSolver | ANY | 16384 | &#10004; | 1 | 1.84521e-12 | 1.19327e+09 | 
| BiCGSTABSolver | HypreDiagScale | 1 | &#10004; | 1 | 2.88499e-12 | 5.68504e+11 | 
| BiCGSTABSolver | HypreDiagScale | 2 | &#10004; | 1 | 1.27984e-11 | 2.39475e+11 | 
| BiCGSTABSolver | HypreDiagScale | 32 | &#10004; | 1 | 1.99794e-12 | 3.46609e+10 | 
| BiCGSTABSolver | HypreDiagScale | 64 | &#10004; | 2 | 9.97519e-13 | 2.22414e+10 | 
| BiCGSTABSolver | HypreDiagScale | 128 | &#10004; | 1 | 1.92607e-12 | 6.32674e+09 | 
| BiCGSTABSolver | HypreDiagScale | 256 | &#10004; | 1 | 5.41864e-12 | 2.11334e+09 | 
| BiCGSTABSolver | HypreDiagScale | 512 | &#10004; | 1 | 1.8446e-12 | 1.24263e+09 | 
| BiCGSTABSolver | HypreDiagScale | 1024 | &#10004; | 2 | 9.879e-13 | 8.83811e+08 | 
| BiCGSTABSolver | HypreDiagScale | 2048 | &#10004; | 1 | 2.17724e-12 | 1.00523e+09 | 
| BiCGSTABSolver | HypreDiagScale | 4096 | &#10004; | 1 | 1.1898e-11 | 5.67741e+08 | 
| BiCGSTABSolver | HypreDiagScale | 8192 | &#10004; | 1 | 2.26667e-12 | 8.75524e+08 | 
| BiCGSTABSolver | HypreDiagScale | 16384 | &#10004; | 1 | 5.92715e-12 | 2.23452e+09 | 
| BiCGSTABSolver | HypreEuclid | 1 | &#10004; | 1 | 1.15013e-12 | 3.71974e+11 | 
| BiCGSTABSolver | HypreEuclid | 2 | &#10004; | 1 | 1.15075e-12 | 4.58015e+11 | 
| BiCGSTABSolver | HypreEuclid | 32 | &#10004; | 1 | 8.77026e-13 | 7.28717e+11 | 
| BiCGSTABSolver | HypreEuclid | 64 | &#10004; | 1 | 8.88153e-13 | 3.10133e+11 | 
| BiCGSTABSolver | HypreEuclid | 128 | &#10004; | 1 | 6.07182e-13 | 1.89732e+11 | 
| BiCGSTABSolver | HypreEuclid | 256 | &#10004; | 1 | 1.00664e-12 | 1.10819e+11 | 
| BiCGSTABSolver | HypreEuclid | 512 | &#10004; | 1 | 1.13245e-12 | 6.54819e+10 | 
| BiCGSTABSolver | HypreEuclid | 1024 | &#10004; | 1 | 1.68911e-12 | 4.04569e+10 | 
| BiCGSTABSolver | HypreEuclid | 2048 | &#10004; | 1 | 8.30399e-13 | 2.42104e+10 | 
| BiCGSTABSolver | HypreEuclid | 4096 | &#10004; | 1 | 9.57949e-13 | 1.38746e+10 | 
| BiCGSTABSolver | HypreEuclid | 8192 | &#10004; | 1 | 1.02619e-12 | 9.0551e+09 | 
| BiCGSTABSolver | HypreEuclid | 16384 | &#10004; | 1 | 6.81028e-13 | 8.71127e+09 | 
| BiCGSTABSolver | HypreBoomerAMG | 1 | &#10004; | 1 | 9.86676e-13 | 2.59346e+11 | 
| BiCGSTABSolver | HypreBoomerAMG | 2 | &#10004; | 1 | 1.01599e-12 | 1.51063e+11 | 
| BiCGSTABSolver | HypreBoomerAMG | 32 | &#10004; | 1 | 8.77936e-13 | 2.57424e+10 | 
| BiCGSTABSolver | HypreBoomerAMG | 64 | &#10004; | 1 | 9.38865e-13 | 1.32965e+10 | 
| BiCGSTABSolver | HypreBoomerAMG | 128 | &#10004; | 1 | 7.72631e-13 | 6.20509e+09 | 
| BiCGSTABSolver | HypreBoomerAMG | 256 | &#10004; | 1 | 1.25515e-12 | 2.14726e+09 | 
| BiCGSTABSolver | HypreBoomerAMG | 512 | &#10004; | 1 | 8.76885e-13 | 1.89128e+09 | 
| BiCGSTABSolver | HypreBoomerAMG | 1024 | &#10004; | 1 | 8.05915e-13 | 1.60674e+09 | 
| BiCGSTABSolver | HypreBoomerAMG | 2048 | &#10004; | 1 | 1.1055e-12 | 1.6647e+09 | 
| BiCGSTABSolver | HypreBoomerAMG | 4096 | &#10004; | 1 | 2.62762e-12 | 1.79875e+09 | 
| BiCGSTABSolver | HypreBoomerAMG | 8192 | &#10004; | 1 | 1.00322e-12 | 1.8871e+09 | 
| BiCGSTABSolver | HypreBoomerAMG | 16384 | &#10004; | 1 | 1.164e-12 | 2.28082e+09 | 
| BiCGSTABSolver | HypreILU | 1 | &#10004; | 1 | 1.08101e-12 | 1.89735e+11 | 
| BiCGSTABSolver | HypreILU | 2 | &#10004; | 1 | 1.1117e-12 | 1.18495e+11 | 
| BiCGSTABSolver | HypreILU | 32 | &#10004; | 1 | 7.19619e-12 | 2.25176e+10 | 
| BiCGSTABSolver | HypreILU | 64 | &#10004; | 1 | 7.79477e-13 | 1.23988e+10 | 
| BiCGSTABSolver | HypreILU | 128 | &#10004; | 1 | 7.96874e-13 | 6.55254e+09 | 
| BiCGSTABSolver | HypreILU | 256 | &#10004; | 1 | 1.2684e-12 | 2.95689e+09 | 
| BiCGSTABSolver | HypreILU | 512 | &#10004; | 1 | 1.06823e-12 | 8.97774e+08 | 
| BiCGSTABSolver | HypreILU | 1024 | &#10004; | 1 | 9.31115e-13 | 6.05788e+08 | 
| BiCGSTABSolver | HypreILU | 2048 | &#10004; | 1 | 1.17833e-12 | 5.03326e+08 | 
| BiCGSTABSolver | HypreILU | 4096 | &#10004; | 1 | 1.45589e-12 | 5.03303e+08 | 
| BiCGSTABSolver | HypreILU | 8192 | &#10004; | 1 | 1.24544e-12 | 5.36014e+08 | 
| BiCGSTABSolver | HypreILU | 16384 | &#10004; | 1 | 2.13015e-12 | 9.18826e+08 | 
| BiCGSTABSolver | ANY | 1 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | ANY | 2 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | ANY | 32 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | ANY | 64 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | ANY | 128 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | ANY | 256 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | ANY | 512 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | ANY | 1024 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | ANY | 2048 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | ANY | 4096 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | ANY | 8192 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | ANY | 16384 | &cross; | &cross; | &cross; | &cross; | 
| MINRESSolver | HypreDiagScale | 1 | &#10004; | 2 | 5.42146e-12 | 1.62611e+11 | 
| MINRESSolver | HypreDiagScale | 2 | &#10004; | 1 | 1.13595e-11 | 8.39732e+10 | 
| MINRESSolver | HypreDiagScale | 32 | &#10004; | 2 | 5.34845e-12 | 1.50716e+10 | 
| MINRESSolver | HypreDiagScale | 64 | &#10004; | 2 | 5.40212e-12 | 7.6117e+09 | 
| MINRESSolver | HypreDiagScale | 128 | &#10004; | 1 | 1.34318e-11 | 2.75783e+09 | 
| MINRESSolver | HypreDiagScale | 256 | &#10004; | 2 | 5.387e-12 | 9.53422e+08 | 
| MINRESSolver | HypreDiagScale | 512 | &#10004; | 2 | 5.39938e-12 | 5.1694e+08 | 
| MINRESSolver | HypreDiagScale | 1024 | &#10004; | 2 | 5.38256e-12 | 3.522e+08 | 
| MINRESSolver | HypreDiagScale | 2048 | &#10004; | 2 | 5.42127e-12 | 2.93124e+08 | 
| MINRESSolver | HypreDiagScale | 4096 | &#10004; | 2 | 5.34643e-12 | 2.8144e+08 | 
| MINRESSolver | HypreDiagScale | 8192 | &#10004; | 2 | 5.39062e-12 | 2.29868e+08 | 
| MINRESSolver | HypreDiagScale | 16384 | &#10004; | 2 | 5.32297e-12 | 2.92963e+08 | 
| MINRESSolver | HypreEuclid | 1 | &#10004; | 1 | 2.84563e-12 | 2.38906e+11 | 
| MINRESSolver | HypreEuclid | 2 | &#10004; | 1 | 2.08134e-12 | 4.03854e+11 | 
| MINRESSolver | HypreEuclid | 32 | &#10004; | 1 | 2.59859e-12 | 7.0531e+11 | 
| MINRESSolver | HypreEuclid | 64 | &#10004; | 1 | 2.34448e-12 | 3.01307e+11 | 
| MINRESSolver | HypreEuclid | 128 | &#10004; | 1 | 1.85058e-12 | 1.85005e+11 | 
| MINRESSolver | HypreEuclid | 256 | &#10004; | 1 | 2.36094e-12 | 1.09198e+11 | 
| MINRESSolver | HypreEuclid | 512 | &#10004; | 1 | 2.5897e-12 | 6.4782e+10 | 
| MINRESSolver | HypreEuclid | 1024 | &#10004; | 1 | 2.17107e-12 | 4.02919e+10 | 
| MINRESSolver | HypreEuclid | 2048 | &#10004; | 1 | 2.53513e-12 | 2.3883e+10 | 
| MINRESSolver | HypreEuclid | 4096 | &#10004; | 1 | 2.39092e-12 | 1.33964e+10 | 
| MINRESSolver | HypreEuclid | 8192 | &#10004; | 1 | 2.57328e-12 | 8.78146e+09 | 
| MINRESSolver | HypreEuclid | 16384 | &#10004; | 1 | 2.67081e-12 | 8.42404e+09 | 
| MINRESSolver | HypreBoomerAMG | 1 | &#10004; | 1 | 3.01518e-12 | 1.86742e+11 | 
| MINRESSolver | HypreBoomerAMG | 2 | &#10004; | 1 | 2.9744e-12 | 1.07441e+11 | 
| MINRESSolver | HypreBoomerAMG | 32 | &#10004; | 1 | 3.04616e-12 | 1.93084e+10 | 
| MINRESSolver | HypreBoomerAMG | 64 | &#10004; | 1 | 3.05328e-12 | 9.91053e+09 | 
| MINRESSolver | HypreBoomerAMG | 128 | &#10004; | 1 | 3.04966e-12 | 4.08936e+09 | 
| MINRESSolver | HypreBoomerAMG | 256 | &#10004; | 1 | 3.22016e-12 | 1.62496e+09 | 
| MINRESSolver | HypreBoomerAMG | 512 | &#10004; | 1 | 3.18081e-12 | 1.31444e+09 | 
| MINRESSolver | HypreBoomerAMG | 1024 | &#10004; | 1 | 3.44884e-12 | 1.2496e+09 | 
| MINRESSolver | HypreBoomerAMG | 2048 | &#10004; | 1 | 3.3948e-12 | 1.20737e+09 | 
| MINRESSolver | HypreBoomerAMG | 4096 | &#10004; | 1 | 3.46967e-12 | 1.37414e+09 | 
| MINRESSolver | HypreBoomerAMG | 8192 | &#10004; | 1 | 3.70236e-12 | 1.29334e+09 | 
| MINRESSolver | HypreBoomerAMG | 16384 | &#10004; | 1 | 3.67913e-12 | 1.29745e+09 | 
| MINRESSolver | HypreILU | 1 | &#10004; | 1 | 1.65034e-12 | 1.24847e+11 | 
| MINRESSolver | HypreILU | 2 | &#10004; | 1 | 3.03756e-12 | 8.39336e+10 | 
| MINRESSolver | HypreILU | 32 | &#10004; | 1 | 3.83505e-12 | 1.45459e+10 | 
| MINRESSolver | HypreILU | 64 | &#10004; | 1 | 4.43504e-12 | 8.21975e+09 | 
| MINRESSolver | HypreILU | 128 | &#10004; | 1 | 4.15845e-12 | 3.98065e+09 | 
| MINRESSolver | HypreILU | 256 | &#10004; | 1 | 4.52914e-12 | 1.66037e+09 | 
| MINRESSolver | HypreILU | 512 | &#10004; | 1 | 5.2658e-12 | 5.64828e+08 | 
| MINRESSolver | HypreILU | 1024 | &#10004; | 1 | 5.99956e-12 | 3.39024e+08 | 
| MINRESSolver | HypreILU | 2048 | &#10004; | 1 | 6.42509e-12 | 2.34299e+08 | 
| MINRESSolver | HypreILU | 4096 | &#10004; | 1 | 6.57518e-12 | 2.12006e+08 | 
| MINRESSolver | HypreILU | 8192 | &#10004; | 1 | 8.18843e-12 | 1.81407e+08 | 
| MINRESSolver | HypreILU | 16384 | &#10004; | 1 | 9.2075e-12 | 2.58607e+08 | 
| MINRESSolver | ANY | 1 | &#10004; | 2 | 9.79573e-13 | 2.98232e+11 | 
| MINRESSolver | ANY | 2 | &#10004; | 2 | 9.96936e-13 | 1.69774e+11 | 
| MINRESSolver | ANY | 32 | &#10004; | 2 | 9.76179e-13 | 2.7601e+10 | 
| MINRESSolver | ANY | 64 | &#10004; | 2 | 9.94973e-13 | 1.47696e+10 | 
| MINRESSolver | ANY | 128 | &#10004; | 2 | 9.89211e-13 | 5.58842e+09 | 
| MINRESSolver | ANY | 256 | &#10004; | 2 | 9.79986e-13 | 1.72975e+09 | 
| MINRESSolver | ANY | 512 | &#10004; | 2 | 9.75355e-13 | 9.92546e+08 | 
| MINRESSolver | ANY | 1024 | &#10004; | 2 | 9.90457e-13 | 6.91275e+08 | 
| MINRESSolver | ANY | 2048 | &#10004; | 2 | 9.76513e-13 | 5.67617e+08 | 
| MINRESSolver | ANY | 4096 | &#10004; | 2 | 9.84411e-13 | 5.54827e+08 | 
| MINRESSolver | ANY | 8192 | &#10004; | 2 | 9.80479e-13 | 6.2937e+08 | 
| MINRESSolver | ANY | 16384 | &#10004; | 2 | 9.9762e-13 | 6.46841e+08 | 



Refinement level 2

| solver | preconditionner | number of mpi | converged | iterations | residu | time |
| :-------------|:---------------|:---------------:|:---------------:|:---------------:|:---------------:|---------------:|
| HypreFGMRES | HypreDiagScale | 128 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 512 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 1024 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 2048 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 4096 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 8192 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreDiagScale | 16384 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreBoomerAMG | 128 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreBoomerAMG | 512 | &#10004; | 8 | 3.50589e-12 | 1.22156e+12 | 
| HypreFGMRES | HypreBoomerAMG | 1024 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreBoomerAMG | 2048 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreBoomerAMG | 4096 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreBoomerAMG | 8192 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreBoomerAMG | 16384 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreILU | 128 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreILU | 512 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreILU | 1024 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreILU | 2048 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreILU | 4096 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreILU | 8192 | &cross; | &cross; | &cross; | &cross; | 
| HypreFGMRES | HypreILU | 16384 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 128 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 512 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 1024 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 2048 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 4096 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 8192 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreDiagScale | 16384 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreBoomerAMG | 128 | &#10004; | 10 | 2.17626e-12 | 6.59336e+12 | 
| HypreGMRES | HypreBoomerAMG | 512 | &#10004; | 10 | 7.21225e-13 | 1.52134e+12 | 
| HypreGMRES | HypreBoomerAMG | 1024 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreBoomerAMG | 2048 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreBoomerAMG | 4096 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreBoomerAMG | 8192 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreBoomerAMG | 16384 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreILU | 128 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreILU | 512 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreILU | 1024 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreILU | 2048 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreILU | 4096 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreILU | 8192 | &cross; | &cross; | &cross; | &cross; | 
| HypreGMRES | HypreILU | 16384 | &cross; | &cross; | &cross; | &cross; | 
| HyprePCG | HypreDiagScale | 128 | &#10004; | 2 | 4.66279e-14 | 1.67521e+11 | 
| HyprePCG | HypreDiagScale | 512 | &#10004; | 2 | 4.65242e-14 | 3.14688e+10 | 
| HyprePCG | HypreDiagScale | 1024 | &#10004; | 2 | 4.77407e-14 | 1.12259e+10 | 
| HyprePCG | HypreDiagScale | 2048 | &#10004; | 2 | 4.64894e-14 | 3.81999e+09 | 
| HyprePCG | HypreDiagScale | 4096 | &#10004; | 2 | 4.75072e-14 | 2.40976e+09 | 
| HyprePCG | HypreDiagScale | 8192 | &#10004; | 2 | 4.65072e-14 | 1.99514e+09 | 
| HyprePCG | HypreDiagScale | 16384 | &#10004; | 2 | 4.64841e-14 | 3.18626e+09 | 
| HyprePCG | HypreBoomerAMG | 128 | &#10004; | 2 | 4.66762e-14 | 1.74074e+11 | 
| HyprePCG | HypreBoomerAMG | 512 | &#10004; | 2 | 4.6517e-14 | 4.17199e+10 | 
| HyprePCG | HypreBoomerAMG | 1024 | &#10004; | 2 | 4.77139e-14 | 1.71938e+10 | 
| HyprePCG | HypreBoomerAMG | 2048 | &#10004; | 2 | 4.64793e-14 | 8.4098e+09 | 
| HyprePCG | HypreBoomerAMG | 4096 | &#10004; | 2 | 4.74872e-14 | 7.191e+09 | 
| HyprePCG | HypreBoomerAMG | 8192 | &#10004; | 2 | 4.65201e-14 | 6.74059e+09 | 
| HyprePCG | HypreBoomerAMG | 16384 | &#10004; | 2 | 4.64184e-14 | 9.06286e+09 | 
| HyprePCG | HypreILU | 128 | &#10004; | 2 | 4.67021e-14 | 1.39792e+11 | 
| HyprePCG | HypreILU | 512 | &#10004; | 2 | 4.64731e-14 | 3.318e+10 | 
| HyprePCG | HypreILU | 1024 | &#10004; | 2 | 4.7618e-14 | 1.58491e+10 | 
| HyprePCG | HypreILU | 2048 | &#10004; | 2 | 4.65057e-14 | 6.26883e+09 | 
| HyprePCG | HypreILU | 4096 | &#10004; | 2 | 4.75119e-14 | 2.42314e+09 | 
| HyprePCG | HypreILU | 8192 | &#10004; | 2 | 4.64937e-14 | 2.1162e+09 | 
| HyprePCG | HypreILU | 16384 | &#10004; | 2 | 4.65279e-14 | 2.41204e+09 | 
| CGSolver | HypreDiagScale | 128 | &#10004; | 1 | 4.82871e-12 | 8.33267e+10 | 
| CGSolver | HypreDiagScale | 512 | &#10004; | 1 | 4.77591e-12 | 1.58043e+10 | 
| CGSolver | HypreDiagScale | 1024 | &#10004; | 1 | 5.80629e-12 | 5.77684e+09 | 
| CGSolver | HypreDiagScale | 2048 | &#10004; | 1 | 5.5351e-12 | 1.90185e+09 | 
| CGSolver | HypreDiagScale | 4096 | &#10004; | 1 | 4.5837e-12 | 1.17885e+09 | 
| CGSolver | HypreDiagScale | 8192 | &#10004; | 1 | 4.41958e-12 | 9.60211e+08 | 
| CGSolver | HypreDiagScale | 16384 | &#10004; | 1 | 4.78575e-12 | 2.16078e+09 | 
| BiCGSTABSolver | HypreDiagScale | 128 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | HypreDiagScale | 512 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | HypreDiagScale | 1024 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | HypreDiagScale | 2048 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | HypreDiagScale | 4096 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | HypreDiagScale | 8192 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | HypreDiagScale | 16384 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | HypreEuclid | 128 | &#10004; | 2 | 5.05833e-13 | 1.83076e+12 | 
| BiCGSTABSolver | HypreEuclid | 512 | &#10004; | 1 | 1.35085e-12 | 3.33505e+11 | 
| BiCGSTABSolver | HypreEuclid | 1024 | &#10004; | 1 | 1.38577e-12 | 1.98211e+11 | 
| BiCGSTABSolver | HypreEuclid | 2048 | &#10004; | 1 | 2.13428e-12 | 1.25239e+11 | 
| BiCGSTABSolver | HypreEuclid | 4096 | &#10004; | 1 | 1.33095e-12 | 7.39624e+10 | 
| BiCGSTABSolver | HypreEuclid | 8192 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | HypreEuclid | 16384 | &#10004; | 1 | 6.85457e-12 | 4.06115e+10 | 
| BiCGSTABSolver | HypreBoomerAMG | 128 | &#10004; | 1 | 1.2778e-12 | 1.28692e+11 | 
| BiCGSTABSolver | HypreBoomerAMG | 512 | &#10004; | 1 | 1.32959e-12 | 3.02628e+10 | 
| BiCGSTABSolver | HypreBoomerAMG | 1024 | &#10004; | 1 | 4.28677e-12 | 1.40287e+10 | 
| BiCGSTABSolver | HypreBoomerAMG | 2048 | &#10004; | 1 | 1.41841e-12 | 5.6194e+09 | 
| BiCGSTABSolver | HypreBoomerAMG | 4096 | &#10004; | 1 | 1.70773e-12 | 6.0374e+09 | 
| BiCGSTABSolver | HypreBoomerAMG | 8192 | &#10004; | 1 | 1.31118e-12 | 4.5316e+09 | 
| BiCGSTABSolver | HypreBoomerAMG | 16384 | &#10004; | 1 | 1.64552e-12 | 6.36326e+09 | 
| BiCGSTABSolver | HypreILU | 128 | &#10004; | 2 | 9.97418e-13 | 1.34514e+11 | 
| BiCGSTABSolver | HypreILU | 512 | &#10004; | 1 | 3.1389e-12 | 3.26538e+10 | 
| BiCGSTABSolver | HypreILU | 1024 | &#10004; | 1 | 1.78469e-12 | 1.53952e+10 | 
| BiCGSTABSolver | HypreILU | 2048 | &#10004; | 1 | 1.68905e-12 | 7.03549e+09 | 
| BiCGSTABSolver | HypreILU | 4096 | &#10004; | 2 | 9.13438e-13 | 2.39929e+09 | 
| BiCGSTABSolver | HypreILU | 8192 | &#10004; | 1 | 2.05545e-12 | 1.70735e+09 | 
| BiCGSTABSolver | HypreILU | 16384 | &#10004; | 1 | 2.38409e-12 | 1.90885e+09 | 
| MINRESSolver | HypreDiagScale | 128 | &#10004; | 2 | 3.86984e-12 | 9.8413e+10 | 
| MINRESSolver | HypreDiagScale | 512 | &#10004; | 2 | 3.86429e-12 | 1.71968e+10 | 
| MINRESSolver | HypreDiagScale | 1024 | &#10004; | 2 | 3.92209e-12 | 6.34614e+09 | 
| MINRESSolver | HypreDiagScale | 2048 | &#10004; | 2 | 3.84437e-12 | 2.02862e+09 | 
| MINRESSolver | HypreDiagScale | 4096 | &#10004; | 2 | 3.90232e-12 | 1.30548e+09 | 
| MINRESSolver | HypreDiagScale | 8192 | &#10004; | 2 | 3.82029e-12 | 1.00985e+09 | 
| MINRESSolver | HypreDiagScale | 16384 | &#10004; | 2 | 3.8292e-12 | 1.82812e+09 | 
| MINRESSolver | HypreEuclid | 128 | &cross; | &cross; | &cross; | &cross; | 
| MINRESSolver | HypreEuclid | 512 | &#10004; | 2 | 1.23541e-12 | 5.92752e+11 | 
| MINRESSolver | HypreEuclid | 1024 | &#10004; | 2 | 1.14049e-12 | 3.56569e+11 | 
| MINRESSolver | HypreEuclid | 2048 | &#10004; | 2 | 1.19802e-12 | 2.28731e+11 | 
| MINRESSolver | HypreEuclid | 4096 | &#10004; | 2 | 1.21281e-12 | 1.39055e+11 | 
| MINRESSolver | HypreEuclid | 8192 | &#10004; | 2 | 1.24416e-12 | 9.09293e+10 | 
| MINRESSolver | HypreEuclid | 16384 | &#10004; | 2 | 1.32794e-12 | 7.44786e+10 | 
| MINRESSolver | HypreBoomerAMG | 128 | &#10004; | 1 | 4.54395e-12 | 8.69332e+10 | 
| MINRESSolver | HypreBoomerAMG | 512 | &#10004; | 1 | 2.46229e-12 | 2.15607e+10 | 
| MINRESSolver | HypreBoomerAMG | 1024 | &#10004; | 1 | 4.7549e-12 | 8.9393e+09 | 
| MINRESSolver | HypreBoomerAMG | 2048 | &#10004; | 1 | 2.55511e-12 | 4.31768e+09 | 
| MINRESSolver | HypreBoomerAMG | 4096 | &#10004; | 1 | 3.53415e-12 | 3.47602e+09 | 
| MINRESSolver | HypreBoomerAMG | 8192 | &#10004; | 1 | 2.65217e-12 | 3.46678e+09 | 
| MINRESSolver | HypreBoomerAMG | 16384 | &#10004; | 1 | 2.77066e-12 | 3.731e+09 | 
| MINRESSolver | HypreILU | 128 | &#10004; | 2 | 2.04345e-12 | 7.52408e+10 | 
| MINRESSolver | HypreILU | 512 | &#10004; | 2 | 2.30936e-12 | 1.78062e+10 | 
| MINRESSolver | HypreILU | 1024 | &#10004; | 2 | 2.37093e-12 | 8.41829e+09 | 
| MINRESSolver | HypreILU | 2048 | &#10004; | 2 | 2.4203e-12 | 3.6989e+09 | 
| MINRESSolver | HypreILU | 4096 | &#10004; | 2 | 2.47327e-12 | 1.28305e+09 | 
| MINRESSolver | HypreILU | 8192 | &#10004; | 2 | 2.78188e-12 | 8.27313e+08 | 
| MINRESSolver | HypreILU | 16384 | &#10004; | 2 | 2.82591e-12 | 8.9135e+08 | 


Refinement level : 3

 
| solver | preconditionner | number of mpi | converged | iterations | residu | time |
| :-------------|:---------------|:---------------:|:---------------:|:---------------:|:---------------:|---------------:|
| HyprePCG | HypreDiagScale | 256 | &#10004; | 3 | 5.99451e-14 | 1.80764e+12 | 
| HyprePCG | HypreDiagScale | 512 | &#10004; | 3 | 5.99392e-14 | 8.47542e+11 | 
| HyprePCG | HypreDiagScale | 1024 | &#10004; | 3 | 6.17399e-14 | 3.66342e+11 | 
| HyprePCG | HypreDiagScale | 2048 | &#10004; | 3 | 5.99354e-14 | 1.74938e+11 | 
| HyprePCG | HypreDiagScale | 4096 | &#10004; | 3 | 6.14796e-14 | 7.67449e+10 | 
| HyprePCG | HypreDiagScale | 8192 | &#10004; | 3 | 5.99556e-14 | 2.91296e+10 | 
| HyprePCG | HypreDiagScale | 16384 | &#10004; | 3 | 5.99843e-14 | 1.33259e+10 | 
| HyprePCG | HypreDiagScale | 32768 | &#10004; | 3 | 5.96553e-14 | 1.08015e+10 | 
| HyprePCG | HypreDiagScale | 65536 | &#10004; | 3 | 6.16933e-14 | 1.3964e+10 | 
| HyprePCG | HypreBoomerAMG | 256 | &#10004; | 2 | 5.99598e-14 | 1.54919e+12 | 
| HyprePCG | HypreBoomerAMG | 512 | &#10004; | 2 | 5.99723e-14 | 7.5246e+11 | 
| HyprePCG | HypreBoomerAMG | 1024 | &#10004; | 2 | 6.17408e-14 | 3.6757e+11 | 
| HyprePCG | HypreBoomerAMG | 2048 | &#10004; | 2 | 5.99671e-14 | 1.93026e+11 | 
| HyprePCG | HypreBoomerAMG | 4096 | &#10004; | 2 | 6.14726e-14 | 9.77932e+10 | 
| HyprePCG | HypreBoomerAMG | 8192 | &#10004; | 2 | 5.99531e-14 | 4.99741e+10 | 
| HyprePCG | HypreBoomerAMG | 16384 | &#10004; | 2 | 5.99578e-14 | 3.93017e+10 | 
| HyprePCG | HypreBoomerAMG | 32768 | &#10004; | 2 | 5.91947e-14 | 3.55983e+10 | 
| HyprePCG | HypreBoomerAMG | 65536 | &#10004; | 2 | 6.17071e-14 | 4.69211e+10 | 
| HyprePCG | HypreILU | 256 | &#10004; | 2 | 5.99416e-14 | 1.21189e+12 | 
| HyprePCG | HypreILU | 512 | &#10004; | 2 | 5.99207e-14 | 6.12919e+11 | 
| HyprePCG | HypreILU | 1024 | &#10004; | 2 | 6.16857e-14 | 3.0916e+11 | 
| HyprePCG | HypreILU | 2048 | &#10004; | 2 | 5.99917e-14 | 1.60395e+11 | 
| HyprePCG | HypreILU | 4096 | &#10004; | 2 | 6.14315e-14 | 9.73826e+10 | 
| HyprePCG | HypreILU | 8192 | &#10004; | 2 | 5.9959e-14 | 3.75567e+10 | 
| HyprePCG | HypreILU | 16384 | &#10004; | 2 | 5.99755e-14 | 2.29695e+10 | 
| HyprePCG | HypreILU | 32768 | &#10004; | 2 | 5.9215e-14 | 9.46843e+09 | 
| HyprePCG | HypreILU | 65536 | &#10004; | 2 | 6.17121e-14 | 1.8007e+10 | 
| CGSolver | HypreDiagScale | 256 | &cross; | &cross; | &cross; | &cross; | 
| CGSolver | HypreDiagScale | 512 | &cross; | &cross; | &cross; | &cross; | 
| CGSolver | HypreDiagScale | 1024 | &cross; | &cross; | &cross; | &cross; | 
| CGSolver | HypreDiagScale | 2048 | &cross; | &cross; | &cross; | &cross; | 
| CGSolver | HypreDiagScale | 4096 | &cross; | &cross; | &cross; | &cross; | 
| CGSolver | HypreDiagScale | 8192 | &cross; | &cross; | &cross; | &cross; | 
| CGSolver | HypreDiagScale | 16384 | &cross; | &cross; | &cross; | &cross; | 
| CGSolver | HypreDiagScale | 32768 | &cross; | &cross; | &cross; | &cross; | 
| CGSolver | HypreDiagScale | 65536 | &cross; | &cross; | &cross; | &cross; | 
| BiCGSTABSolver | HypreBoomerAMG | 256 | &#10004; | 1 | 3.47102e-12 | 1.10947e+12 | 
| BiCGSTABSolver | HypreBoomerAMG | 512 | &#10004; | 1 | 3.03433e-12 | 5.28839e+11 | 
| BiCGSTABSolver | HypreBoomerAMG | 1024 | &#10004; | 2 | 9.73126e-13 | 3.0943e+11 | 
| BiCGSTABSolver | HypreBoomerAMG | 2048 | &#10004; | 1 | 2.77819e-12 | 1.42461e+11 | 
| BiCGSTABSolver | HypreBoomerAMG | 4096 | &#10004; | 1 | 3.60123e-12 | 6.68014e+10 | 
| BiCGSTABSolver | HypreBoomerAMG | 8192 | &#10004; | 1 | 3.88785e-12 | 3.58955e+10 | 
| BiCGSTABSolver | HypreBoomerAMG | 16384 | &#10004; | 2 | 8.61773e-13 | 2.90985e+10 | 
| BiCGSTABSolver | HypreBoomerAMG | 32768 | &#10004; | 2 | 9.83403e-13 | 2.52815e+10 | 
| BiCGSTABSolver | HypreBoomerAMG | 65536 | &#10004; | 2 | 9.87234e-13 | 3.01919e+10 | 
| MINRESSolver | HypreDiagScale | 256 | &cross; | &cross; | &cross; | &cross; | 
| MINRESSolver | HypreDiagScale | 512 | &cross; | &cross; | &cross; | &cross; | 
| MINRESSolver | HypreDiagScale | 1024 | &cross; | &cross; | &cross; | &cross; | 
| MINRESSolver | HypreDiagScale | 2048 | &cross; | &cross; | &cross; | &cross; | 
| MINRESSolver | HypreDiagScale | 4096 | &cross; | &cross; | &cross; | &cross; | 
| MINRESSolver | HypreDiagScale | 8192 | &cross; | &cross; | &cross; | &cross; | 
| MINRESSolver | HypreDiagScale | 16384 | &cross; | &cross; | &cross; | &cross; | 
| MINRESSolver | HypreDiagScale | 32768 | &cross; | &cross; | &cross; | &cross; | 
| MINRESSolver | HypreDiagScale | 65536 | &cross; | &cross; | &cross; | &cross; | 
| MINRESSolver | HypreBoomerAMG | 256 | &#10004; | 2 | 1.37001e-12 | 7.80847e+11 | 
| MINRESSolver | HypreBoomerAMG | 512 | &#10004; | 2 | 1.39195e-12 | 3.778e+11 | 
| MINRESSolver | HypreBoomerAMG | 1024 | &#10004; | 2 | 1.43441e-12 | 1.82843e+11 | 
| MINRESSolver | HypreBoomerAMG | 2048 | &#10004; | 2 | 1.39478e-12 | 9.5925e+10 | 
| MINRESSolver | HypreBoomerAMG | 4096 | &#10004; | 2 | 1.4522e-12 | 4.96601e+10 | 
| MINRESSolver | HypreBoomerAMG | 8192 | &#10004; | 2 | 1.45188e-12 | 2.31648e+10 | 
| MINRESSolver | HypreBoomerAMG | 16384 | &#10004; | 2 | 1.48741e-12 | 1.33982e+10 | 
| MINRESSolver | HypreBoomerAMG | 32768 | &#10004; | 2 | 1.52153e-12 | 1.2014e+10 | 
| MINRESSolver | HypreBoomerAMG | 65536 | &#10004; | 2 | 1.56258e-12 | 1.25803e+10 | 
| MINRESSolver | HypreILU | 256 | &#10004; | 2 | 1.39464e-12 | 6.3744e+11 | 
| MINRESSolver | HypreILU | 512 | &#10004; | 2 | 1.45556e-12 | 3.15677e+11 | 
| MINRESSolver | HypreILU | 1024 | &#10004; | 2 | 1.52908e-12 | 1.56383e+11 | 
| MINRESSolver | HypreILU | 2048 | &#10004; | 2 | 1.57832e-12 | 8.42539e+10 | 
| MINRESSolver | HypreILU | 4096 | &#10004; | 2 | 1.61329e-12 | 5.07457e+10 | 
| MINRESSolver | HypreILU | 8192 | &#10004; | 2 | 1.75378e-12 | 1.92632e+10 | 
| MINRESSolver | HypreILU | 16384 | &#10004; | 2 | 1.73013e-12 | 1.01591e+10 | 
| MINRESSolver | HypreILU | 32768 | &#10004; | 2 | 1.83231e-12 | 7.67382e+09 | 
| MINRESSolver | HypreILU | 65536 | &#10004; | 2 | 1.8579e-12 | 4.71681e+09 | 