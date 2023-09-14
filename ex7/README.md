# Representative Volume Element of Combustible Mixed Oxides for Nuclear Applications

This simulation represents an RVE of MOx (Mixed Oxide) material under uniform macroscopic
deformation. The aim of this simulation is to reproduce and compare the results obtained by
(Fauque et al., 2021; Masson et al., 2020) who used an FFT method.

## Problem solved

```text
    Problem : RVE MOx 2 phases with elasto-viscoplastic behavior laws

    Parameters : 

    start time = 0
    end time = 5s
    number of time step = 40

    Strain Gradient matrix : val = 0.012
    [ - val / 2 ,         0 ,   0 ]
    [ 0         , - val / 2 ,   0 ] 
    [ 0         ,          0, val ]
    
    Solver : HyprePCG
    Preconditionner : HypreBoomerAMG

    Behavior law parameters : ImplicitNortonThreshold
    [ parameters       , matrix   , inclusions ]    
    [ Young Modulus    , 8.182e9  , 2*8.182e9  ];
    [ Poisson Ratio    , 0.364    , 0.364      ];
    [ Stress Threshold , 100.0e6  , 100.0e12   ];
    [ Norton Exponent  , 3.333333 , 3.333333   ];
    [ Temperature      , 293.15   , 293.15     ];

    Element :

    Familly H1
    Order 2
```

![Illustration of a RVE with 634 spheres after 5 seconds.](./results/order2.png){width=50%}

## How to run the simulation "RVE MOX"

## Build the mesh

The mesh is generated with MEROPE and GMSH through the following steps:

- First step, use MEROPE to generate a `.geo` file using the RSA algorithm. Scripts are in directory `script_merope`. Command line:

```bash
# generate .geo file with MEROPE
python3 script_17percent_minimal.py
```

- Second step, use GMSH to mesh the geometry. Files `.geo` are in the directory `file_geo`. Command line:

```bash
# generate the .msh file with GMSH
gmsh -3 OneSphere.geo 
```

## Run the simulation

### Run the minimal version of this example

Run the simulation with the command line:
```bash
# run the simulation by specifying the mesh with --mesh option
./mox2 --mesh OneSphere.msh
```
Run the simulation with parralel solver with the command line:
```bash
# run the simulation by specifying the mesh with --mesh option
mpirun -n 12 ./mox2 --mesh 634Spheres.msh
```

### Available options

Command line | Descritption
---|---
--mesh or -m | specify the mesh ".msh" usedn default = inclusion.msh
--refinement or -r | refinement level of the mesh, default = 0
--order or -o | Finite element order (polynomial degree), default = 2
--verbosity-level or -v | choose the verbosity level, default = 0
--post-processing or -p | run post processing step, default = 1

### CCRT/Topaze

```bash
ccc_mprun -n 8 -c 1 -p milan ./mox2 -r 0 -o 3 --mesh OneSphere.msh
ccc_mprun -n 2048 -c 1 -p milan ./mox2 -r 2 -o 1 --mesh 634Sphere.msh
```

## Display results

Comparison by seeing the value of SZZ in function of the time.

### Extract data

```bash
awk '{if(NR>13) print $1 " " 0.83*$4+0.17*$10}' avgStress > res-mfem-mgis.txt
```

Reference values can be found in the directory `results`, file res-fft.txt. Experimental values are available: `results/res-mfem-mgis-onesphere-o3.txt` and `results/res-mfem-mgis-634sphere-o2.txt`.

### Display results with gnuplot

```bash
plot "res-fft.txt" u 1:10 w l title "fft"
gnuplot> replot "res-mfem-mgis-onesphere-o3.txt" u 1:2 w l title "mfem-mgis"
```

## Illusations

Illustrations are available in the directory `gif`
