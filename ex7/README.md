# How to run the simulation "RVE MOX"

## How to build the mesh

The mesh is generated with MEROPE and GMSH. First step, use MEROPE to generate a file `.geo` using the RSA algorithm. Scripts are in `script_merope`. Second step, use GMSH to mesh the geometry. Files `.geo` are in the directory `file_geo`.

Command lines:

```
python3 script_17percent_minimal.py
gmsh -3 OneSphere.geo 
```

## Run the simulation

### Run the minimal version of this example

```
./mox2 -o 3 --mesh OneSphere.msh
mpirun -n 12 ./mox2 -o 1 --mesh 634Spheres.msh
```

### CCRT/Topaze

```
ccc_mprun -n 512 -c 1 -p milan mox2 -r 0 -o 1 --mesh OneSphere.msh
```

## Display results

Comparison by seeing the value of SZZ in function of the time.

## Extract data

```
awk '{if(NR>13) print $1 " " 0.83*$4+0.17*$10}' avgStress > res-mfem-mgis.txt
```

Reference values are in the directory `results`, file fft.txt. Experimental values are available: `results/res-mfem-mgis-onesphere-o3.txt`.

## Display results with gnuplot

```
plot "res-fft.txt" u 1:10 w l title "fft"
gnuplot> replot "res-mfem-mgis-onesphere-o3.txt" u 1:2 w l title "mfem-mgis"
```
