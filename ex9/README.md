## Usage of third-medium test cases
1. MFEM-MGIS and dependencies loaded
2. Compile using the `compile_thirdmedium.sh` script
3. Move to the build directory and run the executable (`BendingTest`,`BendingTest3D`,`CShapeContact`) with some options. Exemple : 
```sh
cd thirdmedium/build
OMP_NUM_THREADS=1 mpirun -n 4 ./BendingTest -pp 1 -m ../third_medium_quad.sh -g 5e7 -a 1e11 -o output_dir/
```

Notes : Option `-pp 1` to enable post-processing paraview output 

./exec -h gives all the options
