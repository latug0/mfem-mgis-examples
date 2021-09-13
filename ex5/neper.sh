NBG=8
#neper -T -domain "cube(1,1,1)" -n ${NBG} -morpho graingrowth -periodicity all
neper -M -format msh:ascii -rcl 0.3 n${NBG}-id1.tess 
