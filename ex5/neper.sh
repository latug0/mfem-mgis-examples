NBG=8
neper -T -domain "cube(1,1,1)" -n ${NBG} -morpho graingrowth -periodicity all
neper -M -format msh:binary -rcl 0.8 n${NBG}-id1.tess 
