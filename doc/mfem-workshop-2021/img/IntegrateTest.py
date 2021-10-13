import mgis.behaviour as mgis_bv

# strain increment per time step
de = 5.e-5
# time step
d.dt = 180

# setting the temperature
b = mgis_bv.load('src/libBehaviour.so', 'Norton',
                 mgis_bv.Hypothesis.Tridimensional)
d = mgis_bv.BehaviourData(b)
o = mgis_bv.getVariableOffset(b.isvs, 'EquivalentViscoplasticStrain',
                              b.hypothesis)
mgis_bv.setExternalStateVariable(d.s1, 'Temperature', 293.15)

# copy d.s1 in d.s0
mgis_bv.update(d)
d.s1.gradients[0] = de

# integrate the behaviour
for i in range(0, 20):
    mgis_bv.integrate(d, b)
    mgis_bv.update(d)
    d.s1.gradients[0] += de
