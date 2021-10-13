const auto b = load("libBehaviour.so", "Norton", Hypothesis::TRIDIMENSIONAL);
auto d = BehaviourData{b};
const auto o = getVariableOffset(b.isvs, "EquivalentViscoplasticStrain", b.hypothesis);
const auto de = 5.e-5;
d.dt = 180;  
// initialize the states
setExternalStateVariable(d.s1, "Temperature", 293.15);
// copy d.s1 in d.s0
update(d);
d.s1.gradients[0] = de;
// integration
auto p = std::array<real, 21>{};
p[0] = d.s0.internal_state_variables[o];
for (size_type i = 0; i != 20; ++i) {
  auto v = make_view(d);
  v.rdt = 1;
  integrate(v, b);
  update(d);
  d.s1.gradients[0] += de;
}
