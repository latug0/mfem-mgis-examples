>>> import mgis.behaviour as mgis_bv
>>> o = mgis_bv.FiniteStrainBehaviourOptions()
>>> o.stress_measure = mgis_bv.FiniteStrainBehaviourOptionsStressMeasure.PK1
>>> o.tangent_operator = mgis_bv.FiniteStrainBehaviourOptionsTangentOperator.DPK1_DF
>>> b = mgis_bv.load(o, 'src/libBehaviour.so', 'LogarithmicStrainPlasticity',
...                  mgis_bv.Hypothesis.Tridimensional)
>>> for v in b.internal_state_variables:
...     print('-' + v.name + '(' + v.getType() + ')')
- ElasticStrain(Stensor)
- EquivalentPlasticStrain(Scalar)