const mgis_real de = 5.e-5;
int r;  // behaviour integration result
check_status(mgis_bv_load_behaviour(&b, "libBehaviour.so", "Norton", "Tridimensional"));
check_status(mgis_bv_allocate_behaviour_data(&d, b));
// getting the offset of the equivalent plastic strain
check_status(mgis_bv_behaviour_get_internal_state_variable_offset(
    &o, b, "EquivalentViscoplasticStrain"));
// setting the time increment
check_status(mgis_bv_behaviour_data_set_time_increment(d, 180));
// state at the beginning of the time step
check_status(mgis_bv_behaviour_data_get_state_0(&s0, d));
// state at the end of the time step
check_status(mgis_bv_behaviour_data_get_state_0(&s1, d));
// setting the temperature
check_status(mgis_bv_state_set_external_state_variable_by_name(s0, "Temperature", 293.15));
check_status(mgis_bv_state_set_external_state_variable_by_name(s1, "Temperature", 293.15));
// initializing the view
check_status(mgis_bv_make_behaviour_data_view(&v, d));
// strain at the end of the time step
check_status(mgis_bv_state_get_gradient_by_name(&e, s1, "Strain"));
*e += de;
// loop over the time steps
for (i = 0; i != 20; ++i) {
  check_status(mgis_bv_integrate(&r, &v, b));
  check_status(mgis_bv_update_behaviour_data(d));
  *e += de;
}
