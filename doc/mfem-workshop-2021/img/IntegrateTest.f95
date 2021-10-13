use mgis
use mgis_behaviour
!
...
! strain increment
de = 5.d-5
call check_status(load_behaviour(b, 'libBehaviour.so', 'Norton', 'Tridimensional'))
call check_status(allocate_behaviour_data(d,b))
call check_status(behaviour_data_get_state_0(s0, d))
call check_status(behaviour_data_get_state_1(s1, d))
call check_status(state_set_external_state_variable_by_name(s0, 'Temperature', 293.15d0))
call check_status(behaviour_get_internal_state_variable_offset(o, b, 'EquivalentViscoplasticStrain'))
call check_status(behaviour_data_set_time_increment(d, 180d0))
e = (/ de, 0d0, 0d0, 0d0, 0d0, 0d0 /)
call check_status(state_set_gradient_by_name(s1, 'Strain', e))
do i = 1, 20
   call check_status(integrate(ri, d, b))
   r = check(ri.eq.1, 'integration failed')
   call check_status(update_behaviour_data(d))
   e(1) = e(1) + de
   call check_status(state_set_gradient_by_name(s1, 'Strain', e))
end do
