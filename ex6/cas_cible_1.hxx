#pragma once

#include "timer.hpp"
#include "test_parameters.hpp"
#include "solver_name.hxx"
#include "precond_name.hxx"
#include "config_solver.hxx"
#include "data_gathering.hxx"


namespace cas_cible_1
{
	using namespace configuration;
	using solver=solver_name;
	using pc=precond_name;
	using regular_solver_list=std::vector<solver>;
	using regular_pc_list=std::vector<pc>;
	regular_solver_list 	build_solvers_list();
	regular_pc_list 	build_pc_list();
	bool match(solver, pc);

	using petsc_solver=petsc_ksp_type;
	using petsc_pc=petsc_pc_type;
	using petsc_solver_list=std::vector<petsc_solver>;
	using petsc_pc_list=std::vector<petsc_pc>;

	petsc_solver_list 	build_solvers_list_with_petsc();
	petsc_pc_list 		build_pc_list_with_petsc();
	bool match_with_petsc(petsc_solver, petsc_pc);
	int kernel(const TestParameters& p, const bool use_post_processing, const solver_name a_solv, const precond_name a_precond, gather_information& a_info);
};
