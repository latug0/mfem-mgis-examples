#include <parameters/catch_functions.hxx>

std::vector<solver_name> get_solvers(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::build_solvers_list();
		case 2: return fissuration::build_solvers_list();
		default: return cas_cible_1::build_solvers_list();
	}
}

std::vector<precond_name> get_pc(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::build_pc_list();
		case 2: return fissuration::build_pc_list();
		default: return cas_cible_1::build_pc_list();
	}
}

std::function<bool(solver_name,precond_name,int)> get_match(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::match;
		case 2: return fissuration::match;
		default: return cas_cible_1::match;
	}
}

std::vector<petsc_ksp_type> get_solvers_with_petsc(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::build_solvers_list_with_petsc();
		case 2: return fissuration::build_solvers_list_with_petsc();
		default: return cas_cible_1::build_solvers_list_with_petsc();
	}
}

std::vector<petsc_pc_type> get_pc_with_petsc(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::build_pc_list_with_petsc();
		case 2: return fissuration::build_pc_list_with_petsc();
		default: return cas_cible_1::build_pc_list_with_petsc();
	}
}

std::function<bool(petsc_ksp_type,petsc_pc_type,int)> get_match_with_petsc(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::match_with_petsc;
		case 2: return fissuration::match_with_petsc;
		default: return cas_cible_1::match_with_petsc;
	}
}
