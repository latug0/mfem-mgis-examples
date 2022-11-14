#pragma once

#include <solver/solver_name.hxx>
#include <solver/precond_name.hxx>

// get test case (kernel, solver, preconditionner and match function)
#include <test-cases/cas_cible_1.hxx>
#include <test-cases/fissuration.hxx>



template<typename Solver, typename Pc>
std::function<void(const TestParameters&, const bool, const Solver, const Pc, gather_information&)> get_kernel(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::kernel<Solver,Pc>;
		case 2: return fissuration::kernel<Solver,Pc>;
		default: return cas_cible_1::kernel<Solver,Pc>;
	}
}

std::function<int(const TestParameters&, const bool, const char*, gather_information&)> get_petsc_kernel_with_file(int a_case);
std::vector<solver_name> get_solvers(int a_case);
std::vector<petsc_ksp_type> get_solvers_with_petsc(int a_case);
std::vector<precond_name> get_pc(int a_case);
std::vector<petsc_pc_type> get_pc_with_petsc(int a_case);
std::function<bool(solver_name,precond_name,int)> get_match(int a_case);
std::function<bool(petsc_ksp_type,petsc_pc_type,int)> get_match_with_petsc(int a_case);
