#pragma once
#include <mfem_mgis_headers.hxx>
#include <timer.hpp>
#include <solver_name.hxx>
#include <common.hxx>
#include <vector>

using namespace configuration;

// We need this class for test case sources
constexpr double xmax = 1.;
struct TestParameters {
	const char* mesh_file = "cas-cible-1.mesh";
	const char* behaviour = "Elasticity";
	const char* library = "src/libBehaviour.so";
	const char* reference_file = "Elasticity.ref";
	int order = 1;
	int tcase = 1;
	double xmax = 1.;
	double ymax = 1.;
	double zmax = 1.;
	int linearsolver = 1;
	int preconditioner = 1;
	bool parallel = true;
};

TestParameters parseCommandLineOptions(int& argc, char* argv[]);

// get test case (kernel, solver, preconditionner and match function)
#include <cas_cible_1.hxx>
#include <fissuration.hxx>



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

std::vector<solver_name> get_solvers(int a_case);
std::vector<petsc_ksp_type> get_solvers_with_petsc(int a_case);
std::vector<precond_name> get_pc(int a_case);
std::vector<petsc_pc_type> get_pc_with_petsc(int a_case);
std::function<bool(solver_name,precond_name)> get_match(int a_case);
std::function<bool(petsc_ksp_type,petsc_pc_type)> get_match_with_petsc(int a_case);
