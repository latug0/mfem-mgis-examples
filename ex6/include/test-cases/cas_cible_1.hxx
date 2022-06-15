#pragma once

#include <common/timer.hpp>
#include <parameters/test_parameters.hpp>
#include <solver/solver_name.hxx>
#include <solver/precond_name.hxx>
#include <solver/config_solver.hxx>
#include <common/data_gathering.hxx>
#include <common/common.hxx>

namespace cas_cible_1
{
	using namespace configuration;
	using solver=solver_name;
	using pc=precond_name;
	using regular_solver_list=std::vector<solver>;
	using regular_pc_list=std::vector<pc>;
	regular_solver_list 	build_solvers_list();
	regular_pc_list 	build_pc_list();
	bool match(solver s, pc p, int r);

	using petsc_solver=petsc_ksp_type;
	using petsc_pc=petsc_pc_type;
	using petsc_solver_list=std::vector<petsc_solver>;
	using petsc_pc_list=std::vector<petsc_pc>;

	petsc_solver_list 	build_solvers_list_with_petsc();
	petsc_pc_list 		build_pc_list_with_petsc();
	bool match_with_petsc(petsc_solver, petsc_pc, int r);
	
	void setup_properties(const TestParameters& p, mfem_mgis::PeriodicNonLinearEvolutionProblem& problem);

	template<typename Solver, typename Pc>
	int kernel(const TestParameters& p, const bool use_post_processing, const Solver a_solv, const Pc a_precond, gather_information& a_info)
	{
		constexpr const auto dim = mfem_mgis::size_type{3};

		std::string string_solver       = getName(a_solv);
		std::string string_precond      = getName(a_precond);
		std::string timer_name          = "run_" + string_solver + "_" + string_precond;
		START_TIMER(timer_name);

		// creating the finite element workspace
		auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
				mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
				{"FiniteElementFamily", "H1"},
				{"FiniteElementOrder", p.order},
				{"UnknownsSize", dim},
				{"NumberOfUniformRefinements", p.parallel ? p.refinement : 0},
				{"Parallel", p.parallel}});

		mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed);
		//common::print_mesh_information(problem);
		common::print_mesh_information(problem.getImplementation<true>());
		setup_properties(p, problem);

		setLinearSolver(problem, a_solv, a_precond);
		setSolverParameters(problem);

		if(use_post_processing) common::add_post_processings(problem, "OutputFile-cas_cible_1");


		mfem_mgis::NonLinearResolutionOutput solver_statistics;
		double measure = 0.;

		// main function here
		measure = common::solve(problem, 0, 1, solver_statistics, string_solver, string_precond);


		if(use_post_processing)	common::execute_post_processings(problem, 0,1);
		common::fill_statistics(a_info, a_solv, a_precond, solver_statistics, measure);
		common::print_statistics(string_solver, string_precond, measure);

		return(EXIT_SUCCESS);
	}
};
