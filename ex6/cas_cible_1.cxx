#include "test_parameters.hpp"
#include "solver_name.hxx"
#include "timer.hpp"
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

	using petsc_solver=petsc_ksp_type;
	using petsc_pc=petsc_pc_type;
	using petsc_solver_list=std::vector<petsc_solver>;
	using petsc_pc_list=std::vector<petsc_pc>;
	
	/*
	 * brief return a list of solvers that work or do not run slowly for the test case 1
	 */
	regular_solver_list build_solvers_list()
	{
		regular_solver_list res = {
			solver_name::HypreFGMRES,
			solver_name::HypreGMRES,
			solver_name::HyprePCG,
			solver_name::CGSolver,
			solver_name::BiCGSTABSolver,
			solver_name::MINRESSolver
		};
		return res;
	}

	/*
	 * brief return a list of preconditioners that work or do not run slowly for the test case 1
	 */
	regular_pc_list build_pc_list()
	{

		regular_pc_list res = {
			precond_name::HypreBoomerAMG,
			precond_name::HypreILU,
			precond_name::HypreEuclid,
			precond_name::HypreDiagScale
		};
		return res;	
	}
	
	bool match(solver s, pc p)
	{
		if(p == pc::ANY) return true;
		if(isMumps(s)) return false;

		if(p == pc::HypreBoomerAMG)
		{
			if(s==CGSolver) return false; // setup error
			return true; // no mumps
		}
		if(p == pc::HypreDiagScale) 
		{
			return true; // NU MUMPS
		}
		if(p == pc::HypreParaSails) 
		{
			if(s==solver::HyprePCG) return false; // seg fault
			return true; // NU MUMPS
		}
		if(p == pc::HypreEuclid) 
		{
			if(s==solver::CGSolver) return false;
			if(s==solver::HypreFGMRES) return false; // seg fault
			if(s==solver::HypreGMRES) return false; // seg fault
			if(s==solver::HyprePCG) return false; // seg fault
			return true; // NU MUMPS
		}
		if(p == pc::HypreILU) 
		{
			if(s==solver::CGSolver) return false; // setup error
			return true; // NU MUMPS
		}
		return false;
	}

	/*
	 * brief return a list of petsc solvers that work or do not run slowly for the test case 1
	 */
	petsc_solver_list build_solvers_list_with_petsc()
	{
		petsc_solver_list res = {
			petsc_solver::richardson,	petsc_solver::chebyshev,	petsc_solver::cg,
			petsc_solver::gmres,		petsc_solver::tcqmr,		petsc_solver::bcgs,
			petsc_solver::cgs,		petsc_solver::tfqmr,		petsc_solver::cr,
			petsc_solver::lsqr,		petsc_solver::bicg,		petsc_solver::preonly
		};
		return res;

	}

	petsc_pc_list build_pc_list_with_petsc()
	{
		petsc_pc_list res = {
			petsc_pc::jacobi,	petsc_pc::bjacobi,	petsc_pc::sor,
			petsc_pc::eisenstat,	petsc_pc::icc,		petsc_pc::ilu,
			petsc_pc::gasm,		petsc_pc::gamg,		petsc_pc::bddc,
			petsc_pc::ksp,		petsc_pc::lu,		petsc_pc::cholesky,
			petsc_pc::none
		};
		return res;
	}
	
	bool match_with_petsc(petsc_solver s, petsc_pc p)
	{
		return true;
	}
	

	int kernel(const TestParameters& p, const bool use_post_processing, const solver_name a_solv, const precond_name a_precond, gather_information& a_info) 
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
				{"NumberOfUniformRefinements", p.parallel ? 2 : 0},
				{"Parallel", p.parallel}});

		mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed);

		//getMesh
		auto mesh = problem.getImplementation<true>().getFiniteElementSpace().GetMesh();
		//get the number of vertices
		int numbers_of_vertices = mesh->GetNV();
		//get the number of elements
		int numbers_of_elements = mesh->GetNE();
		//get the element size
		double h = mesh->GetElementSize(0);

		profiling::output::printMessage("INFO: number of vertices -> ", numbers_of_vertices);
		profiling::output::printMessage("INFO: number of elements -> ", numbers_of_elements);
		profiling::output::printMessage("INFO: element size -> ", h);

		{
			START_TIMER("set_mgis_stuff");
			problem.addBehaviourIntegrator("Mechanics", 1, p.library, "Elasticity");
			problem.addBehaviourIntegrator("Mechanics", 2, p.library, "Elasticity");
			// materials
			auto& m1 = problem.getMaterial(1);
			auto& m2 = problem.getMaterial(2);
			// setting the material properties
			auto set_properties = [](auto& m, const double l, const double mu) {
				mgis::behaviour::setMaterialProperty(m.s0, "FirstLameCoefficient", l);
				mgis::behaviour::setMaterialProperty(m.s0, "ShearModulus", mu);
				mgis::behaviour::setMaterialProperty(m.s1, "FirstLameCoefficient", l);
				mgis::behaviour::setMaterialProperty(m.s1, "ShearModulus", mu);
			};

			std::array<mfem_mgis::real,2> lambda({100, 200});
			std::array<mfem_mgis::real,2>     mu({75 , 150});
			set_properties(m1, lambda[0], mu[0]);
			set_properties(m2, lambda[1], mu[1]);
			//
			auto set_temperature = [](auto& m) {
				mgis::behaviour::setExternalStateVariable(m.s0, "Temperature", 293.15);
				mgis::behaviour::setExternalStateVariable(m.s1, "Temperature", 293.15);
			};
			set_temperature(m1);
			set_temperature(m2);

			// macroscopic strain
			std::vector<mfem_mgis::real> e(6, mfem_mgis::real{});
			if (p.tcase < 3) {
				e[p.tcase] = 1;
			} else {
				e[p.tcase] = 1.41421356237309504880 / 2;
			}
			problem.setMacroscopicGradientsEvolution([e](const double) { return e; });
		} // end timer set_mgis_stuff

		setLinearSolver(problem, a_solv, a_precond);
		setSolverParameters(problem);

		if(use_post_processing)
		{
			START_TIMER("add_postprocessing_and_outputs");
			problem.addPostProcessing(
					"ParaviewExportResults",
					{{"OutputFileName", "PeriodicTestOutput-" + std::to_string(p.tcase)}}
					);
		} // end timer add_postprocessing_and_outputs

		mfem_mgis::NonLinearResolutionOutput solverStatistics;
		double measure = 0.;

		{
			START_TIMER("Solve");
			// solving the problem
			measure = profiling::timers::chrono_section( [&](){
					solverStatistics = problem.solve(0, 1);
					});

			// check status
			if (!solverStatistics) {
				profiling::output::printMessage("INFO: ", string_solver,"+",string_precond," FAILED");
			}
		} // end timer Solve

		if(use_post_processing)
		{
			START_TIMER("Postprocessing_step");
			problem.executePostProcessings(0, 1);
		} // end timer Postprocessing_step

		{
			START_TIMER("get_statistics");
			measure = profiling::output::reduce_max(measure);
			// fill info data.
			a_info.add(
					info{
					a_solv,
					a_precond,
					solverStatistics.status,
					solverStatistics.iterations,
					solverStatistics.final_residual_norm,
					measure
					}
				  );
			// print information
			profiling::output::printMessage(
					"INFO: ",
					string_solver,
					"+",
					string_precond,
					" works correctly, elapsed time: ",
					measure
					);

		} // end timer get_statistics

		profiling::output::printMessage("");
		return(EXIT_SUCCESS);
	}
};
