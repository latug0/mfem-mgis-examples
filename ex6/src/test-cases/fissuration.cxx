#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/ImposedDirichletBoundaryConditionAtClosestNode.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include <test-cases/fissuration.hxx>

namespace fissuration
{
	using namespace configuration;
	
	using solver=solver_name;
	using pc=precond_name;

	using petsc_solver=petsc_ksp_type;
	using petsc_pc=petsc_pc_type;

	using regular_solver_list=std::vector<solver>;
	using petsc_solver_list=std::vector<petsc_solver>;

	using regular_pc_list=std::vector<pc>;
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
			solver_name::MUMPSSolver,
			solver_name::GMRESSolver,
			solver_name::CGSolver,
			solver_name::BiCGSTABSolver,
			solver_name::UMFPackSolver,
			solver_name::MINRESSolver,
			solver_name::SLISolver,
		};
		return res;
	}

	/*
	 * brief return a list of preconditioners that work or do not run slowly for the test case 1
	 */
	regular_pc_list build_pc_list()
	{

		regular_pc_list res = {
			precond_name::HypreParaSails,
			precond_name::HypreBoomerAMG,
			precond_name::HypreILU,
			precond_name::HypreEuclid,
			precond_name::HypreDiagScale
		};
		return res;	
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

	bool match(solver s, pc p)
	{
		return true;
	}

	bool match_with_petsc(petsc_solver s, petsc_pc p)
	{
		return true;
	}
}
