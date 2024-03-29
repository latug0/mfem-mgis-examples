#include <test-cases/check_manta_rev_elastic.hxx>

namespace check_manta_rev_elastic
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

	/*
	 * 
	 */

	bool match(solver s, pc p, int r)
	{
		if(r>=2)
		{
			// these possibilities work but are too slow
			if(s==solver::HypreFGMRES) return false;
			if(s==solver::HypreGMRES) return false;
			if(s==solver::BiCGSTABSolver && p==pc::HypreDiagScale) return false;
		}

		// too slow 
		if(r >= 3)
		{
			// these possibilities work but are too slow
			if(s==solver::HypreFGMRES && p==pc::HypreBoomerAMG) return false;
			if(s==solver::HypreFGMRES && p==pc::HypreILU) return false;
			if(s==solver::HypreFGMRES && p==pc::HypreDiagScale) return false;
			if(s==solver::HypreGMRES && p==pc::HypreBoomerAMG) return false;
			if(s==solver::HypreGMRES && p==pc::HypreILU) return false;
			if(s==solver::HypreGMRES && p==pc::HypreDiagScale) return false;
			if(s==solver::BiCGSTABSolver && p==pc::HypreILU) return false;
			if(s==solver::BiCGSTABSolver && p==pc::HypreEuclid) return false;
			if(s==solver::BiCGSTABSolver && p==pc::HypreDiagScale) return false;
			if(s==solver::MINRESSolver && p==pc::HypreEuclid) return false;
		}

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
	
	bool match_with_petsc(petsc_solver s, petsc_pc p, int r)
	{
		return true;
	}
	
	void setup_properties(const TestParameters& p, mfem_mgis::PeriodicNonLinearEvolutionProblem& problem)
	{
		using namespace mgis::behaviour;
		using real=mfem_mgis::real;

		CatchTimeSection("set_mgis_stuff");
		problem.addBehaviourIntegrator("Mechanics", 1, p.library, "SaintVenantKirchhoffElasticity");
		problem.addBehaviourIntegrator("Mechanics", 2, p.library, "SaintVenantKirchhoffElasticity");
		// materials
		auto& m1 = problem.getMaterial(1);
		auto& m2 = problem.getMaterial(2);
		auto set_properties = [](auto& m, const double yo, const double po) {
			setMaterialProperty(m.s0, "YoungModulus", yo);
			setMaterialProperty(m.s0, "PoissonRatio", po);
			setMaterialProperty(m.s1, "YoungModulus", yo);
			setMaterialProperty(m.s1, "PoissonRatio", po);
		};

		set_properties(m1, 2.0e11       , 0.3);
		set_properties(m2, 8.0e11       , 0.3);

		//
		auto set_temperature = [](auto& m) {
			setExternalStateVariable(m.s0, "Temperature", 293.15);
			setExternalStateVariable(m.s1, "Temperature", 293.15);
		};
		set_temperature(m1);
		set_temperature(m2);



		// macroscopic strain
		std::vector<real> e(9, real{0});
		e[0] = 1.1;
		e[1] = 1.0;
		e[2] = 1.0;
		problem.setMacroscopicGradientsEvolution([e](const double) { return e; });
	} 
};
