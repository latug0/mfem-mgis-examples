#pragma once

namespace configuration
{


	template<typename Problem>		
		static void setLinearSolver(Problem& p,
				const solver_name a_solver_name, 
				const precond_name a_precond_name,
				const int v,
				const mfem_mgis::real a_abs_tol = 1e-12,
				const mfem_mgis::real a_rel_tol = 1e-12
				)
		{
			CatchTimeSection("set_linear_solver");
			// pilote
			const int verbosity 			= v; 			// VerbosityLevel
			constexpr int defaultMaxNumOfIt	 	= 5000; 		// MaximumNumberOfIterations
			constexpr int adjustMaxNumOfIt 		= 500000; 		// MaximumNumberOfIterations
			const mfem_mgis::real absTol		= a_abs_tol; 		// AbsoluteTolerance
			const mfem_mgis::real relTol		= a_rel_tol;		// RelativeTolerance
			//constexpr mfem_mgis::real relTol	= 1e-12; 		// RelativeTolerance
			const mfem_mgis::real Tol		= a_rel_tol;  		// Tolerance



			auto solverParameters = mfem_mgis::Parameters{};

			if(!isMumps(a_solver_name))
			{
				solverParameters.insert(mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
			};

			if(!isMumps(a_solver_name))
			{
				solverParameters.insert(mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
			}

			if(!isHypre(a_solver_name) && !isMumps(a_solver_name))
			{
				solverParameters.insert(mfem_mgis::Parameters{{"AbsoluteTolerance", absTol}});
				solverParameters.insert(mfem_mgis::Parameters{{"RelativeTolerance", relTol}});
			}

			if(isHypre(a_solver_name))
			{
				solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", Tol}});
			}
#ifdef MFEM_USE_MUMPS
			if(isMumps(a_solver_name))
			{
				solverParameters.insert(mfem_mgis::Parameters{{"Symmetric", true}});
				solverParameters.insert(mfem_mgis::Parameters{{"PositiveDefinite", true}});	
			}
#endif
			if(a_precond_name != precond_name::ANY)
			{
				auto precond = buildPrecond(a_precond_name, verbosity);
				solverParameters.insert(mfem_mgis::Parameters{{"Preconditioner",precond}});
			}


			auto current_solver_name = getName(a_solver_name);
			auto current_precond_name = getName(a_precond_name);
			Profiler::Utils::Message("INFO: solver: ", current_solver_name, " and preconditionner: ",current_precond_name);
			p.setLinearSolver(current_solver_name, solverParameters);
		}

	static void setSolverParameters(
			mfem_mgis::AbstractNonLinearEvolutionProblem& problem) {
		CatchTimeSection("set_solve_parameters");
		problem.setSolverParameters({{"VerbosityLevel", 0},
				{"RelativeTolerance", 1e-12},
				{"AbsoluteTolerance", 1e-12},
				{"MaximumNumberOfIterations", 10}});
	}  // end of setSolverParmeters

};
