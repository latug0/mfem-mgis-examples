#pragma once


bool match(solver_name s, precond_name p)
{
	if(p == precond_name::ANY) return true;
	if(isMumps(s)) return false;

	if(p == precond_name::HypreBoomerAMG)
	{
		if(s==CGSolver) return false; // setup error
		return true; // no mumps
	}
	if(p == precond_name::HypreDiagScale) 
	{
		return true; // NU MUMPS
	}
	if(p == precond_name::HypreParaSails) 
	{
		if(s==solver_name::HyprePCG) return false; // seg fault
		return true; // NU MUMPS
	}
	if(p == precond_name::HypreEuclid) 
	{
		if(s==solver_name::CGSolver) return false;
		if(s==solver_name::HypreFGMRES) return false; // seg fault
		if(s==solver_name::HypreGMRES) return false; // seg fault
		return true; // NU MUMPS
	}
	if(p == precond_name::HypreILU) 
	{
		if(s==solver_name::CGSolver) return false; // setup error
		return true; // NU MUMPS
	}
	return false;
}

static void setLinearSolver(mfem_mgis::AbstractNonLinearEvolutionProblem& p,
                            const solver_name a_solver_name, 
			    const precond_name a_precond_name)
{
	START_TIMER("set_linear_solver");
	// pilote
	constexpr int verbosity 		= 0; 		// VerbosityLevel
	constexpr int defaultMaxNumOfIt	 	= 5000; 	//MaximumNumberOfIterations
	constexpr int adjustMaxNumOfIt 		= 500000; 	//MaximumNumberOfIterations
	constexpr mfem_mgis::real absTol	= 1e-12; 	// AbsoluteTolerance
	constexpr mfem_mgis::real relTol	= 0; 		// RelativeTolerance
	//constexpr mfem_mgis::real relTol	= 1e-12; 	// RelativeTolerance
	constexpr mfem_mgis::real Tol		= 1e-12;        // Tolerance


  
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
	profiling::output::printMessage("INFO: solver: ", current_solver_name, " and preconditionner: ",current_precond_name);
	p.setLinearSolver(current_solver_name, solverParameters);
}

static void setSolverParameters(
	mfem_mgis::AbstractNonLinearEvolutionProblem& problem) {
	START_TIMER("set_solve_parameters");
	problem.setSolverParameters({{"VerbosityLevel", 0},
				       {"RelativeTolerance", 1e-12},
				       {"AbsoluteTolerance", 1e-12},
				       {"MaximumNumberOfIterations", 10}});
}  // end of setSolverParmeters


