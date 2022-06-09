#pragma once

#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"

namespace configuration
{
	enum precond_name
	{
		HypreDiagScale, //JACOBI
		HypreParaSails, 
		HypreEuclid, // Euclid implements the Parallel Incomplete LU factorization technique
		HypreBoomerAMG,
		HypreILU,
		ANY,
		precond_name_count
	};


	enum class petsc_pc_type
	{
		jacobi,
		bjacobi,
		sor,
		eisenstat,
		icc,
		ilu,
		gasm,
		gamg,
		bddc,
		ksp,
		lu,
		cholesky,
		none,
		petsc_pc_type_count
	};

	std::string getName(const precond_name a_name);
	std::string getName(petsc_pc_type in);
	mfem_mgis::Parameters buildPrecond(const precond_name a_name, const int a_verbosity);
};
