#pragma once

#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"

enum precond_name
{
	HypreDiagScale, //JACOBI
	HypreParaSails, 
	HypreEuclid, // Euclid implements the Parallel Incomplete LU factorization technique
	HypreBoomerAMG,
	HypreILU,
	ANY
};

std::string getName(const precond_name a_name)
{
  switch(a_name) 
  {
	  case precond_name::HypreBoomerAMG: return "HypreBoomerAMG";
	  case precond_name::HypreILU: return "HypreILU";
	  case precond_name::HypreDiagScale: return "HypreDiagScale"; // JACOBI
	  case precond_name::HypreParaSails: return "HypreParaSails";
	  case precond_name::HypreEuclid: return "HypreEuclid";
	  case precond_name::ANY: return "ANY";
	  default: return "ANY";
  }
}

mfem_mgis::Parameters buildPrecond(const precond_name a_name, const int a_verbosity)
{
	return mfem_mgis::Parameters{{"Name", getName(a_name)},
                                {"Options",
                                        mfem_mgis::Parameters{
                                                {"VerbosityLevel", a_verbosity}
                                        }
                                }
                        };
}
