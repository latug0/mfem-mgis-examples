#include "precond_name.hxx"

namespace configuration
{
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
	
	std::string getName(petsc_pc_type in)
	{
		switch(in)
		{
			case petsc_pc_type::jacobi : return "jacobi";
			case petsc_pc_type::bjacobi : return "bjacobi";
			case petsc_pc_type::sor : return "sor";
			case petsc_pc_type::eisenstat : return "eisenstat";
			case petsc_pc_type::icc : return "icc";
			case petsc_pc_type::ilu : return "ilu";
			case petsc_pc_type::gasm : return "gasm";
			case petsc_pc_type::gamg : return "gamg";
			case petsc_pc_type::bddc : return "bddc";
			case petsc_pc_type::ksp : return "ksp";
			case petsc_pc_type::lu : return "lu";
			case petsc_pc_type::cholesky : return "cholesky";
			case petsc_pc_type::none : return "none";
			default :return "none";
		}
	}

	mfem_mgis::Parameters buildPrecond(const precond_name a_name, const int a_verbosity)
	{
		auto options = mfem_mgis::Parameters{{"VerbosityLevel", a_verbosity}};
		/*	
			if(a_name == precond_name::HypreILU)
			{
			options.insert(	mfem_mgis::Parameters{
			{"LevelOfFill", 2}
			});
			}*/
		return mfem_mgis::Parameters{{"Name", getName(a_name)}, {"Options",options}};
	}
};
