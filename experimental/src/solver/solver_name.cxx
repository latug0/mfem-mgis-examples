#include <solver/solver_name.hxx>

namespace configuration
{
	bool isHypre(solver_name name)
	{
		if(name == HypreFGMRES || name == HypreGMRES || name == HyprePCG) return true;
		else return false;
	}

	bool isMumps(solver_name name)
	{
		if(name == MUMPSSolver) return true;
		else return false;
	}

	std::string getName(const solver_name a_name)
	{
		switch(a_name) 
		{
			case solver_name::GMRESSolver:		return "GMRESSolver";
			case solver_name::CGSolver:		return "CGSolver";
#ifdef MFEM_USE_SUITESPARSE
			case solver_name::UMFPackSolver:	return "UMFPackSolver";
#endif
#ifdef MFEM_USE_MUMPS
			case solver_name::MUMPSSolver:		return "MUMPSSolver";
#endif
			case solver_name::BiCGSTABSolver:	return "BiCGSTABSolver";
			case solver_name::HypreFGMRES:		return "HypreFGMRES";
			case solver_name::HypreGMRES:		return "HypreGMRES";
			case solver_name::HyprePCG:		return "HyprePCG";
			case solver_name::MINRESSolver:		return "MINRESSolver";
			case solver_name::SLISolver: 		return "SLISolver"; 
			case solver_name::Default: 		return "Default";
			default: 				return "error";
		} 
	}

	std::string getName(const petsc_ksp_type in)
	{
		switch(in)
		{
			case petsc_ksp_type::richardson : return "richardson";
			case petsc_ksp_type::chebyshev 	: return "chebyshev";
			case petsc_ksp_type::cg 	: return "cg";
			case petsc_ksp_type::gmres	: return "gmres";
			case petsc_ksp_type::tcqmr 	: return "tcqmr";
			case petsc_ksp_type::bcgs 	: return "bcgs";
			case petsc_ksp_type::cgs 	: return "cgs";
			case petsc_ksp_type::tfqmr 	: return "tfqmr";
			case petsc_ksp_type::cr 	: return "cr";
			case petsc_ksp_type::lsqr 	: return "lsqr";
			case petsc_ksp_type::bicg 	: return "bicg";
			case petsc_ksp_type::preonly 	: return "preonly";
			default 			: return "error";
		}
	}
};
