#pragma once

enum solver_name
{
	Default,
	HypreFGMRES,
	HypreGMRES,
	HyprePCG,
	MUMPSSolver,
	GMRESSolver,
	CGSolver,
	BiCGSTABSolver,
	UMFPackSolver,
	MINRESSolver,
	SLISolver
};

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
	  case solver_name::GMRESSolver:	return "GMRESSolver";
	  case solver_name::CGSolver:		return "CGSolver";
#ifdef MFEM_USE_SUITESPARSE
	  case solver_name::UMFPackSolver:	return "UMFPackSolver";
#endif
#ifdef MFEM_USE_MUMPS
	  case solver_name::MUMPSSolver:	return "MUMPSSolver";
#endif
	  case solver_name::BiCGSTABSolver:	return "BiCGSTABSolver";
	  case solver_name::HypreFGMRES:	return "HypreFGMRES";
	  case solver_name::HypreGMRES:		return "HypreGMRES";
	  case solver_name::HyprePCG:		return "HyprePCG";
	  case solver_name::MINRESSolver:	return "MINRESSolver";
	  case solver_name::SLISolver: 		return "SLISolver"; 
	  case solver_name::Default: 		return "Default";
	  default: 				return "error";
  } 
}


