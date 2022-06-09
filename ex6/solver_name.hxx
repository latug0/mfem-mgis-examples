#pragma once

#include<string>

namespace configuration
{

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
		SLISolver,
		solver_name_count
	};


	enum class petsc_ksp_type
	{
		richardson, 
		chebyshev, 
		cg, 
		gmres, 
		tcqmr, 
		bcgs, 
		cgs, 
		tfqmr, 
		cr, 
		lsqr, 
		bicg, 
		preonly,
		petsc_ksp_type_count
	};

	bool isHypre(solver_name name);
	bool isMumps(solver_name name);

	std::string getName(const solver_name a_name);
	std::string getName(const petsc_ksp_type in);
};
