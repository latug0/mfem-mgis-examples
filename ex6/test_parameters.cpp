#include "test_parameters.hpp"
#include <cas_cible_1.hxx>
#include <fissuration.hxx>

void common_parameters(mfem::OptionsParser& args, TestParameters& p)
{
	args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
	args.AddOption(&p.library, "-l", "--library", "Material library.");
	args.AddOption(&p.order, "-o", "--order",
			"Finite element order (polynomial degree).");
	args.AddOption(&p.tcase, "-t", "--test-case",
			"identifier of the case : cas_cible_1, fissuration : default = cas_cible_1");
	args.AddOption(&p.refinement, "-r", "--refinement",
			"refinement level of the mesh, default = 0");

	args.Parse();
	
	if (!args.Good()) {
		if (mfem_mgis::getMPIrank() == 0)
			args.PrintUsage(std::cout);
		mfem_mgis::finalize();
		exit(0);
	}
	if (p.mesh_file == nullptr) {
		if (mfem_mgis::getMPIrank() == 0)
			std::cout << "ERROR: Mesh file missing" << std::endl;
		args.PrintUsage(std::cout);
		mfem_mgis::abort(EXIT_FAILURE);
	}
	if (mfem_mgis::getMPIrank() == 0)
		args.PrintOptions(std::cout);
}


TestParameters parseCommandLineOptions(int& argc, char* argv[]) {
	START_TIMER("parse_command_line_options");
	TestParameters p;

	// options treatment
	mfem::OptionsParser args(argc, argv);
	common_parameters(args, p);
	
	if ((p.tcase < 0) || (p.tcase > 5)) {
		std::cerr << "Invalid test case\n";
		mfem_mgis::abort(EXIT_FAILURE);
	}
	return p;
}

// legacy
TestParameters parseCommandLineOptions_one_test(int& argc, char* argv[]) {
	START_TIMER("parse_command_line_options");
	TestParameters p;

	// options treatment
	mfem::OptionsParser args(argc, argv);

	args.AddOption(
			&p.linearsolver, "-ls", "--linearsolver",
			"identifier of the linear solver: 0 -> Default, 1 -> HypreFGMRES, 2 -> HypreGMRES, 3-> HyprePCG, 4 -> MUMPSSolver, 5 -> GMRESSolver, 6 -> CGSolver, 7 -> BiCGSTABSolver, 8 -> UMFPackSolver, 9 ->	MINRESSolver, 10 -> SLISolver");

	args.AddOption(
			&p.preconditioner, "-lp", "--preconditioner",
			" 0 -> HypreDiagScale, 1 -> HypreParaSails, 2 -> HypreEuclid, 3 -> HypreBoomerAMG, 4 -> HypreILU, 5 -> ANY");

	common_parameters(args, p);
	if (!args.Good()) {
		if (mfem_mgis::getMPIrank() == 0)
			args.PrintUsage(std::cout);
		mfem_mgis::finalize();
		exit(0);
	}

	if ((p.tcase < 0) || (p.tcase > 5)) {
		std::cerr << "Invalid test case\n";
		mfem_mgis::abort(EXIT_FAILURE);
	}
	return p;
}

std::vector<solver_name> get_solvers(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::build_solvers_list();
		case 2: return fissuration::build_solvers_list();
		default: return cas_cible_1::build_solvers_list();
	}
}

std::vector<precond_name> get_pc(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::build_pc_list();
		case 2: return fissuration::build_pc_list();
		default: return cas_cible_1::build_pc_list();
	}
}


std::function<bool(solver_name,precond_name)> get_match(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::match;
		case 2: return fissuration::match;
		default: return cas_cible_1::match;
	}
}

std::vector<petsc_ksp_type> get_solvers_with_petsc(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::build_solvers_list_with_petsc();
		case 2: return fissuration::build_solvers_list_with_petsc();
		default: return cas_cible_1::build_solvers_list_with_petsc();
	}
}

std::vector<petsc_pc_type> get_pc_with_petsc(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::build_pc_list_with_petsc();
		case 2: return fissuration::build_pc_list_with_petsc();
		default: return cas_cible_1::build_pc_list_with_petsc();
	}
}

std::function<bool(petsc_ksp_type,petsc_pc_type)> get_match_with_petsc(int a_case)
{
	switch(a_case)
	{
		case 1: return cas_cible_1::match_with_petsc;
		case 2: return fissuration::match_with_petsc;
		default: return cas_cible_1::match_with_petsc;
	}
}
