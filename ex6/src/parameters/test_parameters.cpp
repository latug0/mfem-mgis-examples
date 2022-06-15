#include <parameters/test_parameters.hpp>
#include <test-cases/cas_cible_1.hxx>
#include <test-cases/fissuration.hxx>

void common_parameters(mfem::OptionsParser& args, TestParameters& p)
{
	args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
	args.AddOption(&p.library, "-l", "--library", "Material library.");
	args.AddOption(&p.order, "-o", "--order", "Finite element order (polynomial degree).");
	args.AddOption(&p.tcase, "-t", "--test-case", "identifier of the case : cas_cible_1, fissuration : default = cas_cible_1");
	args.AddOption(&p.refinement, "-r", "--refinement", "refinement level of the mesh, default = 0");
	args.AddOption(&p.post_processing, "-p", "--post-processing", "run post processing step");
	
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

markdown_reader_parameters markdown_reader_parameters_with_parse(int& argc, char* argv[])
{
	START_TIMER("reader_parameters_with_parse");

	mfem::OptionsParser args(argc, argv);

	// default value are start = 1 to last = 1024 
	markdown_reader_parameters obj;

	// options
	args.AddOption(&obj.start, "-s", "--start", "number of mpi processes of the first collect file");
	args.AddOption(&obj.last, "-l", "--last", "number of mpi processes of the last collect file");

	args.Parse();
	if (!args.Good()) {
		if (mfem_mgis::getMPIrank() == 0)
			args.PrintUsage(std::cout);
		mfem_mgis::finalize();
		exit(0);
	}
	return obj;
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




