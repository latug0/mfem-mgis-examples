#include "test_parameters.hpp"

TestParameters parseCommandLineOptions(int& argc, char* argv[]) {
	START_TIMER("parse_command_line_options");
	TestParameters p;

	// options treatment
	mfem::OptionsParser args(argc, argv);
	args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
	args.AddOption(&p.library, "-l", "--library", "Material library.");
	args.AddOption(&p.order, "-o", "--order",
			"Finite element order (polynomial degree).");
	args.AddOption(&p.tcase, "-t", "--test-case",
			"identifier of the case : cas_cible_1, fissuration");
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
	if ((p.tcase < 0) || (p.tcase > 5)) {
		std::cerr << "Invalid test case\n";
		mfem_mgis::abort(EXIT_FAILURE);
	}
	return p;
}
