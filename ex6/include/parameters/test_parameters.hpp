#pragma once
#include <common/mfem_mgis_headers.hxx>
#include <common/timer.hpp>
#include <solver/solver_name.hxx>
#include <common/common.hxx>
#include <vector>

using namespace configuration;

// We need this class for test case sources
constexpr double xmax = 1.;
struct TestParameters {
	const char* mesh_file = "cas-cible-1.mesh";
	const char* behaviour = "Elasticity";
	const char* library = "src/libBehaviour.so";
	const char* reference_file = "Elasticity.ref";
	int order = 1;
	int tcase = 1;
	double xmax = 1.;
	double ymax = 1.;
	double zmax = 1.;
	int linearsolver = 1;
	int preconditioner = 1;
	bool parallel = true;
	int refinement = 0;
	int post_processing = 0; // default value : disabled
	int verbosity_level = 0; // default value : lower level
};

TestParameters parseCommandLineOptions(int& argc, char* argv[]);
TestParameters parseCommandLineOptions_one_test(int& argc, char* argv[]);


struct markdown_reader_parameters
{
	int start = 1;
	int last = 1024;
};
markdown_reader_parameters markdown_reader_parameters_with_parse(int& argc, char* argv[]);

