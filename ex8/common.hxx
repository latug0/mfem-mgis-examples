#pragma once
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/ImposedDirichletBoundaryConditionAtClosestNode.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

/**
 * @brief Parses common command-line parameters and populates the TestParameters structure.
 *
 * This function parses common command-line parameters using the provided `args` object
 * and populates the `TestParameters` structure `p` with the parsed values. The common
 * parameters include options for configuring various aspects of the test.
 *
 * @param args An instance of the `mfem::OptionsParser` class for parsing command-line arguments.
 * @param p    A reference to the `TestParameters` structure where parsed values will be stored.
 *
 * @note The `TestParameters` structure should be defined prior to calling this function,
 *       as it will be modified to contain the parsed values.
 *
 * @see TestParameters
 * @see mfem::OptionsParser
 */
void common_parameters(mfem::OptionsParser& args, TestParameters& p)
{
	args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
	args.AddOption(&p.library, "-l", "--library", "Material library.");
	args.AddOption(&p.order, "-o", "--order", "Finite element order (polynomial degree).");
	args.AddOption(&p.refinement, "-r", "--refinement", "refinement level of the mesh, default = 0");
	args.AddOption(&p.post_processing, "-p", "--post-processing", "run post processing step");
	args.AddOption(&p.verbosity_level, "-v", "--verbosity-level", "choose the verbosity level");

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

/**
 * @brief Adds post-processing options to a given Problem instance.
 *
 * This template function allows you to add post-processing options
 * to a Problem instance of the specified type. Post-processing can include tasks such
 * as generating and saving visualizations, stress analysis, and more.
 *
 * @tparam Problem              The type of the Problem instance to which post-processing will be added.
 *
 * @param p                     A reference to the Problem instance to which post-processing will be added.
 * @param msg                   Paraview file name (results).
 * @param stress_output_name    The name or identifier for stress-related output data.
 * @param paraview_output_name  Paraview file with integration point value.
 *
 * @note This function allows you to customize and configure various post-processing options
 *       specific to the given Problem type. Please refer to the documentation of the Problem
 *       class for details on available post-processing options and their meanings.
 */
	template<typename Problem>
void add_post_processings(Problem& p, std::string msg, std::string stress_output_name, std::string paraview_output_name)
{
	p.addPostProcessing(
			"ParaviewExportResults",
			{{"OutputFileName", msg}}
			);
	p.addPostProcessing(
			"MeanThermodynamicForces",
			{{"OutputFileName", stress_output_name}});
	p.addPostProcessing(
			"ParaviewExportIntegrationPointResultsAtNodes",
			{
			{"OutputFileName", paraview_output_name},
			{"Results", {"Stress"}}
			}
			);
}

/**
 * @brief Adds post-processing options to a given Problem instance.
 *
 * This template function allows you to add post-processing options
 * to a Problem instance of the specified type. Post-processing can include tasks such
 * as generating and saving visualizations, stress analysis, and more.
 *
 * @tparam Problem              The type of the Problem instance to which post-processing will be added.
 *
 * @param p                     A reference to the Problem instance to which post-processing will be added.
 * @param msg                   Paraview file name (results).
 * @param stress_output_name    The name or identifier for stress-related output data.
 * @param paraview_output_name  Paraview file with integration point value.
 *
 * @note This function allows you to customize and configure various post-processing options
 *       specific to the given Problem type. Please refer to the documentation of the Problem
 *       class for details on available post-processing options and their meanings.
 */
	template<typename Problem>
void add_post_processings_2(Problem& p, std::string msg, std::string stress_output_name)
{
	p.addPostProcessing(
			"ParaviewExportResults",
			{{"OutputFileName", msg}}
			);
	p.addPostProcessing(
			"MeanThermodynamicForces",
			{{"OutputFileName", stress_output_name}});
	p.addPostProcessing(
			"ParaviewExportIntegrationPointResultsAtNodes",
			{{"OutputFileName", msg + "IntegrationPointOutput"},
			{"Results", {"FirstPiolaKirchhoffStress"}}});
} // end timer add_postprocessing_and_outputs
/**
 * @brief Executes post-processing steps for a given Problem instance within a specified time range.
 *
 * @tparam Problem  The type of the Problem instance for which post-processing will be executed.
 *
 * @param p         A reference to the Problem instance for which post-processing will be executed.
 * @param start     The starting time of the post-processing range.
 * @param end       The ending time of the post-processing range.
 *
 * @note This function is responsible for executing the defined post-processing steps within the specified time range.
 *       The specific post-processing steps and their configurations should have been previously added to the Problem
 *       instance using the `add_post_processings` function or other appropriate methods.
 */
	template<typename Problem>
void execute_post_processings(Problem& p, double start, double end)
{
	CatchTimeSection("common::post_processing_step");
	p.executePostProcessings(start, end);
}


/**
 * @brief Solves a Problem instance within a specified time range and time step.
 * @param p         A reference to the Problem instance to be solved.
 * @param start     The starting time for the solution.
 * @param dt        The time step for the solution.
 *
 * @note This function is responsible for solving the Problem instance within the specified
 *       time range using the provided time step. It also checks the status of the solver
 *       and aborts if the solve operation fails.
 */
	template<typename Problem>
void run_solve(Problem& p, double start, double dt)
{
	CatchTimeSection("Solve");
	mfem_mgis::Profiler::Utils::Message("Solving problem from: ",start," to: ",start+dt);
	// solving the problem
	auto statistics = p.solve(start, dt);
	// check status
	if (!statistics.status) {
		mfem_mgis::Profiler::Utils::Message("Info solver: solving step has failed");
		std::abort();
	}
}


	template<typename Implementation>
void print_mesh_information(Implementation& impl)
{
	CatchTimeSection("common::print_mesh_information");

	using mfem_mgis::Profiler::Utils::sum;
	using mfem_mgis::Profiler::Utils::Message;

	//getMesh
	auto mesh = impl.getFiniteElementSpace().GetMesh();

	//get the number of vertices
	int64_t numbers_of_vertices_local = mesh->GetNV();
	int64_t  numbers_of_vertices = sum(numbers_of_vertices_local);

	//get the number of elements
	int64_t numbers_of_elements_local = mesh->GetNE();
	int64_t numbers_of_elements = sum(numbers_of_elements_local);

	//get the element size
	double h = mesh->GetElementSize(0);

	// get n dofs
	auto& fespace = impl.getFiniteElementSpace();
	int64_t unknowns_local = fespace.GetTrueVSize(); 
	int64_t unknowns = sum(unknowns_local);

	Message("Info problem: number of vertices -> ", numbers_of_vertices);
	Message("Info problem: number of elements -> ", numbers_of_elements);
	Message("Info prolbem: element size -> ", h);
	Message("Info porblem: Number of finite element unknowns: " , unknowns);
}
