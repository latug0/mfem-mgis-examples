#include<string>
#include<iostream>
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#ifdef MFEM_USE_PETSC
#include "mfem/linalg/petsc.hpp"
#endif /* MFEM_USE_PETSC */

#include<timer.hpp>
#include<solver_name.hxx>
#include<precond_name.hxx>

// make a petsc file
//
using namespace configuration;

void create_petsc_file(const petsc_ksp_type a_solv_name, const petsc_pc_type a_pc_name)
{
	START_TIMER("petsc_stuff::create_petsc_file");
	// options used
	std::string solver 	= "-ksp_type ";
	std::string pc 		= "-pc_type ";
	std::string rtol 	= "-ksp_rtol ";
	std::string atol 	= "-ksp_atol ";
	std::string max_it 	= "-ksp_max_it ";
	std::string name 	= "petscrc_file_" + getName(a_solv_name) + "_" + getName(a_pc_name);

	solver 	+= getName(a_solv_name);
	pc 	+= getName(a_pc_name);
	rtol 	+= "1e-10";
	atol 	+= "1e-10";
	max_it 	+= "10e1";

	profiling::output::printMessage(" creating ... ", name);
	std::ofstream p(name, std::ofstream::out);

	p 	<< 	solver 	<< std::endl;
	p 	<< 	pc 	<< std::endl;
	p 	<< 	rtol 	<< std::endl;
	p 	<< 	atol 	<< std::endl;
	p 	<< 	max_it 	<< std::endl;
		
	p.close();
}

void generate_petsc_files()
{
	START_TIMER("petsc_stuff::generate_petsc_files");
	size_t petsc_ksp_type_size = (size_t)petsc_ksp_type::petsc_ksp_type_count;
	size_t petsc_pc_type_size  = (size_t)petsc_pc_type::petsc_pc_type_count;
	profiling::output::printMessage("number of solver: ", petsc_ksp_type_size);
	profiling::output::printMessage("number of preconditioner: ", petsc_pc_type_size);
	for(int ksp_type = 0 ; ksp_type < petsc_ksp_type_size ; ksp_type++)
	{
		for(int pc_type = 0; pc_type < petsc_pc_type_size ; pc_type++)
		{
			create_petsc_file(
					(petsc_ksp_type)ksp_type, 
					(petsc_pc_type)pc_type
					);
		}
	}
}


int main(int argc, char* argv[]) 
{
	mfem_mgis::initialize(argc, argv);
	profiling::timers::init_timers();

	generate_petsc_files();

	profiling::timers::print_and_write_timers();
	return EXIT_SUCCESS;
}












