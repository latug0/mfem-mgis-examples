#pragma once

#include <common/mfem_mgis_headers.hxx>
#include <common/data_gathering.hxx>

namespace common
{
	template<typename Implementation>
	void print_mesh_information(Implementation& impl)
	{
		START_TIMER("common::print_mesh_information");
		
		//getMesh
		auto mesh = impl.getFiniteElementSpace().GetMesh();
		//get the number of vertices
		int numbers_of_vertices = mesh->GetNV();
		//get the number of elements
		int numbers_of_elements = mesh->GetNE();
		//get the element size
		double h = mesh->GetElementSize(0);

		using profiling::output::printMessage;

		printMessage("INFO: number of vertices -> ", numbers_of_vertices);
		printMessage("INFO: number of elements -> ", numbers_of_elements);
		printMessage("INFO: element size -> ", h);
	}

	template<typename Problem, typename T>
		double solve(Problem& p, double start, double end, T& statistics, std::string string_solver="", std::string string_pc="")
		{
			START_TIMER("common::solve");
			// solving the problem
			double res = profiling::timers::chrono_section( [&](){
					statistics = p.solve(0, 1);
					});

			// check status
			if (!statistics.status) {
				profiling::output::printMessage("INFO: ", string_solver,"+",string_pc," FAILED");
			}
			return res;
		}

	template<typename Problem>
		double add_post_processings(Problem& p, std::string msg)
		{
	//		START_TIMER("common::add_postprocessing_and_outputs");
			using profiling::output::printMessage;
			printMessage("before AddPost");
			p.addPostProcessing(
					"ParaviewExportResults",
					{{"OutputFileName", msg}}
					);
	//		printMessage("after AddPost");
		} // end timer add_postprocessing_and_outputs

	template<typename Problem>
		double execute_post_processings(Problem& p, double start, double end)
		{
			START_TIMER("common::post_processing_step");
			p.executePostProcessings(start, end);
		}

	template<typename Solver, typename Pc, typename Statistics>
		void fill_statistics(gather_information& data, Solver solver, Pc pc, Statistics stats, double time)
		{
			START_TIMER("common::fill_statistics");
			time = profiling::output::reduce_max(time);
			// fill info data.
			data.add(
					info{	solver,		pc,	
					stats.status,	stats.iterations,	
					stats.final_residual_norm,	time
					}
				);
		}

	void print_statistics(std::string string_solver, std::string string_pc, double time);
};
