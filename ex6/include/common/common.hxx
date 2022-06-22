#pragma once

#include <common/mfem_mgis_headers.hxx>
#include <common/data_gathering.hxx>

namespace common
{
	template<typename Implementation>
	void print_mesh_information(Implementation& impl)
	{
		CatchTimeSection("common::print_mesh_information");
		
		using Profiler::Utils::sum;
		using Profiler::Utils::Message;
		
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

		Message("INFO: number of vertices -> ", numbers_of_vertices);
		Message("INFO: number of elements -> ", numbers_of_elements);
		Message("INFO: element size -> ", h);
		Message("INFO: Number of finite element unknowns: " , unknowns);
	}

	template<typename Problem, typename T>
		double solve(Problem& p, double start, double end, T& statistics, std::string string_solver="", std::string string_pc="")
		{
			CatchTimeSection("common::solve");
			// solving the problem
			double res = Profiler::timers::chrono_section( [&](){
					statistics = p.solve(0, 1);
					});

			// check status
			if (!statistics.status) {
				Profiler::Utils::Message("INFO: ", string_solver,"+",string_pc," FAILED");
			}
			return res;
		}

	template<typename Problem>
		void add_post_processings(Problem& p, std::string msg)
		{
	//		CatchTimeSection("common::add_postprocessing_and_outputs");
			using Profiler::Utils::Message;
			Message("before AddPost");
			p.addPostProcessing(
					"ParaviewExportResults",
					{{"OutputFileName", msg}}
					);
			Message("after AddPost");
		} // end timer add_postprocessing_and_outputs

	template<typename Problem>
		void execute_post_processings(Problem& p, double start, double end)
		{
			CatchTimeSection("common::post_processing_step");
			p.executePostProcessings(start, end);
		}

	template<typename Solver, typename Pc, typename Statistics>
		void fill_statistics(gather_information& data, Solver solver, Pc pc, Statistics stats, double time)
		{
			CatchTimeSection("common::fill_statistics");
			time = Profiler::Utils::reduce_max(time);
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
