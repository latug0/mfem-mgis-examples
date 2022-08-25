#include <common/common.hxx>
namespace common 
{
	void print_statistics(bool converged, std::string string_solver, std::string string_pc, double time)
	{
		if(converged)
		{
			Profiler::Utils::Message(
					"INFO: ",
					string_solver,
					"+",
					string_pc,
					" works correctly, elapsed time: ",
					time
					);
		}
		else
		{
			Profiler::Utils::Message(
					"if(s == solver_name::", 
					string_solver, 
					" && p == precond_name::", 
					string_pc,
					") return false;"
					);
		}
		Profiler::Utils::Message("");
	}
}
