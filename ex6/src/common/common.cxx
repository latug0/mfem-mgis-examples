#include <common/common.hxx>
namespace common 
{
	void print_statistics(std::string string_solver, std::string string_pc, double time)
	{
		Profiler::Utils::Message(
				"INFO: ",
				string_solver,
				"+",
				string_pc,
				" works correctly, elapsed time: ",
				time
				);

		Profiler::Utils::Message("");
	}
}
