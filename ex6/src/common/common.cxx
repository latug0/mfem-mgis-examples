#include <common/common.hxx>
namespace common 
{
	void print_statistics(std::string string_solver, std::string string_pc, double time)
	{
		profiling::output::printMessage(
				"INFO: ",
				string_solver,
				"+",
				string_pc,
				" works correctly, elapsed time: ",
				time
				);

		profiling::output::printMessage("");
	}
}
