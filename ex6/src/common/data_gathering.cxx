#include <common/data_gathering.hxx>

using namespace configuration;

gather_information read(const std::string a_input_file_name, bool& a_exist)
{
	CatchTimeSection("gather_information::read");

	assert(a_input_file_name != "" && "this file does not exist");

	gather_information res;
	std::ifstream input_file (a_input_file_name, std::ifstream::in);
	std::string line;


	size_t number_of_lines = 0;
	while(std::getline(input_file, line))
	{
		info data;
		std::istringstream iss(line);
		int solv_id=-1;
		int precond_id=-1;
		if(!(iss>> solv_id 
					>> precond_id 
					>> data.m_converged 
					>> data.m_iteration 
					>> data.m_residu 
					>> data.m_time))
		{
			std::abort();
		}

		data.m_solver = solver_name(solv_id);
		data.m_precond = precond_name(precond_id);
		res.add(data);
		number_of_lines++;
	}

	if(number_of_lines != 0)
	{	
		a_exist = true;
		Profiler::Utils::Message("file ",a_input_file_name,"has been correctly read");
	}
	else 
	{
		a_exist = false;
		Profiler::Utils::Message("file ",a_input_file_name," is empty");
	}
	return res;
}

