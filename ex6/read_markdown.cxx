#include<timer.hpp>
#include<data_gathering.hxx>



void writeMD(std::vector<gather_information> a_vec, std::string a_name)
{
	bool print_banner = true;
	std::ofstream tmp;
	tmp.open(a_name, std::ofstream::out);
	tmp << " markdown - mfem-mgis - analysis" << std::endl;
	tmp.close();

	for(auto it : a_vec)
	{
		it.writeMD(a_name, print_banner);
		print_banner = false;
	}
	profiling::output::printMessage("write: ", a_name);
}


int main(int argc, char* argv[]) 
{
	mfem_mgis::initialize(argc, argv);
	profiling::timers::init_timers();
	const int first = 128;
	const int last = 512;
	const std::string base_name = "collect_";
	//const int last = 4096;

	std::vector<gather_information> all; // easier to manage it with a vector


	// get data
	for(unsigned int nb_proc = first ; nb_proc <= last ; nb_proc*=2)
	{
		const std::string name = base_name + std::to_string(nb_proc);
		auto sub_all = read(name);
		//all.insert(sub_all); 
		all.push_back(std::move(sub_all));
	}

	assert(all.size()>0);
	const auto number_of_simu = all.size();
	const auto number_of_item = all[0].size(); // same size for each elem

	std::vector<std::pair<int,info>> try_classification;

	for(int item_id = 0 ; item_id < number_of_item ; item_id++)
	{
		int proc_id = first;
		for(int simu_id = 0 ; simu_id < number_of_simu ; simu_id++)
		{
			auto& local_cont_data = all[simu_id];
			auto& local_data = local_cont_data.get_data();
			auto& local_info = local_data[item_id]; 

			auto elem = std::make_pair(proc_id,local_info);
			try_classification.push_back(elem);
			proc_id *= 2;
		}
	}	

	//writeMD(all, "all.md");

	// other
	std::sort (try_classification.begin(), try_classification.end(), [](const std::pair<int,info> a, const std::pair<int,info> b)
			{
				if(a.second.m_solver < b.second.m_solver) 
				{
					return true;
				}
				else if (a.second.m_solver == b.second.m_solver)
				{
					if(a.second.m_precond < b.second.m_precond)
					{
						return true;
					}
					else if(a.second.m_precond == b.second.m_precond)
					{
						if(a.first < b.first) return true;
						else return false;
					}
					else
					{
						return false;
					}
				}
				else 
				{
					return false;
				}
			});


	std::ofstream out("all.md", std::ofstream::out);

	out << "| solver | preconditionner | number of mpi | converged | iterations | residu | time |" << std::endl;
	out << "| :-------------|:---------------|:---------------:|:---------------:|:---------------:|:---------------:|---------------:|" << std::endl;
	for(auto it: try_classification)
	{
		info& inf = it.second;
		std::string solv = getName(inf.m_solver);
		std::string prec = getName(inf.m_precond);
		std::string conv = inf.m_converged ? "&#10004;" : "&cross;";
		std::string ite  = inf.m_converged ? std::to_string(inf.m_iteration) : "&cross;";

		std::ostringstream streamRes; // manage precision
		streamRes << inf.m_residu;

		std::ostringstream streamTime; // manage precision
		streamTime << inf.m_time;

		std::string res  = inf.m_converged ? streamRes.str()  : "&cross;";
		std::string time = inf.m_converged ? streamTime.str() : "&cross;";
		std::string proc = std::to_string(it.first);
		std::string line = "| " + solv + " | " + prec + " | " + proc + " | " + conv + " | " + ite + " | " + res + " | " + time + " | ";
		out << line << std::endl;
	}	

	profiling::timers::print_and_write_timers();
	return EXIT_SUCCESS;
}
