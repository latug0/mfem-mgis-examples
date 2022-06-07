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

using simu_name = std::pair<solver_name,precond_name>;
using simu_val  = std::pair<size_t, double>;
using list_val  = std::vector<simu_val>;

class store_speed_up_values : public std::map<simu_name, list_val>
{
	public:

	store_speed_up_values() {}
	
	void fill(const solver_name a_solv, const precond_name a_precond, const size_t a_proc, const double a_time)
	{
		simu_name name {a_solv,a_precond};
		simu_val value{a_proc,a_time};
		auto elem = this->find(name);
		if(elem != this->end())
		{
			auto& list = elem->second;
			list.push_back(value);	
		}
		else
		{
			list_val list;
			list.push_back(value);
			(*this)[name] = list;
		}
	}

	std::string get_line(const solver_name a_solv, const precond_name a_precond, list_val a_list)
	{
		std::string line = getName(a_solv) + "-" + getName(a_precond);
		for (auto val : a_list) line += " (" + std::to_string(val.first) + " " +  std::to_string(val.second) + ")";
		return line;
	}
};


void build_speed_up_for_python_partial(store_speed_up_values& a_in)
{
	for(auto it : a_in)
	{
		auto key = it.first;
		auto values = it.second;
		auto line = a_in.get_line(key.first, key.second, values);
		std::cout << line << std::endl;
	}
}

template<typename T>
void build_speed_up_for_python_full(const T a_in, const int a_first, const int a_last)
{
	START_TIMER("build_speed_up_for_python");
	const std::string name = "speedup.log";
	std::ofstream out(name, std::ofstream::out);

	std::string header = "solver preconditionner";
		
	
	int nb_elem = 0;
	for(int i=a_first; i<= a_last ; i*=2) 
	{
		header += " " + std::to_string(i);
		nb_elem++;
	}	

	out << header << std::endl;
	std::cout << header << std::endl;
	double res[nb_elem]; 
	auto iterator = a_in.begin();

	while (iterator != a_in.end())
	{
		auto solv = getName(iterator->second.m_solver);
		auto prec = getName(iterator->second.m_precond);
		std::string baseline = solv + "-" + prec;
		for(int i = 0 ; i < nb_elem ; i++)
		{
			if(solv != getName(iterator->second.m_solver)) std::cout << "solver is different" << std::endl;
			if(prec != getName(iterator->second.m_precond)) std::cout << "preconditionner is different" << std::endl;
			res[i] = iterator->second.m_time;
			//std::next(iterator);
			iterator++;
		}
		auto base = res[0];
                for(int i = 0 ; i < nb_elem ; i++)
                {
			assert(res[i]>0);
                        res[i] = base / res[i];
                }

                for(int i = 0 ; i < nb_elem ; i++)
		{
			baseline += " " + std::to_string(res[i]);
		}
//		std::cout << baseline << std::endl;
		out << baseline << std::endl;
	}


}

int main(int argc, char* argv[]) 
{
	mfem_mgis::initialize(argc, argv);
	profiling::timers::init_timers();
	const int first = 1;
	const int last = 16384;
	const std::string base_name = "collect_";
	//const int last = 4096;

	std::vector<gather_information> all; // easier to manage it with a vector
	std::vector<bool> if_file; 

	store_speed_up_values storage; // just a map[sovler,precond]={(proc1,val1), .... , (procn,valn)}

	// get data
	for(unsigned int nb_proc = first ; nb_proc <= last ; nb_proc*=2)
	{
		const std::string name = base_name + std::to_string(nb_proc);
		bool exist = false;
		auto sub_all = read(name, exist);
		//all.insert(sub_all); 
		if(exist)
		{
			for(auto it : sub_all.get_data())
			{
				storage.fill(it.m_solver, it.m_precond, nb_proc, it.m_time);
			}
		}
		if_file.push_back(exist);
		all.push_back(std::move(sub_all));
	}
	build_speed_up_for_python_partial(storage);

	assert(all.size()>0);
	const auto number_of_simu = all.size();
	const auto number_of_item = all[0].size(); // same size for each elem

	std::vector<std::pair<int,info>> try_classification;

	for(int item_id = 0 ; item_id < number_of_item ; item_id++)
	{
		int proc_id = first;
		for(int simu_id = 0 ; simu_id < number_of_simu ; simu_id++)
		{
			if(if_file[simu_id] == true)
			{
				auto& local_cont_data = all[simu_id];
				auto& local_data = local_cont_data.get_data();
				auto& local_info = local_data[item_id]; 

				auto elem = std::make_pair(proc_id,local_info);
				try_classification.push_back(elem);
			}
			proc_id *= 2;
		}
	}	

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
	
	//build_speed_up_for_python_full(try_classification, first, last);

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
