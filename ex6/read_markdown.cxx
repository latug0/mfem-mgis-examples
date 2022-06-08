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
using simu_val  = std::tuple<size_t, double, bool>;
using list_val  = std::vector<simu_val>;

class store_speed_up_values : public std::map<simu_name, list_val>
{
	public:

	store_speed_up_values() {}
	
	void fill(const solver_name a_solv, const precond_name a_precond, const size_t a_proc, const double a_time, const bool a_conv)
	{
		simu_name name {a_solv,a_precond};
		simu_val value{a_proc,a_time,a_conv};
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

	std::string get_line(const solver_name a_solv, const precond_name a_precond, list_val& a_list)
	{
		// get the first value
		double coeff = 0.0;
		// a_list is sorted by construction
		for(auto it : a_list)
		{
			if(std::get<2>(it))
			{
				coeff = std::get<0>(it) * std::get<1>(it);
				break;
			}
		}
		
		std::string line = getName(a_solv) + "-" + getName(a_precond);
		
		for (auto val : a_list)
		{
			if(std::get<2>(val)) // converged
			{
				auto speed_up = coeff / std::get<1>(val);
				line += " (" + std::to_string(std::get<0>(val)) + " " +  std::to_string(speed_up) + ")";
			}
		}
		
		return line;
	}
	
	std::string get_label(const solver_name a_solv, const precond_name a_precond)
	{
		std::string label = "label=\"" + getName(a_solv) + "-" + getName(a_precond) + "\"";
		return label;
	}
	
	std::string get_procs(const solver_name a_solv, const precond_name a_precond, list_val& a_list)
	{
		std::string line = "np.array([";
		bool first = true;
		for (auto val : a_list)
		{
			if(std::get<2>(val)) // converged
			{
				if(!first) line += ",";
				line += std::to_string(std::get<0>(val));
				first = false;
			}
		}
		
		line += "])";

		return line;
	}
	
	std::string get_values(const solver_name a_solv, const precond_name a_precond, list_val a_list)
	{
		// get the first value
		double coeff = 0.0;
		// a_list is sorted by construction
		for(auto it : a_list)
		{
			if(std::get<2>(it))
			{
				coeff = std::get<0>(it) * std::get<1>(it);
				break;
			}
		}
		
		std::string line = "np.array([";
		bool first = true;
		for (auto val : a_list)
		{
			if(std::get<2>(val)) // converged
			{
				if(!first) line += ",";
				auto speed_up = coeff / std::get<1>(val) ;
				line += std::to_string(speed_up);
				first = false;
			}
		}
		
		line += "])";

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

void build_speed_up_for_python_plot(store_speed_up_values& a_in)
{
	std::cout << "import numpy as np" << std::endl;
	std::cout << "import matplotlib.pyplot as plt" << std::endl;
	std::cout << "from matplotlib.ticker import FormatStrFormatter" << std::endl;
	std::cout << "import matplotlib.ticker as ticker" << std::endl;
	std::cout << "def my_function():" << std::endl;	
	for(auto it : a_in)
	{
		auto key = it.first;
		auto values = it.second;

		// test if at least one simulation has converged 
		size_t do_print = 0; // we need at least two points
		for(auto it : values)
		{
			if(std::get<2>(it))
			{
				do_print++;
				if(do_print >= 2)
				{
					break;
				}
			}
		}
		
		if(do_print<2) continue;
		
		std::string line = "	plt.plot(";
		line += a_in.get_procs(key.first, key.second, values);
		line += ",";
		line += a_in.get_values(key.first, key.second, values);
		line += ",";
		line += a_in.get_label(key.first, key.second);
		line += ")";
		std::cout << line << std::endl;
	}
	std::cout << "	pass" 	<< std::endl<< std::endl;
	std::cout << "my_function()"			<< std::endl;
	std::cout << "legend_outside = plt.legend(loc=\'center left\',bbox_to_anchor=(1, 0.5))" << std::endl;
	std::cout << "namepdf=\'basic-version.pdf\'" 	<< std::endl;
	std::cout << "plt.savefig(namepdf, bbox_extra_artists=(legend_outside,), bbox_inches=\'tight\')" << std::endl;
	std::cout << "fig, ax = plt.subplots()" 	<< std::endl;
	std::cout << "plt.xscale(\"log\",basex=2)" 	<< std::endl;
	std::cout << "plt.yscale(\"log\",basey=2)" 	<< std::endl;
	std::cout << "ax.xaxis.set_major_formatter(ticker.ScalarFormatter())" 	<< std::endl;
	std::cout << "ax.yaxis.set_major_formatter(ticker.ScalarFormatter())" 	<< std::endl;
	std::cout << "plt.ticklabel_format(axis=\'x\', style=\'plain\')" 	<< std::endl;
	std::cout << "plt.ticklabel_format(axis=\'y\', style=\'plain\')" 	<< std::endl;
	std::cout << "my_function()"			<< std::endl;
	std::cout << "legend_outside = plt.legend(loc=\'center left\',bbox_to_anchor=(1, 0.5))" << std::endl; 
	std::cout << "namepdf=\'log-version.pdf\'" 	<< std::endl;
	std::cout << "plt.savefig(namepdf, bbox_extra_artists=(legend_outside,), bbox_inches=\'tight\')" << std::endl;

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
	std::vector<std::pair<int,info>> try_classification;

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
				storage.fill(it.m_solver, it.m_precond, nb_proc, it.m_time, it.m_converged);
			}
		}
		//all.push_back(std::move(sub_all));
		if(exist)
		{
			auto& local_data = sub_all.get_data();
			for(auto it : local_data)
			{
				auto elem = std::make_pair(nb_proc, it);
				try_classification.push_back(elem);
			}
		}
	}
	
	//build_speed_up_for_python_partial(storage);
	build_speed_up_for_python_plot(storage);


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
