#pragma once

#include<solver/precond_name.hxx>
#include<solver/solver_name.hxx>
#include<MFEMMGIS/Profiler.hxx>

using namespace configuration;

struct info
{
	solver_name 	m_solver;
	precond_name 	m_precond;
	bool 		m_converged;
	int 		m_iteration;
	double 		m_residu;
	double 		m_time;
};

class gather_information
{
	public:
		gather_information() : m_data() {};
		~gather_information() {}

		void add(info& a_info)
		{
			CatchTimeSection("gather_information::add&");
			m_data.push_back(std::move(a_info));
		}

		void add(info&& a_info)
		{
			CatchTimeSection("gather_information::add");
			m_data.push_back(a_info);
		}

		std::string build_name()
		{
			std::string base_name = "collect";
			int size;
			MPI_Comm_size(MPI_COMM_WORLD,&size);
			std::string name = base_name + "_" + std::to_string(size);
			return name;	
		}

		void write()
		{
			CatchTimeSection("gather_information::write");

			if(mfem_mgis::Profiler::Utils::is_master())
			{
				auto name = build_name();
				std::ofstream file (name, std::ofstream::out);
				for(auto it : m_data)
				{
					file << it.m_solver << " "
						<< it.m_precond << " "
						<< it.m_converged << " "
						<< it.m_iteration << " "
						<< it.m_residu << " "
						<< it.m_time << std::endl;	
				}
			}	
		}

		std::vector<info>& get_data()
		{
			std::vector<info>& ref = m_data;
			return ref;
		}

		void insert(gather_information& a_in)
		{
			auto& in_data = a_in.get_data();
			m_data.resize(m_data.size() + in_data.size());
			std::move(in_data.begin(), in_data.end(), std::back_inserter(m_data));
		}

		void writeMD(std::string a_name, bool add_banner = true)
		{
			CatchTimeSection("gather_information::writeMD");
			if(mfem_mgis::Profiler::Utils::is_master())
			{
				std::ofstream file (a_name, std::ofstream::in | std::ofstream::out | std::ofstream::ate);
				if(add_banner)
				{
					file << "| solver | preconditionner | converged | iterations | residu | time |" << std::endl;
				}

				for(auto it : m_data)
				{
					std::string solv = getName(it.m_solver);
					std::string prec = getName(it.m_precond);
					std::string conv = it.m_converged ? "&#10004;" : "&cross;";
					std::string ite  = it.m_converged ? std::to_string(it.m_iteration) : "&cross;";

					std::ostringstream streamRes; // manage precision
					streamRes << it.m_residu;

					std::ostringstream streamTime; // manage precision
					streamTime << it.m_time;

					std::string res  = it.m_converged ? streamRes.str()  : "&cross;";
					std::string time = it.m_converged ? streamTime.str() : "&cross;";
					std::string line = "| " + solv + " | " + prec + " | " + conv + " | " + ite + " | " + res + " | " + time + " | ";
					file << line << std::endl;
				}
			}	
		}

		void writeMD()
		{
			auto name = "markdown_" + build_name();
			writeMD(name);
		}

		void print()
		{
			CatchTimeSection("gather_information::print");
			mfem_mgis::Profiler::Utils::Message("| solver | ", " preconditionner |"," converged |"," iterations |"," residu |", " time |");
			// improve readability
			std::string old_solv =" ";
			std::string old_prec =" ";


			for(auto it : m_data)
			{
				std::string solv = getName(it.m_solver);

				if(solv == old_solv)
				{
					solv = " ";
				}
				else
				{
					old_solv = solv;
				}

				std::string prec = getName(it.m_precond);
				if(prec == old_prec)
				{
					prec = " ";
				}
				else
				{
					old_prec = prec;
				} 	

				std::string conv = it.m_converged ? "&#10004;" : "&cross;";
				std::string ite  = it.m_converged ? std::to_string(it.m_iteration) : "&cross;";

				std::ostringstream streamRes; // manage precision
				streamRes << it.m_residu;

				std::ostringstream streamTime; // manage precision
				streamRes << it.m_time;

				std::string res  = it.m_converged ? streamRes.str()  : "&cross;";
				std::string time = it.m_converged ? streamTime.str() : "&cross;";


				std::string line = "| " + getName(it.m_solver) + " | " + prec + " | " + conv + " | " + ite + " | " + res + " | " + time + " | ";
				mfem_mgis::Profiler::Utils::Message(line);
			}
		}

		void reset()
		{
			m_data.clear();
		}

		size_t size()
		{
			const size_t _size = m_data.size();
			return _size;
		}

	private:
		std::vector<info> m_data;
};

gather_information read(const std::string a_input_file_name, bool& a_exist);
