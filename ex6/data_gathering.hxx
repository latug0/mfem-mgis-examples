#pragma once


struct info
{
	solver_name 	m_solver;
	precond_name 	m_precond;
	bool 		m_converged;
	int 		m_iteration;
	double 		m_residu;
};

class gather_information
{
	public:
	gather_information() : m_data() {};
	~gather_information() {}
	
	void add(info&& a_info)
	{
		m_data.push_back(std::move(a_info));
	}
	
	void print()
	{
		profiling::output::printMessage("solver", "preconditionner","converged","iterations","residu");
		for(auto it : m_data)
		{
			profiling::output::printMessage(getName(it.m_solver),
					getName(it.m_precond),
					it.m_converged,
					it.m_iteration,
					it.m_residu);
		}
	}

	void reset()
	{
		m_data.clear();
	}
	private:
	std::vector<info> m_data;
};


