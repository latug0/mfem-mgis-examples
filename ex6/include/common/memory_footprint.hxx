#pragma once

#include <sys/time.h>
#include <sys/resource.h>
#include <MFEMMGIS/Profiler.hxx>

namespace common
{
	namespace memory
	{
		rusage make_memory_checkpoint()
		{
			rusage obj;
			int who;
			auto res = getrusage(who, &obj);
			assert(res = 0 && "error: getrusage has failed");
			return obj;
		};

		class footprint
		{
			public:
				footprint() : m_data() {}

				void add_memory_checkpoint()
				{
					auto obj = make_memory_checkpoint();
					m_data.push_back(obj);
				}

				using max_use=long;
				std::vector<max_use> reduce()
				{
					// data are reduce on the mpi process 0
					std::vector<max_use> obj;
					size_t nb_points = m_data.size();
					obj.resize(nb_points);

					for(int id = 0 ; id < nb_points ; id++)
					{
						obj[id] = 0;
						MPI_Reduce(&(m_data[id].ru_maxrss), &(obj[id]), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
					}
					return obj;
				}

			private:
				std::vector<rusage> m_data;
		};

		void print_checkpoints(footprint& f)
		{
			auto obj = f.reduce();
			if(mfem_mgis::Profiler::Utils::is_master())
			{
				std::cout << " List (maximum resident size): ";
				for(auto it : obj)
				{
					std::cout << " " << it;
				}
				std::cout << std::endl;
			}
		};

		void print_memory_footprint()
		{
			footprint f;
			f.add_memory_checkpoint();
			auto obj = f.reduce();
			auto last = obj.back() * 1e-6; // conversion kb to Gb
			mfem_mgis::Profiler::Utils::Message(" memory footprint: ", last, " GB");
		};
	};
};
