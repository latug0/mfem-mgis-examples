#include<common/timer.hpp>
#ifdef __MPI
#include"mpi.h"
#endif
#include<numeric>

namespace profiling
{
	namespace output
	{
		constexpr int master=0;

		bool is_master()
		{
#ifdef __MPI
			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			return (rank == master);
#else
			return true;
#endif
		}

		double reduce_max(double a_duration)
		{
#ifdef __MPI
			int size = -1;

			double global = 0.0;
			MPI_Reduce(&a_duration, &global, 1, MPI_DOUBLE, MPI_MAX, master, MPI_COMM_WORLD); // master rank is 0
			return global;
#else
			return a_duration;
#endif
		}

	};

	namespace timers
	{
		using duration = std::chrono::duration<double>;

		// constructor
		profilingTimer::profilingTimer() : m_daughter() // only used for root
		{
			m_name 		= "root";
			m_iteration 	= 1;
			m_level 	= 0;
			m_mother 	= nullptr;

		}

		profilingTimer::profilingTimer(std::string name, profilingTimer* mother): m_daughter(), m_duration(0)
		{
			m_name 		= name;
			m_iteration 	= 1;
			m_level 	= mother->m_level + 1;
			m_mother 	= mother;
		}

		profilingTimer* 
			profilingTimer::find(std::string name)
			{
				assert(this != nullptr);
				for(auto it = m_daughter.begin() ; it < m_daughter.end() ; it++)
				{
					if((*it)->m_name == name)
					{
						(*it)->m_iteration++;
						return (*it);
					}
				}
				profilingTimer* myTmp = new profilingTimer(name, this);
				m_daughter.push_back(myTmp);
				return myTmp;
			}

		void 
			profilingTimer::printReplicate(size_t begin, size_t end, std::string motif)
			{
				for(size_t i = begin ; i < end ; i++) std::cout << motif;
			}


		void 
			profilingTimer::space()
			{
				std::cout << " "; 
			}

		void 
			profilingTimer::column()
			{
				std::cout << "|"; 
			}

		void 
			profilingTimer::endline()
			{
				std::cout << std::endl; 
			}

		void 
			profilingTimer::printBanner(size_t shift)
			{
				if(m_name == "root")
				{
#ifndef __MPI
					profiling::output::printMessage(" MPI feature is disable for timers, if you use MPI please add -D__MPI ");
#else
					if(profiling::output::is_master()) {
						profiling::output::printMessage(" MPI feature activated, rank 0:");
#endif
						std::string start_name = " |-- start timetable "; 
						std::cout << start_name;
						size_t end = shift+ nColumns*(cWidth+1) + 1;	
						printReplicate(start_name.size(), end,"-");
						column(); endline();
						std::string name = " |    name";
						std::cout << name;
						printReplicate(name.size(),shift + 1," ");
						for(size_t i =  0 ; i < nColumns ; i++)
						{
							column();
							int size = cName[i].size();
							printReplicate(0,(int(cWidth)-size - 1), " ");
							std::cout << cName[i];
							space();
						}
						column(); endline();
						space(); column();
						printReplicate(2, end,"-");
						column(); endline();
#ifdef __MPI
					}
#endif
				}
			}

		void 
			profilingTimer::printEnding(size_t shift)
			{
				if(m_name == "root")
				{
					if(profiling::output::is_master()) 
					{
						shift+= nColumns*(cWidth+1) + 1; // +1 for "|";
						std::string end_name = " |-- end timetable " ;
						std::cout << end_name;
						printReplicate(end_name.size(),shift,"-");
						column(); endline();
					}
				}
			}

		duration* 
			profilingTimer::get_ptr_duration()
			{
				return &m_duration;
			}

		void 
			profilingTimer::print(size_t shift, double total_time)
			{
				assert(total_time >= 0);
				std::string cValue[nColumns];
				if(profiling::output::is_master()) 
				{
					size_t realShift = shift;
					space(); column(); space();
					size_t currentShift = 3;
					for(int i = 0 ; i < int(m_level) - 1; i++) 
					{
						int spaceSize = 3;
						for(int j = 0 ; j < spaceSize ; j++) space();
						currentShift += spaceSize;
					}
					if(m_level>0) {
						std::cout << "|--";
						currentShift += 3;
					}
					std::cout << "> "<< m_name;
					currentShift += m_name.size() + 1;
					printReplicate(currentShift, realShift, " ");

					cValue[0] = std::to_string(m_iteration);	
				}
#ifdef __MPI
				double local  = m_duration.count();	
				int size = -1;
				MPI_Comm_size(MPI_COMM_WORLD, &size);

				assert(size > 0);
				std::vector<double> list;

				if(profiling::output::is_master()) list.resize(size);

				MPI_Gather(&local,1,MPI_DOUBLE, list.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // master rank is 0


				if(profiling::output::is_master())
				{
					const auto [min,max]	= std::minmax_element(list.begin(), list.end());
					auto global_max 	= *max;
					auto global_min 	= *min; 
					auto sum 		= std::accumulate(list.begin(), list.end(), double(0.));
					auto global_mean 	= sum / double(size);
					auto part_time 		= (global_max / total_time ) * 100;

					assert(global_mean >= 0);
					assert(global_min >= 0);
					assert(global_max >= 0);

					assert(global_max >= global_mean);
					assert(global_mean >= global_min);

					cValue[1] = std::to_string( global_min);	
					cValue[2] = std::to_string( global_mean);	
					cValue[3] = std::to_string( global_max);	
					cValue[4] = std::to_string( part_time) + "%";	
					cValue[5] = std::to_string( (global_max/global_mean)-1) + "%";
				}
#else
				cValue[1] = std::to_string( m_duration.count());	
				cValue[2] = std::to_string( (m_duration.count()/total_time)*100 );	
#endif
				if(profiling::output::is_master())
				{
					for(size_t i =  0 ; i < nColumns ; i++)
					{
						column();
						int size = cValue[i].size();
						printReplicate(0,(int(cWidth)-size - 1), " ");
						std::cout << cValue[i];
						space();
					}
					column();endline();
				}
			}

		std::string 
			profilingTimer::getName()
			{
				return m_name;
			}

		double 
			profilingTimer::get_duration()
			{
				return m_duration.count();
			}

		std::size_t 
			profilingTimer::get_iteration()
			{
				return m_iteration;
			}

		std::size_t 
			profilingTimer::get_level()
			{
				return m_level;
			}

		std::vector<profilingTimer*>& 
			profilingTimer::get_daughter()
			{
				return m_daughter;
			}

		profilingTimer* 
			profilingTimer::get_mother()
			{
				return m_mother;
			}

		void init_timers()
		{
			profiling::output::printMessage(" Init timers ");
			profilingTimer*& root_timer_ptr 	= profiling::timers::get_timer<ROOT>() ;
			assert(root_timer_ptr == nullptr);	
			root_timer_ptr 			= new profilingTimer(); 
			profilingTimer*& current 	= profiling::timers::get_timer<CURRENT>(); 
			current 			= root_timer_ptr;
			assert(current != nullptr);	
			profiling::timer::start_global_timer<ROOT>();
		}

		void print_and_write_timers()
		{
			profiling::timer::end_global_timer<ROOT>(); 
			profiling::outputManager::writeFile(); 
			profiling::outputManager::printTimeTable();
		}


	};

	namespace timer
	{
		Timer::Timer(duration * acc) 
		{
			m_duration = acc;
			start();
		}

		inline 
			void Timer::start()
			{
				m_start = steady_clock::now();
			}

		inline 
			void Timer::end()
			{
				assert(m_duration != nullptr && "duration has to be initialised");
				m_stop = steady_clock::now();
				*m_duration += m_stop - m_start;
				assert(m_duration->count() >= 0);
			}

		Timer::~Timer() 
		{
			end();
			auto& current_timer = profiling::timers::get_timer<CURRENT>();
			current_timer = current_timer->get_mother();
		}
	};

	namespace outputManager
	{
		std::string build_name()
		{
			std::size_t nthreads=0;
#pragma omp parallel
			{
				nthreads = omp_get_num_threads();
			}
			std::string base_name = "mfem-mgis";
#ifdef __MPI
			int mpiSize;
			MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);
			std::string file_name = base_name + "." + std::to_string(mpiSize) + ".perf";
#else
			std::string file_name = base_name + "." + std::to_string(nthreads) + ".perf";
#endif
			return file_name;
		}

		void printTimeTable()
		{
			profilingTimer* root_timer = profiling::timers::get_timer<ROOT>();
			double runtime = root_timer->get_duration();
			runtime = profiling::output::reduce_max(runtime); // if MPI, else return runtime

			auto my_print = [](profilingTimer* a_ptr, size_t a_shift, double a_runtime)
			{
				a_ptr->print(a_shift, a_runtime);
			};

			auto sort_comp = [] (profilingTimer* a_ptr, profilingTimer* b_ptr)
			{
				return a_ptr->get_duration() > b_ptr->get_duration() ;
			};


			auto max_length = [](profilingTimer* a_ptr, size_t& a_count, size_t& a_nbElem)
			{
				size_t length = a_ptr->get_level()*3 + a_ptr->getName().size();
				a_count = std::max(a_count, length);
				a_nbElem++;
			};
			size_t count(0), nbElem(0);

			recursive_call(max_length, root_timer, count, nbElem);
			count += 6;
			root_timer->printBanner(count);
			recursive_sorted_call(my_print, sort_comp, root_timer, count, runtime);
			root_timer->printEnding(count);
		}

		void writeFile()
		{
			std::string name = build_name();
			writeFile(name);
		}

		void writeFile(std::string a_name)
		{
			using namespace profiling::output;
			//using profiling::output::reduce_max;

			std::ofstream myFile (a_name, std::ofstream::out);	
			profilingTimer* root_timer = profiling::timers::get_timer<ROOT>();
			auto rootTime = root_timer->get_duration();
			rootTime = profiling::output::reduce_max(rootTime);
			auto my_write = [rootTime](profilingTimer* a_ptr, std::ofstream& a_file)
			{
				std::string space;
				std::string motif = "   ";

				for(std::size_t i = 0 ; i < a_ptr->get_level() ; i++) space +=motif;

				const auto max_time = reduce_max(a_ptr->get_duration());

				if(is_master())
				{
					a_file << space << a_ptr->getName() 
						<< " " << a_ptr->get_iteration()
						<< " " << max_time
						<< " " <<(max_time/rootTime)*100
						<< std::endl;
				}
			};

			recursive_call(my_write,root_timer,myFile);
		}
	}
};

