#pragma once

#include<iostream>
#include<chrono>
#include<vector>
#include<omp.h>
#include <fstream>
#include <algorithm>
#include<cassert>

// variables 
#ifdef __MPI
const size_t cWidth =20;
const size_t nColumns=6;
const std::string cName[nColumns]={"number Of Calls","min(s)", "mean(s)", "max(s)" ,"part(%)", "imb(%)"}; // [1-Imax/Imean]% 
#else
const size_t cWidth =20;
const size_t nColumns=3;
const std::string cName[nColumns]={"number Of Calls", "max(s)," ,"part(%)"};
#endif

enum enumTimer
{
	CURRENT,
	ROOT
};

namespace profiling
{
	namespace output
	{
		bool is_master();
		double reduce_max(double a_duration);

		template<typename Arg>
		void printMessage(Arg a_msg)
		{
			if(is_master())
			{
				std::cout << a_msg << std::endl;
			}
		}

		template<typename Arg, typename... Args>
		void printMessage(Arg a_msg, Args... a_msgs)
		{
			if(is_master())
			{
				std::cout << a_msg << " ";
				printMessage(a_msgs...);
			}
		}
	};

	namespace timers
	{
		class profilingTimer
		{
			using duration = std::chrono::duration<double>;

			public:

			profilingTimer();
			profilingTimer(std::string name, profilingTimer* mother);
			profilingTimer* find(std::string name);

			// printer functions
			//
			void printReplicate(size_t begin, size_t end, std::string motif);
			void space();
			void column();
			void endline();
			void printBanner(size_t shift);
			void printEnding(size_t shift);
			duration* get_ptr_duration();
			void print(size_t shift, double runtime);

			// accessor
			//
			std::string getName();
			std::size_t get_iteration();
			std::size_t get_level();
			std::vector<profilingTimer*>& get_daughter();
			profilingTimer* get_mother();
			double get_duration();
		

			private:

			std::string m_name;
			std::size_t m_iteration;
			std::size_t m_level;
			std::vector<profilingTimer*> m_daughter;
			profilingTimer* m_mother;
			duration m_duration;
		};

		void init_timers();
		void print_and_write_timers();

		template<enumTimer T>
		profilingTimer*& get_timer()
		{
			static profilingTimer* __current;
			return __current;
		}

		template<typename Lambda>
		double chrono_section(Lambda&& lambda)
		{
			using duration = std::chrono::duration<double>;
			using steady_clock = std::chrono::steady_clock;
			using time_point = std::chrono::time_point<steady_clock>;
			time_point tic, toc;
			tic = steady_clock::now();
			lambda();
			toc = steady_clock::now();
			auto measure = toc - tic;
			return measure.count();	
		}
	};

	namespace timer
	{
		using duration = std::chrono::duration<double>;
		using steady_clock = std::chrono::steady_clock;
		using time_point = std::chrono::time_point<steady_clock>;
		class Timer
		{
			public:
			Timer(duration * acc);
			void start();
			void end();
			~Timer(); 

			private:
			time_point m_start;
			time_point m_stop;
			duration * m_duration; 
		};

		template<enumTimer T>
		Timer*& get_timer()
		{
			static Timer* __timer;
			return __timer;
		}

		template<enumTimer T>
		void start_global_timer()
		{
			assert(T == enumTimer::ROOT);
			auto& timer = get_timer<T>(); 
			timer = new Timer(profiling::timers::get_timer<T>()->get_ptr_duration());
			timer->start(); // reset start
		}

		template<enumTimer T>
		void end_global_timer()
		{
			assert(T == enumTimer::ROOT);
			auto timer = get_timer<T>();
			timer->end();
		}
	}

	namespace outputManager
	{
		using profiling::timers::profilingTimer;

		std::string build_name();
		void printTimeTable();
		void writeFile();
		void writeFile(std::string a_name);

		template<typename Func, typename... Args>
		void recursive_call(Func& func, profilingTimer* ptr, Args&... arg)
		{
			func(ptr, arg...);
			auto& daughters = ptr->get_daughter();
			for(auto& it: daughters)
				recursive_call(func,it,arg...);
		}

		template<typename Func, typename Sort, typename... Args>
		void recursive_sorted_call(Func& func, Sort mySort, profilingTimer* ptr, Args&... arg)
		{
			func(ptr, arg...);
			auto& daughters = ptr->get_daughter();
			std::sort(daughters.begin(),daughters.end(),mySort);
			for(auto& it:daughters)
				recursive_sorted_call(func, mySort, it,arg...);
		}


	};
};


#define START_TIMER(X) profiling::timers::profilingTimer*& current = profiling::timers::get_timer<CURRENT>();\
	assert(current != nullptr && "do not use an undefined profilingTimer"); \
	current = current->find(X); \
        profiling::timer::Timer tttt(current->get_ptr_duration());

