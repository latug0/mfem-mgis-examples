#include<timer.hpp>
#include<data_gathering.hxx>

int main()
{

	const int first = 1;
	const int last = 4096;

	std::vector<data_gathering> all; // easier to manage it with a vector

	// get data
	for(unsigned int nb_proc = first ; nb_proc <= last ; nb_proc)
	{
		const std::string name = bas_name + std::to_string(nb_proc);
		auto sub_all = read(name);
		//all.insert(sub_all); 
		all.push_back(std::move(sub_all));
	}

	assert(all.size()>0);
	const auto number_of_simu = all.size()
	const auto number_of_item = all[0].size(); // same size for each elem

	for(int item_id = 0 ; item_id < number_of_item ; item_id++)
	{
		for(int simu_id = 0 ; data_id < number_of_simu ; data_id++)
		{
			const auto proc_id = std::pow(2, first + simu_id);
			auto& local_cont_data = all[simu_id];
			auto& local_data = local_cont_data.get_data();
			auto& local_info = local_data[item_id]; 
			const auto get_markdown_line();
		}
	}	

	all.writeMD("all.md")


	return 0;
}
