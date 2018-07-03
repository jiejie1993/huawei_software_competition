/*************************************************************************
    > File Name: distserver.h
    > Author: geeker
    > Mail: 932834897@qq.com 
    > Created Time: Tue 20 Mar 2018 06:07:39 AM PDT
 ************************************************************************/

#include<iostream>
using namespace std;





#include <string>
#include <vector>

#define CPU true
#define MEM false




//Flavor类，描述虚拟机信息：名称 内存 CPU
struct Flavor {
	int id;
	int mem;
	int cpu;
	Flavor(int _id,int _mem, int _cpu):
		id(_id),mem(_mem),cpu(_cpu){}
	Flavor() {

	}
};

//Server类，描述物理服务器信息：总内存 总CPU 可用内存 可用CPU 已存放的虚拟机列表
class Server {
	public:
		std::vector<Flavor> flavors;  //物理服务器已存放虚拟机列表

		vector<int> flavors_vec;

		///Server构造函数，参数为内存大小和CPU大小
		Server(int mem, int cpu);
		///放置虚拟机函数，参数为虚拟机对象，返回值为是否放置成功
		bool put_flavor(Flavor flavor);
		///获取服务器CPU使用率
		double get_cpu_usage_rate();
		///获取服务器内存使用率
		double get_mem_usage_rate();

		double get_my_usage_rate(vector<int>&);

	private:
		int total_mem;  // 物理服务器总内存
		int total_cpu;  // 物理服务器总CPU
		int free_mem;   // 物理服务器剩余可用内存
		int free_cpu;   // 物理服务器剩余可用CPU
};


///Server构造函数，参数为内存大小和CPU大小
Server::Server(int mem, int cpu) :total_mem(mem), total_cpu(cpu) {
	free_cpu = cpu;  //初始化时剩余CPU等于总CPU
	free_mem = mem;  //初始化时剩余内存等于总内存
	flavors_vec = vector<int>(25, 0);
}
///放置虚拟机函数，参数为虚拟机对象，返回值为是否放置成功
///首先检查剩余CPU和内存是否足够放置该虚拟机
///如果能够放下虚拟机，则将虚拟机放入服务器，并更新服务器可用内存和可用CPU，并返回true
///如果剩余内存和CPU不足以放下该虚拟机，则返回false
bool Server::put_flavor(Flavor flavor) {
	if (free_cpu >= flavor.cpu && free_mem >= flavor.mem) {
		free_cpu -= flavor.cpu;
		free_mem -= flavor.mem;
		flavors.push_back(flavor);
		flavors_vec[flavor.id]++;
		return true;
	}
	return false;
}

///获取服务器CPU使用率
double Server::get_cpu_usage_rate() {
	return 1 - free_cpu / static_cast<double>(total_cpu);
}
///获取服务器内存使用率
double Server::get_mem_usage_rate() {
	return 1 - free_mem / static_cast<double>(total_mem);
}

double Server::get_my_usage_rate(vector<int>& pred_vec) {
	//Tools::print_vec(flavors_vec);
	double sum = 0.0;
	int count = 0;
	for(int i = 0; i < pred_vec.size(); i++){
		if(pred_vec[i] >= 0){
			if(pred_vec[i] != 0)
				sum += (double)flavors_vec[i] / pred_vec[i];
			count++;
		}
	}
	//cout << "count" << count << endl;
	//cout << "count" << count << endl;
	return sum/count-0.0001;
}


vector<int> unfold(vector<int>& pred_vec, vector<vector<int> >& Flavor_Info){

	vector<int> res;
	for(int i = 0; i < Flavor_Info.size(); i++){
		int flavId = Flavor_Info[i][0];		
		int flavnum = pred_vec[flavId];
		for(int j = 0; j < flavnum; j++)
			res.push_back(flavId);
	}
	return res;
}

vector<Server> put_flavors_to_servers(vector<int>& pred_vec, vector<vector<int> >& Flavor_Info, int& Max_Server_Cpu, int& Max_Server_Mem, bool CPUorMEM){
	
	Tools::print_vec(pred_vec);


	vector<int> unfold_vec = unfold(pred_vec, Flavor_Info);
	
	vector<Flavor> vec_flavors(unfold_vec.size());
	
	for(int i = 0; i < vec_flavors.size(); i++){
		int temp_flavor_cpu, temp_flavor_mem;
		for(int j = 0; j < Flavor_Info.size(); j++){
			if(Flavor_Info[j][0] == unfold_vec[i]){
				temp_flavor_cpu = Flavor_Info[j][2];
				temp_flavor_mem = Flavor_Info[j][1];
			}
		}

		Flavor temp_flavor(unfold_vec[i], temp_flavor_cpu, temp_flavor_mem);
		vec_flavors[i] = temp_flavor;
	}
	
	//cout << "vec_flavors.size==>" << vec_flavors.size() << endl;

	//for(int k = 0; k < vec_flavors.size(); k++){
	//	cout << vec_flavors[k].id << ", "
	//		<< vec_flavors[k].cpu << ", "
	//		<< vec_flavors[k].mem << endl;
	
	//}

	cout << endl;

	int server_mem = Max_Server_Mem * 1024;
	int server_cpu = Max_Server_Cpu;

	//=========================================================================
	//模拟退火算法找最优解
	double min_server = vec_flavors.size() + 1;
	vector<Server> res_servers;  //用于存放最好结果（服务器使用数量最少）
	double T = 100.0;  //模拟退火初始温度
	double Tmin = 1;   //模拟退火终止温度
	double r = 0.9999; //温度下降系数
	vector<int> dice;  //骰子，每次随机投掷，取vector前两个变量作为每次退火需要交换顺序的虚拟机
	for (int i = 0; i < vec_flavors.size(); i++) {
		dice.push_back(i);
	}
	double server_num;
	while (T > Tmin) {
		//cout << T << endl;
		//投掷骰子，如vector前两个数为3和9，则把vec_flavors[3]和vec_flavors[9]进行交换作为新的flavors顺序
		std::random_shuffle(dice.begin(), dice.end());
		auto new_vec_flavors = vec_flavors;
		std::swap(new_vec_flavors[dice[0]], new_vec_flavors[dice[1]]);

		//把上一步计算出来的虚拟机尝试加入到服务器中

		//先使用一个服务器用于放置虚拟机
		vector<Server> servers;
		Server firstserver(server_mem, server_cpu);
		servers.push_back(firstserver);  

		//放置虚拟机主要逻辑
		//如果当前所有服务器都放不下虚拟机，就新建一个服务器用于存放
		for (auto element : new_vec_flavors) {
			auto iter = servers.begin();
			for (; iter != servers.end(); ++iter) {
				if (iter->put_flavor(element)) {
					break;
				}
			}
			if (iter == servers.end()) {
				Server newserver(server_mem, server_cpu);
				newserver.put_flavor(element);
				servers.push_back(newserver);
			}
		}

		//cout << "servers.size->" << servers.size() << endl;

		//计算本次放置虚拟机耗费服务器评价分数(double型)
		//如果使用了N个服务器，则前N-1个服务器贡献分数为1，第N个服务器分数为资源利用率
		//模拟退火就是得到取得分数最小的放置方式
		//double server_num;

		//对于题目关心CPU还是MEM，需要分开讨论，资源利用率计算方法不同
		/*if (CPUorMEM == CPU)
			server_num = servers.size() - 1 + servers.rbegin()->get_cpu_usage_rate();
		else
			server_num = servers.size() - 1 + servers.rbegin()->get_mem_usage_rate();*/

		server_num = servers.size()-1 + servers.rbegin()->get_my_usage_rate(pred_vec);
		
		//server_num = servers.size() - 1 + servers.rbegin()->get_cpu_usage_rate();
	
		//cout << server_num << endl;

		//如果分数更低，则保存结果
		if (server_num <= min_server) {
			//cout << "dsafdsa" << server_num << endl;
			min_server = server_num;
			res_servers = servers;
			vec_flavors = new_vec_flavors;
		}
		//如果分数更高，则以一定概率保存结果，防止优化陷入局部最优解
		else {
			if (exp((min_server - server_num) / T) > rand() / RAND_MAX) {
				//cout << "dsafdsafdafds" << server_num << endl;
				min_server = server_num;
				res_servers = servers;
				vec_flavors = new_vec_flavors;
			}
		}
		T = r * T;  //一次循环结束，温度降低
	}

	if(server_num - (int)server_num < 0.1 && res_servers.size() > 1){
		//cout << "nani???" << endl;
		res_servers.pop_back();
	}

	//Tools::print_vec(res_servers[0].flavors_vec);

	

	//cout << "begin" << endl;
	
	auto iter = res_servers.begin();
	for(; iter != res_servers.end(); ++iter){
		bool flag = true;
		while(flag){
			flag = false;
			for(int i = Flavor_Info.size()-1; i >= 0; i--){
				/*
				for (; iter != servers.end(); ++iter) {
					if (iter->put_flavor(element)) {
						break;
					}
				}*/

				if(iter->put_flavor(Flavor( Flavor_Info[i][0], Flavor_Info[i][2], Flavor_Info[i][1] ))){
				
					//cout << "hello" << i << endl;	
					flag = true;
				}
			}
		}

	}

	//cout << server_num << endl;
	return res_servers;

}

void dist_server(vector<vector<int> > &distri_plan, vector<int> &pred_vec, int &Max_Server_Cpu, int &Max_Server_Mem, int &Flavor_Num, vector<vector<int> > &Flavor_Info, bool& opt_target){

	//vector<int> unfold_vec = unfold(pred_vec, Flavor_Info);

	//Tools::print_vec(pred_vec);

	vector<Server> server_vec = put_flavors_to_servers(pred_vec, Flavor_Info, Max_Server_Cpu, Max_Server_Mem, opt_target);

	/*for(auto element : server_vec){
		vector<int> temp_dist(25, 0);
		for(auto flavor : element.flavors){
			temp_dist[flavor.id]++;
		}
		distri_plan.push_back(temp_dist);
	}*/

	for(auto element : server_vec){
		distri_plan.push_back(element.flavors_vec);
	}

	if(distri_plan == vector<vector<int> >{})
		distri_plan.push_back(vector<int>(25, 0));
	
	for(int i = 0; i < pred_vec.size(); i++){
		if(pred_vec[i] >= 0){
			int sum = 0;
			for(int j = 0; j < distri_plan.size(); j++){
				sum += distri_plan[j][i];
			}
			pred_vec[i] = sum;
		}
	}

}


/***************************************************************************************
Function: dist_server
Description: Based on predicted data, the virtual machine is allocated to the server
			to achive the highest utilization of resources
Input: distri_plan //a two-dimention array with each row representing a physical server
		pred_vec //prediction array
		Max_Server_Cpu
		Max_Server_Mem
		Flavor_Num
		Flavor_Info
Output: void
****************************************************************************************/


void dist_server1(vector<vector<int>> &distri_plan, vector<int> &pred_vec, int &Max_Server_Cpu, int &Max_Server_Mem, int &Flavor_Num, vector<vector<int>> &Flavor_Info){

	vector<int> my_server(25, 0);
	vector<int> curr_server = my_server;
	int server_cpu = Max_Server_Cpu, server_mem = Max_Server_Mem * 1024;

	for(int i = 0; i < Flavor_Num; i++){
		int curr_flavor = Flavor_Info[i][0];

		int curr_num = pred_vec[curr_flavor];
		int flavor_cpu = Flavor_Info[i][1];
		int flavor_mem = Flavor_Info[i][2];
		while(1){
			if(server_cpu < flavor_cpu || server_mem < flavor_mem){
				distri_plan.push_back(curr_server);
				curr_server = my_server;
				server_cpu = Max_Server_Cpu; 
				server_mem = Max_Server_Mem * 1024;
				continue;
			}
			if(curr_num == 0)
				break;
			server_cpu -= flavor_cpu;
			server_mem -= flavor_mem;
			curr_num--;
			curr_server[curr_flavor]++;
		}
	}
	if(curr_server != my_server || distri_plan == vector<vector<int>>{})
		distri_plan.push_back(curr_server);
		
}
