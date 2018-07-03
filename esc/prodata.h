/*************************************************************************
    > File Name: prodata.h
    > Author: geeker
    > Mail: 932834897@qq.com 
    > Created Time: Thu 15 Mar 2018 07:09:28 AM PDT
 ************************************************************************/

#include<iostream>
using namespace std;

class DealWithFileData{
	private:
	typedef struct Data{
		string vm;
		string date;
		int vm_num;
		int date_num;
		Data(char data[MAX_DATA_NUM]){
			string temp = data;
			int begin = 0, count = 0;
			for(size_t j = 0; j < temp.size(); j++){
				if(temp[j] == '\t'){
					count++;
					if(count == 1)
						begin = j+1;
					if(count == 2){
						vm = temp.substr(begin, j-begin);
						begin = j+1;
					}
				}
				if(temp[j] == ' '){
					date = temp.substr(begin, j-begin);
					break;
				}
			}
			vm_num = de_vm(vm);
			date_num = de_date(date);
		}
	}*pData;

	public:
	static int de_vm(string myvm){
		int temp;
		string temp_sub = myvm.substr(6, myvm.size()-6);
		str2int(temp, temp_sub);
		return temp;
	}

	static string de_vm(int myvm){
		return "flavor" + to_string(myvm);	
	}


	static int de_date(string mydate){
		int temp = 0;
		temp += de_year(mydate.substr(0, 4));
		temp += de_month(mydate.substr(5, 2));
		int temp_num;
		str2int(temp_num, mydate.substr(8, 2));
		temp += temp_num;
		return temp;
	}


	static void str2int(int &int_temp, string string_temp){
		stringstream stream(string_temp);
		stream>>int_temp;
	}
	
	static int de_month(string mon){
		int temp, ans = 0;
		int lemp[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
		str2int(temp, mon);
		for(int i = 0; i < temp-1; i++)
			ans += lemp[i];
		return ans;
	}

	static int de_year(string year){
		int temp_num = 0;
		int ans = 0;
		str2int(temp_num, year);
		ans += temp_num-2015;
		return ans * 365;
	}
	
	static void decomp_input(char * info[MAX_INFO_NUM], int& max_server_cpu, int &max_server_mem, int &max_server_disk, int &my_vm_num, vector<vector<int>> &myvec, int &dura_time, bool& opt_target){
		
		string temp = info[0];
		int count = 0, begin;
		for(size_t i = 0; i < temp.size(); i++){
			if(count == 0 && temp[i] == ' '){
				str2int(max_server_cpu, temp.substr(0, i));
				begin = i+1;
				count++;
				continue;
			}
			if(count == 1 && temp[i] == ' '){	
				str2int(max_server_mem, temp.substr(begin, i));
				str2int(max_server_disk, temp.substr(i+1, temp.size()));
				break;
			}
		}
		
		temp = info[2];
		str2int(my_vm_num, temp);

		count = 0; begin = 0;
		for(int i = 3; i < 3+my_vm_num; i++){
			temp = info[i]; 
			vector<int> temp_vec;
			for(size_t j = 0; j < temp.size(); j++){
				if(count == 0 && temp[j] == ' '){
					temp_vec.push_back(de_vm(temp.substr(begin, j) ));
					begin = j+1;
					count++;
					continue;
				}
				if(count == 1 && temp[j] == ' '){
					int a;
					str2int(a, temp.substr(begin, j));
					temp_vec.push_back(a);
					str2int(a, temp.substr(j+1, temp.size()));
					temp_vec.push_back(a);
					break;
				}
			}
			myvec.push_back(temp_vec);
			count = 0;
			begin = 0;
		}
	
		temp = info[4+my_vm_num];
		//cout << temp << endl;
		if(temp[0] == 'C'){
			//cout << "cpu" << endl;
			opt_target = true;
		}else{
			//cout << "mem" << endl;
			opt_target = false;
		}
			


		temp = info[6+my_vm_num];
		int start = de_date(temp.substr(0, 11));
		temp = info[7+my_vm_num];
		int end = de_date(temp.substr(0, 11));
		dura_time = end-start;

	}




	static void decomp_data(char* data[MAX_DATA_NUM], vector<vector<int>>& my_data, int &data_time_length){

		int start_date = Data(data[0]).date_num;
		int data_length;
		for(data_length = 0; data[data_length] != NULL; data_length++){}
		int end_date = Data(data[data_length-1]).date_num;
		data_time_length = end_date - start_date + 1;
		
		vector<vector<int >> temp_data(25, vector<int>(data_time_length, 0));
		my_data = temp_data;
		
		for(int i = 0; i < data_length; i++){
			pData curr_data = new Data(data[i]);
			int curr_flavor = curr_data->vm_num;
			int curr_date = curr_data->date_num - start_date;
			my_data[curr_flavor][curr_date]++;
		}
	}

	static void write_to_res(vector<int> &pred_vec, vector<vector<int>> &distri_plan, char* &result_file){

		int sum_flavor = 0;
		for(size_t i = 0; i < pred_vec.size(); i++){
			if(pred_vec[i] >= 0)
				sum_flavor += pred_vec[i];
		}
		vector<string> my_res_file;
		my_res_file.push_back(to_string(sum_flavor)+"\n");
		for(size_t i = 0; i < pred_vec.size(); i++){
			if(pred_vec[i] >= 0){
				string temp = DealWithFileData::de_vm(i) + " " + to_string(pred_vec[i]) + "\n";
				my_res_file.push_back(temp);
			}
		}
		my_res_file.push_back("\n");
		my_res_file.push_back(to_string(distri_plan.size()) + "\n");
		for(size_t i = 0; i < distri_plan.size(); i++){
			string temp = "";
			if(i == 0){
				for(size_t j = 0; j < pred_vec.size(); j++){
					if(pred_vec[j] == 0)
						temp += " " + DealWithFileData::de_vm(j) + " 0";
				}
			}
			for(size_t j = 0; j < distri_plan[i].size(); j++){
				if(distri_plan[i][j] != 0){
					temp += " " + DealWithFileData::de_vm(j) + " " + to_string(distri_plan[i][j]);
				}
			}
			if(i == distri_plan.size()-1)
				my_res_file.push_back(to_string(i+1) + temp);
			else
				my_res_file.push_back(to_string(i+1) + temp + "\n");
		}
		
		string temp_str = "";
		for(size_t i = 0; i < my_res_file.size(); i++){
			temp_str += my_res_file[i];
		}
		result_file = (char*)malloc(sizeof(char)*(temp_str.length()+1));
		temp_str.copy(result_file, temp_str.length(), 0);
		*(result_file + temp_str.length()) = '\0';
		string temp = result_file;

	}
};
