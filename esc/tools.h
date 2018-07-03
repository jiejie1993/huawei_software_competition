/*************************************************************************
    > File Name: tools.h
    > Author: geeker
    > Mail: 932834897@qq.com 
    > Created Time: Tue 20 Mar 2018 06:05:29 AM PDT
 ************************************************************************/

#include<iostream>
#include <fstream>
using namespace std;


class Tools{

	public:
	template<typename T>
	static void print_vec(vector<T> a){
		for(size_t i = 0; i < a.size(); i++)
			cout << i << "-->" << a[i] << " ";
		cout << endl;
	}

	template<typename T>
	static void print_vec(vector<vector<T>> a){
		for(size_t i = 0; i < a.size(); i++){
			for(size_t j = 0; j < a[i].size(); j++)
				cout << a[i][j] << " ";
			cout << endl;
		}
	}

	template<typename T>
	static void print_cir(vector<T> a){
		for(size_t i = 0; i < a.size(); i++){
			cout << "*";
			for(size_t j = 0; j < a[i]*5; j++){
				cout << "*";
			}
			cout << endl;
		}	
	}

	template<typename T>
	static vector<T> integral(vector<T> a){
		vector<T> b(a.size(), 0);
		T sum = 0;
		for(int i = 0; i < a.size(); i++){
			sum += a[i];
			b[i] = sum;
		}
		return b;
	}

	template<typename T>
	static double pre_vec(vector<T> a){
		int n = a.size(), count = 0;
		//return (double)accumulate(a.begin(), a.end(), 0) / n;
		for(int i = 0; i < n; i++)
			if(a[i] > 0)
				count++;
		return (double)accumulate(a.begin(), a.end(), 0) / count;
	}

	template<typename T>
	static void out2file(const char * filename, vector<vector<T>> data){
		ofstream ofile;
		ofile.open(filename, ios::out | ios::trunc);
		for(int j = 0; j < data[0].size(); j++){
			ofile << j;
			for(int i = 0; i < data.size(); i++){
				ofile << ", " << data[i][j];
			}
			ofile << endl;
		}
	}

	template<typename T>
	static vector<T> norm(vector<T>& curr_data, T& dMaxValue, T& dMinValue){
		double ymax = 1; //归一i化数据范围  
		double ymin = 0;   
		std::vector<T> features; //临时特征向量  
		for (int d = 0; d < curr_data.size(); ++d)   
			  features.push_back(curr_data[d]);  
		
		for (int f = 0; f < features.size(); ++f) {  
			features[f] = (ymax-ymin)*(features[f]-dMinValue)/(dMaxValue-dMinValue+1e-8)+ymin;      
		}
		return features;
	}

	template<typename T>
	static vector<T> denorm(vector<T> data, T& dMaxValue, T& dMinValue){
		double ymax = 1; //归一i化数据范围  
		double ymin = 0;   
		
		for (int f = 0; f < data.size(); ++f) {  
			data[f] = (dMaxValue - dMinValue) * (data[f] - ymin) / (ymax - ymin) + dMinValue;
		} 
		return data;
	}



	//template<typename T>
	static vector<double> linearSmooth(vector<double> in)
	{
		int i, n = in.size();
		vector<double> out(n, 0.0);
		if (n < 5)
		{
			for ( i = 0; i <= n - 1; i++ )
			{
				out[i] = in[i];
			}
		}
		else
		{
			out[0] = ( 3.0 * in[0] + 2.0 * in[1] + in[2] - in[4] ) / 5.0;
			out[1] = ( 4.0 * in[0] + 3.0 * in[1] + 2 * in[2] + in[3] ) / 10.0;
			for ( i = 2; i <= n - 3; i++ )
			{
				out[i] = ( in[i - 2] + in[i - 1] + in[i] + in[i + 1] + in[i + 2] ) / 5.0;
			}
			out[n - 2] = ( 4.0 * in[n - 1] + 3.0 * in[n - 2] + 2 * in[n - 3] + in[n - 4] ) / 10.0;
			out[n - 1] = ( 3.0 * in[n - 1] + 2.0 * in[n - 2] + in[n - 3] - in[n - 5] ) / 5.0;
		return out;
		}
	}

	template<typename T>
	static vector<double> my_diffdata(vector<T> data){
		vector<double> mydata(data.begin(), data.end());
		vector<double> res(mydata.size(), 0.0);
		for(int i = 0; i < mydata.size(); i++){
			if(i != 0)
				res[i] += mydata[i] - mydata[i-1];
			if(i != data.size()-1)
				//res[i] = mydata[i] - mydata[i-1];
				res[i] += mydata[i] - mydata[i+1];
			if(res[i] < 0)
				res[i] = 0;
		}
		return res;
	}

	template<typename T>
	static vector<double> diffdata(vector<T> data){
		vector<double> mydata(data.begin(), data.end());
		vector<double> res(mydata.size(), 0.0);
		res[0] = mydata[0];
		for(int i = 1; i < mydata.size(); i++){
			res[i] = mydata[i] - mydata[i-1];
		}
		return res;
	}

	template<typename T>
	static vector<double> undiffdata(vector<T> data){
		vector<double> mydata(data.begin(), data.end());
		vector<double> res(mydata.size(), 0.0);
		res[0] = mydata[0];
		for(int i = 1; i < mydata.size(); i++){
			res[i] = mydata[i] + res[i-1];
		}
		return res;
	}

	/*template<typename T>
	static vector<double> denoise(vector<T> data){	
		vector<double> mydata(data.begin(), data.end());
		
		vector<double> res = mydata;
		vector<double> mydiff = Tools::my_diffdata(data);
		double pre = Tools::pre_vec(mydata);
		for(int i = 0; i < 0.05 * data.size(); i++){
			vector<double>::iterator biggest = max_element(mydiff.begin(), mydiff.end());
			int posi = distance(mydiff.begin(), biggest);
			res[posi] = pre;
		}
		return res;
	}*/

	static bool denoise(vector<double>& data)
	{
		bool flag_res = false;
		vector<double> my_data;
		for (auto element : data) {
			if (element > 0)
				my_data.push_back(element);
		}
		double mean = accumulate(my_data.begin(), my_data.end(), 0.0) / my_data.size();
		double var = 0;
		for (auto element : my_data) {
			var += pow(element - mean, 2);
		}
		var = sqrt(var / my_data.size());
		//cout << mean << ", " << var << endl;
		double low_thre = mean - var;
		double high_thre = mean + var;
		//cout << low_thre << ", " << high_thre << endl;
		for (int i = 0; i < data.size(); i++) {
			if (data[i] > high_thre) {
				flag_res = true;
				data[i] = high_thre;
			}
		}

		return flag_res;
	}
	








	template<typename T>
	static vector<double> data_trans(vector<T> data, int n){
	     int len = data.size();
	     vector<T> res(len-n+1, 0);                                              	 for(int i = 0; i < len-n+1; i++){
	         res[i] = accumulate(data.begin()+i, data.begin()+i+n, 0.0);
	     }
	     return vector<double>(res.begin(), res.end());

	 }

};
