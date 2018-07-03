#include "predict.h"
#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <string.h>
#include <vector>
#include <sstream>
#include <numeric>
#include "prodata.h"
#include <time.h>

#include <fstream>
#include <cmath>
#include <math.h>
#include "tools.h"
#include "distserver.h"


#define innode  9       //输入结点数，将输入2个加数
#define hidenode 40       //隐藏结点数，存储“携带位”
#define outnode  1      //输出结点数，将输出一个预测数字
#define alpha  0.0015      //学习速率
//#define binary_dim 80    //训练序列的长度
static int binary_dim;

#define randval(high) ( (double)rand() / RAND_MAX * high )
#define uniform_plus_minus_one ( (double)( 2.0 * rand() ) / ((double)RAND_MAX + 1.0) - 1.0 )  //均匀随机分布

using namespace std;


double sigmoid(double x)
{
	return 1.0 / (1.0 + exp(-x));
}

//激活函数的导数，y为激活函数值
double dsigmoid(double y)
{
	return y * (1.0 - y);
}

//tanh的导数，y为tanh值
double dtanh(double y)
{
	y = tanh(y);
	return 1.0 - y * y;
}

double dtanh2(double y)
{
	return 1.0 - y * y;
}


class RNN
{
public:
	RNN();
	virtual ~RNN();
	void train(vector<vector<double>>& , vector<double>&);
	vector<double> predict(vector<vector<double>>&, int);

public:
	double W_I[innode][hidenode];     //连接输入与隐含层单元中输入门的权值矩阵
	double U_I[hidenode][hidenode];   //连接上一隐层输出与本隐含层单元中输入门的权值矩阵
	double W_F[innode][hidenode];     //连接输入与隐含层单元中遗忘门的权值矩阵
	double U_F[hidenode][hidenode];   //连接上一隐含层与本隐含层单元中遗忘门的权值矩阵
	double W_O[innode][hidenode];     //连接输入与隐含层单元中遗忘门的权值矩阵
	double U_O[hidenode][hidenode];   //连接上一隐含层与现在时刻的隐含层的权值矩阵
	double W_G[innode][hidenode];     //用于产生新记忆的权值矩阵
	double U_G[hidenode][hidenode];   //用于产生新记忆的权值矩阵
	double W_out[hidenode][outnode];  //连接隐层与输出层的权值矩阵

	double *x;             //layer 0 输出值，由输入向量直接设定
						   //double *layer_1;     //layer 1 输出值
	double *y;             //layer 2 输出值
};

void winit(double w[], int n) //权值初始化
{
	for (int i = 0; i<n; i++)
		w[i] = uniform_plus_minus_one;  //均匀随机分布
}

RNN::RNN()
{
	x = new double[innode];
	y = new double[outnode];
	winit((double*)W_I, innode * hidenode);
	winit((double*)U_I, hidenode * hidenode);
	winit((double*)W_F, innode * hidenode);
	winit((double*)U_F, hidenode * hidenode);
	winit((double*)W_O, innode * hidenode);
	winit((double*)U_O, hidenode * hidenode);
	winit((double*)W_G, innode * hidenode);
	winit((double*)U_G, hidenode * hidenode);
	winit((double*)W_out, hidenode * outnode);
}

RNN::~RNN()
{
	delete x;
	delete y;
}


void RNN::train(vector<vector<double>>& data, vector<double>& lable)
{
	cout << "trainning..." << endl;
	int epoch, i, j, k, m, p;
	vector<double*> I_vector;      //输入门
	vector<double*> F_vector;      //遗忘门
	vector<double*> O_vector;      //输出门
	vector<double*> G_vector;      //新记忆
	vector<double*> S_vector;      //状态值
	vector<double*> h_vector;      //输出值
	vector<double> y_delta;        //保存误差关于输出层的偏导
	binary_dim = lable.size();

	for (epoch = 0; epoch<20000; epoch++)  //训练次数
	{
		double e = 0.0;  //误差
		
		//vector<double> a(binary_dim, 0);
		//vector<double> b(binary_dim, 0);
		//vector<double> c(binary_dim, 0);
		vector<double> predict(binary_dim, 0);
	
		vector<double*> S_vector;      //状态值
		vector<double*> h_vector;      //输出值
									   //在0时刻是没有之前的隐含层的，所以初始化一个全为0的
		double *S = new double[hidenode];     //状态值
		double *h = new double[hidenode];     //输出值

		for (int i = 0; i<hidenode; i++)
		{
			S[i] = 0;
			h[i] = 0;
		}
		S_vector.push_back(S);
		h_vector.push_back(h);


		/*
		for (int g = 0; g < binary_dim; g++)
		{
			a[g] = data[0][g];
		}

		for (int g = 0; g < binary_dim; g++)
		{
			b[g] = 1.0;
		}

		for (int g = 0; g < binary_dim; g++)
		{
			c[g] = lable[g];
		}*/
		
		//正向传播
		for (p = 0; p<binary_dim; p++)           //循环遍历二进制数组，从最低位开始
		{
			//靠靠
			for(int i = 0; i < innode-1; i++)
				x[i] = data[i][p];
			x[innode-1] = 1.0;
			
			double t = (double)lable[p];          //实际值
			double *in_gate = new double[hidenode];     //输入门
			double *out_gate = new double[hidenode];    //输出门
			double *forget_gate = new double[hidenode]; //遗忘门
			double *g_gate = new double[hidenode];      //新记忆
			double *state = new double[hidenode];       //状态值
			double *h = new double[hidenode];           //隐层输出值

			for (j = 0; j<hidenode; j++)
			{
				//输入层转播到隐层
				double inGate = 0.0;
				double outGate = 0.0;
				double forgetGate = 0.0;
				double gGate = 0.0;
				double s = 0.0;

				for (m = 0; m<innode; m++)
				{
					inGate += x[m] * W_I[m][j];
					outGate += x[m] * W_O[m][j];
					forgetGate += x[m] * W_F[m][j];
					gGate += x[m] * W_G[m][j];
				}

				double *h_pre = h_vector.back();
				double *state_pre = S_vector.back();
				for (m = 0; m<hidenode; m++)
				{
					inGate += h_pre[m] * U_I[m][j];
					outGate += h_pre[m] * U_O[m][j];
					forgetGate += h_pre[m] * U_F[m][j];
					gGate += h_pre[m] * U_G[m][j];
				}

				in_gate[j] = sigmoid(inGate);
				out_gate[j] = sigmoid(outGate);
				forget_gate[j] = sigmoid(forgetGate);
				g_gate[j] = tanh(gGate);

				double s_pre = state_pre[j];
				state[j] = forget_gate[j] * s_pre + g_gate[j] * in_gate[j];
				h[j] = out_gate[j] * tanh(state[j]);
			}


			for (k = 0; k<outnode; k++)
			{
				//隐藏层传播到输出层
				double out = 0.0;
				for (j = 0; j<hidenode; j++)
					out += h[j] * W_out[j][k];
				y[k] = sigmoid(out);               //输出层各单元输出
			}

			predict[p] = y[0];   //记录预测值

												   //保存隐藏层，以便下次计算
			I_vector.push_back(in_gate);
			F_vector.push_back(forget_gate);
			O_vector.push_back(out_gate);
			S_vector.push_back(state);
			G_vector.push_back(g_gate);
			h_vector.push_back(h);

			//保存标准误差关于输出层的偏导
			y_delta.push_back((t - y[0]) * dsigmoid(y[0]));
			e += fabs(t - y[0]);          //误差
		}


		//误差反向传播

		//隐含层偏差，通过当前之后一个时间点的隐含层误差和当前输出层的误差计算
		double h_delta[hidenode];
		double *O_delta = new double[hidenode];
		double *I_delta = new double[hidenode];
		double *F_delta = new double[hidenode];
		double *G_delta = new double[hidenode];
		double *state_delta = new double[hidenode];
		//当前时间之后的一个隐藏层误差
		double *O_future_delta = new double[hidenode];
		double *I_future_delta = new double[hidenode];
		double *F_future_delta = new double[hidenode];
		double *G_future_delta = new double[hidenode];
		double *state_future_delta = new double[hidenode];
		double *forget_gate_future = new double[hidenode];
		for (j = 0; j<hidenode; j++)
		{
			O_future_delta[j] = 0;
			I_future_delta[j] = 0;
			F_future_delta[j] = 0;
			G_future_delta[j] = 0;
			state_future_delta[j] = 0;
			forget_gate_future[j] = 0;
		}
		for (p = binary_dim - 1; p >= 0; p--)
		{
			//x[0] = data[0][p];
			//x[1] = 1.0;

			for(int i = 0; i < innode-1; i++)
				x[i] = data[i][p];
			x[innode-1] = 1.0;
			
			//当前隐藏层
			double *in_gate = I_vector[p];     //输入门
			double *out_gate = O_vector[p];    //输出门
			double *forget_gate = F_vector[p]; //遗忘门
			double *g_gate = G_vector[p];      //新记忆
			double *state = S_vector[p + 1];     //状态值
			double *h = h_vector[p + 1];         //隐层输出值

												 //前一个隐藏层
			double *h_pre = h_vector[p];
			double *state_pre = S_vector[p];

			for (k = 0; k < outnode; k++)  //对于网络中每个输出单元，更新权值
			{
				//更新隐含层和输出层之间的连接权
				for (j = 0; j<hidenode; j++)
					W_out[j][k] += alpha * y_delta[p] * h[j];
			}

			//对于网络中每个隐藏单元，计算误差项，并更新权值
			for (j = 0; j<hidenode; j++)
			{
				h_delta[j] = 0.0;
				for (k = 0; k<outnode; k++)
				{
					h_delta[j] += y_delta[p] * W_out[j][k];
				}
				for (k = 0; k<hidenode; k++)
				{
					h_delta[j] += I_future_delta[k] * U_I[j][k];
					h_delta[j] += F_future_delta[k] * U_F[j][k];
					h_delta[j] += O_future_delta[k] * U_O[j][k];
					h_delta[j] += G_future_delta[k] * U_G[j][k];
				}

				O_delta[j] = 0.0;
				I_delta[j] = 0.0;
				F_delta[j] = 0.0;
				G_delta[j] = 0.0;
				state_delta[j] = 0.0;

				//隐含层的校正误差
				O_delta[j] = h_delta[j] * tanh(state[j]) * dsigmoid(out_gate[j]);
				state_delta[j] = h_delta[j] * out_gate[j] * dtanh(state[j]);// + state_future_delta[j] * forget_gate_future[j];
				F_delta[j] = state_delta[j] * state_pre[j] * dsigmoid(forget_gate[j]);
				I_delta[j] = state_delta[j] * g_gate[j] * dsigmoid(in_gate[j]);
				G_delta[j] = state_delta[j] * in_gate[j] * dtanh2(g_gate[j]);

				//更新前一个隐含层和现在隐含层之间的权值
				for (k = 0; k<hidenode; k++)
				{
					U_I[k][j] += alpha * I_delta[j] * h_pre[k];
					U_F[k][j] += alpha * F_delta[j] * h_pre[k];
					U_O[k][j] += alpha * O_delta[j] * h_pre[k];
					U_G[k][j] += alpha * G_delta[j] * h_pre[k];
				}

				//更新输入层和隐含层之间的连接权
				for (k = 0; k<innode; k++)
				{
					W_I[k][j] += alpha * I_delta[j] * x[k];
					W_F[k][j] += alpha * F_delta[j] * x[k];
					W_O[k][j] += alpha * O_delta[j] * x[k];
					W_G[k][j] += alpha * G_delta[j] * x[k];
				}

			}

			if (p == binary_dim - 1)
			{
				delete  O_future_delta;
				delete  F_future_delta;
				delete  I_future_delta;
				delete  G_future_delta;
				delete  state_future_delta;
				delete  forget_gate_future;
			}

			O_future_delta = O_delta;
			F_future_delta = F_delta;
			I_future_delta = I_delta;
			G_future_delta = G_delta;
			state_future_delta = state_delta;
			forget_gate_future = forget_gate;
		}
		delete  O_future_delta;
		delete  F_future_delta;
		delete  I_future_delta;
		delete  G_future_delta;
		delete  state_future_delta;

		
		//if(epoch == 4800){
		//	Tools::out2file("myresult.csv", vector<vector<double>>{predict, lable});
		//}
		
		if(epoch % 2000 == 0){
			cout << "erro-->" << e << endl;
		}

		//cout << epoch << endl;


		for (i = 0; i<I_vector.size(); i++)
			delete I_vector[i];
		for (i = 0; i<F_vector.size(); i++)
			delete F_vector[i];
		for (i = 0; i<O_vector.size(); i++)
			delete O_vector[i];
		for (i = 0; i<G_vector.size(); i++)
			delete G_vector[i];

		I_vector.clear();
		F_vector.clear();
		O_vector.clear();
		G_vector.clear();
		S_vector.clear();
		h_vector.clear();
		y_delta.clear();
	}
}

vector<double> RNN::predict(vector<vector<double>>& data, int Dura_Time)
{
	int i, j, k, m, p;
	double currval = data[0][0];
	binary_dim = data[0].size();
	vector<double> predict(binary_dim, 0);

	//vector<double> a(binary_dim, 0);
	//vector<double> b(binary_dim, 0);

	vector<double*> S_vector;      //状态值
	vector<double*> h_vector;      //输出值
								   //在0时刻是没有之前的隐含层的，所以初始化一个全为0的
	//靠縎 H
	double *S = new double[hidenode];     //状态值
	double *h = new double[hidenode];     //输出值

	for (int i = 0; i<hidenode; i++)
	{
		S[i] = 0;
		h[i] = 0;
	}
	S_vector.push_back(S);
	h_vector.push_back(h);

	/*
	for (int g = 0; g < binary_dim; g++)
	{
		a[g] = data[0][g];
		b[g] = 1.0;
	}*/

	//靠縫redict[0]
	//predict[0] = 0;

	//正向传播
	for (p = 0; p < binary_dim; p++)           //循环遍历二进制数组，从最低位开始
	{
		//x[0] = data[0][p];
		//x[1] = 1.0;
		for(int i = 0; i < innode-1; i++)
			x[i] = data[i][p];
		x[innode-1] = 1.0;
		double *in_gate = new double[hidenode];     //输入门
		double *out_gate = new double[hidenode];    //输出门
		double *forget_gate = new double[hidenode]; //遗忘门
		double *g_gate = new double[hidenode];      //新记忆
		double *state = new double[hidenode];       //状态值
		double *h = new double[hidenode];           //隐层输出值

		for (j = 0; j < hidenode; j++)
		{
			//输入层转播到隐层
			double inGate = 0.0;
			double outGate = 0.0;
			double forgetGate = 0.0;
			double gGate = 0.0;
			double s = 0.0;

			for (m = 0; m < innode; m++)
			{
				inGate += x[m] * W_I[m][j];
				outGate += x[m] * W_O[m][j];
				forgetGate += x[m] * W_F[m][j];
				gGate += x[m] * W_G[m][j];
			}

			double *h_pre = h_vector.back();
			double *state_pre = S_vector.back();
			for (m = 0; m < hidenode; m++)
			{
				inGate += h_pre[m] * U_I[m][j];
				outGate += h_pre[m] * U_O[m][j];
				forgetGate += h_pre[m] * U_F[m][j];
				gGate += h_pre[m] * U_G[m][j];
			}

			in_gate[j] = sigmoid(inGate);
			out_gate[j] = sigmoid(outGate);
			forget_gate[j] = sigmoid(forgetGate);
			g_gate[j] = tanh(gGate);

			double s_pre = state_pre[j];
			state[j] = forget_gate[j] * s_pre + g_gate[j] * in_gate[j];
			h[j] = out_gate[j] * tanh(state[j]);
		}


		for (k = 0; k < outnode; k++)
		{
			double out = 0.0;
			for (j = 0; j < hidenode; j++)
				out += h[j] * W_out[j][k];
			y[k] = sigmoid(out);               //输出层各单元输出
		}

		predict[p] = y[0];
		
		S_vector.push_back(state);
		h_vector.push_back(h);
	}

	return predict;
}

//你要完成的功能总入口

/**********************predict single flavor function**************************************/
vector<double> pred_si_flavor(vector<double> in_curr_data, int Dura_Time, int flavor){
	
	static int nextstep = 1;
	static int weeknum = 2;	
	static int daynum = 6;


	vector<double> features, predict, temp_data, temp_lable; //临时特征向量 
	vector<vector<double>> my_data1(daynum, vector<double>{});



	//cout << "data3-->" << endl;
	//Tools::print_vec(in_curr_data);
	//cout << endl;
	
	//vector<double> in_curr_data_norm;

	vector<double> my_temp;	
	my_temp = Tools::data_trans(in_curr_data, Dura_Time);
	
	double dMaxValue = *max_element(my_temp.begin(),my_temp.end());  //求最大值  
	double dMinValue = *min_element(my_temp.begin(),my_temp.end());  //求最小值 

	my_temp = Tools::norm(my_temp, dMaxValue, dMinValue); //normalize features

	vector<double> oridata = Tools::norm(in_curr_data, dMaxValue, dMinValue);
	int lablen = 0, index = 0;
	
	for(index = my_temp.size()-1; index > daynum; index -= Dura_Time){
		lablen++;
	}
	index += Dura_Time;

	//cout << "hello" << endl;

	for(; index <= my_temp.size()-1; index += Dura_Time){
		//cout << index << endl;
		temp_lable.push_back(my_temp[index]);
		for(int i = index-1; i >= index-daynum; i--){
			my_data1[i - index + daynum].push_back(oridata[i]);
		}
	}

	vector<vector<double>> my_data2(weeknum, vector<double>{});
	
	for(int i = 0; i < weeknum; i++){
		for(int j = 0; j < temp_lable.size()-weeknum; j++){
			my_data2[i].push_back(temp_lable[i+j]);
		}
	}

	vector<double> my_lable(temp_lable.size()-weeknum, 0.0);
	vector<vector<double>> my_data(daynum+weeknum, vector<double>(temp_lable.size()-weeknum, 0.0));	
	for(int i = 0; i < daynum; i++){
		for(int j = 0; j < temp_lable.size()-weeknum; j++){
			my_data[i][j] = my_data1[i][j + weeknum];
		}
	}


	for(int i = daynum; i < daynum + weeknum; i++){
		for(int j = 0; j < temp_lable.size()-weeknum; j++){
			my_data[i][j] = my_data2[i-daynum][j];
		}
	}


	for(int j = 0; j < temp_lable.size()-weeknum; j++){
		my_lable[j] = temp_lable[j+weeknum];
	}

/*
	for(int i = 0; i < my_data.size(); i++){
		my_data[i] = vector<double>(my_data[i].begin() + weeknum, my_data[i].end());
	}

	for(int i = 0; i < my_data[0].size(); i++){
		for(int j = 0; j < weeknum; j++){
			my_data[daynum + j].push_back(temp_lable[i+j]);
		}
	}*/


	/*
	cout << my_data.size() << endl;
	for(int i = 0; i < my_data.size(); i++)
		cout << my_data[i].size() << ", ";
	cout << endl;
	*/
		
	//cout << oridata.size() << endl;
	//cout << my_temp.size() << endl;

	/*
	vector<double> temp;

	for(int n = my_temp.size()-1; n >= 0; n -= Dura_Time){
		temp.push_back(my_temp[n]);
	}
	reverse(temp.begin(), temp.end());		


	//vector<double> curr_data = temp;
	//vector<double> my_curr_data = temp;
	vector<double> data = temp;
	
	for(int i = 0; i < data.size()-nextstep-innode+2; i++){
		//temp_data.push_back(data[i]);
		temp_lable.push_back(data[i+nextstep+innode-2]);
	}//init data, lables

	
	for(int j = 0; j < innode-1; j++){
		for(int i = 0; i < data.size()-nextstep-innode+2; i++){
			temp_data.push_back(data[i+j]);		
		}
		my_data.push_back(temp_data);
		temp_data.clear();
	}
	*/


	RNN rnn;

	//rnn.train(my_data, temp_lable);
	
	rnn.train(my_data, my_lable);

	//my_data[0].push_back(temp_lable.back());

	//my_data[0].push_back(my_data[1].back());
	//my_data[1].push_back(my_data[2].back());	
	//my_data[2].push_back(temp_lable.back());

	for(int i = daynum-1; i >= 0; i--){
		my_data[i].push_back(oridata[oridata.size() - 1 + (i - daynum + 1)]);
	}
	for(int i = daynum; i < weeknum+daynum; i++){
		if(i == weeknum + daynum -1)
			my_data[i].push_back(my_lable.back());
		else
			my_data[i].push_back(my_data[i+1].back());
	}
	
	//cout << "mydata.size" << my_data.size() << endl; 

	predict = rnn.predict(my_data, nextstep);


	//cout << "predict1-------------->" << endl;
	//Tools::print_vec(predict);
	//cout << endl;

	vector<double> my_temp_smooth(my_lable.begin(), my_lable.begin() + predict.size());
	
	//predict = Tools::undiffdata(predict);

	string temp_s = string("result") + 	to_string(flavor) + string(".csv");
	Tools::out2file(temp_s.c_str(), vector<vector<double>>{my_temp_smooth, predict });
	
	
	//data.push_back(predict.back());
	//predict = data;

	predict = Tools::denorm(predict, dMaxValue, dMinValue);
	
	//cout << "predict2-------------->" << endl;
	//Tools::print_vec(predict);
	//cout << endl;

	//predict = Tools::undiffdata(predict);
	
	cout << "predict3-------------->" << endl;
	Tools::print_vec(predict);
	cout << endl;
	
	return predict;
}


/***************************predict all flavor function*****************************/
void pred_flavor(vector<int> &pred_vec, vector<vector<int>> &Flavor_Info, vector<vector<int>> &My_Data, int &Dura_Time){

	bool flag = false;
	for(size_t i = 0, j; i < My_Data.size(); i++){
		for(j = 0; j < Flavor_Info.size(); j++){
			if((int)i == Flavor_Info[j][0]){
				flag = true;
				break;
			}
		}
		if(flag){

			/*pred_vec[i] = 0;	
			
			// using the latest n days
			
			int sum = 0;
			for(size_t k = My_Data[i].size()-1; k >= My_Data[i].size()-Dura_Time; k--){
				sum += My_Data[i][k];
			}
			pred_vec[i] = sum;*/
			
			
		
			
			// using my predict single flavor function

			
			vector<double> my_temp(My_Data[i].begin(), My_Data[i].end());
	
			//Tools::denoise(my_temp);		
	

			vector<double> temp_pred = pred_si_flavor(my_temp, Dura_Time, i);

			if(temp_pred.back() > 0)
				pred_vec[i] = round(temp_pred.back());

			//pred_vec[i] = i*10;
			
			//cout << "pred_vec--------------------->";
			//cout << pred_vec[i] << endl;
			
				
		}
		flag = false;
	}
	
/*	
	pred_vec[1] = 3;
	pred_vec[2] = 11;
	pred_vec[3] = 0;
	pred_vec[4] = 1;
	pred_vec[5] = 10;
*/

}


void predict_server(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, char * filename)
{
	long beginTime = clock();
	static vector<vector<int> > Flavor_Info;
	static int Max_Server_Cpu;
	static int Max_Server_Mem;
	static int Max_Server_Disk;
	static int Flavor_Num;
	static int Dura_Time;
	static int Data_Time_Length;
	static vector<vector<int>> My_Data;
	static bool Opt_Target;

	DealWithFileData::decomp_input(info, Max_Server_Cpu, Max_Server_Mem, Max_Server_Disk, Flavor_Num, Flavor_Info, Dura_Time, Opt_Target);
	DealWithFileData::decomp_data(data, My_Data, Data_Time_Length);
	
	cout << Opt_Target << endl;

	vector<int> pred_vec(25, -1);
	vector<vector<int>> distri_plan;

	pred_flavor(pred_vec, Flavor_Info, My_Data, Dura_Time);
	dist_server(distri_plan, pred_vec, Max_Server_Cpu, Max_Server_Mem, Flavor_Num, Flavor_Info, Opt_Target);

	Tools::print_vec(distri_plan);


	char * result_file;
	DealWithFileData::write_to_res(pred_vec, distri_plan, result_file);

	write_result(result_file, filename);

	long endTime = clock();
	cout << "endTime-beginTime:" << endTime-beginTime << "ms" << endl;
}

