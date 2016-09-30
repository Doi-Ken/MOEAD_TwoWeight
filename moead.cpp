#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<algorithm>
#include<vector>
#include<list>
#include<queue>
#include<stack>
#include<fstream>
#include<iomanip>
#include<ctime>
#include <functional>
#include"moead.h"


void t_neighborhood(int **B, double **lambda, int N, int T,int ob){

	int **index = new int*[N];
	double **len = new double*[N];
	int *including = new int [N];
	for (int i = 0; i < N; i++){
		including[i] = 0;
	}
	
	for (int n = 0; n < N; n++){

		//std::priority_queue<double, std::vector<double>, std::greater<double> > lenmin;
		index[n] = new int[N];
		len[n] = new double[N];
		for (int i = 0; i < N; i++){
			index[n][i] = i;
		}

		for (int i = 0; i < N; i++){
			len[n][i] = 0.0;
			for (int o = 0; o < ob; o++){
				len[n][i] += (lambda[n][o] - lambda[i][o]) * (lambda[n][o] - lambda[i][o]);
			}
			len[n][i] = sqrt(len[n][i]);
			//lenmin.push(len[n][i]);
		}
		/*
		int t = 0;
		while (t < T){
			for (int i = 0; i < N; i++){
				if (lenmin.top() >= len[n][i] && including[i] == 0){
					including[i] = 1;
					B[n][t++] = i;

				}
			}
			lenmin.pop();
		}*/

		double min = 0.0;
		double nextmin = 0.0;
		int t = 0;
		while (t < T){
			min = nextmin;
			nextmin = 100000000.0;
			for (int i = 0; i < N; i++){
				if (min >= len[n][i]){
					if (including[i] == 0){
						including[i] = 1;
						B[n][t++] = i;
					}
				}
				else if (nextmin > len[n][i] && min < len[n][i]){
					nextmin = len[n][i];
				}
				if (t == T){
					break;
				}
			}
		}

		for (int i = 0; i < N; i++){
			including[i] = 0;
		}

	}



	//sort

	//int t = 0;
	//for (int n = 0; n < N; n++){
	//	for (int i = 0; i < N - 1; i++) {
	//		for (int j = N - 1; j > i; j--) {
	//			if (len[n][index[n][j - 1]] > len[n][index[n][j]]) {
	//				t = index[n][j];
	//				index[n][j] = index[n][j - 1];
	//				index[n][j - 1] = t;
	//			}
	//		}
	//	}
	//}

	//

	//for (int i = 0; i < N; i++){
	//	for (int t = 0; t < T; t++){
	//		B[i][t] = index[i][t];
	//	}
	//}
	/*		for (int i = 0; i < N; i++){
	for (int t = 0; t < T; t++){
	cout << len[i][B[i][t]] << " ";
	}
	cout << endl;
	}

	*/
	delete[] including;

}

double a_g_t_e(double *fit, double *lambda, double *z,int ob){
	double min = lambda[0] * (double)fabs(z[0] - fit[0]);
	for (int o = 0; o < ob; o++){
		double g = lambda[o] * (double)fabs(z[o] - fit[o]);
		if (min > g){
			min = g;
		}
	}
	return min;
}

double g_t_e(double *fit, double *lambda, double *z, int ob){
	double max = lambda[0] * (double)fabs(z[0] - fit[0]);
	for (int o = 0; o < ob; o++){
		double g = lambda[o] * (double)fabs(z[o] - fit[o]);
		if (max < g){
			max = g;
		}
	}
	return max;
}

double p_b_i_min_d2(double *fit, double *lambda, double *z, int ob, double PENALTY){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = fit[o] - z[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return d_n;
}


double p_b_i_d2(double *fit, double *lambda, double *z, int ob, double PENALTY){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = z[o] - fit[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return d_n;
}

double p_b_i(double *fit, double *lambda, double *z,int ob,double PENALTY){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = z[o] - fit[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return (d_t + PENALTY * d_n);
}

double weight_sum(double *fit, double *lambda,int ob){

	double sum = 0.0;

	for (int o = 0; o < ob; o++){
		sum += lambda[o] * fit[o];
	}


	return sum;
}

void max_z(double *fit, double *z,int ob,double ALPHA, double *fitmax){
	for (int o = 0; o < ob; o++){
		if (ALPHA * fit[o] > z[o]){
			z[o] = ALPHA * fit[o];
		}
	}

	for (int o = 0; o < ob; o++){
		if (fit[o] > fitmax[o]){
			fitmax[o] = fit[o];
		}
	}
}


void min_z(double *fit, double *z, int ob, double ALPHA, double *fitmin){

	for (int o = 0; o < ob; o++){
		if (ALPHA * fit[o] < z[o]){
			z[o] = ALPHA * fit[o];
		}
	}
	
	for (int o = 0; o < ob; o++){
		if (fit[o] < fitmin[o]){
			fitmin[o] = fit[o];
		}
	}

}

double p_b_i_min(double *fit, double *lambda, double *z, int ob, double PENALTY){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = fit[o] - z[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return (d_t + PENALTY * d_n);
}

double quadratic_p_b_i(double *fit, double *lambda, double *z, int ob, double PENALTY, double shita){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = z[o] - fit[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return (d_t + PENALTY * d_n * d_n/ shita);
}

double quadratic_p_b_i_min(double *fit, double *lambda, double *z, int ob, double PENALTY, double shita){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = fit[o] - z[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;
	//std::cout <<d_n << " " << shita << std::endl;

	delete[] z_minus_fit;
	return (d_t + PENALTY * d_n * d_n/shita);
}

double cubic_p_b_i(double *fit, double *lambda, double *z, int ob, double PENALTY){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = z[o] - fit[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return (d_t + PENALTY * d_n * d_n * d_n);
}

double cubic_p_b_i_min(double *fit, double *lambda, double *z, int ob, double PENALTY){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = fit[o] - z[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return (d_t + PENALTY * d_n * d_n * d_n);
}

double p_b_i_min_double(double *fit, double *lambda, double *z, int ob, double PENALTY, double PENALTY2,double shita){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double d_n_limit = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = fit[o] - z[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	d_n_limit = abs(d_t) * tan(shita);

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	if (d_n < d_n_limit){
		return (d_t + PENALTY * d_n);
	}
	else{
	/*	std::cout << "lambda: ";
		for (int o = 0; o < ob; o++){
			std::cout << lambda[o] << " ";
		}
		std::cout << std::endl;
		std::cout << "f: ";
		for (int o = 0; o < ob; o++){
			std::cout << fit[o] << " ";
		}
		std::cout << std::endl;
		std::cout << "z: ";
		for (int o = 0; o < ob; o++){
			std::cout << z[o] << " ";
		}
		std::cout << std::endl << std::endl;
*/

		return (d_t + PENALTY * d_n_limit + PENALTY2 * (d_n - d_n_limit));
	}
}

double p_b_i_min_double_normalize_fix(double *fit, double *lambda, double *z, int ob, double PENALTY, double PENALTY2, double shita){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double d_n_limit = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = fit[o] - z[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	d_n_limit = shita;

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;
	//std::cout << shita << std::endl;


	delete[] z_minus_fit;
	if (d_n < d_n_limit){
		return (d_t + PENALTY * d_n);
	}
	else{
		/*	std::cout << "lambda: ";
		for (int o = 0; o < ob; o++){
		std::cout << lambda[o] << " ";
		}
		std::cout << std::endl;
		std::cout << "f: ";
		for (int o = 0; o < ob; o++){
		std::cout << fit[o] << " ";
		}
		std::cout << std::endl;
		std::cout << "z: ";
		for (int o = 0; o < ob; o++){
		std::cout << z[o] << " ";
		}
		std::cout << std::endl << std::endl;
		*/

		return (d_t + PENALTY * d_n_limit + PENALTY2 * (d_n - d_n_limit));
	}
}

double p_b_i_double(double *fit, double *lambda, double *z, int ob, double PENALTY, double PENALTY2, double shita){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double d_n_limit = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = z[o] - fit[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	d_n_limit = abs(d_t) * tan(shita);

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;

	delete[] z_minus_fit;
	if (d_n < d_n_limit){
		return (d_t + PENALTY * d_n);
	}
	else{
		return (d_t + PENALTY * d_n_limit + PENALTY2 * (d_n - d_n_limit));
	}
}

double p_b_i_double_normalize_fix(double *fit, double *lambda, double *z, int ob, double PENALTY, double PENALTY2, double shita){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double d_n_limit = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = z[o] - fit[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	d_n_limit = shita;

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;
	//std::cout << shita << std::endl;

	delete[] z_minus_fit;
	if (d_n < d_n_limit){
		return (d_t + PENALTY * d_n);
	}
	else{
		return (d_t + PENALTY * d_n_limit + PENALTY2 * (d_n - d_n_limit));
	}
}
double p_b_i_min_para0(double *fit, double *lambda, double *z, int ob, double shita){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double d1 = 0.0;
	double d2 = 0.0;

	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = fit[o] - z[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	d1 = d_t + d_n * tan(shita);

	d2 = d_n / cos(shita);
	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return (d1 + d2);
}


double p_b_i_para0(double *fit, double *lambda, double *z, int ob, double shita){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double d1 = 0.0;
	double d2 = 0.0;

	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = z[o] - fit[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}
	d_n = sqrt(d_n);

	d1 = d_t + d_n * tan(shita);

	d2 = d_n / cos(shita);
	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return (d1 + d2);
}



double inverted_p_b_i(double *fit, double *lambda, double *z, int ob, double PENALTY){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = z[o] - fit[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;
	

	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}

	d_n = sqrt(d_n);
	d_t = abs(d_t);
	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return (d_t - PENALTY * d_n);
}




double inverted_p_b_i_max(double *fit, double *lambda, double *z, int ob, double PENALTY){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = fit[o] - z[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;


	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}

	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return (d_t - PENALTY * d_n);
}

double inverted_p_b_i_double_normalize_fix(double *fit, double *lambda, double *z, int ob, double PENALTY, double PENALTY2, double shita){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double d_n_limit = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = fit[o] - z[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	d_n_limit = shita;
	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}

	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;
	//std::cout << shita << std::endl;

	delete[] z_minus_fit;
	if (d_n < d_n_limit){
		return (d_t - PENALTY * d_n);
	}
	else{
		/*	std::cout << "lambda: ";
		for (int o = 0; o < ob; o++){
		std::cout << lambda[o] << " ";
		}
		std::cout << std::endl;
		std::cout << "f: ";
		for (int o = 0; o < ob; o++){
		std::cout << fit[o] << " ";
		}
		std::cout << std::endl;
		std::cout << "z: ";
		for (int o = 0; o < ob; o++){
		std::cout << z[o] << " ";
		}
		std::cout << std::endl << std::endl;
		*/

		return (d_t - PENALTY * d_n_limit - PENALTY2 * (d_n - d_n_limit));
	}
}

double quadratic_inverted_p_b_i_max(double *fit, double *lambda, double *z, int ob, double PENALTY,double shita){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = fit[o] - z[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;


	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}

	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return (d_t - PENALTY * d_n * d_n/shita);
}






double inverted_p_b_i_min(double *fit, double *lambda, double *z, int ob, double PENALTY){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = z[o] - fit[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;


	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}

	d_n = sqrt(d_n);
	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return (d_t - PENALTY * d_n);
}

double inverted_p_b_i_min_double_normalize_fix(double *fit, double *lambda, double *z, int ob, double PENALTY, double PENALTY2, double shita){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double d_n_limit = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = z[o] - fit[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;

	d_n_limit = shita;
	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}

	d_n = sqrt(d_n);

	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;
	//std::cout << shita << std::endl;

	delete[] z_minus_fit;
	if (d_n < d_n_limit){
		return (d_t - PENALTY * d_n);
	}
	else{
		/*	std::cout << "lambda: ";
		for (int o = 0; o < ob; o++){
		std::cout << lambda[o] << " ";
		}
		std::cout << std::endl;
		std::cout << "f: ";
		for (int o = 0; o < ob; o++){
		std::cout << fit[o] << " ";
		}
		std::cout << std::endl;
		std::cout << "z: ";
		for (int o = 0; o < ob; o++){
		std::cout << z[o] << " ";
		}
		std::cout << std::endl << std::endl;
		*/

		return (d_t - PENALTY * d_n_limit - PENALTY2 * (d_n - d_n_limit));
	}
}

double quadratic_inverted_p_b_i_min(double *fit, double *lambda, double *z, int ob, double PENALTY, double shita){
	double lambda_norm = 0.0;
	double d_t = 0.0;
	double d_n = 0.0;
	double *z_minus_fit = new double[ob];

	for (int o = 0; o < ob; o++){
		lambda_norm += lambda[o] * lambda[o];
	}
	lambda_norm = sqrt(lambda_norm);

	for (int o = 0; o < ob; o++){
		z_minus_fit[o] = z[o] - fit[o];
	}
	//d_t
	for (int o = 0; o < ob; o++){
		d_t += z_minus_fit[o] * lambda[o];
	}
	d_t /= lambda_norm;


	//d_n
	for (int o = 0; o < ob; o++){
		d_n += (z_minus_fit[o] - d_t * lambda[o] / lambda_norm) * (z_minus_fit[o] - d_t * lambda[o] / lambda_norm);
	}

	d_n = sqrt(d_n);
	//double check = 0.0;
	//for(int o = 0;o < ob;o++){
	//	check += z_minus_fit[o] * z_minus_fit[o];
	//}
	//cout << check << " " << d_n * d_n + d_t * d_t << endl;


	delete[] z_minus_fit;
	return (d_t - PENALTY * d_n * d_n/ shita);
}


bool function_valuation(char *problem, char *scalarizing_function,
	double *original_objective, double *children_objective,
	bool improve_lambda2, double *lambda, double *lambda2,
	double *maxz, double *minz,
	int ob,
	double PENALTY, double PENALTY2, double shita, char *max_min){

	if (strcmp(max_min, "min") == 0){
		double y_g_t_e = 0.0;
		double cure = 0.0;

		if (improve_lambda2 == false){
			if (strcmp(scalarizing_function, "tchebycheff") == 0){
				y_g_t_e = g_t_e(children_objective, lambda, minz, ob);
				cure = g_t_e(original_objective, lambda, minz, ob);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "weightedsum") == 0){
				y_g_t_e = weight_sum(children_objective, lambda, ob);
				cure = weight_sum(original_objective, lambda, ob);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "pbi") == 0){
				y_g_t_e = p_b_i_min(children_objective, lambda, minz, ob, PENALTY);
				cure = p_b_i_min(original_objective, lambda, minz, ob, PENALTY);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "quadratic_pbi") == 0){
				y_g_t_e = quadratic_p_b_i_min(children_objective, lambda, minz, ob, PENALTY, shita);
				cure = quadratic_p_b_i_min(original_objective, lambda, minz, ob, PENALTY,shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "cubic_pbi") == 0){
				y_g_t_e = cubic_p_b_i_min(children_objective, lambda, minz, ob, PENALTY);
				cure = cubic_p_b_i_min(original_objective, lambda, minz, ob, PENALTY);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "double_pbi") == 0){
				y_g_t_e = p_b_i_min_double(children_objective, lambda, minz, ob, PENALTY, PENALTY2, shita);
				cure = p_b_i_min_double(original_objective, lambda, minz, ob, PENALTY, PENALTY2, shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "fix_double_pbi") == 0){
				y_g_t_e = p_b_i_min_double_normalize_fix(children_objective, lambda, minz, ob, PENALTY, PENALTY2, shita);
				cure = p_b_i_min_double_normalize_fix(original_objective, lambda, minz, ob, PENALTY, PENALTY2, shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "no_penalty_pbi") == 0){
				y_g_t_e = p_b_i_min_para0(children_objective, lambda, minz, ob, shita);
				cure = p_b_i_min_para0(original_objective, lambda, minz, ob, shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "inverted_pbi") == 0){
				y_g_t_e = inverted_p_b_i_min(children_objective, lambda, maxz, ob, PENALTY);
				cure = inverted_p_b_i_min(original_objective, lambda, maxz, ob, PENALTY);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "quadratic_inverted_pbi") == 0){
				y_g_t_e = quadratic_inverted_p_b_i_min(children_objective, lambda, maxz, ob, PENALTY,shita);
				cure = quadratic_inverted_p_b_i_min(original_objective, lambda, maxz, ob, PENALTY,shita);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "fix_double_inverted_pbi") == 0){
				y_g_t_e = inverted_p_b_i_min_double_normalize_fix(children_objective, lambda, maxz, ob, PENALTY, PENALTY2, shita);
				cure = inverted_p_b_i_min_double_normalize_fix(original_objective, lambda, maxz, ob, PENALTY, PENALTY2, shita);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "asf") == 0){
				y_g_t_e = a_g_t_e(children_objective, lambda, maxz, ob);
				cure = a_g_t_e(original_objective, lambda, maxz, ob);
				return y_g_t_e >= cure;
			}
			else{
				std::cout << "function is not available" << std::endl;
				exit(1);
			}
		}
		else{
			if (strcmp(scalarizing_function, "tchebycheff") == 0){
				y_g_t_e = g_t_e(children_objective, lambda2, minz, ob);
				cure = g_t_e(original_objective, lambda2, minz, ob);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "weightedsum") == 0){
				y_g_t_e = weight_sum(children_objective, lambda2, ob);
				cure = weight_sum(original_objective, lambda2, ob);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "pbi") == 0){
				y_g_t_e = p_b_i_min(children_objective, lambda2, minz, ob, PENALTY);
				cure = p_b_i_min(original_objective, lambda2, minz, ob, PENALTY);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "quadratic_pbi") == 0){
				y_g_t_e = quadratic_p_b_i_min(children_objective, lambda2, minz, ob, PENALTY,shita);
				cure = quadratic_p_b_i_min(original_objective, lambda2, minz, ob, PENALTY,shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "cubic_pbi") == 0){
				y_g_t_e = cubic_p_b_i_min(children_objective, lambda2, minz, ob, PENALTY);
				cure = cubic_p_b_i_min(original_objective, lambda2, minz, ob, PENALTY);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "double_pbi") == 0){
				y_g_t_e = p_b_i_min_double(children_objective, lambda2, minz, ob, PENALTY, PENALTY2, shita);
				cure = p_b_i_min_double(original_objective, lambda2, minz, ob, PENALTY, PENALTY2, shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "fix_double_pbi") == 0){
				y_g_t_e = p_b_i_min_double_normalize_fix(children_objective, lambda2, minz, ob, PENALTY, PENALTY2, shita);
				cure = p_b_i_min_double_normalize_fix(original_objective, lambda2, minz, ob, PENALTY, PENALTY2, shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "no_penalty_pbi") == 0){
				y_g_t_e = p_b_i_min_para0(children_objective, lambda2, minz, ob, shita);
				cure = p_b_i_min_para0(original_objective, lambda2, minz, ob, shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "inverted_pbi") == 0){
				y_g_t_e = inverted_p_b_i_min(children_objective, lambda2, maxz, ob, PENALTY);
				cure = inverted_p_b_i_min(original_objective, lambda2, maxz, ob, PENALTY);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "quadratic_inverted_pbi") == 0){
				y_g_t_e = quadratic_inverted_p_b_i_min(children_objective, lambda2, maxz, ob, PENALTY,shita);
				cure = quadratic_inverted_p_b_i_min(original_objective, lambda2, maxz, ob, PENALTY,shita);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "fix_double_inverted_pbi") == 0){
				y_g_t_e = inverted_p_b_i_min_double_normalize_fix(children_objective, lambda2, maxz, ob, PENALTY, PENALTY2, shita);
				cure = inverted_p_b_i_min_double_normalize_fix(original_objective, lambda2, maxz, ob, PENALTY, PENALTY2, shita);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "asf") == 0){
				y_g_t_e = a_g_t_e(children_objective, lambda2, maxz, ob);
				cure = a_g_t_e(original_objective, lambda2, maxz, ob);
				return y_g_t_e >= cure;
			}
			else{
				std::cout << "function is not available" << std::endl;
				exit(1);
			}
		}

	}
	else if (strcmp(max_min, "max") == 0){
		double y_g_t_e = 0.0;
		double cure = 0.0;
		if (improve_lambda2 == false){
			if (strcmp(scalarizing_function, "tchebycheff") == 0){
				y_g_t_e = g_t_e(children_objective, lambda, maxz, ob);
				cure = g_t_e(original_objective, lambda, maxz, ob);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "weightedsum") == 0){
				y_g_t_e = weight_sum(children_objective, lambda, ob);
				cure = weight_sum(original_objective, lambda, ob);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "pbi") == 0){
				y_g_t_e = p_b_i(children_objective, lambda, maxz, ob, PENALTY);
				cure = p_b_i(original_objective, lambda, maxz, ob, PENALTY);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "quadratic_pbi") == 0){
				y_g_t_e = quadratic_p_b_i(children_objective, lambda, maxz, ob, PENALTY,shita);
				cure = quadratic_p_b_i(original_objective, lambda, maxz, ob, PENALTY,shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "cubic_pbi") == 0){
				y_g_t_e = cubic_p_b_i(children_objective, lambda, maxz, ob, PENALTY);
				cure = cubic_p_b_i(original_objective, lambda, maxz, ob, PENALTY);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "double_pbi") == 0){
				y_g_t_e = p_b_i_double(children_objective, lambda, maxz, ob, PENALTY, PENALTY2, shita);
				cure = p_b_i_double(original_objective, lambda, maxz, ob, PENALTY, PENALTY2, shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "fix_double_pbi") == 0){
				y_g_t_e = p_b_i_double_normalize_fix(children_objective, lambda, maxz, ob, PENALTY, PENALTY2, shita);
				cure = p_b_i_double_normalize_fix(original_objective, lambda, maxz, ob, PENALTY, PENALTY2, shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "no_penalty_pbi") == 0){
				y_g_t_e = p_b_i_para0(children_objective, lambda, maxz, ob, shita);
				cure = p_b_i_para0(original_objective, lambda, maxz, ob, shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "inverted_pbi") == 0){
				y_g_t_e = inverted_p_b_i_max(children_objective, lambda, minz, ob, PENALTY);
				cure = inverted_p_b_i_max(original_objective, lambda, minz, ob, PENALTY);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "quadratic_inverted_pbi") == 0){
				y_g_t_e = quadratic_inverted_p_b_i_max(children_objective, lambda, minz, ob, PENALTY,shita);
				cure = quadratic_inverted_p_b_i_max(original_objective, lambda, minz, ob, PENALTY,shita);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "fix_double_inverted_pbi") == 0){
				y_g_t_e = inverted_p_b_i_double_normalize_fix(children_objective, lambda, minz, ob, PENALTY, PENALTY2, shita);
				cure = inverted_p_b_i_double_normalize_fix(original_objective, lambda, minz, ob, PENALTY, PENALTY2, shita);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "asf") == 0){
				y_g_t_e = a_g_t_e(children_objective, lambda, minz, ob);
				cure = a_g_t_e(original_objective, lambda, minz, ob);
				return y_g_t_e >= cure;
			}
			else{
				std::cout << "function is not available" << std::endl;
				exit(1);
			}
		}
		else{
			if (strcmp(scalarizing_function, "tchebycheff") == 0){
				y_g_t_e = g_t_e(children_objective, lambda2, maxz, ob);
				cure = g_t_e(original_objective, lambda2, maxz, ob);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "weightedsum") == 0){
				y_g_t_e = weight_sum(children_objective, lambda2, ob);
				cure = weight_sum(original_objective, lambda2, ob);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "pbi") == 0){
				y_g_t_e = p_b_i(children_objective, lambda2, maxz, ob, PENALTY);
				cure = p_b_i(original_objective, lambda2, maxz, ob, PENALTY);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "quadratic_pbi") == 0){
				y_g_t_e = quadratic_p_b_i(children_objective, lambda2, maxz, ob, PENALTY,shita);
				cure = quadratic_p_b_i(original_objective, lambda2, maxz, ob, PENALTY,shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "cubic_pbi") == 0){
				y_g_t_e = cubic_p_b_i(children_objective, lambda2, maxz, ob, PENALTY);
				cure = cubic_p_b_i(original_objective, lambda2, maxz, ob, PENALTY);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "double_pbi") == 0){
				y_g_t_e = p_b_i_double(children_objective, lambda2, maxz, ob, PENALTY, PENALTY2, shita);
				cure = p_b_i_double(original_objective, lambda2, maxz, ob, PENALTY, PENALTY2, shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "fix_double_pbi") == 0){
				y_g_t_e = p_b_i_double_normalize_fix(children_objective, lambda2, maxz, ob, PENALTY, PENALTY2, shita);
				cure = p_b_i_double_normalize_fix(original_objective, lambda2, maxz, ob, PENALTY, PENALTY2, shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "no_penalty_pbi") == 0){
				y_g_t_e = p_b_i_para0(children_objective, lambda2, maxz, ob, shita);
				cure = p_b_i_para0(original_objective, lambda2, maxz, ob, shita);
				return y_g_t_e <= cure;
			}
			else if (strcmp(scalarizing_function, "inverted_pbi") == 0){
				y_g_t_e = inverted_p_b_i_max(children_objective, lambda2, minz, ob, PENALTY);
				cure = inverted_p_b_i_max(original_objective, lambda2, minz, ob, PENALTY);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "quadratic_inverted_pbi") == 0){
				y_g_t_e = quadratic_inverted_p_b_i_max(children_objective, lambda2, minz, ob, PENALTY,shita);
				cure = quadratic_inverted_p_b_i_max(original_objective, lambda2, minz, ob, PENALTY,shita);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "fix_double_inverted_pbi") == 0){
				y_g_t_e = inverted_p_b_i_double_normalize_fix(children_objective, lambda2, minz, ob, PENALTY, PENALTY2, shita);
				cure = inverted_p_b_i_double_normalize_fix(original_objective, lambda2, minz, ob, PENALTY, PENALTY2, shita);
				return y_g_t_e >= cure;
			}
			else if (strcmp(scalarizing_function, "asf") == 0){
				y_g_t_e = a_g_t_e(children_objective, lambda2, minz, ob);
				cure = a_g_t_e(original_objective, lambda2, minz, ob);
				return y_g_t_e >= cure;
			}
			else{
				std::cout << "function is not available" << std::endl;
				exit(1);
			}
		}
	}
	else{
		std::cout << "problem is not available" << std::endl;
		exit(1);
	}

}