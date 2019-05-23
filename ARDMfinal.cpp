#include<bits/stdc++.h>
using namespace std;
#define numRange 1

void transpose(double** A, double** B, int n);
void matrixTimesVector( double** A, double* V, double* ans, int n);
void grad(double* v,double **a,double*b, double* ans, int n );  
long double matrixnorm(double ** A, int n);
long double arrnorm(double* x, int n);
double innerProduct( double* U, double* V, int n );
int main(int argc, char** argv)
{
    if(argc < 4)
	{
		cout << "Format is ./testNest matrixSize _numItrer_ _error_" << endl;
		exit(1);
	}

	int i, j, k;
	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	double error = atof(argv[3]);
	clock_t start, end;
	cout << "Size of n x n Square matrix is : " << n; 
	srand((unsigned int)time(NULL));
	double** a = new double*[n];
	double** a1 = new double*[n];
	double** a2 = new double*[n];
	for(i=0; i<n; i++){
		a[i]= new double[n];
		a1[i]= new double[n];
		a2[i]= new double[n];
	}
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			a1[i][j] = ((double) rand() / (RAND_MAX));
		}
	}
	transpose(a1, a2, n);
	for( i=0; i<n; i++){
		for(j=0; j< n; j++){
			a[i][j] = (a1[i][j]+a2[i][j])/2.0;
			if(i==j)a[i][j]+=n;
		}
	}
	double* b = new double[n];
	for (i = 0; i < n; i++){
		b[i] = ((double) rand() / (RAND_MAX));
	}
	cout << "Generated Matrix a is " << endl;
    for(i=0; i<n; ++i){
        for(j=0; j<n; ++j){
			cout << a[i][j] << " ";
        }
        cout << endl;
	}

	cout << "Generated Matrix b is " << endl;
    for(i=0; i<n; ++i){
		cout << b[i] << " ";
	}
    cout << endl;
    start = clock();
    //================================================
	long double alpha,beta;
	alpha = matrixnorm(a, n);
	cout << "alpha--" << alpha << endl;
	if(alpha==0.0){
		cout << "Zero Matrix " << endl;
		return 0;
	}
	alpha = 1/(alpha);
	// cout << "alpha = " << alpha << endl;
	//=================================================
	double *uc = new double[n]();
	double *up = new double[n];
	double *v = new double[n];
	double *ans = new double[n];	
	double *graduc = new double[n];
	long double gc=0, gp=0;
	int cnt=0;
	bool flag=false;
	for(i=0; i<n; i++){
		graduc[i]=-b[i];
	}
	gc = arrnorm(graduc, n);
    for(k=0;cnt<m;k++){	
		cnt++;
	//	cout << "k--" << k << endl;
		if(k==0){
			for(i=0;i<n;i++){
				v[i] = uc[i] - alpha*graduc[i];
			}
			for(i=0; i<n; i++){
				up[i]=uc[i];
			}
			grad(v, a, b, ans, n);
			for(i=0; i<n; i++){
				uc[i]= v[i]-alpha*ans[i];
			}
			//gc = arrnorm(graduc, n);
		}else{
			// gp = gc;
			// gc = arrnorm(graduc, n);
			beta = gc/gp;
			// cout << beta << " this is beta" << endl;
			for(i=0;i<n;i++){
				v[i] = uc[i] + beta*(uc[i]-up[i]) - alpha*(1+beta)*graduc[i];
			}
			for(i=0; i<n; i++){
				up[i]=uc[i];
			}
			grad(v, a, b, ans, n);
			for(i=0; i<n; i++){
				uc[i]= v[i]-alpha*ans[i];
			}

		}
		grad(uc, a, b, graduc, n);
		gp = gc;
		gc = arrnorm(graduc, n);
		if(gc>gp&&k!=0){	
		//	cout << "restarting k--" << k << " gc gp " << gc << " "<<gp << endl;
			for(i=0; i<n; i++){
				uc[i] = up[i];
			}
			k=-1;	
		}
		if(gc<error&&k!=0){
			//cout<< "exiting error --" << error << " gc--" << gc << endl;
			flag=true;
			break;
		}
		// cout << "uc---->" << endl;
		// for(i=0; i<n; i++){
		// 	cout << uc[i] << " ";
		// }
		// cout << endl;

	}
	end=clock();
	if(flag){
		cout << "Converged in given number of iterations " <<  cnt <<endl;
		long double time_take = (long double)(end-start)/(long double)(CLOCKS_PER_SEC);
		cout << "Time taken " << time_take << " sec"<<endl;
	}else{
		cout << "not Converged" << endl;
	}
	cout << "X" << endl;
	for(i=0; i<n; i++){
		cout << uc[i] << " ";
	}
	cout << endl;
	cout << "AX" << endl;
	matrixTimesVector(a, uc, v, n);
	for(i=0; i<n; i++){
		cout << v[i] << " ";
	}
	cout << endl;
    return 0;
}

void transpose(double** A, double** B, int n){
	int i, j;
	for (i = 0; i <n ; i++){
		for (j = 0; j < n; j++){
			B[i][j] = A[j][i];
		}
	}
}
double innerProduct( double* U, double* V, int n ){
	double ans=0.0;
	int i;
	for(i=0; i<n; i++){
		ans+=U[i]*V[i];
	}
	return ans;
}

void matrixTimesVector( double** A, double* V, double* ans, int n){
	int i;
	for ( i = 0; i < n; i++ ) ans[i] = innerProduct( A[i], V, n );
}

void grad(double* v,double **a,double*b, double* ans, int n ){
	matrixTimesVector(a, v, ans, n);
	for(int i=0; i<n; i++){
		ans[i]-=b[i];
	}
}

long double arrnorm(double* x, int n){
	long double ans=0.0;
	int i;
	for(i=0; i<n; i++){
		ans+=(long double)(x[i])*(x[i]);
	}
	ans = sqrt(ans);
	return ans;
}
long double matrixnorm(double ** A, int n){
	long double ans=0.0;
	int i, j;
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			ans+=(long double)A[i][j]*A[i][j];
		}
	}
	return sqrt(ans);
}