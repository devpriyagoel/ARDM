#include<bits/stdc++.h>
using namespace std;
#define numRange 1
const double NEARZERO = 1.0e-10;  
void transpose(double** A, double** B, int n);
void matrixTimesVector( double** A, double* V, double* ans, int n);
void conjugateGradientSolver( double** A, double* B, double* X, int n );
void grad(double* v,double **a,double*b, double* ans, int n );  
long double matrixnorm(double ** A, int n);
long double arrnorm(double* x, int n);
double innerProduct( double* U, double* V, int n );
void vectorCombination( double a, double* U, double b, double* V, double* W, int n );
int main()
{	
	clock_t start, end; 
	
    int i, j, k;
	int n=22, m=10000;
	long double error = (long double)(sqrt(n));
	double** a = new double*[n];
	for(i=0; i<n; i++){
		a[i]= new double[n];
	}
	double* b = new double[n];
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			cin >> a[i][j];
		}
	}
	for(i=0; i<n; i++){
		cin >> b[i];
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
	alpha = 1/(alpha);
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
		graduc[i] = -b[i];
	}
	gc = arrnorm(graduc, n);
    for(k=0;cnt<m;k++){	
		cnt++;
		/*cout << "k--" << k << endl;*/
		if(k==0){
			// alpha = 0.2;
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
		}else{
			// alpha = (2.0)*(1/(gc*gc));
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
			//cout << "restarting k--" << k << " gc gp " << gc << " "<<gp << endl;
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
	double* B = new double[n];
	matrixTimesVector(a, uc, B, n);
	for(i=0; i<n; i++){
		B[i] = b[i]-B[i];
	}
	conjugateGradientSolver(a, B, uc, n);
	end=clock();
	cout << "iter " << cnt << endl;
	long double time_take = (long double)(end-start)/(long double)(CLOCKS_PER_SEC);
	cout << "Time taken " << time_take << " sec"<<endl;
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
void vectorCombination( double a, double* U, double b, double* V, double* W, int n ){        // Linear combination of vectors
	int j;
	for ( j = 0; j < n; j++ ) W[j] = a * U[j] + b * V[j];
}



void conjugateGradientSolver( double** A, double* B, double* X, int n ){
	double TOLERANCE = 1;
	double* P = new double[n];
	double* R = new double[n];
	double* Rold = new double[n];
	double* AP = new double[n];
	int i, k=0;
	for(i=0; i<n; i++){
		P[i]=B[i];
		R[i]=B[i];
	}

	while (k<n){
		for(i=0; i<n; i++){
			Rold[i] = R[i];
		}
		// for(i=0; i<n; i++){
		// 	cout << X[i] << " ";
		// }
		//cout << endl;
		matrixTimesVector(A, P, AP, n);
		double alpha = innerProduct( R, R, n ) / max( innerProduct( P, AP, n ), NEARZERO );
		//cout<<"alpha   "<<alpha<<endl;
		vectorCombination( 1.0, X, alpha, P, X, n );            // Next estimate of solution
		vectorCombination( 1.0, R, -alpha, AP, R, n );          // Residual

		if ( arrnorm( R, n ) < TOLERANCE ) break;             // Convergence test

		double beta = innerProduct( R, R, n ) / max( innerProduct( Rold, Rold, n ), NEARZERO );
		vectorCombination( 1.0, R, beta, P,  P, n);             // Next gradient
		k++;
	}
	cout << "iter --conju  " << k << endl;
	delete [] P;
	delete [] R;
	delete [] Rold;
	delete [] AP;

}
