#include<bits/stdc++.h>
using namespace std;
#define numRange 1
const double NEARZERO = 1.0e-10;  
void matrixTimesVector( double* A, double* V, double* ans, int n,int*rowind,int*colind,int nnz);
void conjugateGradientSolver( double* A, double* B, double* X, int n,int*rowind,int*colind,int nnz );
void grad(double* v,double *A,double*b, double* ans, int n,int*rowind,int*colind,int nnz );  
long double matrixnorm(double * A, int nnz);
long double arrnorm(double* x, int n);
double innerProduct( double* U, double* V, int n );
void vectorCombination( double a, double* U, double b, double* V, double* W, int n );
int main()
{	
	clock_t start, end; 
	
    int i, j, k;
	int n, m=200;
	long double error = (long double)n;
	int nnz;
	cin>>n>>nnz;
	double*a=new double[nnz];
	double*b=new double[n];
	for(i=0;i<nnz;i++)
	{
		cin>>a[i];
	}
	int*rowind=new int[nnz];
	for(i=0;i<nnz;i++)
	{
		cin>>rowind[i];
	}
	int*colind=new int[n+1];
	for(i=0;i<n+1;i++)
	{
		cin>>colind[i];
	}

	for(i=0; i<n; i++){
		cin >> b[i];
	}
	cout << "Generated Matrix a is " << endl;
   

	cout << "Generated Matrix b is " << endl;
    for(i=0; i<n; ++i){
		cout << b[i] << " ";
	}
    cout << endl;
    start = clock();
    double*uc = new double[n]();
    double*v = new double[n];
   
    //================================================
	
	conjugateGradientSolver(a,b,uc,n,rowind,colind,nnz );
	end=clock();
	// cout << "iter " << cnt << endl;
	long double time_take = (long double)(end-start)/(long double)(CLOCKS_PER_SEC);
	cout << "Time taken " << time_take << " sec"<<endl;
	cout << "X" << endl;
	for(i=0; i<n; i++){
		cout << uc[i] << " ";
	}
	cout << endl;
	cout << "AX" << endl;
	matrixTimesVector( a, uc, v,n,rowind,colind,nnz);
	for(i=0; i<n; i++){
		cout << v[i] << " ";
	}
	cout << endl;
    return 0;
}


double innerProduct( double* U, double* V, int n ){
	double ans=0.0;
	int i;
	for(i=0; i<n; i++){
		ans+=U[i]*V[i];
	}
	return ans;
}

void matrixTimesVector( double* A, double* V, double* ans, int n,int*rowind,int*colind,int nnz){
	int i, j=0;
	int k;
	int *temp = new int[nnz];
	for(i=0; i<n; i++){
		ans[i]=0;
	}
	for(i=0; i<n; i++){
		k = colind[i+1]-colind[i];
		while(k--){
			temp[j] = i;
			j++;
		}
	}
	for(i=0;i<nnz;i++)
	{
		ans[rowind[i]]+= A[i]*V[temp[i]];
	}
	delete [] temp;
}

void grad(double* v,double *A,double*b, double* ans, int n,int*rowind,int*colind,int nnz ){
	matrixTimesVector( A, v, ans,  n,rowind,colind, nnz);
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
long double matrixnorm(double * A, int nnz){
	long double ans=0.0;
	int j;

	for(j=0; j<nnz; j++){
		ans+=(long double)A[j]*A[j];
		
	}
	return sqrt(ans);
}
void vectorCombination( double a, double* U, double b, double* V, double* W, int n ){        // Linear combination of vectors
	int j;
	for ( j = 0; j < n; j++ ) W[j] = a * U[j] + b * V[j];
}



void conjugateGradientSolver( double* A, double* B, double* X, int n,int*rowind,int*colind,int nnz ){
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
		matrixTimesVector( A, P, AP,  n,rowind,colind, nnz);
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
