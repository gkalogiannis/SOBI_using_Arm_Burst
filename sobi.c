#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "defs_and_types.h"
#include "lapacke.h"

int i,j,k,l,o,q,r,s; //indices
float ***X,**v,*w, **X2; //Source Signal x(t)=y(t)+n(t)=A*s(t)+n(t)
float *M; //various arrays used throughout the programm
int m=2, N=3, ntrials=2; //X dimensions
void allocate(); //Allocate memory for used arrays
void print(float ***m, int N, int M, int K);
void mean();	 
//svd functions
static double PYTHAG(double a, double b);
int dsvd(float **a, int m, int n, float *w, float **v);
//inverse functions
void minor(float **b,float **a,int i,int n);
float det(float **a,int n);
void transpose(float **c,float **d,float n,float m, float det);
void cofactor(float **a,float **d,float n,float determinte);
void inverse(float **a,float **d,int n,float det);
//Matrix Multiplication
void multiply(float **a, float **b,float **c, int Na, int Ma, int Mb, int k, int ntrials);
//Quicksort
void swap(float* a, float* b);
int partition (float *arr, float *ind, int low, int high);
void quickSort(float *arr, float *ind, int low, int high);

int main(){
	
	//Trial X=[2,4,3];
	allocate();
	//Test values for X
	for (i=0; i<m; i++){
		for (j=0; j<N; j++){
			for (k=0; k<ntrials; k++){
				X[i][j][k]=(float)rand() / (float)RAND_MAX;
			}
		}
	}
	//Print X
	print(X,m,N,ntrials);
	
	/*//Make the data zero mean
	mean();
	
	//kron
	for (i=0; i<m; i++){
		for (j=0; j<N; j++){
			for (k=0; k<ntrials; k++){
				X[i][j][k]-=M[i];
			}
		}
	}
	print(X,m,N,ntrials);*/
	
	//3D to 2D
	for (i=0; i<m; i++){
		for (j=0; j<N; j++){
			for (k=0; k<ntrials; k++){
				X2[j+k*N][i] = X[i][j][k];
			}
		}
	}
	//Print X2
	printf("\n\n X2=\n");
	for (i=0; i<N * ntrials; i++){
		for (j=0; j<m; j++){
						printf("%03f ", X2[i][j]);
		}
			printf("\n");
	}
	//Singular Value Decomposition (SVD) of X2 (2D version of X)
	/*dsvd(X2, N * ntrials, m, w, v);
	
	//Print V
	printf("\n\n v=\n");
	for (i=0; i<m ; i++){
		for (j=0; j<m; j++){
						printf("%03f ", v[i][j]);
		}
			printf("\n");
	}
	
	//Print w
	printf("\n\n w=\n");
	for (i=0; i<m ; i++){
		printf("%03f ", w[i]);
	}
	printf("\n");
	
	//Turn w 2D (w2)
	//Memory allocation
	float ** w2 = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		w2[i] = (float *) malloc(m * sizeof(float ));
	}
	//Turn w into a diagonal 2d matrix
	for (i=0; i<m ; i++){
		for (j=0; j<m; j++){
						if (i==j){
							w2[i][j]=w[i];
						}
						else{
							w2[i][j]=0;
						}
		}
	}
	//Find the inverse of w
	float ** Winv = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++){
		Winv[i] = (float *) malloc(m * sizeof(float ));
	}
	float deter = (float) det(w2,m);
	inverse(w2,Winv,m,deter);
	//Print winv
	printf("\n\n w^-1=\n");
	for (i=0; i<m ; i++){
		for (j=0; j<m; j++){
						printf("%03f ", Winv[i][j]);
		}
			printf("\n");
	}
	
	//Transpose v
	transpose(v,v,m,m,1);
	
	//Q=Sinv*v'
	float ** Q = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		Q[i] = (float *) malloc(m * sizeof(float ));
	}
			//Multiplication
	multiply(Winv, v, Q, m, m, m, 0, 0);*/
	
	//Transpose X2
	float **X2T = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		X2T[i] = (float *) malloc(N*ntrials * sizeof(float ));
	}
	transpose(X2,X2T,N*ntrials,m, 1);
	
	//X2=Q*X2
	/*float **X2T2 = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		X2T2[i] = (float *) malloc(N*ntrials * sizeof(float ));
	}
	
	multiply(Q, X2T, X2T2, m, m, N*ntrials, 0, 0);*/
	
	
	//Print X2T
	printf("\n\n X2T=\n");
	for (i=0; i<m; i++){
		for (j=0; j<N * ntrials; j++){
						printf("%03f ", X2T[i][j]);
		}
			printf("\n");
	}
	
	
	
	printf("\n\n Finding Corellation matrixes\n");
	//Estimate the correlation matrices
	
	int p=1; //number of correlation matrices to be diagonalised (placeholder)
	//Find p = number of correlation matrices to be diagonalised
	//Rxp stores the correlation matrices
	float **Rxp = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		Rxp[i] = (float *) malloc(N*ntrials * sizeof(float ));
	}
	float **Rxp2 = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		Rxp2[i] = (float *) malloc(N*ntrials * sizeof(float ));
	}

	int c=0; //k in matlab code
	float frobenius, sq_sum=0; //Used for the calculation of the Frobenius Norm of Rxp
	
	int pm =p*m;
	float **MM = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		MM[i] = (float *) malloc(pm * sizeof(float ));
	}
	//Imaginary part of MM
	float **MM_im = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		MM_im[i] = (float *) malloc(pm * sizeof(float ));
	}
	float **MM2 = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		MM2[i] = (float *) malloc(pm * sizeof(float ));
	}
	float **MM2_im = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		MM2_im[i] = (float *) malloc(pm * sizeof(float ));
	}
	for (i=0; i<pm; i+=m){
		c+=1;
		//div=/(N-c+2)/ntrials;
		for (j=0; j<ntrials; j++){
			if (j==0){
				multiply(X2T, X2, Rxp, m, N, m, c, j);	
				//Divide by (N-c+1)/ntrials
				for (k=0; k<m; k++){
					for (l=0; l<m; l++){
						Rxp[k][l]=Rxp[k][l] /(N-c)/ntrials;
					}
				}
				continue;
			}
			multiply(X2T, X2, Rxp2, m, N, m, c, j);
			//Add Rxp and Rxp2 element wise
			for (k=0; k<m; k++){
				for (l=0; l<m; l++){
					Rxp[k][l]+=Rxp2[k][l] / (N-c)/ntrials;
				}
			}
		}
	
		//Frobenius norm of Rxp
			//Find square sum of all it's elements
		for (k=0; k<m; k++){
			for (l=0; l<m; l++){
				sq_sum+=Rxp[k][l]*Rxp[k][l];
			}
		}
		//Find the square root of sq_sum
		frobenius=sqrt(sq_sum);
		printf("\n\nfrobemius=%f",frobenius);
		for (k=0; k<m; k++){
			for (l=i; l<(i+m); l++){
				MM[k][l]=Rxp[k][l-i]*frobenius;
			}
		}
	}	

	//Print Rxp
	printf("\n\n Rxp=\n");
	for (i=0; i<m; i++){
		for (j=0; j<m; j++){
						printf("%03f ", Rxp[i][j]);
		}
			printf("\n");
	}
	//Print M
	printf("\n\n M=\n");
	for (i=0; i<m; i++){
		for (j=0; j<pm; j++){
						printf("%03f ", MM[i][j]);
		}
			printf("\n");
	}
	//free oloi oi pinakes ektos apo MM
	
	//Perform Joint Diagonalization
	float epsil=1 / sqrt(N) / 100;
	int encore=1;
	
	//V is an mxm !!complex!! identity matrix
	float **V = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		V[i] = (float *) malloc(m * sizeof(float ));
	}
	float **V_im = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		V_im[i] = (float *) malloc(m * sizeof(float ));
	}
	float **V2 = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		V2[i] = (float *) malloc(m * sizeof(float ));
	}
	float **V2_im = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		V2_im[i] = (float *) malloc(m * sizeof(float ));
	}
	
	for (i=0; i<m; i++){
		for (j=0; j<m; j++){
			if (i==j) V[i][j]=1;
		}
	}
	//g is used to store the result of the givens rotation
	float **g = (float **) malloc(3 * sizeof(float *));
	for(i=0; i<3; i++)
	{
		g[i] = (float *) malloc(p * sizeof(float ));
	}
	//Various matrices and variables for use inside the while loop
	float *DD = (float *) malloc(3 * sizeof(float));
	float *angles = (float *) malloc(3 * sizeof(float));
	float cc, sr_re, sr_im, sc_re, sc_im;
	int oui;
	
	while (encore){
		encore=0;
		for (i=0; i<m-1; i++){
			for (j=p; j<m; j++){
				//Perform Givens Rotation (g)
				//Fist Row
				for (k=0; k<p; k++){
						//Real and imaginary
						g[0][k]=MM[i][i+k*m] - MM[j][j+k*m];
						g[0][k]=MM_im[i][i+k*m] - MM_im[j][j+k*m];
					
						g[1][k]=MM[i][j+k*m] + MM[j][i+k*m];
						g[1][k]=MM_im[i][j+k*m] + MM_im[j][i+k*m];
					
						//!!!!!!!!!!!!!!!!!!!!!!!!
						g[2][k]=MM[j][i+k*m] - MM[i][j+k*m];				
						g[2][k]=MM_im[j][i+k*m] - MM_im[i][j+k*m];
						//!!!!!!!!!!!!!!!!!!!!!!!!!!
				}	
				//[vcp,D] = eig(real(g*g'))
				//Test values
				//vcp = the right eigenvector
				float vcp[3][3]={
											{-0.9834, 0, 0.1812},
											{-0.1812, 0, 0.9834},
											{0, 1, 1}
										};
				//D=diagonal matrix of eigenvalues
				float D[3][3]={
											{0, 0, 0},
											{0, 0.00036, 0},
											{0,0, 1.6525}
										};
				//Keep only the diagonal of D
				for (l=0; l<3; l++){
					DD[l]=D[l][l];
				}
				float ind[3]={0,1,2}; //Keeps the indexes of DD after DD is sorted
				//Sort Eigenvalues using quicksort
				quickSort(DD, ind, 0, 2);
				
				for (l=0; l<3; l++){
					angles[l]=vcp[l][(int)ind[2]];
				}
				//Sign
				if (angles[0]<0){
					for (l=0; l<3; l++){
						angles[l]*=-1;
					}
				}
				cc = sqrt(0.5+angles[0]/2);
				//sr=0.5*(angles(2)-j*angles(3))/c;
				sr_re=0.5*(angles[1])/cc;
				sr_im=0.5*(-angles[2])/cc;
				//sc=conjugate of sr
				sc_re = sr_re;
				sc_im = (-1)*sr_im;
				
				oui = sqrt(sr_re*sr_re + sr_im * sr_im) > epsil;
				encore = encore || oui;
				if (oui){  //Update the V and M matrices
					//Store previous M
					for (l=0; l<m; l++){
						for (k=0; k<pm; k++){
							MM2[l][k] = MM[l][k];
							MM2_im[l][k] = MM_im[l][k];
						}
					}
					for (l=0; l<m; l++){
						for (k=0; k<p; k++){
								//Real part
								MM[l][i+k*m] = cc * MM2[l][i+k*m] + sr_re * MM2[l][j+k*m] - sr_im * MM2_im[l][j+k*m]; 
								//Imaginary part
								MM_im[l][i+k*m] = cc * MM2_im[l][i+k*m] + sr_re * MM2_im[l][j+k*m] + sr_im * MM2[l][j+k*m];
						}
					}
					for (l=0; l<m; l++){
						for (k=0; k<p; k++){
							//Real part
							MM[l][j+k*m] = cc * MM2[l][j+k*m] - sc_re * MM2[l][i+k*m] + sc_im * MM2_im[l][i+k*m];
							//Imaginary part
							MM_im[l][j+k*m] = cc * MM2_im[l][j+k*m] - sc_re * MM2_im[l][i+k*m] - sc_im * MM2[l][i+k*m];
						}
					}
					//Print M
					printf("\n\n M1=\n");
					for (l=0; l<m; l++){
						for (k=0; k<pm; k++){
										printf("|%03f %03fi| ", MM[l][k], MM_im[l][k]);
						}
							printf("\n");
					}
					//Store previous M (re and im parts) (again)
					for (l=0; l<m; l++){
						for (k=0; k<pm; k++){
							MM2[l][k] = MM[l][k];
							MM2_im[l][k] = MM_im[l][k];
						}
					}
						for (k=0; k<pm; k++){
							//Real part
							MM[i][k] = cc * MM2[i][k] + sc_re * MM2[j][k] - sc_im * MM2_im[j][k];
							//Imaginary part
							MM_im[i][k] = cc * MM2_im[i][k] + sc_re * MM2_im[j][k] + sc_im * MM2[j][k];
						}
						
						for (k=0; k<pm; k++){
							//Real part
							MM[j][k] = cc * MM2[j][k] - sr_re * MM2[i][k] + sr_im * MM2_im[i][k];
							//Imaginary part
							MM_im[j][k] = cc * MM2_im[j][k] - sr_re * MM2_im[i][k] - sr_im * MM2[i][k];
						}
						
					//Print M
					printf("\n\n M2=\n");
					for (l=0; l<m; l++){
						for (k=0; k<pm; k++){
										printf("|%03f %03fi| ", MM[l][k], MM_im[l][k]);
						}
							printf("\n");
					}
					
					//V
					//Vtemp
					for (l=0; l<m; l++){
						for (k=0; k<pm; k++){
							V2[l][k] = V[l][k];
							V2_im[l][k] = V_im[l][k];
						}
					}
					for (k=0; k<m; k++){
						//Real part
						V[k][i] = cc * V[k][i] + sr_re * V[k][j] - sr_im * V_im[k][j];
						//Imaginary part
						V_im[k][i] = cc * V_im[k][i] + sr_re * V_im[k][j] + sr_im * V[k][j];
					}
					for (k=0; k<pm; k++){
						//Real part
						V[k][j] = cc * V[k][j] - sc_re * V2[k][i] + sc_im * V2_im[k][i];
						//Imaginary part
						V_im[k][j] = cc * V_im[k][j] - sc_re * V2_im[k][i] - sc_im * V2[k][i];
					}
					//Print V
					printf("\n\n V=\n");
					for (l=0; l<m; l++){
						for (k=0; k<pm; k++){
										printf("|%03f %03fi| ", V[l][k], V_im[l][k]);
						}
							printf("\n");
					}

					printf("\n\n sr = %03f j%03f", sr_re, sr_im);
					printf("\n\n sc = %03f j%03f", sc_re, sc_im);
					printf("\n\n Abs(sr) = %03f", sqrt(sr_re*sr_re + sr_im * sr_im));
					printf("\n\n epsil = %03f", epsil);
					printf("\n\n oui = %d", oui);
					printf("\n\n END");
					return 0;
				}
			}
		}
	}
		
		

	printf("\n\n sr = %03f j%03f", sr_re, sr_im);
	printf("\n\n sc = %03f j%03f", sc_re, sc_im);
	printf("\n\n Abs(sr) = %03f", sqrt(sr_re*sr_re + sr_im * sr_im));
	printf("\n\n epsil = %03f", epsil);
	printf("\n\n oui = %d", oui);
	printf("\n\n END");
}

void allocate(){
	printf("Allocating...\n");
	//3D array allocation
	X = (float ***) malloc(m * sizeof(float **));
	for(i=0; i<m; i++)
	{
		X[i] = (float **) malloc(N * sizeof(float *));
		for (j=0; j<N;j++)
		{
			X[i][j] = (float *) malloc(ntrials* sizeof(float));
		}
	}
	//2d version of X
	X2 = (float **) malloc(N * ntrials * sizeof(float *));
	for(i=0; i<N*ntrials; i++)
	{
		X2[i] = (float *) malloc(m * sizeof(float ));
	}
	//M is used in the "mean" function
	M = (float *) malloc(m * sizeof(float));
	//v is the output of svd for pre whitening
	v = (float **) malloc(m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		v[i] = (float *) malloc(m * sizeof(float ));
	}
	w = (float *) malloc(m * sizeof(float ));
	
}

void print(float ***P, int N, int M, int K){
	for (i=0; i<K; i++){
		for (j=0; j<N; j++){
			for (k=0; k<M; k++){
						printf("%03f ", P[j][k][i]);
			}
			printf("\n");
		}
		printf("\n\n");
    
	}
}

void mean(){
	float sum=0;
	
	for (i=0; i<m; i++){
		for (k=0; k<ntrials; k++){
			for (j=0; j<N; j++){
				sum+=X[i][j][k];
			}
		}
		M[i]=sum/(N*ntrials);
		sum=0;
	}
}


/*--------------SVD Code--------------*/
static double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}


int dsvd(float **a, int m, int n, float *w, float **v)
{
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) 
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }
  
    rv1 = (double *)malloc((unsigned int) n*sizeof(double));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) 
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
            for (k = i; k < m; k++) 
                scale += fabs((double)a[k][i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k][i] = (float)((double)a[k][i]/scale);
                    s += ((double)a[k][i] * (double)a[k][i]);
                }
                f = (double)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (float)(f - g);
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += ((double)a[k][i] * (double)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k][j] += (float)(f * (double)a[k][i]);
                    }
                }
                for (k = i; k < m; k++) 
                    a[k][i] = (float)((double)a[k][i]*scale);
            }
        }
        w[i] = (float)(scale * g);
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs((double)a[i][k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i][k] = (float)((double)a[i][k]/scale);
                    s += ((double)a[i][k] * (double)a[i][k]);
                }
                f = (double)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (float)(f - g);
                for (k = l; k < n; k++) 
                    rv1[k] = (double)a[i][k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += ((double)a[j][k] * (double)a[i][k]);
                        for (k = l; k < n; k++) 
                            a[j][k] += (float)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++) 
                    a[i][k] = (float)((double)a[i][k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }
  
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++)
                    v[j][i] = (float)(((double)a[i][j] / (double)a[i][l]) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += ((double)a[i][k] * (double)v[k][j]);
                    for (k = l; k < n; k++) 
                        v[k][j] += (float)(s * (double)v[k][i]);
                }
            }
            for (j = l; j < n; j++) 
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
  
    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = (double)w[i];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i][j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += ((double)a[k][i] * (double)a[k][j]);
                    f = (s / (double)a[i][i]) * g;
                    for (k = i; k < m; k++) 
                        a[k][j] += (float)(f * (double)a[k][i]);
                }
            }
            for (j = i; j < m; j++) 
                a[j][i] = (float)((double)a[j][i]*g);
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm]) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = (double)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (float)h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                        {
                            y = (double)a[j][nm];
                            z = (double)a[j][i];
                            a[j][nm] = (float)(y * c + z * s);
                            a[j][i] = (float)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k];
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {              /* make singular value nonnegative */
                    w[k] = (float)(-z);
                    for (j = 0; j < n; j++) 
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                fprintf(stderr, "No convergence after 30,000! iterations \n");
                return(0);
            }
    
            /* shift from bottom 2 x 2 minor */
            x = (double)w[l];
            nm = k - 1;
            y = (double)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) 
                {
                    x = (double)v[jj][j];
                    z = (double)v[jj][i];
                    v[jj][j] = (float)(x * c + z * s);
                    v[jj][i] = (float)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (float)z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = (double)a[jj][j];
                    z = (double)a[jj][i];
                    a[jj][j] = (float)(y * c + z * s);
                    a[jj][i] = (float)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (float)x;
        }
    }
    free((void*) rv1);
    return(1);
}
/*xxxxxxxxxxxxx End of SVD Code xxxxxxxxxxxxx*/

/*--------------Matrix Inversion Code--------------*/
//---------------------------------------------------
//	calculate minor of matrix OR build new matrix : k-had = minor
void minor(float **b,float **a,int i,int n){
	int j,l,h=0,k=0;
	for(l=1;l<n;l++)
		for( j=0;j<n;j++){
			if(j == i)
				continue;
			b[h][k] = a[l][j];
			k++;
			if(k == (n-1)){
				h++;
				k=0;
			}
		}
}// end function

//---------------------------------------------------
//	calculate determinte of matrix
float det(float **a,int n){
	int i;
	float **b,sum=0; 
	//Allocate memory for B
	b = (float **) malloc(n * sizeof(float *));
	for(i=0; i<n; i++)
	{
		b[i] = (float *) malloc(n * sizeof(float ));
	}
	if (n == 1)
return a[0][0];
	else if(n == 2)
return (a[0][0]*a[1][1]-a[0][1]*a[1][0]);
	else
		for(i=0;i<n;i++){
			minor(b,a,i,n);	// read function
			sum = (float) (sum+a[0][i]*pow(-1,i)*det(b,(n-1)));	// read function	// sum = determinte matrix
		}
return sum;
}// end function

//---------------------------------------------------
//	calculate transpose of matrix
void transpose(float **c,float **d,float n, float m, float det){
	int i,j;
	
	//Allocate space for b
	float **b = (float **) malloc( m * sizeof(float *));
	for(i=0; i<m; i++)
	{
		b[i] = (float *) malloc(n * sizeof(float ));
	}
	
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
			b[i][j] = c[j][i];
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
			d[i][j] = b[i][j]/det;	// array d[][] = inverse matrix
	free(b);
}// end function

//---------------------------------------------------
//	calculate cofactor of matrix
void cofactor(float **a,float **d,float n,float determinte){
	float **b,**c;
	//Allocate memory for b and c
	b = (float **) malloc(n * sizeof(float *));
	for(i=0; i<n; i++)
	{
		b[i] = (float *) malloc(n * sizeof(float ));
	}
	c = (float **) malloc(n * sizeof(float *));
	for(i=0; i<n; i++)
	{
		c[i] = (float *) malloc(n * sizeof(float ));
	}
	int l,h,m,k,i,j;
	for (h=0;h<n;h++)
		for (l=0;l<n;l++){
			m=0;
			k=0;
			for (i=0;i<n;i++)
				for (j=0;j<n;j++)
					if (i != h && j != l){
						b[m][k]=a[i][j];
						if (k<(n-2))
							k++;
						else{
							k=0;
							m++;
						}
					}
			c[h][l] = pow(-1,(h+l))*det(b,(n-1));	// c = cofactor Matrix
		}
	transpose(c,d,n,n,determinte);	// read function
}// end function

//---------------------------------------------------
//	calculate inverse of matrix
void inverse(float **a,float **d,int n,float det){
	if(det == 0)
		printf("\nInverse of Entered Matrix is not possible\n");
	else if(n == 1)
		d[0][0] = 1;
	else
		cofactor(a,d,n,det);	// read function
}// end function

//Matrix multiplication
void multiply(float **a, float **b,float **c, int Na, int Ma, int Mb, int kk, int nt){
	//performs Matrix Multiplication
	//c=a*b
	int u,p,l;
	
	float sum=0;
	//if kk=0 and nt=0 o1=0 and o2=Ma
	int o1=kk+(Ma*nt);
	int o2=Ma*(nt+1);
	
	//Print a
	/*
	printf("\n\n a=\n");
	for (i=0; i<m; i++){
		for (j=o1; j<o2; j++){
						printf("%03f ", a[i][j]);
		}
			printf("\n");
	}
	
	//Print b
	printf("\n\n b=\n");
	for (i=o1; i<o2; i++){
		for (j=0; j<m; j++){
						printf("%03f ", b[i-kk][j]);
		}
			printf("\n");
	}
	*/
	for (u=0; u<Na; u++)
	{
		for (p=0; p<Mb; p++)
		{
			sum=0;	
			for (l=o1; l<o2; l++)
			{

				sum+=a[u][l]*b[l-kk][p];

			}
			c[u][p]=sum;
		}
	}
}

//Quicksort functions
// A utility function to swap two elements
void swap(float* a, float* b)
{
    float t = *a;
    *a = *b;
    *b = t;
}
/* This function takes last element as pivot, places
   the pivot element at its correct position in sorted
    array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot */
int partition (float *arr, float *ind, int low, int high)
{
    int pivot = arr[high];    // pivot
    int i = (low - 1);  // Index of smaller element
 
    for (int j = low; j <= high- 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (arr[j] <= pivot)
        {
            i++;    // increment index of smaller element
            swap(&arr[i], &arr[j]);
						swap(&ind[i], &ind[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
		swap(&ind[i + 1], &ind[high]);
		
    return (i + 1);
}
 
/* The main function that implements QuickSort
 arr[] --> Array to be sorted,
  low  --> Starting index,
  high  --> Ending index */
void quickSort(float *arr, float *ind, int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
           at right place */
        int pi = partition(arr, ind, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quickSort(arr, ind, low, pi - 1);
        quickSort(arr, ind,  pi + 1, high);
    }
}








