#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

#include "matrix.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define sign(x)   ( ((x) > 0) ? 1 : (((x) < 0) ? (-1) : 0) )

///< Minimum positive floating point number x such that 1.0 + x != 1.0
///< in STM32 MCUs this value is 1.19209290e-7F
#ifndef FLT_EPSILON
#define FLT_EPSILON  1.19209290e-7F
#endif

float *mat(int m, int n)
{
    float *p;
    
    if (n <= 0 || m <= 0) return NULL;
    if (!(p = (float *) malloc(sizeof(float) * m * n))) {
        debug_printf("Matrix memory allocation error: m=%d, n=%d\n", m, n);
    }
	
    return p;
}

int matcpy(float *A, const float *B, int m, int n)
{
	if (n <= 0 || m <= 0) 
		return 1;
	else
		memcpy(A, B, sizeof(float) * m * n);  // ?? TODO: Is true?!?!?
}

float *zeros(int m, int n)
{
	float *p;
    
    if (n <= 0 || m <= 0) return NULL;
    if (!(p = (float *) calloc(sizeof(float), m * n))) {
        debug_printf("Matrix memory allocation error: m=%d, n=%d\n", m, n);
    }

    return p;	
}

float *eye(int n)
{
    float *p;
    int i;
    
    if ((p=zeros(n,n))) for (i=0; i<n; i++) p[i + i * n] = 1.0;
	
    return p;
}

void matadd(int m, int n, const float *A, const float *B, float *C)
{
    int i, j;
    
    for (i=0; i<m; i++) 
		for (j=0; j<n; j++)
			C[i+j*m] = A[i+j*m] + B[i+j*m];
}

void matmul(int m, int k, int n, const float *A, const float *B, float *C)
{
    float d;
    int i, j, x;
	
    for (i=0; i<m; i++) 
		for (j=0; j<n; j++) {
			d = 0.0;
			for (x=0; x<k; x++)
				d += A[i+x*m] * B[x+j*k];
			C[i+j*m] = d;
        }
}

void matsmul(int m, int n, const float s, const float *A, float *B)
{
    int i, j;
    
    for (i=0; i<m; i++) 
		for (j=0; j<n; j++)
			B[i+j*m] = s * A[i+j*m];
}

void transpose(int m, int n, const float *A, float *T)
{
	int i, j;
	
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			T[j+i*n] = A[i+j*m];
}

float norm2(int n, const float *a)
{
	float c = 0.0;
    
    while (--n >= 0) c += a[n] * a[n];
	
    return sqrt(c);	
}

void matCol(const float *A, int m, int j, float *v)
{
	int k;
	
	for (k = 0; k < m; k++)
		v[k] = A[k + j * m];   // A[:, j]
}

void permuteCols(float *A, int m, int j1, int j2)
{
	float *v = mat(m, 1);
	int k;
	
	for (k = 0; k < m; k++)
	{
		v[k] = A[k + j2 * m];
		A[k + j2 * m] = A[k + j1 * m];
	}

	for (k = 0; k < m; k++)
		A[k + j1 * m] = v[k];
	
	FREE(v);
}

float dotProduct(const float *u, const float *v, const int n)
{
	int k;
	float result = 0.0f;
	
	for (k = 0; k < n; k++)
		result += u[k] * v[k];
	
	return result;
}

int svn_jacobi(int m, int n, const float *A, float *U, int *r, float *S, float *V)
{
	int cnt, i, j, k;
	const int MAX_STEPS = 40;
	float TOL = 1.e-2;
	float converge;
	bool done = true;
	double aii, ajj, aij, tau;
	double Anorm;    // Frobenius norm of A matrix
	float t, c, s;
	float *ai = mat(m, 1);
	float *aj = mat(m, 1);
	float tmp;
	
	double delta = FLT_EPSILON * norm2(m * n, A);   // FLT_EPSILON = 1.19209290e-7F
	double delta2 = delta * delta;
	
	float *tmpV = eye(n);
	float *tmpA = mat(m, n);
	matcpy(tmpA, A, m, n);
	
	for (cnt = 1; cnt < MAX_STEPS; cnt++)
	{
		converge = 0.0f;
		for (j = 1; j < n; j++)
		{   // Sweep from columns (1, 2), ..., (1, n), (2, 3), ..., (n-1, n) : (i<j)
			for (i = 0; i < j; i++)
			{
				matCol(tmpA, m, i, ai);
				matCol(tmpA, m, j, aj);

				// Compute [aii aij; aij ajj] = the (i,j) submatrix of (A' * A)
				aii = dotProduct(ai, ai, m);
				ajj = dotProduct(aj, aj, m);
				aij = dotProduct(ai, aj, m);
				
				if (ajj < delta2)    // or sqrt(ajj) < delta?
					break; // ???
				else if (aii < delta2)          // means norm2(ai) < norm2(aj)
				{
					permuteCols(tmpA, m, i, j); // Permute ai, aj
					permuteCols(tmpV, n, i, j); // Permute vi, vj
				}
				else if (aij * aij / aii * ajj < TOL * TOL)
				{
					if (aii < ajj)                  // means norm2(ai) < norm2(aj)
					{
						permuteCols(tmpA, m, i, j); // Permute ai, aj
						permuteCols(tmpV, n, i, j); // Permute vi, vj
					}
				}
				else    // Solve 2x2 symmetric eigenvalue problem J'*B*J=D
				{
					converge = max(converge, fabs(aij) / sqrt(aii) / sqrt(ajj));
					/* Compute Jacobi rotation that diagonalizes following matrix:
						 [aii  aij
						  aij  ajj]
					*/
					tau = (ajj - aii) / (2.0f * aij);
					t = sign(tau) / (fabs(tau) + sqrt(1.0f + tau * tau));
					c = 1.0f  / sqrt(1 + t * t);
					s = c * t;
					// J = [c  s; -s  c]
					
					/* Update column i, j of tmpA: A=A*J 
					   [A_ki  A_kj]_new = [A_ki  A_kj]_old * [c  s; -s  c]
					*/
					for (k = 0; k < m; k++)
					{
						tmp = tmpA[k + i * m];    // Column i of tmpA matrix
						tmpA[k + i * m] = c * tmp - s * tmpA[k + j * m];
						tmpA[k + j * m] = s * tmp + c * tmpA[k + j * m];
					}
					
					/* Update the matrix V of right singular vectors: V=V*J 
					   [V_ki  V_kj]_new = [V_ki  V_kj]_old * [c  s; -s  c]
					*/
					for (k = 0; k < n; k++)
					{
						tmp = tmpV[k + i * n];   // Column i of tmpV matrix
						tmpV[k + i * n] = c * tmp - s * tmpV[k + j * n];
						tmpV[k + j * n] = s * tmp + c * tmpV[k + j * n];
					}
				}
			}
		}
	
		if (converge < TOL) break;
	} //of cnt
	
	if (cnt >= MAX_STEPS)
	{
		return 1;
	}
	
	*r = 0;
	for (j = 0; j < n; j++)
	{
		matCol(tmpA, m, j, aj);
		s = norm2(m, aj);
		
		if (s > FLT_EPSILON)
		{
			*r += 1;
			for (i = 0; i < m; i++)
			{
				U[i + m * j] = tmpA[i + m * j] / s;
				if (i == j) S[i + m * j] = s;
			}
		}
	}
	
	matcpy(V, tmpV, n, n);
	
	FREE(ai);
	FREE(aj);
	FREE(tmpA);
	FREE(tmpV);
	
	return 0;
}

int inv(const int n, const float *src, float *dst)
{
    float det;

    if (n==3)
	{
		/* Compute adjoint: */

		dst[0] = + src[4] * src[8] - src[5] * src[7];
		dst[1] = - src[1] * src[8] + src[2] * src[7];
		dst[2] = + src[1] * src[5] - src[2] * src[4];
		dst[3] = - src[3] * src[8] + src[5] * src[6];
		dst[4] = + src[0] * src[8] - src[2] * src[6];
		dst[5] = - src[0] * src[5] + src[2] * src[3];
		dst[6] = + src[3] * src[7] - src[4] * src[6];
		dst[7] = - src[0] * src[7] + src[1] * src[6];
		dst[8] = + src[0] * src[4] - src[1] * src[3];

		/* Compute determinant: */

		det = src[0] * dst[0] + src[1] * dst[3] + src[2] * dst[6];
		if (det == 0.0f)
			return 1;

		/* Multiply adjoint with reciprocal of determinant: */

		det = 1.0f / det;

		dst[0] *= det;
		dst[1] *= det;
		dst[2] *= det;
		dst[3] *= det;
		dst[4] *= det;
		dst[5] *= det;
		dst[6] *= det;
		dst[7] *= det;
		dst[8] *= det;
	}
	
    if (n==2)
	{
		/* Compute adjoint: */

		dst[0] = + src[3];
		dst[1] = - src[1];
		dst[2] = - src[2];
		dst[3] = + src[0];

		/* Compute determinant: */

		det = src[0] * dst[0] + src[1] * dst[2];
		if (det == 0.0f)
			return 1;

		/* Multiply adjoint with reciprocal of determinant: */

		det = 1.0f / det;

		dst[0] *= det;
		dst[1] *= det;
		dst[2] *= det;
		dst[3] *= det;		
	}
	
	return 0;
}
