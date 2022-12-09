/*!
 * @file   matrix.h
 * @author Davood Rezaei
 * @date Created on Sept. 24, 2022, 11:17 AM
 */
 
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdint.h>
#include <stdlib.h>

#define FREE(X)    free_and_null((void *) (X))

inline void free_and_null(void *var)
{
	if (var != NULL)
	{
		free(var);
		var = NULL;
	}
}

/*! 
* @brief New matrix. Allocate memory to matrix.
*
* A is a mxn matrix: A_ij = p[i+j*m]: Column-ward Counting {i=0..(m-1) , j=0..(n-1)}
*
*	  A11(p[0+0*m])        A12(p[0+1*m])       ...     A1n(p[0+(n-1)*m])
*	  A21(p[1+0*m])        A22(p[1+1*m])       ...     A2n(p[1+(n-1)*m])
*	  ...                  ...                 ...     ...
*	  Am1(p[m-1+0*m])      Am2(p[m-1+1*m])     ...     Amn(p[m-1+(n-1)*m])
*
* @param m  The number of rows in matrix
* @param n  The number of columns in matrix
* @return  Matrix pointer (if n<=0 or m<=0, return NULL)
*/
float *mat(int m, int n);

/*!
* @brief  Copy matrix 
* @param  A The address of destination matrix A (m x n)
* @param  B The address of source matrix B (m x n)
* @param  m The number of rows in matrix
* @param  n The number of columns in matrix
* @return success/fail flag
*/
int matcpy(float *A, const float *B, int m, int n);

/*! 
* @brief Zero matrix
* @param m The number of rows in matrix
* @param  n The number of columns in matrix
* @return matrix pointer (if m<=0 or n<=0, return NULL)
*/
float *zeros(int m, int n);

/*!
* @brief  Identity matrix
* @param  n The number of rows and columns of matrix
* @return matrix pointer (if n<=0, return NULL)
*/
float *eye(int n);

/*!
* @brief Matrix addition
* @param m The number of rows in matrixs
* @param n The number of columns in matrixs
* @param A The address of first matrix
* @param B The address of second matrix
* @param C The address of result matrix
* @return void
*/
void matadd(int m, int n, const float *A, const float *B, float *C);

/*!
* @brief Matrix scalar multiplication
* @param m The number of rows in matrixs
* @param n The number of columns in matrixs
* @param s multiplier
* @param A The address of input matrix
* @param B The address of result matrix
* @return void
*/
void matsmul(int m, int n, const float s, const float *A, float *B);

/*!
* @brief  Matrix multiplication
*
*   C_mn = A_mk x B_kn
*   C_ij = SUM(A_ix * B_xj) ---> C[i+j*m] = SUM_x(A[i+x*m] * B[x+j*k]) , x=(0..k-1)
*  
* @param m The number of rows in A matrix
* @param k The number rows in A matrix and columns in B matrix
* @param n The number of columns in B matrix
* @param A The address of A matrix
* @param B The address of B matrix
* @param C The address of C(result) matrix
* @return void
*/
void matmul(int m, int k, int n, const float *A, const float *B, float *C);

/*!
* @brief Compute SVN of a matrix using one-sided jacobi method
*
*  A = U * S * V'
*
* @param m, row of matrix
* @param n, column of matrix
* @param A, mxn input matrix
* @param U, mxm matrix
* @param S, mxn matrix
* @param V, nxn matrix
* @return success/fail flag
*/
int svn_jacobi(int m, int n, const float *A, float *U, int *r, float *S, float *V);

/*!
* @brief  Matrix inversion for 2x2 and 3x3 matrixs
* @param n The dimension of input matrix
* @param A The address of A matrix
* @param B The address of B(inverted) matrix
* @return success/fail flag
*/
int inv(const int n, const float *A, float *B);

/*! 
* @brief  Matrix transpose
*
*   A: mxn
*   T: nxm
*   T_ji = A_ij  ==> t[j+i*n] = a[i+j*m]
*
* @param m The number of rows in input matrix
* @param n The number of columns in input matrix
* @param A The address of input matrix
* @param T The address of T(result) matrix
* @return void
*/
void transpose(int m, int n, const float *A, float *T);

/*!
* @brief  Euclidian norm of a vector
* @param n The dimension of input vector
* @param a The address of input vector
* @return vector norm
*/
float norm2(int n, const float *a);

/*!
* @brief  Fill in the rotation matrix
* @param axis The rotation axis
* @param ang The angle of rotation
* @param R The address of R(result) matrix
* @return void
*/

/*!
* @brief  Return a specific column of a matrix
* @param A The input matrix
* @param m The number of row in matrix
* @param j The column number that will be returned
* @param j The returned column
* @return void
*/
void matCol(const float *A, int m, int j, float *v);

/*!
* @brief  Return dot product of two vectors
* @param u, v Two input vectors
* @param n The length of row vectors
* @return dot product of vectors
*/
float dotProduct(const float *u, const float *v, const int n);

#endif    // _MATRIX_H_