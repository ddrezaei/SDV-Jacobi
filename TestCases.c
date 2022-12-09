#include <stdio.h>
#include "TestCases.h"
#include "matrix.h"

void print_matrix(const char *desc, int m, int n, const float *A)
{
	int i, j;

	printf(" %s = [  \r\n", desc);
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
			printf("    %f", A[i + m * j]); 
		
		if (i < m - 1)
			printf("\r\n"); 
		else
			printf("]\r\n"); 
	}
}

void Test_matcpy(void)
{
	float *A = zeros(3, 3);
	float *B = mat(3, 3);
	
	A[0] = 1.0f; A[1] = -3.2f; A[2] = 0.34f;
	A[4] = 5.7f; A[5] = 0.36f; A[6] = 4.06f;
	A[8] = -2.5f;
	
	matcpy(B, A, 3, 3);
	print_matrix("A Matrix", 3, 3, A);
	print_matrix("B Matrix", 3, 3, B);
	
	FREE(A);
	FREE(B);
}

void Test_SVD(int m, int n, const float *A)
{
	int rank;
	float *U = mat(m, m);
	float *V = mat(n, n);
	float *S = zeros(m, n);
	
	print_matrix("A Matrix", m, n, A);
	
	if (svn_jacobi(m, n, A, U, &rank, S, V) == 0)
	{
		printf("  Rank = %d \r\n", rank);
		print_matrix("U ", m, m, U);
		print_matrix("V ", n, n, V);
		print_matrix("S ", m, n, S);
		
		///< -------------------------- Test U'xU = Imxm
		float *Ut = mat(m, m);
		float *tmp1 = mat(m, m);
		
		transpose(m, m, U, Ut);
		matmul(m, m, m, Ut, U, tmp1);
		print_matrix("U'xU ", m, m, tmp1);
		
		FREE(Ut);
		FREE(tmp1);
		
		///< -------------------------- Test V'xV = Inxn
		float *Vt = mat(n, n);
		float *tmp2 = mat(n, n);
		
		transpose(n, n, V, Vt);
		matmul(n, n, n, Vt, V, tmp2);
		print_matrix("V'xV ", n, n, tmp2);
		
		FREE(Vt);
		FREE(tmp2);
		
		///< --------------------------  Test UxSxV' = A
		float *Vt2 = mat(n, n);
		float *tmp3 = mat(m, n);
		float *tmp4 = mat(m, n);
		
		transpose(n, n, V, Vt2);
		matmul(m, m, n, U, S, tmp3);
		matmul(m, n, n, tmp3, Vt2, tmp4);
		print_matrix("A-SVD ", m, n, tmp4);
		
		FREE(Vt2);
		FREE(tmp3);
		FREE(tmp4);
	}
	else
		printf("Error: Jacobi SVD failed to converge!' \r\n");
	
	FREE(U);
	FREE(S);
	FREE(V);
}