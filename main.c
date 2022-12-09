#include <stdio.h>
#include "matrix.h"
#include "TestCases.h"

int main(void)
{
	float *A3x3 = zeros(3, 3);
	float *A8x3 = zeros(8, 3);
	float *A8x8 = zeros(6, 4);
	
	A3x3[0] = 100.0f; A3x3[1] = -31.2f; A3x3[2] = 0.34f;
	A3x3[4] = 5.7f; A3x3[5] = 0.36f; A3x3[6] = 40.6f;
	A3x3[8] = -2.5f;
	
	
	A8x3[0] = 3.67; A8x3[4] = -50.5; A8x3[6] = 10.456;
	A8x3[9] = 0.45; A8x3[10] = 112; A8x3[15] = 1.2;
	A8x3[16] = -171.7; A8x3[22] = 1000.9;

	A8x8[0] = 3.67; A8x8[4] = -50.5; A8x8[6] = 10.456;
	A8x8[9] = 0.45; A8x8[10] = 112; A8x8[15] = 1.2;
	A8x8[16] = -171.7; A8x8[22] = 1000.9;
	
	Test_SVD(8, 8, A8x8);
	
	FREE(A3x3);
	FREE(A8x3);
	FREE(A8x8);
}
