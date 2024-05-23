#ifndef POLY_MUL_H
#define POLY_MUL_H

#include "SABER_params.h"
#include <stdint.h>
#include <string.h>
#include <immintrin.h>

#define AVX_N (SABER_N >> 4)
#define small_len_avx (AVX_N >> 2)

#define SCHB_N 16

#define N_SB (SABER_N >> 2)
#define N_SB_RES (2*N_SB-1)

#define N_SB_16 (N_SB >> 2)
#define N_SB_16_RES (2*N_SB_16-1)

#define AVX_N1 16 /*N/16*/ 

#define SCM_SIZE 16

// The dimension of a vector. i.e vector has NUM_POLY elements and Matrix has NUM_POLY X NUM_POLY elements
#define NUM_POLY SABER_K
//int NUM_POLY=2; 


void poly_mul_acc(const uint16_t a[SABER_N], const uint16_t b[SABER_N], uint16_t res[SABER_N]);

/*
// Helper function to print avx vector
void print_avx2(__m256i a)
{
	uint16_t a_array[16];
	_mm256_maskstore_epi32(&a_array, mask, a);
	int16_t i;
	for(i=0; i<16; i++)
	printf("%u, ", a_array[i]);
	printf("\n");	
}
*/

void transpose_n1(__m256i *M)
{
	//int i;
	register __m256i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
	register __m256i temp, temp0, temp1, temp2;

	//for(i=0; i<8; i=i+1)
	//{
		r0 = _mm256_unpacklo_epi16(M[0], M[1]); 
		r1 = _mm256_unpacklo_epi16(M[2], M[3]); 
		r2 = _mm256_unpacklo_epi16(M[4], M[5]); 
		r3 = _mm256_unpacklo_epi16(M[6], M[7]);
		r4 = _mm256_unpacklo_epi16(M[8], M[9]); 
		r5 = _mm256_unpacklo_epi16(M[10], M[11]);
		r6 = _mm256_unpacklo_epi16(M[12], M[13]); 
		r7 = _mm256_unpacklo_epi16(M[14], M[15]); 


		temp = _mm256_unpacklo_epi32(r0, r1); 
		temp0 = _mm256_unpacklo_epi32(r2, r3); 
		temp1 = _mm256_unpacklo_epi32(r4, r5); 
		temp2 = _mm256_unpacklo_epi32(r6, r7); 

		r8 = _mm256_unpackhi_epi32(r0, r1); 
		r9 = _mm256_unpackhi_epi32(r2, r3); 
		r10 = _mm256_unpackhi_epi32(r4, r5); 
		r11 = _mm256_unpackhi_epi32(r6, r7);

		r0 = _mm256_unpacklo_epi64(temp, temp0); 
		r2 = _mm256_unpackhi_epi64(temp, temp0); 

		r1 = _mm256_unpacklo_epi64(temp1, temp2); 
		r3 = _mm256_unpackhi_epi64(temp1, temp2);

		temp = _mm256_unpackhi_epi16(M[0], M[1]); 
		temp0 = _mm256_unpackhi_epi16(M[2], M[3]); 
		temp1 = _mm256_unpackhi_epi16(M[4], M[5]); 
		temp2 = _mm256_unpackhi_epi16(M[6], M[7]); 
		r4 = _mm256_unpackhi_epi16(M[8], M[9]); 

		M[0] = _mm256_permute2f128_si256(r0, r1, 0x20);
		M[8] = _mm256_permute2f128_si256(r0, r1, 0x31);
		M[1] = _mm256_permute2f128_si256(r2, r3, 0x20);
		M[9] = _mm256_permute2f128_si256(r2, r3, 0x31);


		r5 = _mm256_unpackhi_epi16(M[10], M[11]); 
		r6 = _mm256_unpackhi_epi16(M[12], M[13]); 
		r7 = _mm256_unpackhi_epi16(M[14], M[15]); 



		r0 = _mm256_unpacklo_epi64(r8, r9); 
		r1 = _mm256_unpacklo_epi64(r10, r11); 

		r2 = _mm256_unpackhi_epi64(r8, r9); 
		r3 = _mm256_unpackhi_epi64(r10, r11); 



		M[3] = _mm256_permute2f128_si256(r2, r3, 0x20);
		M[11] = _mm256_permute2f128_si256(r2, r3, 0x31);
		M[2] = _mm256_permute2f128_si256(r0, r1, 0x20);
		M[10] = _mm256_permute2f128_si256(r0, r1, 0x31);


	//for(i=0; i<4; i=i+1)
	//{
		r0 = _mm256_unpacklo_epi32(temp, temp0); 
		r1 = _mm256_unpacklo_epi32(temp1, temp2);
		r2 = _mm256_unpacklo_epi32(r4, r5); 
		r3 = _mm256_unpacklo_epi32(r6, r7); 

	//}


	//for(i=0; i<2; i=i+1)
	//{
		r8 = _mm256_unpacklo_epi64(r0, r1); 
		r10 = _mm256_unpackhi_epi64(r0, r1); 

		r9 = _mm256_unpacklo_epi64(r2, r3); 
		r11 = _mm256_unpackhi_epi64(r2, r3); 

		M[4] = _mm256_permute2f128_si256(r8, r9, 0x20);
		M[12] = _mm256_permute2f128_si256(r8, r9, 0x31);
		M[5] = _mm256_permute2f128_si256(r10, r11, 0x20);
		M[13] = _mm256_permute2f128_si256(r10, r11, 0x31);

		r0 = _mm256_unpackhi_epi32(temp, temp0); 
		r1 = _mm256_unpackhi_epi32(temp1, temp2); 
		r2 = _mm256_unpackhi_epi32(r4, r5); 
		r3 = _mm256_unpackhi_epi32(r6, r7); 

	//}
//	for(i=0; i<2; i=i+1)
//	{
		r4 = _mm256_unpacklo_epi64(r0, r1); 
		r6 = _mm256_unpackhi_epi64(r0, r1); 

		r5 = _mm256_unpacklo_epi64(r2, r3); 
		r7 = _mm256_unpackhi_epi64(r2, r3); 

//	}

	//-------------------------------------------------------

	M[6] = _mm256_permute2f128_si256(r4, r5, 0x20);
	M[14] = _mm256_permute2f128_si256(r4, r5, 0x31);
	M[7] = _mm256_permute2f128_si256(r6, r7, 0x20);
	M[15] = _mm256_permute2f128_si256(r6, r7, 0x31);
}

/*
void transpose_unrolled(__m256i *M)
{
	int i;
	__m256i tL[8], tH[8];
	__m256i bL[4], bH[4], cL[4], cH[4];
	__m256i dL[2], dH[2], eL[2], eH[2], fL[2], fH[2], gL[2], gH[2];

	__m256i r0, r1, r2, r3, r4, r5, r6, r7;

	//for(i=0; i<8; i=i+1)
	//{
		tL[0] = _mm256_unpacklo_epi16(M[0], M[1]); 
		tH[0] = _mm256_unpackhi_epi16(M[0], M[1]); 

		tL[1] = _mm256_unpacklo_epi16(M[2], M[3]); 
		tH[1] = _mm256_unpackhi_epi16(M[2], M[3]); 

		tL[2] = _mm256_unpacklo_epi16(M[4], M[5]); 
		tH[2] = _mm256_unpackhi_epi16(M[4], M[5]); 

		tL[3] = _mm256_unpacklo_epi16(M[6], M[7]); 
		tH[3] = _mm256_unpackhi_epi16(M[6], M[7]); 

		tL[4] = _mm256_unpacklo_epi16(M[8], M[9]); 
		tH[4] = _mm256_unpackhi_epi16(M[8], M[9]); 

		tL[5] = _mm256_unpacklo_epi16(M[10], M[11]); 
		tH[5] = _mm256_unpackhi_epi16(M[10], M[11]); 

		tL[6] = _mm256_unpacklo_epi16(M[12], M[13]); 
		tH[6] = _mm256_unpackhi_epi16(M[12], M[13]); 

		tL[7] = _mm256_unpacklo_epi16(M[14], M[15]); 
		tH[7] = _mm256_unpackhi_epi16(M[14], M[15]); 

	//}

	//-------------------------------------------------------
	//for(i=0; i<4; i=i+1)
	//{
		bL[0] = _mm256_unpacklo_epi32(tL[0], tL[1]); 
		bH[0] = _mm256_unpackhi_epi32(tL[0], tL[1]); 

		bL[1] = _mm256_unpacklo_epi32(tL[2], tL[3]); 
		bH[1] = _mm256_unpackhi_epi32(tL[2], tL[3]); 

		bL[2] = _mm256_unpacklo_epi32(tL[4], tL[5]); 
		bH[2] = _mm256_unpackhi_epi32(tL[4], tL[5]); 

		bL[3] = _mm256_unpacklo_epi32(tL[6], tL[7]); 
		bH[3] = _mm256_unpackhi_epi32(tL[6], tL[7]); 

	//}

	//for(i=0; i<2; i=i+1)
	//{
		dL[0] = _mm256_unpacklo_epi64(bL[0], bL[1]); 
		dH[0] = _mm256_unpackhi_epi64(bL[0], bL[1]); 

		dL[1] = _mm256_unpacklo_epi64(bL[2], bL[3]); 
		dH[1] = _mm256_unpackhi_epi64(bL[2], bL[3]);

		M[0] = _mm256_permute2f128_si256(dL[0], dL[1], 0x20);
		M[8] = _mm256_permute2f128_si256(dL[0], dL[1], 0x31);
		M[1] = _mm256_permute2f128_si256(dH[0], dH[1], 0x20);
		M[9] = _mm256_permute2f128_si256(dH[0], dH[1], 0x31);

	//}
	//for(i=0; i<2; i=i+1)
	//{
		eL[0] = _mm256_unpacklo_epi64(bH[0], bH[1]); 
		eH[0] = _mm256_unpackhi_epi64(bH[0], bH[1]); 

		eL[1] = _mm256_unpacklo_epi64(bH[2], bH[3]); 
		eH[1] = _mm256_unpackhi_epi64(bH[2], bH[3]); 

	//}

	//-------------------------------------------------------

	//-------------------------------------------------------
	for(i=0; i<4; i=i+1)
	{
		cL[i] = _mm256_unpacklo_epi32(tH[2*i], tH[2*i+1]); 
		cH[i] = _mm256_unpackhi_epi32(tH[2*i], tH[2*i+1]); 
	}


	for(i=0; i<2; i=i+1)
	{
		fL[i] = _mm256_unpacklo_epi64(cL[2*i], cL[2*i+1]); 
		fH[i] = _mm256_unpackhi_epi64(cL[2*i], cL[2*i+1]); 
	}
	for(i=0; i<2; i=i+1)
	{
		gL[i] = _mm256_unpacklo_epi64(cH[2*i], cH[2*i+1]); 
		gH[i] = _mm256_unpackhi_epi64(cH[2*i], cH[2*i+1]); 
	}

	//-------------------------------------------------------



	M[2] = _mm256_permute2f128_si256(eL[0], eL[1], 0x20);
	M[10] = _mm256_permute2f128_si256(eL[0], eL[1], 0x31);
	M[3] = _mm256_permute2f128_si256(eH[0], eH[1], 0x20);
	M[11] = _mm256_permute2f128_si256(eH[0], eH[1], 0x31);

	M[4] = _mm256_permute2f128_si256(fL[0], fL[1], 0x20);
	M[12] = _mm256_permute2f128_si256(fL[0], fL[1], 0x31);
	M[5] = _mm256_permute2f128_si256(fH[0], fH[1], 0x20);
	M[13] = _mm256_permute2f128_si256(fH[0], fH[1], 0x31);

	M[6] = _mm256_permute2f128_si256(gL[0], gL[1], 0x20);
	M[14] = _mm256_permute2f128_si256(gL[0], gL[1], 0x31);
	M[7] = _mm256_permute2f128_si256(gH[0], gH[1], 0x20);
	M[15] = _mm256_permute2f128_si256(gH[0], gH[1], 0x31);
}


void transpose1(__m256i *M)
{
	int i;
	__m256i tL[8], tH[8];
	__m256i bL[4], bH[4], cL[4], cH[4];
	__m256i dL[2], dH[2], eL[2], eH[2], fL[2], fH[2], gL[2], gH[2];

	for(i=0; i<8; i=i+1)
	{
		tL[i] = _mm256_unpacklo_epi16(M[2*i], M[2*i+1]); 
		tH[i] = _mm256_unpackhi_epi16(M[2*i], M[2*i+1]); 
	}

	for(i=0; i<4; i=i+1)
	{
		bL[i] = _mm256_unpacklo_epi32(tL[2*i], tL[2*i+1]); 
		bH[i] = _mm256_unpackhi_epi32(tL[2*i], tL[2*i+1]); 
	}
	for(i=0; i<4; i=i+1)
	{
		cL[i] = _mm256_unpacklo_epi32(tH[2*i], tH[2*i+1]); 
		cH[i] = _mm256_unpackhi_epi32(tH[2*i], tH[2*i+1]); 
	}

	for(i=0; i<2; i=i+1)
	{
		dL[i] = _mm256_unpacklo_epi64(bL[2*i], bL[2*i+1]); 
		dH[i] = _mm256_unpackhi_epi64(bL[2*i], bL[2*i+1]); 
	}
	for(i=0; i<2; i=i+1)
	{
		eL[i] = _mm256_unpacklo_epi64(bH[2*i], bH[2*i+1]); 
		eH[i] = _mm256_unpackhi_epi64(bH[2*i], bH[2*i+1]); 
	}

	for(i=0; i<2; i=i+1)
	{
		fL[i] = _mm256_unpacklo_epi64(cL[2*i], cL[2*i+1]); 
		fH[i] = _mm256_unpackhi_epi64(cL[2*i], cL[2*i+1]); 
	}
	for(i=0; i<2; i=i+1)
	{
		gL[i] = _mm256_unpacklo_epi64(cH[2*i], cH[2*i+1]); 
		gH[i] = _mm256_unpackhi_epi64(cH[2*i], cH[2*i+1]); 
	}

	M[0] = _mm256_permute2f128_si256(dL[0], dL[1], 0x20);
	M[8] = _mm256_permute2f128_si256(dL[0], dL[1], 0x31);
	M[1] = _mm256_permute2f128_si256(dH[0], dH[1], 0x20);
	M[9] = _mm256_permute2f128_si256(dH[0], dH[1], 0x31);

	M[2] = _mm256_permute2f128_si256(eL[0], eL[1], 0x20);
	M[10] = _mm256_permute2f128_si256(eL[0], eL[1], 0x31);
	M[3] = _mm256_permute2f128_si256(eH[0], eH[1], 0x20);
	M[11] = _mm256_permute2f128_si256(eH[0], eH[1], 0x31);

	M[4] = _mm256_permute2f128_si256(fL[0], fL[1], 0x20);
	M[12] = _mm256_permute2f128_si256(fL[0], fL[1], 0x31);
	M[5] = _mm256_permute2f128_si256(fH[0], fH[1], 0x20);
	M[13] = _mm256_permute2f128_si256(fH[0], fH[1], 0x31);

	M[6] = _mm256_permute2f128_si256(gL[0], gL[1], 0x20);
	M[14] = _mm256_permute2f128_si256(gL[0], gL[1], 0x31);
	M[7] = _mm256_permute2f128_si256(gH[0], gH[1], 0x20);
	M[15] = _mm256_permute2f128_si256(gH[0], gH[1], 0x31);
}
*/

//#define SCM_SIZE 16

//#pragma STDC FP_CONTRACT ON


inline __m256i mul_add(__m256i a, __m256i b, __m256i c) { 
    return _mm256_add_epi16(_mm256_mullo_epi16(a, b), c);
}


void schoolbook_avx_new3_acc(__m256i* a, __m256i* b, __m256i* c_avx) ////8 coefficients of a and b has been prefetched
									      //the c_avx are added cummulatively
{

	register __m256i a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7;
	register __m256i temp;


	a0=a[0];
	a1=a[1];
	a2=a[2];
	a3=a[3];
	a4=a[4];
	a5=a[5];
	a6=a[6];
	a7=a[7];

	b0=b[0];
	b1=b[1];
	b2=b[2];
	b3=b[3];
	b4=b[4];
	b5=b[5];
	b6=b[6];
	b7=b[7];

	// New Unrolled first triangle

	//otherwise accumulate
	c_avx[0] = mul_add(a0, b0, c_avx[0]);
	

	temp = _mm256_mullo_epi16 (a0, b1);
	temp=mul_add(a1, b0, temp);
	c_avx[1] = _mm256_add_epi16(temp, c_avx[1]);


	temp = _mm256_mullo_epi16 (a0, b2);
	temp = mul_add(a1, b1, temp);
	temp=mul_add(a2, b0, temp);
	c_avx[2] = _mm256_add_epi16(temp, c_avx[2]);
	

	temp = _mm256_mullo_epi16 (a0, b3);
	temp = mul_add(a1, b2, temp);
	temp = mul_add(a2, b1, temp);
	temp=mul_add(a3, b0, temp);
	c_avx[3] = _mm256_add_epi16(temp, c_avx[3]);

	temp = _mm256_mullo_epi16 (a0, b4);
	temp = mul_add(a1, b3, temp);
	temp = mul_add(a3, b1, temp);
	temp = mul_add(a4, b0, temp);
	temp=mul_add(a2, b2, temp);
	c_avx[4] = _mm256_add_epi16(temp, c_avx[4]);


	temp = _mm256_mullo_epi16 (a0, b5);
	temp = mul_add(a1, b4 , temp);
	temp = mul_add(a2, b3, temp);
	temp = mul_add(a3, b2, temp);
	temp = mul_add( a4, b1, temp);
	temp=mul_add(a5, b0, temp);
	c_avx[5] = _mm256_add_epi16(temp, c_avx[5]);
	
	temp = _mm256_mullo_epi16 (a0, b6);
	temp = mul_add(a1, b5, temp);
	temp = mul_add(a5, b1, temp);
	temp = mul_add(a6, b0, temp);
	temp = mul_add(a2, b4, temp);
	temp = mul_add(a3, b3, temp);
	temp=mul_add(a4, b2, temp);
	c_avx[6] = _mm256_add_epi16(temp, c_avx[6]);


	temp = _mm256_mullo_epi16 (a0, b7);
	temp = mul_add(a1, b6, temp);
	temp = mul_add (a6, b1, temp);
	temp = mul_add (a7, b0, temp);
	temp = mul_add(a2, b5, temp);
	temp = mul_add (a3, b4, temp);
	temp = mul_add (a4, b3, temp);
	temp=mul_add(a5, b2, temp);
	c_avx[7] = _mm256_add_epi16(temp, c_avx[7]);

	temp = _mm256_mullo_epi16 (a0, b[8]);
	temp = mul_add (a1, b7, temp);
	temp = mul_add (a7, b1, temp);
	temp = mul_add (a[8], b0, temp);
	temp = mul_add (a2, b6,temp);
	temp = mul_add(a3, b5, temp);
	temp = mul_add (a4, b4,temp);
	temp = mul_add (a5, b3, temp);
	
		temp=mul_add(a6, b2, temp);
		c_avx[8] = _mm256_add_epi16(temp, c_avx[8]);


	temp = _mm256_mullo_epi16 (a0, b[9]);
	temp = mul_add (a1, b[8], temp);
	temp = mul_add (a[8], b1, temp);
	temp = mul_add (a[9], b0, temp);
	temp = mul_add (a2, b7, temp);
	temp = mul_add (a3, b6, temp);
	temp = mul_add (a4, b5, temp);
	temp = mul_add (a5, b4, temp);
	temp = mul_add (a6, b3, temp);
		temp=mul_add(a7, b2, temp);
		c_avx[9] = _mm256_add_epi16(temp, c_avx[9]);


	temp= _mm256_mullo_epi16 (a0, b[10]);
	temp = mul_add (a1, b[9], temp);
	temp = mul_add (a[9], b1, temp);
	temp = mul_add (a[10], b0, temp);
	temp = mul_add (a2, b[8], temp);
	temp = mul_add (a3, b7, temp);
	temp = mul_add (a4, b6, temp);
	temp = mul_add (a5, b5, temp);
	temp = mul_add (a6, b4, temp);
	temp = mul_add (a7, b3, temp);
		temp=mul_add(a[8], b2, temp);
		c_avx[10] = _mm256_add_epi16(temp, c_avx[10]);


	temp = _mm256_mullo_epi16 (a0, b[11]);
	temp = mul_add (a1, b[10], temp );
	temp = mul_add (a[10], b1, temp );
	temp = mul_add (a[11], b0, temp );
	temp = mul_add (a2, b[9], temp );
	temp = mul_add (a3, b[8], temp );
	temp = mul_add (a4, b7, temp );
	temp = mul_add (a5, b6, temp );
	temp = mul_add (a6, b5, temp );
	temp = mul_add (a7, b4, temp );
	temp = mul_add (a[8], b3, temp );
		temp=mul_add(a[9], b2, temp);
		c_avx[11] = _mm256_add_epi16(temp, c_avx[11]);


	temp = _mm256_mullo_epi16 (a0, b[12]);
	temp = mul_add (a1, b[11], temp);
	temp = mul_add (a[11], b1, temp);
	temp = mul_add (a[12], b0, temp);
	temp = mul_add (a2, b[10], temp);
	temp = mul_add (a3, b[9], temp);
	temp = mul_add (a4, b[8], temp);
	temp = mul_add (a5, b7, temp);
	temp = mul_add (a6, b6, temp);
	temp = mul_add (a7, b5, temp);
	temp = mul_add (a[8], b4, temp);
	temp = mul_add (a[9], b3, temp);
		temp=mul_add(a[10], b2, temp);
		c_avx[12] = _mm256_add_epi16(temp, c_avx[12]);


	temp = _mm256_mullo_epi16 (a0, b[13]);
	temp = mul_add (a1, b[12], temp );
	temp = mul_add (a[12], b1, temp );
	temp = mul_add (a[13], b0, temp );
	temp = mul_add (a2, b[11], temp );
	temp = mul_add (a3, b[10], temp );
	temp = mul_add (a4, b[9], temp );
	temp = mul_add (a5, b[8], temp );
	temp = mul_add (a6, b7, temp );
	temp = mul_add (a7, b6, temp );
	temp = mul_add (a[8], b5, temp );
	temp = mul_add (a[9], b4, temp );
	temp = mul_add (a[10], b3, temp );
		temp=mul_add(a[11], b2, temp);
		c_avx[13] = _mm256_add_epi16(temp, c_avx[13]);



	temp = _mm256_mullo_epi16 (a0, b[14]);
	temp = mul_add (a1, b[13], temp );
	temp = mul_add (a[13], b1, temp );
	temp = mul_add (a[14], b0, temp );
	temp = mul_add (a2, b[12], temp );
	temp = mul_add (a3, b[11], temp );
	temp = mul_add (a4, b[10], temp );
	temp = mul_add (a5, b[9], temp );
	temp = mul_add (a6, b[8], temp );
	temp = mul_add (a7, b7, temp );
	temp = mul_add (a[8], b6, temp );
	temp = mul_add (a[9], b5, temp );
	temp = mul_add (a[10], b4, temp );
	temp = mul_add (a[11], b3, temp );
		temp=mul_add(a[12], b2, temp);
		c_avx[14] = _mm256_add_epi16(temp, c_avx[14]);


	temp = _mm256_mullo_epi16 (a0, b[15]);
	temp = mul_add (a1, b[14], temp );
	temp = mul_add (a[14], b1, temp );
	temp = mul_add (a[15], b0, temp );
	temp = mul_add (a2, b[13], temp );
	temp = mul_add (a3, b[12], temp );
	temp = mul_add (a4, b[11], temp );
	temp = mul_add (a5, b[10], temp );
	temp = mul_add (a6, b[9], temp );
	temp = mul_add (a7, b[8], temp );
	temp = mul_add (a[8], b7, temp );
	temp = mul_add (a[9], b6, temp );
	temp = mul_add (a[10], b5, temp );
	temp = mul_add (a[11], b4, temp );
	temp = mul_add (a[12], b3, temp );
		temp=mul_add(a[13], b2, temp);
		c_avx[15] = _mm256_add_epi16(temp, c_avx[15]);


	// unrolled second triangle
	a0=a[14];
	a1=a[15];
	a2=a[13];
	a3=a[12];
	a4=a[11];
	a5=a[10];
	a6=a[9];
	a7=a[8];

	b0=b[14];
	b1=b[15];
	b2=b[13];
	b3=b[12];
	b4=b[11];
	b5=b[10];
	b6=b[9];
	b7=b[8];

	temp = _mm256_mullo_epi16 (a[1], b1);
	temp = mul_add (a[2], b0, temp );
	temp = mul_add (a[3], b2, temp );
	temp = mul_add (a[4], b3, temp );
	temp = mul_add (a[5], b4, temp );
	temp = mul_add (a[6], b5, temp );
	temp = mul_add (a[7], b6, temp );
	temp = mul_add (a7, b7, temp );
	temp = mul_add (a6, b[7], temp );
	temp = mul_add (a5, b[6], temp );
	temp = mul_add (a4, b[5], temp );
	temp = mul_add (a3, b[4], temp );
	temp = mul_add (a2, b[3], temp );
	temp = mul_add (a0, b[2], temp );
		temp=mul_add(a1, b[1], temp);
		c_avx[16] = _mm256_add_epi16(temp, c_avx[16]);


	temp = _mm256_mullo_epi16 (a[2], b1);
	temp = mul_add (a[3], b0, temp );
	temp = mul_add (a[4], b2, temp );
	temp = mul_add (a[5], b3, temp );
	temp = mul_add (a[6], b4, temp );
	temp = mul_add (a[7], b5, temp );
	temp = mul_add (a7, b6, temp );
	temp = mul_add (a6, b7, temp );
	temp = mul_add (a5, b[7], temp );
	temp = mul_add (a4, b[6], temp );
	temp = mul_add (a3, b[5], temp );
	temp = mul_add (a2, b[4], temp );
	temp = mul_add (a0, b[3], temp );
		temp=mul_add(a1, b[2], temp);
		c_avx[17] = _mm256_add_epi16(temp, c_avx[17]);


	temp = _mm256_mullo_epi16 (a[3], b1);
	temp = mul_add (a[4], b0, temp );
	temp = mul_add (a[5], b2, temp );
	temp = mul_add (a[6], b3, temp );
	temp = mul_add (a[7], b4, temp );
	temp = mul_add (a7, b5, temp );
	temp = mul_add (a6, b6, temp );
	temp = mul_add (a5, b7, temp );
	temp = mul_add (a4, b[7], temp );
	temp = mul_add (a3, b[6], temp );
	temp = mul_add (a2, b[5], temp );
	temp = mul_add (a0, b[4], temp );
		temp=mul_add(a1, b[3], temp);
		c_avx[18] = _mm256_add_epi16(temp, c_avx[18]);


	temp = _mm256_mullo_epi16 (a[4], b1);
	temp = mul_add (a[5], b0, temp );
	temp = mul_add (a[6], b2, temp );
	temp = mul_add (a[7], b3, temp );
	temp = mul_add (a7, b4, temp );
	temp = mul_add (a6, b5, temp );
	temp = mul_add (a5, b6, temp );
	temp = mul_add (a4, b7, temp );
	temp = mul_add (a3, b[7], temp );
	temp = mul_add (a2, b[6], temp );
	temp = mul_add (a0, b[5], temp );
		temp=mul_add(a1, b[4], temp);
		c_avx[19] = _mm256_add_epi16(temp, c_avx[19]);


	temp = _mm256_mullo_epi16 (a[5], b1);
	temp = mul_add (a[6], b0, temp );
	temp = mul_add (a[7], b2, temp );
	temp = mul_add (a7, b3, temp );
	temp = mul_add (a6, b4, temp );
	temp = mul_add (a5, b5, temp );
	temp = mul_add (a4, b6, temp );
	temp = mul_add (a3, b7, temp );
	temp = mul_add (a2, b[7], temp );
	temp = mul_add (a0, b[6], temp );
		temp=mul_add(a1, b[5], temp);
		c_avx[20] = _mm256_add_epi16(temp, c_avx[20]);


	temp = _mm256_mullo_epi16 (a[6], b1);
	temp = mul_add (a[7], b0, temp );
	temp = mul_add (a7, b2, temp );
	temp = mul_add (a6, b3, temp );
	temp = mul_add (a5, b4, temp );
	temp = mul_add (a4, b5, temp );
	temp = mul_add (a3, b6, temp );
	temp = mul_add (a2, b7, temp );
	temp = mul_add (a0, b[7], temp );
		temp=mul_add(a1, b[6], temp);
		c_avx[21] = _mm256_add_epi16(temp, c_avx[21]);


	temp = _mm256_mullo_epi16 (a[7], b1);
	temp = mul_add (a7, b0, temp );
	temp = mul_add (a6, b2, temp );
	temp = mul_add (a5, b3, temp );
	temp = mul_add (a4, b4, temp );
	temp = mul_add (a3, b5, temp );
	temp = mul_add (a2, b6, temp );
	temp = mul_add (a0, b7, temp );
		temp=mul_add(a1, b[7], temp);
		c_avx[22] = _mm256_add_epi16(temp, c_avx[22]);


	temp = _mm256_mullo_epi16 (a7, b1);
	temp = mul_add (a6, b0, temp );
	temp = mul_add (a5, b2, temp );
	temp = mul_add (a4, b3, temp );
	temp = mul_add (a3, b4, temp );
	temp = mul_add (a2, b5, temp );
	temp = mul_add (a0, b6, temp );
		temp=mul_add(a1, b7, temp);
		c_avx[23] = _mm256_add_epi16(temp, c_avx[23]);


	temp = _mm256_mullo_epi16 (a6, b1);
	temp = mul_add (a5, b0, temp );
	temp = mul_add (a4, b2, temp );
	temp = mul_add (a3, b3, temp );
	temp = mul_add (a2, b4, temp );
	temp = mul_add (a0, b5, temp );
		temp=mul_add(a1, b6, temp);
		c_avx[24] = _mm256_add_epi16(temp, c_avx[24]);


	temp = _mm256_mullo_epi16 (a5, b1);
	temp = mul_add (a4, b0, temp );
	temp = mul_add (a3, b2, temp );
	temp = mul_add (a2, b3, temp );
	temp = mul_add (a0, b4, temp );
		temp=mul_add(a1, b5, temp);
		c_avx[25] = _mm256_add_epi16(temp, c_avx[25]);


	temp = _mm256_mullo_epi16 (a4, b1);
	temp = mul_add (a3, b0, temp );
	temp = mul_add (a2, b2, temp );
	temp = mul_add (a0, b3, temp );
		temp=mul_add(a1, b4, temp);
		c_avx[26] = _mm256_add_epi16(temp, c_avx[26]);


	temp = _mm256_mullo_epi16 (a3, b1);
	temp = mul_add (a2, b0, temp );
	temp = mul_add (a0, b2, temp );
		temp=mul_add(a1, b3, temp);
		c_avx[27] = _mm256_add_epi16(temp, c_avx[27]);


	temp = _mm256_mullo_epi16 (a2, b1);
	temp = mul_add (a0, b0, temp );
		temp=mul_add(a1, b2, temp);
		c_avx[28] = _mm256_add_epi16(temp, c_avx[28]);


	temp = _mm256_mullo_epi16 (a0, b1);
		temp=mul_add(a1, b0, temp);
		c_avx[29] = _mm256_add_epi16(temp, c_avx[29]);


		c_avx[30] = mul_add(a1, b1, c_avx[30]);



	c_avx[2*SCM_SIZE-1] = _mm256_set_epi64x(0, 0, 0, 0);


}



void schoolbook_avx_new2(__m256i* a, __m256i* b, __m256i* c_avx) ////8 coefficients of a and b has been prefetched
									      //the c_avx are not added cummulatively
{

	__m256i a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7;
	__m256i temp;


	a0=a[0];
	a1=a[1];
	a2=a[2];
	a3=a[3];
	a4=a[4];
	a5=a[5];
	a6=a[6];
	a7=a[7];

	b0=b[0];
	b1=b[1];
	b2=b[2];
	b3=b[3];
	b4=b[4];
	b5=b[5];
	b6=b[6];
	b7=b[7];

	// New Unrolled first triangle
	c_avx[0] = _mm256_mullo_epi16 (a0, b0);

	temp = _mm256_mullo_epi16 (a0, b1);
	c_avx[1]=mul_add(a1, b0, temp);

	temp = _mm256_mullo_epi16 (a0, b2);

	temp = mul_add(a1, b1, temp);
	c_avx[2]= mul_add(a2, b0, temp);

	temp = _mm256_mullo_epi16 (a0, b3);
	temp = mul_add(a1, b2, temp);
	temp = mul_add(a2, b1, temp);
	c_avx[3]= mul_add(a3, b0, temp);

	temp = _mm256_mullo_epi16 (a0, b4);
	temp = mul_add(a1, b3, temp);
	temp = mul_add(a3, b1, temp);
	temp = mul_add(a4, b0, temp);
	c_avx[4]= mul_add(a2, b2, temp);

	temp = _mm256_mullo_epi16 (a0, b5);
	temp = mul_add(a1, b4 , temp);
	temp = mul_add(a2, b3, temp);
	temp = mul_add(a3, b2, temp);
	temp = mul_add( a4, b1, temp);
	c_avx[5] = mul_add(a5, b0, temp);
	
	temp = _mm256_mullo_epi16 (a0, b6);
	temp = mul_add(a1, b5, temp);
	temp = mul_add(a5, b1, temp);
	temp = mul_add(a6, b0, temp);
	temp = mul_add(a2, b4, temp);
	temp = mul_add(a3, b3, temp);
	c_avx[6] = mul_add(a4, b2, temp);

	temp = _mm256_mullo_epi16 (a0, b7);
	temp = mul_add(a1, b6, temp);
	temp = mul_add (a6, b1, temp);
	temp = mul_add (a7, b0, temp);
	temp = mul_add(a2, b5, temp);
	temp = mul_add (a3, b4, temp);
	temp = mul_add (a4, b3, temp);
	c_avx[7] = mul_add (a5, b2, temp);

	temp = _mm256_mullo_epi16 (a0, b[8]);
	temp = mul_add (a1, b7, temp);
	temp = mul_add (a7, b1, temp);
	temp = mul_add (a[8], b0, temp);
	temp = mul_add (a2, b6,temp);
	temp = mul_add(a3, b5, temp);
	temp = mul_add (a4, b4,temp);
	temp = mul_add (a5, b3, temp);
	c_avx[8] = mul_add (a6, b2, temp);

	temp = _mm256_mullo_epi16 (a0, b[9]);
	temp = mul_add (a1, b[8], temp);
	temp = mul_add (a[8], b1, temp);
	temp = mul_add (a[9], b0, temp);
	temp = mul_add (a2, b7, temp);
	temp = mul_add (a3, b6, temp);
	temp = mul_add (a4, b5, temp);
	temp = mul_add (a5, b4, temp);
	temp = mul_add (a6, b3, temp);
	c_avx[9] = mul_add (a7, b2, temp);

	temp= _mm256_mullo_epi16 (a0, b[10]);
	temp = mul_add (a1, b[9], temp);
	temp = mul_add (a[9], b1, temp);
	temp = mul_add (a[10], b0, temp);
	temp = mul_add (a2, b[8], temp);
	temp = mul_add (a3, b7, temp);
	temp = mul_add (a4, b6, temp);
	temp = mul_add (a5, b5, temp);
	temp = mul_add (a6, b4, temp);
	temp = mul_add (a7, b3, temp);
	c_avx[10] = mul_add (a[8], b2, temp);

	temp = _mm256_mullo_epi16 (a0, b[11]);
	temp = mul_add (a1, b[10], temp );
	temp = mul_add (a[10], b1, temp );
	temp = mul_add (a[11], b0, temp );
	temp = mul_add (a2, b[9], temp );
	temp = mul_add (a3, b[8], temp );
	temp = mul_add (a4, b7, temp );
	temp = mul_add (a5, b6, temp );
	temp = mul_add (a6, b5, temp );
	temp = mul_add (a7, b4, temp );
	temp = mul_add (a[8], b3, temp );
	c_avx[11] = mul_add (a[9], b2, temp );

	temp = _mm256_mullo_epi16 (a0, b[12]);
	temp = mul_add (a1, b[11], temp);
	temp = mul_add (a[11], b1, temp);
	temp = mul_add (a[12], b0, temp);
	temp = mul_add (a2, b[10], temp);
	temp = mul_add (a3, b[9], temp);
	temp = mul_add (a4, b[8], temp);
	temp = mul_add (a5, b7, temp);
	temp = mul_add (a6, b6, temp);
	temp = mul_add (a7, b5, temp);
	temp = mul_add (a[8], b4, temp);
	temp = mul_add (a[9], b3, temp);
	c_avx[12] = mul_add (a[10], b2, temp);

	temp = _mm256_mullo_epi16 (a0, b[13]);
	temp = mul_add (a1, b[12], temp );
	temp = mul_add (a[12], b1, temp );
	temp = mul_add (a[13], b0, temp );
	temp = mul_add (a2, b[11], temp );
	temp = mul_add (a3, b[10], temp );
	temp = mul_add (a4, b[9], temp );
	temp = mul_add (a5, b[8], temp );
	temp = mul_add (a6, b7, temp );
	temp = mul_add (a7, b6, temp );
	temp = mul_add (a[8], b5, temp );
	temp = mul_add (a[9], b4, temp );
	temp = mul_add (a[10], b3, temp );
	c_avx[13] = mul_add (a[11], b2, temp );

	temp = _mm256_mullo_epi16 (a0, b[14]);
	temp = mul_add (a1, b[13], temp );
	temp = mul_add (a[13], b1, temp );
	temp = mul_add (a[14], b0, temp );
	temp = mul_add (a2, b[12], temp );
	temp = mul_add (a3, b[11], temp );
	temp = mul_add (a4, b[10], temp );
	temp = mul_add (a5, b[9], temp );
	temp = mul_add (a6, b[8], temp );
	temp = mul_add (a7, b7, temp );
	temp = mul_add (a[8], b6, temp );
	temp = mul_add (a[9], b5, temp );
	temp = mul_add (a[10], b4, temp );
	temp = mul_add (a[11], b3, temp );
	c_avx[14] = mul_add (a[12], b2, temp );

	temp = _mm256_mullo_epi16 (a0, b[15]);
	temp = mul_add (a1, b[14], temp );
	temp = mul_add (a[14], b1, temp );
	temp = mul_add (a[15], b0, temp );
	temp = mul_add (a2, b[13], temp );
	temp = mul_add (a3, b[12], temp );
	temp = mul_add (a4, b[11], temp );
	temp = mul_add (a5, b[10], temp );
	temp = mul_add (a6, b[9], temp );
	temp = mul_add (a7, b[8], temp );
	temp = mul_add (a[8], b7, temp );
	temp = mul_add (a[9], b6, temp );
	temp = mul_add (a[10], b5, temp );
	temp = mul_add (a[11], b4, temp );
	temp = mul_add (a[12], b3, temp );
	c_avx[15] = mul_add (a[13], b2, temp );


	// unrolled second triangle
	a0=a[14];
	a1=a[15];
	a2=a[13];
	a3=a[12];
	a4=a[11];
	a5=a[10];
	a6=a[9];
	a7=a[8];

	b0=b[14];
	b1=b[15];
	b2=b[13];
	b3=b[12];
	b4=b[11];
	b5=b[10];
	b6=b[9];
	b7=b[8];
	

	temp = _mm256_mullo_epi16 (a[1], b1);
	temp = mul_add (a[2], b0, temp );
	temp = mul_add (a[3], b2, temp );
	temp = mul_add (a[4], b3, temp );
	temp = mul_add (a[5], b4, temp );
	temp = mul_add (a[6], b5, temp );
	temp = mul_add (a[7], b6, temp );
	temp = mul_add (a7, b7, temp );
	temp = mul_add (a6, b[7], temp );
	temp = mul_add (a5, b[6], temp );
	temp = mul_add (a4, b[5], temp );
	temp = mul_add (a3, b[4], temp );
	temp = mul_add (a2, b[3], temp );
	temp = mul_add (a0, b[2], temp );
	c_avx[16] = mul_add (a1, b[1], temp );

	temp = _mm256_mullo_epi16 (a[2], b1);
	temp = mul_add (a[3], b0, temp );
	temp = mul_add (a[4], b2, temp );
	temp = mul_add (a[5], b3, temp );
	temp = mul_add (a[6], b4, temp );
	temp = mul_add (a[7], b5, temp );
	temp = mul_add (a7, b6, temp );
	temp = mul_add (a6, b7, temp );
	temp = mul_add (a5, b[7], temp );
	temp = mul_add (a4, b[6], temp );
	temp = mul_add (a3, b[5], temp );
	temp = mul_add (a2, b[4], temp );
	temp = mul_add (a0, b[3], temp );
	c_avx[17] = mul_add (a1, b[2], temp );

	temp = _mm256_mullo_epi16 (a[3], b1);
	temp = mul_add (a[4], b0, temp );
	temp = mul_add (a[5], b2, temp );
	temp = mul_add (a[6], b3, temp );
	temp = mul_add (a[7], b4, temp );
	temp = mul_add (a7, b5, temp );
	temp = mul_add (a6, b6, temp );
	temp = mul_add (a5, b7, temp );
	temp = mul_add (a4, b[7], temp );
	temp = mul_add (a3, b[6], temp );
	temp = mul_add (a2, b[5], temp );
	temp = mul_add (a0, b[4], temp );
	c_avx[18] = mul_add (a1, b[3], temp );

	temp = _mm256_mullo_epi16 (a[4], b1);
	temp = mul_add (a[5], b0, temp );
	temp = mul_add (a[6], b2, temp );
	temp = mul_add (a[7], b3, temp );
	temp = mul_add (a7, b4, temp );
	temp = mul_add (a6, b5, temp );
	temp = mul_add (a5, b6, temp );
	temp = mul_add (a4, b7, temp );
	temp = mul_add (a3, b[7], temp );
	temp = mul_add (a2, b[6], temp );
	temp = mul_add (a0, b[5], temp );
	c_avx[19] = mul_add (a1, b[4], temp );

	temp = _mm256_mullo_epi16 (a[5], b1);
	temp = mul_add (a[6], b0, temp );
	temp = mul_add (a[7], b2, temp );
	temp = mul_add (a7, b3, temp );
	temp = mul_add (a6, b4, temp );
	temp = mul_add (a5, b5, temp );
	temp = mul_add (a4, b6, temp );
	temp = mul_add (a3, b7, temp );
	temp = mul_add (a2, b[7], temp );
	temp = mul_add (a0, b[6], temp );
	c_avx[20] = mul_add (a1, b[5], temp );

	temp = _mm256_mullo_epi16 (a[6], b1);
	temp = mul_add (a[7], b0, temp );
	temp = mul_add (a7, b2, temp );
	temp = mul_add (a6, b3, temp );
	temp = mul_add (a5, b4, temp );
	temp = mul_add (a4, b5, temp );
	temp = mul_add (a3, b6, temp );
	temp = mul_add (a2, b7, temp );
	temp = mul_add (a0, b[7], temp );
	c_avx[21] = mul_add (a1, b[6], temp );

	temp = _mm256_mullo_epi16 (a[7], b1);
	temp = mul_add (a7, b0, temp );
	temp = mul_add (a6, b2, temp );
	temp = mul_add (a5, b3, temp );
	temp = mul_add (a4, b4, temp );
	temp = mul_add (a3, b5, temp );
	temp = mul_add (a2, b6, temp );
	temp = mul_add (a0, b7, temp );
	c_avx[22] = mul_add (a1, b[7], temp );

	temp = _mm256_mullo_epi16 (a7, b1);
	temp = mul_add (a6, b0, temp );
	temp = mul_add (a5, b2, temp );
	temp = mul_add (a4, b3, temp );
	temp = mul_add (a3, b4, temp );
	temp = mul_add (a2, b5, temp );
	temp = mul_add (a0, b6, temp );
	c_avx[23] = mul_add (a1, b7, temp );

	temp = _mm256_mullo_epi16 (a6, b1);
	temp = mul_add (a5, b0, temp );
	temp = mul_add (a4, b2, temp );
	temp = mul_add (a3, b3, temp );
	temp = mul_add (a2, b4, temp );
	temp = mul_add (a0, b5, temp );
	c_avx[24] = mul_add (a1, b6, temp );

	temp = _mm256_mullo_epi16 (a5, b1);
	temp = mul_add (a4, b0, temp );
	temp = mul_add (a3, b2, temp );
	temp = mul_add (a2, b3, temp );
	temp = mul_add (a0, b4, temp );
	c_avx[25] = mul_add (a1, b5, temp );

	temp = _mm256_mullo_epi16 (a4, b1);
	temp = mul_add (a3, b0, temp );
	temp = mul_add (a2, b2, temp );
	temp = mul_add (a0, b3, temp );
	c_avx[26] = mul_add (a1, b4, temp );

	temp = _mm256_mullo_epi16 (a3, b1);
	temp = mul_add (a2, b0, temp );
	temp = mul_add (a0, b2, temp );
	c_avx[27] = mul_add (a1, b3, temp );

	temp = _mm256_mullo_epi16 (a2, b1);
	temp = mul_add (a0, b0, temp );
	c_avx[28] = mul_add (a1, b2, temp );

	temp = _mm256_mullo_epi16 (a0, b1);
	c_avx[29] = mul_add (a1, b0, temp);

	c_avx[30] = _mm256_mullo_epi16 (a1, b1);


	c_avx[2*SCM_SIZE-1] = _mm256_set_epi64x(0, 0, 0, 0);

}
void TC_interpol(__m256i *c_bucket, __m256i* res_avx_output);
void KARA_interpol(__m256i *c_bucket, __m256i* result_final0, __m256i* result_final1, __m256i* result_final2, __m256i* result_final3, __m256i* result_final4, __m256i* result_final5, __m256i* result_final6);
void KARA_eval(__m256i* b, __m256i *b_bucket);
void TC_eval(__m256i* b_avx, __m256i* b_bucket);


void batch_64coefficient_multiplications_new(__m256i* a, __m256i* b_bucket, __m256i* c_bucket, int f)//all 7 Karatsuba evaluation and interpolation are done in AVX.
{
	__m256i a_bucket[SCM_SIZE*4]; //SCM_SIZE = 16; Holds evaluation (a & b) for 7 Karatsuba at a time

	//uint16_t i;

	register __m256i r0_avx, r1_avx, r2_avx, r3_avx;



		//CLOCK1=cpucycles();
		
		//------------------AVX evaluation for 1st poly-----------------------

                    r0_avx=a[0];
                    r1_avx=a[1];
                    r2_avx=a[2];
                    r3_avx=a[3];
		    a_bucket[0]=r0_avx;
		    a_bucket[1]=r1_avx;
		    a_bucket[2]=r2_avx;
		    a_bucket[3]=r3_avx;
		    a_bucket[4]= _mm256_add_epi16(r0_avx, r1_avx);
		    a_bucket[5]= _mm256_add_epi16(r2_avx, r3_avx);
		    a_bucket[6]= _mm256_add_epi16(r0_avx, r2_avx);
		    a_bucket[7]= _mm256_add_epi16(r1_avx, r3_avx);
		    a_bucket[8]= _mm256_add_epi16(a_bucket[6],a_bucket[7]);


		//------------------AVX evaluation for 1st poly ends------------------


		//------------------AVX evaluation for 2nd poly-----------------------
                    r0_avx=a[small_len_avx];
                    r1_avx=a[small_len_avx+1];
                    r2_avx=a[small_len_avx+2];
                    r3_avx=a[small_len_avx+3];
		    a_bucket[0+9]=r0_avx;
		    a_bucket[1+9]=r1_avx;
		    a_bucket[2+9]=r2_avx;
		    a_bucket[3+9]=r3_avx;
		    a_bucket[4+9]= _mm256_add_epi16(r0_avx, r1_avx);
		    a_bucket[5+9]= _mm256_add_epi16(r2_avx, r3_avx);
		    a_bucket[6+9]= _mm256_add_epi16(r0_avx, r2_avx);
		    a_bucket[7+9]= _mm256_add_epi16(r1_avx, r3_avx);
		    a_bucket[8+9]= _mm256_add_epi16(a_bucket[6+9],a_bucket[7+9]);

	
		//------------------AVX evaluation for 2nd poly ends------------------


		//------------------AVX evaluation for 3rd poly-----------------------
                    r0_avx=a[2*small_len_avx];
                    r1_avx=a[2*small_len_avx+1];
                    r2_avx=a[2*small_len_avx+2];
                    r3_avx=a[2*small_len_avx+3];
		    a_bucket[0+18]=r0_avx;
		    a_bucket[1+18]=r1_avx;
		    a_bucket[2+18]=r2_avx;
		    a_bucket[3+18]=r3_avx;
		    a_bucket[4+18]= _mm256_add_epi16(r0_avx, r1_avx);
		    a_bucket[5+18]= _mm256_add_epi16(r2_avx, r3_avx);
		    a_bucket[6+18]= _mm256_add_epi16(r0_avx, r2_avx);
		    a_bucket[7+18]= _mm256_add_epi16(r1_avx, r3_avx);
		    a_bucket[8+18]= _mm256_add_epi16(a_bucket[6+18],a_bucket[7+18]);
		
		//------------------AVX evaluation for 3rd poly ends------------------


		//------------------AVX evaluation for 4th poly-----------------------

                    r0_avx=a[3*small_len_avx];
                    r1_avx=a[3*small_len_avx+1];
                    r2_avx=a[3*small_len_avx+2];
                    r3_avx=a[3*small_len_avx+3];
		    a_bucket[0+27]=r0_avx;
		    a_bucket[1+27]=r1_avx;
		    a_bucket[2+27]=r2_avx;
		    a_bucket[3+27]=r3_avx;
		    a_bucket[4+27]= _mm256_add_epi16(r0_avx, r1_avx);
		    a_bucket[5+27]= _mm256_add_epi16(r2_avx, r3_avx);
		    a_bucket[6+27]= _mm256_add_epi16(r0_avx, r2_avx);
		    a_bucket[7+27]= _mm256_add_epi16(r1_avx, r3_avx);
		    a_bucket[8+27]= _mm256_add_epi16(a_bucket[6+27],a_bucket[7+27]);
		
		//------------------AVX evaluation for 4th poly ends------------------

		//------------------AVX evaluation for 5th poly-----------------------
		
                    r0_avx=a[4*small_len_avx+0];
                    r1_avx=a[4*small_len_avx+1];
                    r2_avx=a[4*small_len_avx+2];
                    r3_avx=a[4*small_len_avx+3];
		    a_bucket[0+36]=r0_avx;
		    a_bucket[1+36]=r1_avx;
		    a_bucket[2+36]=r2_avx;
		    a_bucket[3+36]=r3_avx;
		    a_bucket[4+36]= _mm256_add_epi16(r0_avx, r1_avx);
		    a_bucket[5+36]= _mm256_add_epi16(r2_avx, r3_avx);
		    a_bucket[6+36]= _mm256_add_epi16(r0_avx, r2_avx);
		    a_bucket[7+36]= _mm256_add_epi16(r1_avx, r3_avx);
		    a_bucket[8+36]= _mm256_add_epi16(a_bucket[6+36],a_bucket[7+36]);
		
		//------------------AVX evaluation for 5th poly ends------------------


		//------------------AVX evaluation for 6th poly-----------------------
                    r0_avx=a[5*small_len_avx];
                    r1_avx=a[5*small_len_avx+1];
                    r2_avx=a[5*small_len_avx+2];
                    r3_avx=a[5*small_len_avx+3];
		    a_bucket[0+45]=r0_avx;
		    a_bucket[1+45]=r1_avx;
		    a_bucket[2+45]=r2_avx;
		    a_bucket[3+45]=r3_avx;
		    a_bucket[4+45]= _mm256_add_epi16(r0_avx, r1_avx);
		    a_bucket[5+45]= _mm256_add_epi16(r2_avx, r3_avx);
		    a_bucket[6+45]= _mm256_add_epi16(r0_avx, r2_avx);
		    a_bucket[7+45]= _mm256_add_epi16(r1_avx, r3_avx);
		    a_bucket[8+45]= _mm256_add_epi16(a_bucket[6+45],a_bucket[7+45]);
		
		//------------------AVX evaluation for 6th poly ends------------------

		//------------------AVX evaluation for 7th poly-----------------------

                    r0_avx=a[6*small_len_avx];
                    r1_avx=a[6*small_len_avx+1];
                    r2_avx=a[6*small_len_avx+2];
                    r3_avx=a[6*small_len_avx+3];
		    a_bucket[0+54]=r0_avx;
		    a_bucket[1+54]=r1_avx;
		    a_bucket[2+54]=r2_avx;
		    a_bucket[3+54]=r3_avx;
		    a_bucket[4+54]= _mm256_add_epi16(r0_avx, r1_avx);
		    a_bucket[5+54]= _mm256_add_epi16(r2_avx, r3_avx);
		    a_bucket[6+54]= _mm256_add_epi16(r0_avx, r2_avx);
		    a_bucket[7+54]= _mm256_add_epi16(r1_avx, r3_avx);
		    a_bucket[8+54]= _mm256_add_epi16(a_bucket[6+54],a_bucket[7+54]);

		//------------------AVX evaluation for 7th poly ends------------------
		
	

		//CLOCK2=cpucycles();
		//CLOCK_EVAL=CLOCK_EVAL+(CLOCK2-CLOCK1);
		//printf("\nTime for multiplication : %llu\n", CLOCK2-CLOCK1);


		//CLOCK1=cpucycles();
		//-----------------Forward transposes--------------------------------------
			transpose_n1(a_bucket);
			transpose_n1(a_bucket+16);
			transpose_n1(a_bucket+32);
			transpose_n1(a_bucket+48);

		//-----------------Forwatrd transposes ends---------------------------------

		//----------------------all multiplications---------------------------------
		if(f==0){
			schoolbook_avx_new2(a_bucket, b_bucket, c_bucket);
			schoolbook_avx_new2(a_bucket+16, b_bucket+16, c_bucket+2*SCM_SIZE);
			schoolbook_avx_new2(a_bucket+32, b_bucket+32, c_bucket+4*SCM_SIZE);
			schoolbook_avx_new2(a_bucket+48, b_bucket+48, c_bucket+6*SCM_SIZE);
		}
		else{
			schoolbook_avx_new3_acc(a_bucket, b_bucket, c_bucket);
			schoolbook_avx_new3_acc(a_bucket+16, b_bucket+16, c_bucket+2*SCM_SIZE);
			//schoolbook_avx_new3_acc_fused(a_bucket, b_bucket, c_bucket);
			schoolbook_avx_new3_acc(a_bucket+32, b_bucket+32, c_bucket+4*SCM_SIZE);
			schoolbook_avx_new3_acc(a_bucket+48, b_bucket+48, c_bucket+6*SCM_SIZE);
		}
		/*
		schoolbook_avx_new2_acc(a_bucket, b_bucket, c_bucket, f);
		schoolbook_avx_new2_acc(a_bucket+16, b_bucket+16, c_bucket+2*SCM_SIZE, f);
		schoolbook_avx_new2_acc(a_bucket+32, b_bucket+32, c_bucket+4*SCM_SIZE, f);
		schoolbook_avx_new2_acc(a_bucket+48, b_bucket+48, c_bucket+6*SCM_SIZE, f);
		*/


		//----------------------all multiplications ends-----------------------------


		//-----------------Reverse transposes--------------------------------------

			/*
			transpose(c_bucket);
			transpose(c_bucket+16);

			transpose(c_bucket+2*SCM_SIZE);
			transpose(c_bucket+16+2*SCM_SIZE);

			transpose(c_bucket+4*SCM_SIZE);
			transpose(c_bucket+16+4*SCM_SIZE);

			transpose(c_bucket+6*SCM_SIZE);
			transpose(c_bucket+16+6*SCM_SIZE);
			*/
		//-----------------Reverse transposes ends---------------------------------

		//CLOCK2=cpucycles();
		//CLOCK_MULT=CLOCK_MULT+(CLOCK2-CLOCK1);

		//KARA_interpol(c_bucket, result_final0, result_final1, result_final2, result_final3, result_final4, result_final5, result_final6);
		
}

void KARA_eval(__m256i* b, __m256i *b_bucket){

	__m256i r0_avx, r1_avx, r2_avx, r3_avx;


		//-------1st poly----------------------------------------------------
                    r0_avx=b[0];
                    r1_avx=b[1];
                    r2_avx=b[2];
                    r3_avx=b[3];
		    b_bucket[0]=r0_avx;
		    b_bucket[1]=r1_avx;
		    b_bucket[2]=r2_avx;
		    b_bucket[3]=r3_avx;
		    b_bucket[4]= _mm256_add_epi16(r0_avx, r1_avx);
		    b_bucket[5]= _mm256_add_epi16(r2_avx, r3_avx);
		    b_bucket[6]= _mm256_add_epi16(r0_avx, r2_avx);
		    b_bucket[7]= _mm256_add_epi16(r1_avx, r3_avx);
		    b_bucket[8]= _mm256_add_epi16(b_bucket[6],b_bucket[7]);
		//-------2nd poly----------------------------------------------------

                    r0_avx=b[small_len_avx];
                    r1_avx=b[small_len_avx+1];
                    r2_avx=b[small_len_avx+2];
                    r3_avx=b[small_len_avx+3];
		    b_bucket[0+9]=r0_avx;
		    b_bucket[1+9]=r1_avx;
		    b_bucket[2+9]=r2_avx;
		    b_bucket[3+9]=r3_avx;
		    b_bucket[4+9]= _mm256_add_epi16(r0_avx, r1_avx);
		    b_bucket[5+9]= _mm256_add_epi16(r2_avx, r3_avx);
		    b_bucket[6+9]= _mm256_add_epi16(r0_avx, r2_avx);
		    b_bucket[7+9]= _mm256_add_epi16(r1_avx, r3_avx);
		    b_bucket[8+9]= _mm256_add_epi16(b_bucket[6+9],b_bucket[7+9]);

		//-------3rd poly----------------------------------------------------

                    r0_avx=b[2*small_len_avx+0];
                    r1_avx=b[2*small_len_avx+1];
                    r2_avx=b[2*small_len_avx+2];
                    r3_avx=b[2*small_len_avx+3];
		    b_bucket[0+18]=r0_avx;
		    b_bucket[1+18]=r1_avx;
		    b_bucket[2+18]=r2_avx;
		    b_bucket[3+18]=r3_avx;
		    b_bucket[4+18]= _mm256_add_epi16(r0_avx, r1_avx);
		    b_bucket[5+18]= _mm256_add_epi16(r2_avx, r3_avx);
		    b_bucket[6+18]= _mm256_add_epi16(r0_avx, r2_avx);
		    b_bucket[7+18]= _mm256_add_epi16(r1_avx, r3_avx);
		    b_bucket[8+18]= _mm256_add_epi16(b_bucket[6+18],b_bucket[7+18]);

		//-------4th poly----------------------------------------------------
                    r0_avx=b[3*small_len_avx];
                    r1_avx=b[3*small_len_avx+1];
                    r2_avx=b[3*small_len_avx+2];
                    r3_avx=b[3*small_len_avx+3];
		    b_bucket[0+27]=r0_avx;
		    b_bucket[1+27]=r1_avx;
		    b_bucket[2+27]=r2_avx;
		    b_bucket[3+27]=r3_avx;
		    b_bucket[4+27]= _mm256_add_epi16(r0_avx, r1_avx);
		    b_bucket[5+27]= _mm256_add_epi16(r2_avx, r3_avx);
		    b_bucket[6+27]= _mm256_add_epi16(r0_avx, r2_avx);
		    b_bucket[7+27]= _mm256_add_epi16(r1_avx, r3_avx);
		    b_bucket[8+27]= _mm256_add_epi16(b_bucket[6+27],b_bucket[7+27]);

		//-------5th poly----------------------------------------------------

                    r0_avx=b[4*small_len_avx];
                    r1_avx=b[4*small_len_avx+1];
                    r2_avx=b[4*small_len_avx+2];
                    r3_avx=b[4*small_len_avx+3];
		    b_bucket[0+36]=r0_avx;
		    b_bucket[1+36]=r1_avx;
		    b_bucket[2+36]=r2_avx;
		    b_bucket[3+36]=r3_avx;
		    b_bucket[4+36]= _mm256_add_epi16(r0_avx, r1_avx);
		    b_bucket[5+36]= _mm256_add_epi16(r2_avx, r3_avx);
		    b_bucket[6+36]= _mm256_add_epi16(r0_avx, r2_avx);
		    b_bucket[7+36]= _mm256_add_epi16(r1_avx, r3_avx);
		    b_bucket[8+36]= _mm256_add_epi16(b_bucket[6+36],b_bucket[7+36]);

		//-------6th poly----------------------------------------------------

                    r0_avx=b[5*small_len_avx];
                    r1_avx=b[5*small_len_avx+1];
                    r2_avx=b[5*small_len_avx+2];
                    r3_avx=b[5*small_len_avx+3];
		    b_bucket[0+45]=r0_avx;
		    b_bucket[1+45]=r1_avx;
		    b_bucket[2+45]=r2_avx;
		    b_bucket[3+45]=r3_avx;
		    b_bucket[4+45]= _mm256_add_epi16(r0_avx, r1_avx);
		    b_bucket[5+45]= _mm256_add_epi16(r2_avx, r3_avx);
		    b_bucket[6+45]= _mm256_add_epi16(r0_avx, r2_avx);
		    b_bucket[7+45]= _mm256_add_epi16(r1_avx, r3_avx);
		    b_bucket[8+45]= _mm256_add_epi16(b_bucket[6+45],b_bucket[7+45]);

		//-------7th poly----------------------------------------------------

                    r0_avx=b[6*small_len_avx];
                    r1_avx=b[6*small_len_avx+1];
                    r2_avx=b[6*small_len_avx+2];
                    r3_avx=b[6*small_len_avx+3];
		    b_bucket[0+54]=r0_avx;
		    b_bucket[1+54]=r1_avx;
		    b_bucket[2+54]=r2_avx;
		    b_bucket[3+54]=r3_avx;
		    b_bucket[4+54]= _mm256_add_epi16(r0_avx, r1_avx);
		    b_bucket[5+54]= _mm256_add_epi16(r2_avx, r3_avx);
		    b_bucket[6+54]= _mm256_add_epi16(r0_avx, r2_avx);
		    b_bucket[7+54]= _mm256_add_epi16(r1_avx, r3_avx);
		    b_bucket[8+54]= _mm256_add_epi16(b_bucket[6+54],b_bucket[7+54]);

		//--------------Evaluating B poly ends-------------------------------

			transpose_n1(b_bucket);
			transpose_n1(b_bucket+16);
			transpose_n1(b_bucket+32);
			transpose_n1(b_bucket+48);	
}

void KARA_interpol(__m256i *c_bucket, __m256i* result_final0, __m256i* result_final1, __m256i* result_final2, __m256i* result_final3, __m256i* result_final4, __m256i* result_final5, __m256i* result_final6){

		//int64_t i;
		register __m256i res_avx0, res_avx1, res_avx2, res_avx3, res_avx4, res_avx5, res_avx6, res_avx7; // to hold each 64X64 poly mul results

		__m256i temp, c6_avx, c7_avx, c8_avx, c20_avx, c21_avx, c22_avx, c23_avx, c24_avx;

		//CLOCK1=cpucycles();

		   //------------------------AVX interpolation for 1st poly external-------------------		
			
			//loop1
				res_avx0 = c_bucket[0];
				res_avx2 = c_bucket[1];
				res_avx4 = c_bucket[2];
				res_avx6 = c_bucket[3];

				c6_avx=c_bucket[6];
				c7_avx=c_bucket[7];
		
				c8_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[8], c6_avx), c7_avx);

				res_avx1 = c_bucket[16];
				res_avx3 = c_bucket[17];
				res_avx5 = c_bucket[18];
				res_avx7 = c_bucket[19];

				c22_avx=c_bucket[22];
				c23_avx=c_bucket[23];

				c21_avx=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[21], res_avx5),res_avx7);

				c24_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[24], c22_avx), c23_avx);

				c20_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[20], res_avx1), res_avx3);

				temp=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[5], res_avx4),res_avx6);
				res_avx5 = _mm256_add_epi16(res_avx5, temp);

				temp = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[4], res_avx0), res_avx2);
				res_avx1 = _mm256_add_epi16(res_avx1, temp);

				c22_avx=_mm256_add_epi16(c22_avx, c8_avx);

				res_avx6 = _mm256_add_epi16(res_avx6, c21_avx);

				res_avx2 = _mm256_add_epi16(res_avx2, c20_avx);

				c7_avx=_mm256_add_epi16(c7_avx, c24_avx);


			//loop4

				c6_avx=_mm256_sub_epi16(_mm256_sub_epi16(c6_avx, res_avx0), res_avx4);
				c22_avx=_mm256_sub_epi16(_mm256_sub_epi16(c22_avx, res_avx1), res_avx5);

				c7_avx=_mm256_sub_epi16(_mm256_sub_epi16(c7_avx, res_avx2), res_avx6);
				c23_avx=_mm256_sub_epi16(_mm256_sub_epi16(c23_avx, res_avx3), res_avx7);

			//loop5
				result_final0[0]=res_avx0;
				result_final0[1]=res_avx1;

				result_final0[2]=_mm256_add_epi16(res_avx2, c6_avx);
				result_final0[3]=_mm256_add_epi16(res_avx3, c22_avx);


				result_final0[4]=_mm256_add_epi16(res_avx4, c7_avx);
				result_final0[5]=_mm256_add_epi16(res_avx5, c23_avx);

				result_final0[6]=res_avx6;
				result_final0[7]=res_avx7;


		   //------------------------AVX interpolation for 1st poly ends--------------


		   //------------------------AVX interpolation for 2nd poly external-------------------		
			
			//loop1
				res_avx0 = c_bucket[9]; //c_bucket0
				res_avx2 = c_bucket[10]; //c_bucket1
				res_avx4 = c_bucket[11]; //c_bucket2
				res_avx6 = c_bucket[12]; //c_bucket3

				c6_avx=c_bucket[15]; //c_bucket6
				c7_avx=c_bucket[32]; //c_bucket7
		
				c8_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[33], c6_avx), c7_avx);

				res_avx1 = c_bucket[25]; //c_bucket0
				res_avx3 = c_bucket[26]; //c_bucket1
				res_avx5 = c_bucket[27]; //c_bucket2
				res_avx7 = c_bucket[28]; //c_bucket3

				c22_avx=c_bucket[31];
				c23_avx=c_bucket[48];

				c21_avx=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[30], res_avx5),res_avx7);

				c24_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[49], c22_avx), c23_avx);

				c20_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[29], res_avx1), res_avx3);

				temp=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[14], res_avx4),res_avx6);
				res_avx5 = _mm256_add_epi16(res_avx5, temp);

				temp = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[13], res_avx0), res_avx2);
				res_avx1 = _mm256_add_epi16(res_avx1, temp);

				c22_avx=_mm256_add_epi16(c22_avx, c8_avx);

				res_avx6 = _mm256_add_epi16(res_avx6, c21_avx);

				res_avx2 = _mm256_add_epi16(res_avx2, c20_avx);

				c7_avx=_mm256_add_epi16(c7_avx, c24_avx);


			//loop4

				c6_avx=_mm256_sub_epi16(_mm256_sub_epi16(c6_avx, res_avx0), res_avx4);
				c22_avx=_mm256_sub_epi16(_mm256_sub_epi16(c22_avx, res_avx1), res_avx5);

				c7_avx=_mm256_sub_epi16(_mm256_sub_epi16(c7_avx, res_avx2), res_avx6);
				c23_avx=_mm256_sub_epi16(_mm256_sub_epi16(c23_avx, res_avx3), res_avx7);

			//loop5
				result_final1[0]=res_avx0;
				result_final1[1]=res_avx1;

				result_final1[2]=_mm256_add_epi16(res_avx2, c6_avx);
				result_final1[3]=_mm256_add_epi16(res_avx3, c22_avx);


				result_final1[4]=_mm256_add_epi16(res_avx4, c7_avx);
				result_final1[5]=_mm256_add_epi16(res_avx5, c23_avx);

				result_final1[6]=res_avx6;
				result_final1[7]=res_avx7;


		   //------------------------AVX interpolation for 2nd poly ends--------------

		   //------------------------AVX interpolation for 3rd poly external-------------------		
			
			//loop1
				res_avx0 = c_bucket[34]; //c_bucket0
				res_avx2 = c_bucket[35]; //c_bucket1
				res_avx4 = c_bucket[36];
				res_avx6 = c_bucket[37];

				c6_avx=c_bucket[40];
				c7_avx=c_bucket[41];
		
				c8_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[42], c6_avx), c7_avx);

				res_avx1 = c_bucket[50]; //c_bucket0
				res_avx3 = c_bucket[51]; //c_bucket1
				res_avx5 = c_bucket[52];
				res_avx7 = c_bucket[53];

				c22_avx=c_bucket[56];
				c23_avx=c_bucket[57];

				c21_avx=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[55], res_avx5),res_avx7);

				c24_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[58], c22_avx), c23_avx);

				c20_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[54], res_avx1), res_avx3);

				temp=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[39], res_avx4),res_avx6);
				res_avx5 = _mm256_add_epi16(res_avx5, temp);

				temp = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[38], res_avx0), res_avx2);
				res_avx1 = _mm256_add_epi16(res_avx1, temp);

				c22_avx=_mm256_add_epi16(c22_avx, c8_avx);

				res_avx6 = _mm256_add_epi16(res_avx6, c21_avx);

				res_avx2 = _mm256_add_epi16(res_avx2, c20_avx);

				c7_avx=_mm256_add_epi16(c7_avx, c24_avx);

			//loop4
				c6_avx=_mm256_sub_epi16(_mm256_sub_epi16(c6_avx, res_avx0), res_avx4);
				c22_avx=_mm256_sub_epi16(_mm256_sub_epi16(c22_avx, res_avx1), res_avx5);

				c7_avx=_mm256_sub_epi16(_mm256_sub_epi16(c7_avx, res_avx2), res_avx6);
				c23_avx=_mm256_sub_epi16(_mm256_sub_epi16(c23_avx, res_avx3), res_avx7);
			//loop5
				result_final2[0]=res_avx0;
				result_final2[1]=res_avx1;

				result_final2[2]=_mm256_add_epi16(res_avx2, c6_avx);
				result_final2[3]=_mm256_add_epi16(res_avx3, c22_avx);


				result_final2[4]=_mm256_add_epi16(res_avx4, c7_avx);
				result_final2[5]=_mm256_add_epi16(res_avx5, c23_avx);

				result_final2[6]=res_avx6;
				result_final2[7]=res_avx7;

		   //------------------------AVX interpolation for 3rd poly ends--------------
		
		   //------------------------AVX interpolation for 4th poly external-------------------		
			
			//loop1
				res_avx0 = c_bucket[43];
				res_avx2 = c_bucket[44];
				res_avx4 = c_bucket[45];
				res_avx6 = c_bucket[46];

				c6_avx=c_bucket[65];
				c7_avx=c_bucket[66];
		
				c8_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[67], c6_avx), c7_avx);

				res_avx1 = c_bucket[59];
				res_avx3 = c_bucket[60];
				res_avx5 = c_bucket[61];
				res_avx7 = c_bucket[62];

				c22_avx=c_bucket[81];
				c23_avx=c_bucket[82];

				c21_avx=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[80], res_avx5),res_avx7);

				c24_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[83], c22_avx), c23_avx);

				c20_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[63], res_avx1), res_avx3);

				temp=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[64], res_avx4),res_avx6);
				res_avx5 = _mm256_add_epi16(res_avx5, temp);

				temp = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[47], res_avx0), res_avx2);
				res_avx1 = _mm256_add_epi16(res_avx1, temp);

				c22_avx=_mm256_add_epi16(c22_avx, c8_avx);

				res_avx6 = _mm256_add_epi16(res_avx6, c21_avx);

				res_avx2 = _mm256_add_epi16(res_avx2, c20_avx);

				c7_avx=_mm256_add_epi16(c7_avx, c24_avx);


			//loop4

				c6_avx=_mm256_sub_epi16(_mm256_sub_epi16(c6_avx, res_avx0), res_avx4);
				c22_avx=_mm256_sub_epi16(_mm256_sub_epi16(c22_avx, res_avx1), res_avx5);

				c7_avx=_mm256_sub_epi16(_mm256_sub_epi16(c7_avx, res_avx2), res_avx6);
				c23_avx=_mm256_sub_epi16(_mm256_sub_epi16(c23_avx, res_avx3), res_avx7);

			//loop5
				result_final3[0]=res_avx0;
				result_final3[1]=res_avx1;

				result_final3[2]=_mm256_add_epi16(res_avx2, c6_avx);
				result_final3[3]=_mm256_add_epi16(res_avx3, c22_avx);


				result_final3[4]=_mm256_add_epi16(res_avx4, c7_avx);
				result_final3[5]=_mm256_add_epi16(res_avx5, c23_avx);

				result_final3[6]=res_avx6;
				result_final3[7]=res_avx7;


		   //------------------------AVX interpolation for 4th poly ends--------------

		   //------------------------AVX interpolation for 5th poly external-------------------		
			
			//loop1
				res_avx0 = c_bucket[68];
				res_avx2 = c_bucket[69];
				res_avx4 = c_bucket[70];
				res_avx6 = c_bucket[71];

				c6_avx=c_bucket[74];
				c7_avx=c_bucket[75];
		
				c8_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[76], c6_avx), c7_avx);

				res_avx1 = c_bucket[84];
				res_avx3 = c_bucket[85];
				res_avx5 = c_bucket[86];
				res_avx7 = c_bucket[87];

				c22_avx=c_bucket[90];
				c23_avx=c_bucket[91];

				c21_avx=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[89], res_avx5),res_avx7);

				c24_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[92], c22_avx), c23_avx);

				c20_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[88], res_avx1), res_avx3);

				temp=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[73], res_avx4),res_avx6);
				res_avx5 = _mm256_add_epi16(res_avx5, temp);

				temp = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[72], res_avx0), res_avx2);
				res_avx1 = _mm256_add_epi16(res_avx1, temp);

				c22_avx=_mm256_add_epi16(c22_avx, c8_avx);

				res_avx6 = _mm256_add_epi16(res_avx6, c21_avx);

				res_avx2 = _mm256_add_epi16(res_avx2, c20_avx);

				c7_avx=_mm256_add_epi16(c7_avx, c24_avx);


			//loop4

				c6_avx=_mm256_sub_epi16(_mm256_sub_epi16(c6_avx, res_avx0), res_avx4);
				c22_avx=_mm256_sub_epi16(_mm256_sub_epi16(c22_avx, res_avx1), res_avx5);

				c7_avx=_mm256_sub_epi16(_mm256_sub_epi16(c7_avx, res_avx2), res_avx6);
				c23_avx=_mm256_sub_epi16(_mm256_sub_epi16(c23_avx, res_avx3), res_avx7);

			//loop5
				result_final4[0]=res_avx0;
				result_final4[1]=res_avx1;

				result_final4[2]=_mm256_add_epi16(res_avx2, c6_avx);
				result_final4[3]=_mm256_add_epi16(res_avx3, c22_avx);


				result_final4[4]=_mm256_add_epi16(res_avx4, c7_avx);
				result_final4[5]=_mm256_add_epi16(res_avx5, c23_avx);

				result_final4[6]=res_avx6;
				result_final4[7]=res_avx7;


		   //------------------------AVX interpolation for 5th poly ends--------------

		   //------------------------AVX interpolation for 6th poly external-------------------		
			
			//loop1
				res_avx0 = c_bucket[77];
				res_avx2 = c_bucket[78];
				res_avx4 = c_bucket[79];
				res_avx6 = c_bucket[96];

				c6_avx=c_bucket[99];
				c7_avx=c_bucket[100];
		
				c8_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[101], c6_avx), c7_avx);

				res_avx1 = c_bucket[93];
				res_avx3 = c_bucket[94];
				res_avx5 = c_bucket[95];
				res_avx7 = c_bucket[112];

				c22_avx=c_bucket[115];
				c23_avx=c_bucket[116];

				c21_avx=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[114], res_avx5),res_avx7);

				c24_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[117], c22_avx), c23_avx);

				c20_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[113], res_avx1), res_avx3);

				temp=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[98], res_avx4),res_avx6);
				res_avx5 = _mm256_add_epi16(res_avx5, temp);

				temp = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[97], res_avx0), res_avx2);
				res_avx1 = _mm256_add_epi16(res_avx1, temp);

				c22_avx=_mm256_add_epi16(c22_avx, c8_avx);

				res_avx6 = _mm256_add_epi16(res_avx6, c21_avx);

				res_avx2 = _mm256_add_epi16(res_avx2, c20_avx);

				c7_avx=_mm256_add_epi16(c7_avx, c24_avx);


			//loop4

				c6_avx=_mm256_sub_epi16(_mm256_sub_epi16(c6_avx, res_avx0), res_avx4);
				c22_avx=_mm256_sub_epi16(_mm256_sub_epi16(c22_avx, res_avx1), res_avx5);

				c7_avx=_mm256_sub_epi16(_mm256_sub_epi16(c7_avx, res_avx2), res_avx6);
				c23_avx=_mm256_sub_epi16(_mm256_sub_epi16(c23_avx, res_avx3), res_avx7);

			//loop5
				result_final5[0]=res_avx0;
				result_final5[1]=res_avx1;

				result_final5[2]=_mm256_add_epi16(res_avx2, c6_avx);
				result_final5[3]=_mm256_add_epi16(res_avx3, c22_avx);


				result_final5[4]=_mm256_add_epi16(res_avx4, c7_avx);
				result_final5[5]=_mm256_add_epi16(res_avx5, c23_avx);

				result_final5[6]=res_avx6;
				result_final5[7]=res_avx7;


		   //------------------------AVX interpolation for 6th poly ends--------------

		   //------------------------AVX interpolation for 7th poly external-------------------		
			
			//loop1
				res_avx0 = c_bucket[102];
				res_avx2 = c_bucket[103];
				res_avx4 = c_bucket[104];
				res_avx6 = c_bucket[105];

				c6_avx=c_bucket[108];
				c7_avx=c_bucket[109];
		
				c8_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[110], c6_avx), c7_avx);

				res_avx1 = c_bucket[118];
				res_avx3 = c_bucket[119];
				res_avx5 = c_bucket[120];
				res_avx7 = c_bucket[121];

				c22_avx=c_bucket[124];
				c23_avx=c_bucket[125];

				c21_avx=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[123], res_avx5),res_avx7);

				c24_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[126], c22_avx), c23_avx);

				c20_avx = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[122], res_avx1), res_avx3);

				temp=_mm256_sub_epi16(_mm256_sub_epi16(c_bucket[107], res_avx4),res_avx6);
				res_avx5 = _mm256_add_epi16(res_avx5, temp);

				temp = _mm256_sub_epi16(_mm256_sub_epi16(c_bucket[106], res_avx0), res_avx2);
				res_avx1 = _mm256_add_epi16(res_avx1, temp);

				c22_avx=_mm256_add_epi16(c22_avx, c8_avx);

				res_avx6 = _mm256_add_epi16(res_avx6, c21_avx);

				res_avx2 = _mm256_add_epi16(res_avx2, c20_avx);

				c7_avx=_mm256_add_epi16(c7_avx, c24_avx);


			//loop4

				c6_avx=_mm256_sub_epi16(_mm256_sub_epi16(c6_avx, res_avx0), res_avx4);
				c22_avx=_mm256_sub_epi16(_mm256_sub_epi16(c22_avx, res_avx1), res_avx5);

				c7_avx=_mm256_sub_epi16(_mm256_sub_epi16(c7_avx, res_avx2), res_avx6);
				c23_avx=_mm256_sub_epi16(_mm256_sub_epi16(c23_avx, res_avx3), res_avx7);

			//loop5
				result_final6[0]=res_avx0;
				result_final6[1]=res_avx1;

				result_final6[2]=_mm256_add_epi16(res_avx2, c6_avx);
				result_final6[3]=_mm256_add_epi16(res_avx3, c22_avx);


				result_final6[4]=_mm256_add_epi16(res_avx4, c7_avx);
				result_final6[5]=_mm256_add_epi16(res_avx5, c23_avx);

				result_final6[6]=res_avx6;
				result_final6[7]=res_avx7;


		   //------------------------AVX interpolation for 7th poly ends--------------

		//CLOCK2=cpucycles();
		//CLOCK_INTER=CLOCK_INTER+(CLOCK2-CLOCK1);
		//printf("\nTime for interpolation : %llu\n", CLOCK2-CLOCK1);



}

void toom_cook_4way_avx_n1(__m256i* a_avx,__m256i* b_bucket, __m256i *c_bucket, int f){ 

	int i;

//---------------AVX data-----------------------------

	__m256i r0_avx, r1_avx, r2_avx, r3_avx, r4_avx, r5_avx, r6_avx;
	__m256i aw_avx[7*small_len_avx];

//----------------AVX data----------------------------


// EVALUATION

	//CLOCK1=cpucycles();

	for (i=0; i<small_len_avx; i++){
		r0_avx=a_avx[i];
		r1_avx=a_avx[i + small_len_avx];
		r2_avx=a_avx[i + 2*small_len_avx];
		r3_avx=a_avx[i + 3*small_len_avx];
		r4_avx= _mm256_add_epi16(r0_avx, r2_avx);
		r5_avx= _mm256_add_epi16(r1_avx, r3_avx);
		aw_avx[2*small_len_avx+i]= _mm256_add_epi16(r4_avx, r5_avx);
		aw_avx[3*small_len_avx+i]= _mm256_sub_epi16(r4_avx, r5_avx);
		r4_avx=_mm256_slli_epi16(r0_avx,2);
		r4_avx=_mm256_add_epi16(r4_avx,r2_avx);
		r4_avx=_mm256_slli_epi16(r4_avx,1);
		r5_avx=	_mm256_slli_epi16(r1_avx, 2);
		r5_avx=_mm256_add_epi16(r5_avx, r3_avx);
		aw_avx[4*small_len_avx+i]= _mm256_add_epi16(r4_avx, r5_avx);
		aw_avx[5*small_len_avx+i]= _mm256_sub_epi16(r4_avx, r5_avx);
		r4_avx= _mm256_slli_epi16(r3_avx, 3);
		r6_avx= _mm256_slli_epi16(r2_avx, 2);
		r4_avx= _mm256_add_epi16(r4_avx, r6_avx);
		r6_avx= _mm256_slli_epi16(r1_avx, 1);
		r4_avx= _mm256_add_epi16(r4_avx, r6_avx);
		aw_avx[small_len_avx+i]= _mm256_add_epi16(r4_avx, r0_avx);
		aw_avx[6*small_len_avx+i]= r0_avx; 
		aw_avx[i]= r3_avx;
	}


	//CLOCK2=cpucycles();
	//CLOCK_TC_EVAL=CLOCK_TC_EVAL+(CLOCK2-CLOCK1);

	batch_64coefficient_multiplications_new(aw_avx, b_bucket, c_bucket, f);//New

}

void TC_eval(__m256i* b_avx, __m256i* b_bucket){

	int i;
	__m256i bw_avx[7*small_len_avx];

	__m256i r0_avx, r1_avx, r2_avx, r3_avx, r4_avx, r5_avx, r6_avx;

	for (i=0; i<small_len_avx; i++){
		
		r0_avx=b_avx[i];
		r1_avx=b_avx[i + small_len_avx];
		r2_avx=b_avx[i + 2*small_len_avx];
		r3_avx=b_avx[i + 3*small_len_avx];
		r4_avx= _mm256_add_epi16(r0_avx, r2_avx);
		r5_avx= _mm256_add_epi16(r1_avx, r3_avx);
		bw_avx[2*small_len_avx+i]= _mm256_add_epi16(r4_avx, r5_avx);
		bw_avx[3*small_len_avx+i]= _mm256_sub_epi16(r4_avx, r5_avx);
		r4_avx=_mm256_slli_epi16(r0_avx,2);
		r4_avx=_mm256_add_epi16(r4_avx,r2_avx);
		r4_avx=_mm256_slli_epi16(r4_avx,1);
		r5_avx=	_mm256_slli_epi16(r1_avx, 2);
		r5_avx=_mm256_add_epi16(r5_avx, r3_avx);
		bw_avx[4*small_len_avx+i]= _mm256_add_epi16(r4_avx, r5_avx);
		bw_avx[5*small_len_avx+i]= _mm256_sub_epi16(r4_avx, r5_avx);
		r4_avx= _mm256_slli_epi16(r3_avx, 3);
		r6_avx= _mm256_slli_epi16(r2_avx, 2);
		r4_avx= _mm256_add_epi16(r4_avx, r6_avx);
		r6_avx= _mm256_slli_epi16(r1_avx, 1);
		r4_avx= _mm256_add_epi16(r4_avx, r6_avx);
		bw_avx[small_len_avx+i]= _mm256_add_epi16(r4_avx, r0_avx);
		bw_avx[6*small_len_avx+i]= r0_avx;
		bw_avx[i]= r3_avx;
	}

	KARA_eval(bw_avx, b_bucket);

}


void TC_interpol(__m256i *c_bucket, __m256i* res_avx){

	int i;

	register __m256i r0_avx, r1_avx, r2_avx, r3_avx, r4_avx, r5_avx, r6_avx, temp_avx;

	__m256i w1_avx[2*small_len_avx],w2_avx[2*small_len_avx],w3_avx[2*small_len_avx],w4_avx[2*small_len_avx],w5_avx[2*small_len_avx],w6_avx[2*small_len_avx],w7_avx[2*small_len_avx];

	__m256i res_avx_output[2*AVX_N1];

	//CLOCK1=cpucycles();

	
	transpose_n1(c_bucket);
	transpose_n1(c_bucket+16);

	transpose_n1(c_bucket+2*SCM_SIZE);
	transpose_n1(c_bucket+16+2*SCM_SIZE);

	transpose_n1(c_bucket+4*SCM_SIZE);
	transpose_n1(c_bucket+16+4*SCM_SIZE);

	transpose_n1(c_bucket+6*SCM_SIZE);
	transpose_n1(c_bucket+16+6*SCM_SIZE);
	

	KARA_interpol(c_bucket, w1_avx, w2_avx, w3_avx, w4_avx, w5_avx, w6_avx, w7_avx);

	for (i = 0; i < 2*small_len_avx; i++) {

		r0_avx = w1_avx[i];
		r1_avx = w2_avx[i];
		r2_avx = w3_avx[i];
		r3_avx = w4_avx[i];
		r4_avx = w5_avx[i];
		r5_avx = w6_avx[i];
		r6_avx = w7_avx[i];
		r1_avx = _mm256_add_epi16(r1_avx, r4_avx);
		r5_avx = _mm256_sub_epi16(r5_avx, r4_avx);
		r3_avx = _mm256_sub_epi16(r3_avx, r2_avx);
		r3_avx = _mm256_srli_epi16(r3_avx, 1);
		r4_avx = _mm256_sub_epi16(r4_avx, r0_avx);
		temp_avx = _mm256_slli_epi16(r6_avx, 6);
		r4_avx = _mm256_sub_epi16(r4_avx, temp_avx);
		r4_avx = _mm256_slli_epi16(r4_avx, 1);
		r4_avx = _mm256_add_epi16(r4_avx, r5_avx);
		r2_avx = _mm256_add_epi16(r2_avx, r3_avx);
		temp_avx = _mm256_slli_epi16(r2_avx, 6);
		r1_avx = _mm256_sub_epi16(r1_avx, temp_avx);
		r1_avx = _mm256_sub_epi16(r1_avx, r2_avx);
		r2_avx = _mm256_sub_epi16(r2_avx, r6_avx);
		r2_avx = _mm256_sub_epi16(r2_avx, r0_avx);
		temp_avx = _mm256_mullo_epi16 (r2_avx,int45_avx);
		r1_avx = _mm256_add_epi16(r1_avx, temp_avx);
		temp_avx = _mm256_slli_epi16(r2_avx, 3);
		r4_avx = _mm256_sub_epi16(r4_avx, temp_avx);
		r4_avx = _mm256_mullo_epi16 (r4_avx,inv3_avx);
		r4_avx = _mm256_srli_epi16(r4_avx, 3);
		r5_avx = _mm256_add_epi16(r5_avx, r1_avx);
		temp_avx = _mm256_slli_epi16(r3_avx, 4);
		r1_avx= _mm256_add_epi16(r1_avx, temp_avx);
		r1_avx = _mm256_mullo_epi16 (r1_avx, inv9_avx);
		r1_avx= _mm256_srli_epi16(r1_avx, 1); 	
		r3_avx= _mm256_add_epi16(r1_avx, r3_avx);
		r3_avx= _mm256_sub_epi16(int0_avx, r3_avx);
		temp_avx= _mm256_mullo_epi16 (r1_avx,int30_avx);
		temp_avx= _mm256_sub_epi16(temp_avx, r5_avx);
		temp_avx= _mm256_mullo_epi16 (temp_avx ,inv15_avx);
		r5_avx= _mm256_srli_epi16(temp_avx, 2);
		r2_avx = _mm256_sub_epi16(r2_avx, r4_avx);
		r1_avx = _mm256_sub_epi16(r1_avx, r5_avx);

		if(i<small_len_avx){
			res_avx_output[0*small_len_avx+i]=r6_avx;
			res_avx_output[1*small_len_avx+i]=r5_avx;
			res_avx_output[2*small_len_avx+i]=r4_avx;
			res_avx_output[3*small_len_avx+i]=r3_avx;
			res_avx_output[4*small_len_avx+i]=r2_avx;
			res_avx_output[5*small_len_avx+i]=r1_avx;
			res_avx_output[6*small_len_avx+i]=r0_avx;
		}
		else{
			res_avx_output[0*small_len_avx+i]=_mm256_add_epi16 (res_avx_output[0*small_len_avx+i], r6_avx);
			res_avx_output[1*small_len_avx+i]=_mm256_add_epi16 (res_avx_output[1*small_len_avx+i], r5_avx);
			res_avx_output[2*small_len_avx+i]=_mm256_add_epi16 (res_avx_output[2*small_len_avx+i], r4_avx);
			res_avx_output[3*small_len_avx+i]=_mm256_add_epi16 (res_avx_output[3*small_len_avx+i], r3_avx);
			res_avx_output[4*small_len_avx+i]=_mm256_add_epi16 (res_avx_output[4*small_len_avx+i], r2_avx);
			res_avx_output[5*small_len_avx+i]=_mm256_add_epi16 (res_avx_output[5*small_len_avx+i], r1_avx);
			res_avx_output[6*small_len_avx+i]=r0_avx;
		}
	}

	//CLOCK2=cpucycles();
	//CLOCK_TC_INTER=CLOCK_TC_INTER+(CLOCK2-CLOCK1);

	// Reduction by X^256 + 1
	for(i=0; i<16; i++)
		res_avx[i] = _mm256_sub_epi16(res_avx_output[i], res_avx_output[i+16]);

}

#endif