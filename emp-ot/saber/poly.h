#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "SABER_params.h"
#include <stdio.h>
#include "api.h"
#include "pack_unpack.h"
#include "cbd.h"
#include "fips202.h"
#include "poly_mul.h"
#include <immintrin.h>

typedef struct
{
  uint16_t coeffs[SABER_N];
} poly;

typedef struct{
  poly vec[SABER_K];
} polyvec;

void poly_getnoise(uint16_t *r,const unsigned char *seed, unsigned char nonce);

void poly_getnoise4x(uint16_t *r0, uint16_t *r1, uint16_t *r2, const unsigned char *seed, unsigned char nonce0, unsigned char nonce1, unsigned char nonce2, unsigned char nonce3);

uint16_t Bits(uint16_t x, int i, int j);


void GenMatrix(polyvec *a, const unsigned char seed[SABER_SEEDBYTES]);
// void GenSecret(uint16_t s[SABER_K][SABER_N], const unsigned char seed[SABER_NOISE_SEEDBYTES]);
void GenSecret(uint16_t **s, const unsigned char seed[SABER_NOISE_SEEDBYTES]);

uint16_t Bits(uint16_t x, int i, int j)
{
    return (x >> (i - j)) & ((1 << j) - 1);
}


// from AVX_Implementation_KEM/SABER_indcpa.c

void POL2MSG(uint16_t *message_dec_unpacked, unsigned char *message_dec);

/*--------------------------------------------------------------------------------------
	This routine loads the constant values for Toom-Cook multiplication 
----------------------------------------------------------------------------------------*/
void load_values(){

	int64_t i;

	int64_t inv3=43691;
	int64_t inv9=36409;
	int64_t inv15=61167;

	int64_t int45=45;
	int64_t int30=30;
	int64_t int0=0;


	int16_t inv3_avx_load[16],inv9_avx_load[16],inv15_avx_load[16],int45_avx_load[16],int30_avx_load[16],int0_avx_load[16];

	for(i=0;i<16;i++){
		inv3_avx_load[i]=inv3;
		inv9_avx_load[i]=inv9;
		inv15_avx_load[i]=inv15;
		int45_avx_load[i]=int45;
		int30_avx_load[i]=int30;
		int0_avx_load[i]=int0;
	}

	inv3_avx = _mm256_loadu_si256 ((__m256i const *) (&inv3_avx_load));
	inv9_avx = _mm256_loadu_si256 ((__m256i const *) (&inv9_avx_load));
	inv15_avx = _mm256_loadu_si256 ((__m256i const *) (&inv15_avx_load));
	int45_avx = _mm256_loadu_si256 ((__m256i const *) (&int45_avx_load));
	int30_avx = _mm256_loadu_si256 ((__m256i const *) (&int30_avx_load));
	int0_avx = _mm256_loadu_si256 ((__m256i const *) (&int0_avx_load));
	mask = _mm256_loadu_si256 ((__m256i const *)mask_ar);	
}



/*-----------------------------------------------------------------------------------
	This routine generates a=[Matrix K x K] of 256-coefficient polynomials 
-------------------------------------------------------------------------------------*/



void BS2POLq(const unsigned char *bytes, uint16_t data[SABER_N]){
	
	uint32_t j;
	uint32_t offset_data=0,offset_byte=0;	
	
	offset_byte=0;

		for(j=0;j<SABER_N/8;j++){
			offset_byte=13*j;
			offset_data=8*j;
			data[offset_data + 0]= ( bytes[ offset_byte + 0 ] & (0xff)) | ((bytes[offset_byte + 1] & 0x1f)<<8);
			data[offset_data + 1]= ( bytes[ offset_byte + 1 ]>>5 & (0x07)) | ((bytes[offset_byte + 2] & 0xff)<<3) | ((bytes[offset_byte + 3] & 0x03)<<11);
			data[offset_data + 2]= ( bytes[ offset_byte + 3 ]>>2 & (0x3f)) | ((bytes[offset_byte + 4] & 0x7f)<<6);
			data[offset_data + 3]= ( bytes[ offset_byte + 4 ]>>7 & (0x01)) | ((bytes[offset_byte + 5] & 0xff)<<1) | ((bytes[offset_byte + 6] & 0x0f)<<9);
			data[offset_data + 4]= ( bytes[ offset_byte + 6 ]>>4 & (0x0f)) | ((bytes[offset_byte + 7] & 0xff)<<4) | ((bytes[offset_byte + 8] & 0x01)<<12);
			data[offset_data + 5]= ( bytes[ offset_byte + 8]>>1 & (0x7f)) | ((bytes[offset_byte + 9] & 0x3f)<<7);
			data[offset_data + 6]= ( bytes[ offset_byte + 9]>>6 & (0x03)) | ((bytes[offset_byte + 10] & 0xff)<<2) | ((bytes[offset_byte + 11] & 0x07)<<10);
			data[offset_data + 7]= ( bytes[ offset_byte + 11]>>3 & (0x1f)) | ((bytes[offset_byte + 12] & 0xff)<<5);
		}


}


void GenMatrix(polyvec *a, const unsigned char seed[SABER_SEEDBYTES]) 
{
  unsigned int one_vector=13*SABER_N/8;
  unsigned int byte_bank_length=SABER_K*SABER_K*one_vector;
  unsigned char buf[byte_bank_length];

  uint16_t temp_ar[SABER_N];

  int i,j,k;
  uint16_t mod = (SABER_Q-1);

  shake128(buf,byte_bank_length,seed,SABER_SEEDBYTES);
  
  for(i=0;i<SABER_K;i++)
  {
    for(j=0;j<SABER_K;j++)
    {
	BS2POLq(buf+(i*SABER_K+j)*one_vector,temp_ar);
	for(k=0;k<SABER_N;k++){
		a[i].vec[j].coeffs[k] = (temp_ar[k])& mod ;
	}
    }
  }
}

void GenSecret(uint16_t **r,const unsigned char seed[SABER_NOISE_SEEDBYTES]) {


		uint32_t i;

		int32_t buf_size= SABER_MU*SABER_N*SABER_K/8;

		uint8_t buf[buf_size];

		shake128(buf, buf_size, seed,SABER_NOISE_SEEDBYTES);

		for(i=0;i<SABER_K;i++)
		{
			cbd(r[i],buf+i*SABER_MU*SABER_N/8);
		}
}

//********************************matrix-vector mul routines*****************************************************
void matrix_vector_mul(__m256i a1_avx_combined[NUM_POLY][NUM_POLY][AVX_N1], __m256i b_bucket[NUM_POLY][SCHB_N*4], __m256i res_avx[NUM_POLY][AVX_N1], int isTranspose);
void vector_vector_mul(__m256i a_avx[NUM_POLY][AVX_N1], __m256i b_bucket[NUM_POLY][SCHB_N*4], __m256i res_avx[AVX_N1]);

//********************************matrix-vector mul routines*****************************************************

void POL2MSG(uint16_t *message_dec_unpacked, unsigned char *message_dec){

	int32_t i,j;

	for(j=0; j<SABER_KEYBYTES; j++)
	{
		message_dec[j] = 0;
		for(i=0; i<8; i++)
		message_dec[j] = message_dec[j] | (message_dec_unpacked[j*8 + i] <<i);
	} 

}

void matrix_vector_mul(__m256i a1_avx_combined[NUM_POLY][NUM_POLY][AVX_N1], __m256i b_bucket[NUM_POLY][SCHB_N*4], __m256i res_avx[NUM_POLY][AVX_N1], int isTranspose){



	int64_t i,j;

	__m256i c_bucket[2*SCM_SIZE*4]; //Holds results for 9 Karatsuba at a time

	for(i=0;i<NUM_POLY;i++){
		for(j=0;j<NUM_POLY;j++){

			if(isTranspose==0){
				toom_cook_4way_avx_n1(a1_avx_combined[i][j], b_bucket[j], c_bucket, j);
			}
			else{
				toom_cook_4way_avx_n1(a1_avx_combined[j][i], b_bucket[j], c_bucket, j);
			}
		}

		TC_interpol(c_bucket, res_avx[i]);
	}

}

void vector_vector_mul(__m256i a_avx[NUM_POLY][AVX_N1], __m256i b_bucket[NUM_POLY][SCHB_N*4], __m256i res_avx[AVX_N1]){

	int64_t i;

	__m256i c_bucket[2*SCM_SIZE*4]; //Holds results for 9 Karatsuba at a time

	for(i=0;i<NUM_POLY;i++){
		toom_cook_4way_avx_n1(a_avx[i], b_bucket[i], c_bucket, i);		
	}
	TC_interpol(c_bucket, res_avx);
}


#endif