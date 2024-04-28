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

uint16_t Bits(uint16_t x, int i, int j);

// void MatrixVectorMul(const uint16_t a[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose);
void MatrixVectorMul(uint16_t A[SABER_L][SABER_L][SABER_N], uint16_t **s, uint16_t **res, int16_t transpose);

// void RoundingMul(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose);
void RoundingMul(uint16_t A[SABER_L][SABER_L][SABER_N], uint16_t **s, uint16_t **res, int16_t transpose);

// void InnerProd_plush1(const uint16_t b[SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_N]);
void InnerProd_plush1(uint16_t **b, uint16_t **s, uint16_t *res);


// void InnerProd(const uint16_t b[SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_N]);
void InnerProd(uint16_t **b, uint16_t **s, uint16_t *res);


void GenMatrix(uint16_t a[SABER_L][SABER_L][SABER_N], const uint8_t seed[SABER_SEEDBYTES]);
// void GenSecret(uint16_t s[SABER_L][SABER_N], const uint8_t seed[SABER_NOISE_SEEDBYTES]);
void GenSecret(uint16_t **s, const uint8_t seed[SABER_NOISE_SEEDBYTES]);

uint16_t Bits(uint16_t x, int i, int j)
{
    return (x >> (i - j)) & ((1 << j) - 1);
}


// void MatrixVectorMul(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose)
void MatrixVectorMul( uint16_t A[SABER_L][SABER_L][SABER_N],  uint16_t **s, uint16_t **res, int16_t transpose)
{
	int i, j;
	for (i = 0; i < SABER_L; i++)
	{
		for (j = 0; j < SABER_L; j++)
		{
			if (transpose == 1)
			{
				poly_mul_acc(A[j][i], s[j], res[i]);
			}
			else
			{
				poly_mul_acc(A[i][j], s[j], res[i]);
			}	
		}
	}
}

// void RoundingMul(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose) {
void RoundingMul( uint16_t A[SABER_L][SABER_L][SABER_N],  uint16_t **s, uint16_t **res, int16_t transpose) {
	MatrixVectorMul(A, s, res, transpose);

	for (int i = 0; i < SABER_L; i++) {
		for (int j = 0; j < SABER_N; j++) {
			res[i][j] = Bits(res[i][j] + h1, SABER_EQ, SABER_EP);
		}
	}
	
	return;
}

// void InnerProd(const uint16_t ** b[SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_N])
void InnerProd( uint16_t ** b,  uint16_t **s, uint16_t *res)
{
	int j;
	for (j = 0; j < SABER_L; j++)
	{
		poly_mul_acc(b[j], s[j], res);
	}
}

// void InnerProd_plush1(const uint16_t b[SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_N]) 
void InnerProd_plush1( uint16_t **b,  uint16_t **s, uint16_t *res) {
	InnerProd(b, s, res);
	// res += Bits(h1,SABER_EP,SABER_EP);
	uint16_t h1_bits = Bits(h1, SABER_EP, SABER_EP);
	for (int i = 0; i < SABER_N; i++) {
		res[i] = res[i] + h1_bits;
	}
	return;
}

void GenMatrix(uint16_t A[SABER_L][SABER_L][SABER_N], const uint8_t seed[SABER_SEEDBYTES])
{
	uint8_t buf[SABER_L * SABER_POLYVECBYTES];
	int i;

	shake128(buf, sizeof(buf), seed, SABER_SEEDBYTES);

	for (i = 0; i < SABER_L; i++)
	{
		BS2POLVECq(buf + i * SABER_POLYVECBYTES, A[i]);
	}
}

// void GenSecret(uint16_t s[SABER_L][SABER_N], const uint8_t seed[SABER_NOISE_SEEDBYTES])
void GenSecret(uint16_t **s, const uint8_t seed[SABER_NOISE_SEEDBYTES])
{
	uint8_t buf[SABER_L * SABER_POLYCOINBYTES];
	size_t i;

	shake128(buf, sizeof(buf), seed, SABER_NOISE_SEEDBYTES);

	for (i = 0; i < SABER_L; i++)
	{
		cbd(s[i], buf + i * SABER_POLYCOINBYTES);
	}
}

#endif