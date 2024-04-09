#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "SABER_params.h"

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

#endif