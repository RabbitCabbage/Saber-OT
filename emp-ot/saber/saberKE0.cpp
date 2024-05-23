#include "poly.h"
#include "poly_mul.h"
#include "api.h"
#include "pack_unpack.h"
#include "cbd.h"
#include "fips202.h"
#include "SABER_params.h"
#include "rng.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <emp-tool/emp-tool.h>
#include <immintrin.h>

int main()
{
    // temporary variables as buffer
    uint16_t temp[SABER_K][SABER_N];

    // simulate the sender side, without prime
    unsigned char seed[SABER_SEEDBYTES];
    unsigned char noiseseed[SABER_COINBYTES];
    uint32_t i, j, k;
    polyvec a[SABER_K];
    uint16_t **sk = new uint16_t *[SABER_K];
    for (i = 0; i < SABER_K; i++)
    {
        sk[i] = new uint16_t[SABER_N];
    }
    

    //--------------AVX declaration------------------

    __m256i sk_avx[SABER_K][SABER_N / 16];
    __m256i mod, mod_p;
    __m256i res_avx[SABER_K][SABER_N / 16];
    __m256i a_avx[SABER_K][SABER_K][SABER_N / 16];
    //__m256i acc[2*SABER_N/16];


    mask_ar[0] = ~(0UL);
    mask_ar[1] = ~(0UL);
    mask_ar[2] = ~(0UL);
    mask_ar[3] = ~(0UL);
    mask_load = _mm256_loadu_si256((__m256i const *)mask_ar);

    mod = _mm256_set1_epi16(SABER_Q - 1);
    mod_p = _mm256_set1_epi16(SABER_P - 1);

    floor_round = _mm256_set1_epi16(4);

    H1_avx = _mm256_set1_epi16(h1);

    __m256i b_bucket[NUM_POLY][SCHB_N * 4];

    //--------------AVX declaration ends------------------

    load_values();
    randombytes(seed, SABER_SEEDBYTES);

    shake128(seed, SABER_SEEDBYTES, seed, SABER_SEEDBYTES); // for not revealing system RNG state
    randombytes(noiseseed, SABER_COINBYTES);
    GenMatrix(a, seed);
    GenSecret(sk, noiseseed);

    // ----------- Load sk into avx vectors ----------
    for (i = 0; i < SABER_K; i++)
    {
        for (j = 0; j < SABER_N / 16; j++)
        {
            sk_avx[i][j] = _mm256_loadu_si256((__m256i const *)(&sk[i][j * 16]));
        }
    }

    // ----------- Load sk into avx vectors ----------
    for (i = 0; i < SABER_K; i++)
    {
        for (j = 0; j < SABER_K; j++)
        {
            for (k = 0; k < SABER_N / 16; k++)
            {
                a_avx[i][j][k] = _mm256_loadu_si256((__m256i const *)(&a[i].vec[j].coeffs[k * 16]));
            }
        }
    }
    //-----------------matrix-vector multiplication and rounding

    for (j = 0; j < NUM_POLY; j++)
    {
        TC_eval(sk_avx[j], b_bucket[j]);
    }
    matrix_vector_mul(a_avx, b_bucket, res_avx, 0); // Matrix-vector multiplication; Matrix in normal order

    // Now truncation
    for (i = 0; i < SABER_K; i++)
    { // shift right EQ-EP bits
        for (j = 0; j < SABER_N / 16; j++)
        {
            res_avx[i][j] = _mm256_add_epi16(res_avx[i][j], H1_avx);
            res_avx[i][j] = _mm256_srli_epi16(res_avx[i][j], (SABER_EQ - SABER_EP));
            res_avx[i][j] = _mm256_and_si256(res_avx[i][j], mod);
        }
    }

    // res_avx is round A dot s
    //-----this result should be put in b for later use in server.
    for (i = 0; i < SABER_K; i++)
    { // first store in 16 bit arrays
        for (j = 0; j < SABER_N / 16; j++)
        {
            _mm256_maskstore_epi32((int *)(temp[i] + j * 16), mask_load, res_avx[i][j]);
            // temp is for temponary use of the result b.
        }
    }

    // use POLVEC2BS to pack our vector b, pack b into pack_b
    unsigned char *pack_b = new unsigned char[SABER_CIPHERTEXTBYTES];
    POLVEC2BS(pack_b, temp, SABER_P);
    // then pack_b will be sent to the receiver

    // simulate the receiver side, with prime
    unsigned char seed_prime[SABER_SEEDBYTES];
    unsigned char noiseseed_prime[SABER_COINBYTES];
    polyvec a_prime[SABER_K];
    uint16_t **sk_prime = new uint16_t *[SABER_K];
    for (i = 0; i < SABER_K; i++)
    {
        sk_prime[i] = new uint16_t[SABER_N];
    }

    //--------------AVX declaration------------------
    __m256i sk_avx_prime[SABER_K][SABER_N / 16];
    __m256i res_avx_prime[SABER_K][SABER_N / 16];
    __m256i b_bucket_prime[NUM_POLY][SCHB_N * 4];

    //--------------AVX declaration ends------------------

    load_values();
    randombytes(seed_prime, SABER_SEEDBYTES);

    shake128(seed_prime, SABER_SEEDBYTES, seed_prime, SABER_SEEDBYTES); // for not revealing system RNG state
    randombytes(noiseseed_prime, SABER_COINBYTES);
    GenMatrix(a_prime, seed_prime);
    GenSecret(sk_prime, noiseseed_prime);

    // ----------- Load sk into avx vectors ----------
    for (i = 0; i < SABER_K; i++)
    {
        for (j = 0; j < SABER_N / 16; j++)
        {
            sk_avx_prime[i][j] = _mm256_loadu_si256((__m256i const *)(&sk_prime[i][j * 16]));
        }
    }

    // ----------- Load sk into avx vectors ----------
    for (i = 0; i < SABER_K; i++)
    {
        for (j = 0; j < SABER_K; j++)
        {
            for (k = 0; k < SABER_N / 16; k++)
            {
                a_avx[i][j][k] = _mm256_loadu_si256((__m256i const *)(&a_prime[i].vec[j].coeffs[k * 16]));
            }
        }
    }

    //-----------------matrix-vector multiplication and rounding
    for (j = 0; j < NUM_POLY; j++)
    {
        TC_eval(sk_avx_prime[j], b_bucket_prime[j]);
    }
    matrix_vector_mul(a_avx, b_bucket_prime, res_avx_prime, 1); // Matrix-vector multiplication; Matrix in normal order

    // Now truncation
    for (i = 0; i < SABER_K; i++)
    { // shift right EQ-EP bits
        for (j = 0; j < SABER_N / 16; j++)
        {
            res_avx_prime[i][j] = _mm256_add_epi16(res_avx_prime[i][j], H1_avx);
            res_avx_prime[i][j] = _mm256_srli_epi16(res_avx_prime[i][j], (SABER_EQ - SABER_EP));
            res_avx_prime[i][j] = _mm256_and_si256(res_avx_prime[i][j], mod);
        }
    }

    // res_avx is round A dot s
    //-----this result should be put in b for later use in server.
    for (i = 0; i < SABER_K; i++)
    { // first store in 16 bit arrays
        for (j = 0; j < SABER_N / 16; j++)
        {
            _mm256_maskstore_epi32((int *)(temp[i] + j * 16), mask_load, res_avx_prime[i][j]);
            // temp is for temponary use of the result b.
        }
    }

    // use POLVEC2BS to pack our vector b, pack b into pack_b
    unsigned char *pack_bp = new unsigned char[SABER_CIPHERTEXTBYTES];
    POLVEC2BS(pack_bp, temp, SABER_P);
    // then pack_b will be sent to the receiver

    // simulate the sender side, without prime
    // the sender will compute v and c for the receiver
        //--------------AVX declaration------------------
    __m256i v_avx[SABER_N / 16];
    __m256i bp_avx[SABER_K][SABER_N / 16];
    //--------------AVX declaration ends------------------

    // define bp
    uint16_t bp[SABER_K][SABER_N];
    BS2POLVEC(pack_bp, bp, SABER_P); // will this be wrong after compression? TODO
    for (i = 0; i < SABER_K; i++)
    {
        for (j = 0; j < SABER_N / 16; j++)
        {
            bp_avx[i][j] = _mm256_loadu_si256((__m256i const *)(&bp[i][j * 16]));
        }
    }

    vector_vector_mul(bp_avx, b_bucket, v_avx);
    // Computation of v'+h1
    for (i = 0; i < SABER_N / 16; i++)
    { // adding h1
        v_avx[i] = _mm256_add_epi16(v_avx[i], H1_avx);
    }

   
   // copy a duplicate of v_avx
    __m256i v_avx_copy[SABER_N / 16];
    for (i = 0; i < SABER_N / 16; i++)
    {
        // mod p, so that this is a copy of v+h1 mod p
        v_avx_copy[i]= _mm256_and_si256(v_avx[i], mod_p);
    }


    // SHIFTRIGHT(v+h1-m mod p, EP-ET)
    for (k = 0; k < SABER_N / 16; k++)
    {
        v_avx[k] = _mm256_and_si256(v_avx[k], mod_p);
        v_avx[k] = _mm256_srli_epi16(v_avx[k], (SABER_EP - SABER_ET));
    }

    // Unpack avx
    for (j = 0; j < SABER_N / 16; j++)
    {
        _mm256_maskstore_epi32((int *)(temp[0] + j * 16), mask_load, v_avx[j]);
    }
    unsigned char msk_c[SABER_SCALEBYTES_KEM];

#if Saber_type == 1
    SABER_pack_3bit(msk_c, temp[0]);
#elif Saber_type == 2
    SABER_pack_4bit(msk_c, temp[0]);
#elif Saber_type == 3
    SABER_pack_6bit(msk_c, temp[0]);
#endif
    // so msk_c is the c for the receiver
    
    // simulate the receiver side, with prime
    // the receiver will compute v_prime and use c to compute the shared secret
        //--------------AVX declaration------------------
    __m256i v_prime_avx[SABER_N / 16];
    //--------------AVX declaration ends------------------
    
    // clear v_prime_avx
    for (i = 0; i < SABER_N / 16; i++)
    {
        v_prime_avx[i] = _mm256_xor_si256(v_prime_avx[i], v_prime_avx[i]);
    }

    // InnerProduct(b_prime, s_prime, mod_p)
    vector_vector_mul(bp_avx, b_bucket_prime, v_prime_avx);

    // Unpack avx to temp
    uint16_t message_dec_unpacked[SABER_KEYBYTES * 8];
    for(i=0; i<SABER_N/16; i++){
		_mm256_maskstore_epi32 ((int *)(message_dec_unpacked+i*16), mask_load, v_avx[i]);
	}

uint16_t c[SABER_N];
#if Saber_type == 1
    SABER_pack_3bit(msk_c, c);
#elif Saber_type == 2
    SABER_pack_4bit(msk_c, c);
#elif Saber_type == 3
    SABER_pack_6bit(msk_c, c);
#endif

    // addition of h2
    for(i=0;i<SABER_N;i++){
		message_dec_unpacked[i]= ( ( message_dec_unpacked[i] + h2 - (c[i]<<(SABER_EP-SABER_ET)) ) & (SABER_P-1) ) >> (SABER_EP-1);
	} // the message_dec_unpacked is the shared bits
    
    // the sender's result is from v, round to 1 bit, from the v_avx_copy
    uint16_t message_dec_unpacked_sender[SABER_KEYBYTES * 8];
    // Unpack avx to temp
    for(i=0; i<SABER_N/16; i++){
        _mm256_maskstore_epi32 ((int *)(message_dec_unpacked_sender+i*16), mask_load, v_avx_copy[i]);
    }

    // compute the highest bit of the unpacked message_dec_unpacked_sender
    for(i=0;i<SABER_N;i++){
        message_dec_unpacked_sender[i]= ( ( message_dec_unpacked_sender[i]) & (SABER_P-1) ) >> (SABER_EP-1);
    }

    // compare the two results
    for(i=0;i<SABER_N;i++){
        if(message_dec_unpacked[i]!=message_dec_unpacked_sender[i]){
            std::cout<<"Error: "<<i<<std::endl;
        }
    }
    

    for (i = 0; i < SABER_K; i++)
    {
        delete[] sk[i];
    }
    delete[] sk;
    for (i = 0; i < SABER_K; i++)
    {
        delete[] sk_prime[i];
    }
    delete[] sk_prime;
    delete[] pack_b;
    delete[] pack_bp;
    return 0;
}