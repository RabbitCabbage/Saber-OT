// Naor Pinkas Saber OT with one secrets
// Naor-Pinkas Saber OT with two secrets
#ifndef NPSaber_v2_H__
#define NPSaber_v2_H__
#include "emp-ot/saber/api.h"
#include "emp-ot/saber/cbd.h"
#include "emp-ot/saber/fips202.h"
#include "emp-ot/saber/pack_unpack.h"
#include "emp-ot/saber/poly.h"
#include "emp-ot/saber/poly_mul.h"
#include "emp-ot/saber/rng.h"
#include "emp-ot/saber/SABER_params.h"
#include <emp-tool/emp-tool.h>
#include <immintrin.h>
#include <iostream>
#include <stdint.h>
#include <string.h>


namespace emp {
template<typename IO>
class NPSaber2: public OT<IO> {
    public:
    IO* io;
    unsigned char seed_A[SABER_SEEDBYTES];
    uint16_t r[SABER_K * SABER_N];

    NPSaber2(IO* io, uint8_t* _seed_A = nullptr, uint16_t* _r = nullptr) {
        this->io = io;
        if(_seed_A != nullptr) {
            memcpy(seed_A, _seed_A, SABER_SEEDBYTES);
        }
        if(_r != nullptr) {
            memcpy(r, _r, SABER_K * SABER_N * sizeof(uint16_t));
        }
    }

    ~NPSaber2() {

    }

    void send(const block* data0, const block* data1,int64_t length) override {
        uint32_t i,j,k;
        polyvec a[SABER_K];
        uint16_t skpv1[SABER_K][SABER_N];
        uint16_t temp[SABER_K][SABER_N];
        uint16_t bp[SABER_K][SABER_N];

        //--------------AVX declaration------------------
        
        __m256i sk_avx[SABER_K][SABER_N/16];
        __m256i mod, mod_p;
        __m256i res_avx[SABER_K][SABER_N/16];
        __m256i v_avx[SABER_N/16];
        __m256i a_avx[SABER_K][SABER_K][SABER_N/16];
        //__m256i acc[2*SABER_N/16];

        __m256i bp_avx[SABER_K][SABER_N/16];

        mask_ar[0]=~(0UL);mask_ar[1]=~(0UL);mask_ar[2]=~(0UL);mask_ar[3]=~(0UL);
        mask_load = _mm256_loadu_si256 ((__m256i const *)mask_ar);

        mod=_mm256_set1_epi16(SABER_Q-1);
        mod_p=_mm256_set1_epi16(SABER_P-1);

        

        floor_round=_mm256_set1_epi16(4);

        H1_avx=_mm256_set1_epi16(h1);
    
        __m256i b_bucket[NUM_POLY][SCHB_N*4];

        //--------------AVX declaration ends------------------

        load_values();
        GenMatrix(a, seed_A);
        GenSecret(skpv1,noiseseed);

        // ----------- Load skpv1 into avx vectors ---------- 
        for(i=0;i<SABER_K;i++){ 
            for(j=0; j<SABER_N/16; j++){
                sk_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&skpv1[i][j*16]));
            }
        }

        // ----------- Load skpv1 into avx vectors ---------- 
        for(i=0;i<SABER_K;i++){ 
            for(j=0;j<SABER_K;j++){
                for(k=0;k<SABER_N/16;k++){
                    a_avx[i][j][k]=_mm256_loadu_si256 ((__m256i const *) (&a[i].vec[j].coeffs[k*16]));
                }
            }
        }
        //-----------------matrix-vector multiplication and rounding

        for(j=0;j<NUM_POLY;j++){
            TC_eval(sk_avx[j], b_bucket[j]);
        }
        matrix_vector_mul(a_avx, b_bucket, res_avx, 0);// Matrix-vector multiplication; Matrix in normal order
        
        // Now truncation
        for(i=0;i<SABER_K;i++){ //shift right EQ-EP bits
            for(j=0;j<SABER_N/16;j++){
                res_avx[i][j]=_mm256_add_epi16 (res_avx[i][j], H1_avx);
                res_avx[i][j]=_mm256_srli_epi16 (res_avx[i][j], (SABER_EQ-SABER_EP) );
                res_avx[i][j]=_mm256_and_si256 (res_avx[i][j], mod);			
            }
        }

        // res_avx is round A dot s
        //-----this result should be put in b for later use in server.
        for(i=0;i<SABER_K;i++){ // first store in 16 bit arrays
            for(j=0;j<SABER_N/16;j++){
                _mm256_maskstore_epi32 ((int *)(temp[i]+j*16), mask_load, res_avx[i][j]);
                //temp is for temponary use of the result b.
            }
        }

        // use POLVEC2BS to pack our vector b, pack b into ciphertext
        unsigned char *ciphertext;
        POLVEC2BS(ciphertext, temp, SABER_P);

        // receive pack_bp from the receiver
        io->
        // unpack vector bp using BS2POLVEC
        auto pack_bp;
        BS2POLVEC(bp, pack_bp, SABER_P);
        for(i=0;i<SABER_K;i++){
            for(j=0; j<SABER_N/16; j++){
                bp_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&pack_bp[i][j*16]));
            }
        }

        vector_vector_mul(bp_avx, b_bucket, v_avx);
        // Computation of v'+h1
        for(i=0;i<SABER_N/16;i++){//adding h1
            v_avx[i]=_mm256_add_epi16(v_avx[i], H1_avx);
        }

        // SHIFTRIGHT(v'+h1-m mod p, EP-ET)
        for(k=0;k<SABER_N/16;k++)
        {
            // v_avx[k]=_mm256_sub_epi16(v_avx[k], message_avx[k]); // KE has no message.
            v_avx[k]=_mm256_and_si256(v_avx[k], mod_p);
            v_avx[k]=_mm256_srli_epi16 (v_avx[k], (SABER_EP-SABER_ET) );
        }

        // Unpack avx
        for(j=0;j<SABER_N/16;j++)
        {
                _mm256_maskstore_epi32 ((int *) (temp[0]+j*16), mask_load, v_avx[j]);
        }
        
        #if Saber_type == 1
            SABER_pack_3bit(msk_c, temp[0]);
        #elif Saber_type == 2
            SABER_pack_4bit(msk_c, temp[0]);
        #elif Saber_type == 3
            SABER_pack_6bit(msk_c, temp[0]);
        #endif

        for(j=0;j<SABER_SCALEBYTES_KEM;j++){
            ciphertext[SABER_CIPHERTEXTBYTES + j] = msk_c[j];
        }
	}

    void recv(block* data, const bool* x, int64_t length) override {
        
    }
};

}



#endif