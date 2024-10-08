#include "poly.h"
#include "poly_mul.h"
#include "api.h"
#include "pack_unpack.h"
#include "cbd.h"
#include "fips202.h"
#include "SABER_params.h"
#include "rng.h"
#include "SABER_indcpa.h"
#include <emp-tool/emp-tool.h>
#include <iostream>


int main(){
    int length = 1;
    bool b[length] = {1};
    emp::block data0[length];
    emp::block data1[length];
    randombytes(reinterpret_cast<unsigned char *>(data0), length * sizeof(emp::block));
    randombytes(reinterpret_cast<unsigned char *>(data1), length * sizeof(emp::block));
    emp::block data[length];

    // generate receiver message
    uint8_t *r_data = new uint8_t[length * 2 * SABER_INDCPA_PUBLICKEYBYTES];

    uint8_t ***r = new uint8_t**[length];
        for (int64_t i = 0; i < length; i++) {
            r[i] = new uint8_t*[2];
            r[i][0] = r_data + (i * 2 * SABER_INDCPA_PUBLICKEYBYTES);
            r[i][1] = r_data + (i * 2 * SABER_INDCPA_PUBLICKEYBYTES) + SABER_INDCPA_PUBLICKEYBYTES;
        }
        uint8_t ** sk = new uint8_t*[length];
        for (int64_t i = 0; i < length; i++) {
            sk[i] = new uint8_t[SABER_INDCPA_SECRETKEYBYTES];
        }
        uint8_t seed_A[SABER_SEEDBYTES];
        uint8_t hash_r_1xb[SABER_INDCPA_PUBLICKEYBYTES];
        uint16_t A[SABER_L][SABER_L][SABER_N];
        for (int i = 0; i < length; i++) {
            // for b, get a pair of sk and pk
            uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES];
            indcpa_kem_keypair(pk, sk[i]);

            // print pk as hashed for debug
            std::cout << "pk: ";
            for (int j = 0; j < SABER_POLYVECCOMPRESSEDBYTES; j++) {
                std::cout << (int)pk[j];
            }
            std::cout << std::endl;
            std::cout << "pk: ";
            for (int j = SABER_POLYVECCOMPRESSEDBYTES; j < SABER_INDCPA_PUBLICKEYBYTES; j++) {
                std::cout << (int)pk[j];
            }
            std::cout << std::endl;

            // for not b, sample random pk without sk
            randombytes(seed_A, SABER_SEEDBYTES);
            shake128(seed_A, SABER_SEEDBYTES, seed_A, SABER_SEEDBYTES);
            GenMatrix(A, seed_A);
            POLVECp2BS(r[i][!b[i]], A[0]);
            memcpy(r[i][!b[i]] + SABER_POLYVECCOMPRESSEDBYTES, seed_A, sizeof(seed_A));
            // compute H(r_{not b})
            emp::Hash hash;
            unsigned char hash_seed[SABER_SEEDBYTES];
            hash.put(r[i][!b[i]], SABER_POLYVECCOMPRESSEDBYTES);
            hash.digest(hash_seed);
            GenMatrix(A, hash_seed);
            POLVECp2BS(hash_r_1xb, A[0]);
            memcpy(hash_r_1xb + SABER_POLYVECCOMPRESSEDBYTES, seed_A, sizeof(seed_A));

            // compute r_b = pk-H(r_{not b})
            for (int j = 0; j < SABER_INDCPA_PUBLICKEYBYTES; j++) {
                r[i][b[i]][j] =  (pk[j] - hash_r_1xb[j] + 256) % 256;
            }
        }


        // generate sender message
        uint8_t *m_data = new uint8_t[length * 2 * SABER_BYTES_CCA_DEC];
        uint8_t ***m = new uint8_t**[length];
        for (int64_t i = 0; i < length; i++) {
            m[i] = new uint8_t*[2];
            m[i][0] = m_data + (i * 2 * SABER_BYTES_CCA_DEC);
            m[i][1] = m_data + (i * 2 * SABER_BYTES_CCA_DEC) + SABER_BYTES_CCA_DEC;
        }
        uint8_t seed_s[SABER_NOISE_SEEDBYTES];
        uint8_t hash_seed[SABER_SEEDBYTES];
        uint8_t hash_r[SABER_INDCPA_PUBLICKEYBYTES];
        uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES];
        for (int64_t i = 0; i < length; i++) {
            // generate sender's message
            // use r_1 as hash input for pk_0
            emp::Hash hash;
            hash.put(r[i][1], SABER_POLYVECCOMPRESSEDBYTES);
            hash.digest(hash_seed);
            GenMatrix(A, hash_seed);
            POLVECp2BS(hash_r, A[0]);
            memcpy(hash_r + SABER_POLYVECCOMPRESSEDBYTES, r[i][1], SABER_SEEDBYTES);
            // pk_0 = r_0 + H(r_1)
            for (int j = 0; j < SABER_INDCPA_PUBLICKEYBYTES; j++) {
                pk[j] = ((int)(r[i][0][j] + hash_r[j])) % 256;
            }

            // print pk_0 hashed for debug
            std::cout << "pk_0: ";
            for (int j = 0; j <  SABER_INDCPA_PUBLICKEYBYTES; j++) {
                std::cout << (int)pk[j];
            }
            std::cout << std::endl;
            
            // encrypt to get m_0
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            uint8_t ptxt0[SABER_KEYBYTES];
            memset(ptxt0, 0, sizeof(uint8_t)*SABER_KEYBYTES);
            memcpy(ptxt0, &data0[i], sizeof(emp::block));
            indcpa_kem_enc(ptxt0, seed_s, pk, m[i][0]);

            // use r_0 as hash input for pk_1
            emp::Hash hash_;
            hash_.put(r[i][0], SABER_POLYVECCOMPRESSEDBYTES);
            hash_.digest(hash_seed);
            GenMatrix(A, hash_seed);
            POLVECp2BS(hash_r, A[0]);
            memcpy(hash_r + SABER_POLYVECCOMPRESSEDBYTES, r[i][0] + SABER_POLYVECCOMPRESSEDBYTES, SABER_SEEDBYTES);

            // pk_1 = r_1 + H(r_0)
            for (int j = 0; j < SABER_INDCPA_PUBLICKEYBYTES; j++) {
                pk[j] = ((int)(r[i][1][j] + hash_r[j])) % 256;
            }
            // print pk_1 hashed for debug
            std::cout << "pk_1: ";
            for (int j = 0; j <    SABER_INDCPA_PUBLICKEYBYTES; j++) {
                std::cout << (int)pk[j];
            }
            std::cout << std::endl;

            // encrypt to get m_1
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            uint8_t ptxt1[SABER_KEYBYTES];
            memset(ptxt1, 0, sizeof(uint8_t)*SABER_KEYBYTES);
            memcpy(ptxt1, &data1[i], sizeof(emp::block));
            // print ptxt1 for debug
            std::cout << "ptxt1: ";
            for (int j = 0; j < SABER_KEYBYTES; j++) {
                std::cout << (int)ptxt1[j];
            }
            std::cout << std::endl;
            indcpa_kem_enc(ptxt1, seed_s, pk, m[i][1]);
        }

         // decrypt to get data
        for (int64_t i = 0; i < length; i++) {
            uint8_t ptxt[SABER_KEYBYTES];
            indcpa_kem_dec(sk[i], m[i][b[i]], ptxt);
            // print ptxt for debug
            std::cout << "ptxt: ";
            for (int j = 0; j < SABER_KEYBYTES; j++) {
                std::cout << (int)ptxt[j];
            }
            std::cout << std::endl;
            memcpy(data + i, ptxt, sizeof(emp::block));
        }
        // print data0
        for (int i = 0; i < length; i++) {
            std::cout <<  *(reinterpret_cast<uint64_t *>(&data0[i])) << std::dec << std::endl;
        }
        // print data1
        for (int i = 0; i < length; i++) {
            std::cout <<  *(reinterpret_cast<uint64_t *>(&data1[i])) << std::dec << std::endl;
        }
        // print data   
        for (int i = 0; i < length; i++) {
            std::cout <<  *(reinterpret_cast<uint64_t *>(&data[i])) << std::dec << std::endl;
        }
        // free memory
        delete [] r_data;
        for (int64_t i = 0; i < length; i++) {
            delete [] r[i];
        }
        delete [] r;
        for (int64_t i = 0; i < length; i++) {
            delete [] sk[i];
        }
        delete [] sk;
        delete [] m_data;
        for (int64_t i = 0; i < length; i++) {
            delete [] m[i];
        }
        delete [] m;
    return 0;
}