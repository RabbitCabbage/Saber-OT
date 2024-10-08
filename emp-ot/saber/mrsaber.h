#ifndef MRSABER_H
#define MRSABER_H

#include "emp-ot/saber/SABER_indcpa.h"
#include "emp-ot/saber/SABER_params.h"
#include "emp-ot/saber/pack_unpack.h"

namespace emp {
template<typename IO>
class MRSaber: public OT<IO> {
    public:
    IO* io;

    MRSaber(IO* io) {
        this->io = io;
    }

    ~MRSaber() {
    }

    void send (const block* data0, const block* data1, int64_t length) override {
        uint8_t *r_data = new uint8_t[length * 2 * SABER_INDCPA_PUBLICKEYBYTES];
        io->recv_data(r_data, length * 2 * SABER_INDCPA_PUBLICKEYBYTES);
        uint8_t ***r = new uint8_t**[length];
        for (int64_t i = 0; i < length; i++) {
            r[i] = new uint8_t*[2];
            r[i][0] = r_data + (i * 2 * SABER_INDCPA_PUBLICKEYBYTES);
            r[i][1] = r_data + (i * 2 * SABER_INDCPA_PUBLICKEYBYTES) + SABER_INDCPA_PUBLICKEYBYTES;
        }

        uint8_t *m_data = new uint8_t[length * 2 * SABER_BYTES_CCA_DEC];
        uint8_t ***m = new uint8_t**[length];
        for (int64_t i = 0; i < length; i++) {
            m[i] = new uint8_t*[2];
            m[i][0] = m_data + (i * 2 * SABER_BYTES_CCA_DEC);
            m[i][1] = m_data + (i * 2 * SABER_BYTES_CCA_DEC) + SABER_BYTES_CCA_DEC;
        }
        uint8_t seed_s[SABER_NOISE_SEEDBYTES];
        uint8_t hash_seed[SABER_SEEDBYTES];
        uint16_t A[SABER_L][SABER_L][SABER_N];
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
            memcpy(hash_r + SABER_POLYVECCOMPRESSEDBYTES, r[i][1] + SABER_POLYVECCOMPRESSEDBYTES, SABER_SEEDBYTES);
            // pk_0 = r_0 + H(r_1)
            for (int j = 0; j < SABER_INDCPA_PUBLICKEYBYTES; j++) {
                pk[j] = ((int)(r[i][0][j] + hash_r[j])) % 256;
            }
            // encrypt to get m_0
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            uint8_t ptxt0[SABER_KEYBYTES];
            memset(ptxt0, 0, sizeof(block));
            memcpy(ptxt0, &data0[i], sizeof(block));
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
            // encrypt to get m_1
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            uint8_t ptxt1[SABER_KEYBYTES];
            memset(ptxt1, 0, sizeof(block));
            memcpy(ptxt1, &data1[i], sizeof(block));
            indcpa_kem_enc(ptxt1, seed_s, pk, m[i][1]);
        }      
        // send m
        io->send_data(m_data, length * 2 * SABER_BYTES_CCA_DEC); 

        // free memory
        delete [] r_data;
        for (int64_t i = 0; i < length; i++) {
            delete [] r[i];
        }
        delete [] r;
        delete [] m_data;
        for (int64_t i = 0; i < length; i++) {
            delete [] m[i];
        }
        delete [] m;
    }

    void recv (block* data, const bool* b, int64_t length) override {
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
        // send r
        io->send_data(r_data, length * 2 * SABER_INDCPA_PUBLICKEYBYTES);

        // receive m
        uint8_t *m_data = new uint8_t[length * 2 * SABER_BYTES_CCA_DEC];
        io->recv_data(m_data, length * 2 * SABER_BYTES_CCA_DEC);
        uint8_t ***m = new uint8_t**[length];
        for (int64_t i = 0; i < length; i++) {
            m[i] = new uint8_t*[2];
            m[i][0] = m_data + (i * 2 * SABER_BYTES_CCA_DEC);
            m[i][1] = m_data + (i * 2 * SABER_BYTES_CCA_DEC) + SABER_BYTES_CCA_DEC;
        }
        // decrypt to get data
        for (int64_t i = 0; i < length; i++) {
            uint8_t ptxt[SABER_KEYBYTES];
            indcpa_kem_dec(sk[i], m[i][b[i]], ptxt);
            memcpy(data + i, ptxt, sizeof(block));
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
    }
};
} // namespace emp

#endif // MRSABER_H