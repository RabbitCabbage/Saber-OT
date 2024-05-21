#ifndef SimpleSaber_h
#define SimpleSaber_h
#include "emp-ot/saber/api.h"
#include "emp-ot/saber/cbd.h"
#include "emp-ot/saber/fips202.h"
#include "emp-ot/saber/pack_unpack.h"
#include "emp-ot/saber/poly.h"
#include "emp-ot/saber/poly_mul.h"
#include "emp-ot/saber/rng.h"
#include "emp-ot/saber/SABER_params.h"

namespace emp{
template<typename IO>
class SimpleSaber: public OT<IO>{
   public:
    IO* io;
    uint8_t seed_A[SABER_SEEDBYTES];

    SimpleSaber(IO* io, uint8_t* _seed_A = nullptr) {
        this->io = io;
        if(_seed_A != nullptr) {
            memcpy(seed_A, _seed_A, SABER_SEEDBYTES);
        }
    }

    ~SimpleSaber(){
    }

    void send(const block *data0, const block* data1, int64_t length) override {
        uint16_t A[SABER_L][SABER_L][SABER_N];
        GenMatrix(A, seed_A);

        uint16_t ***s = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            s[i] = new uint16_t*[SABER_L];
            for (int j = 0; j < SABER_L; j++) {
                s[i][j] = new uint16_t[SABER_N];
            }
        }

        uint16_t ***b = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            b[i] = new uint16_t*[SABER_L];
            *b[i] = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; j++) {
                b[i][j] = b[i][j - 1] + SABER_N;
            }
        }

        uint16_t ***bp = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            bp[i] = new uint16_t*[SABER_L];
            *bp[i] = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; j++) {
                bp[i][j] = bp[i][j - 1] + SABER_N;
            }
        }

        uint8_t seed_s[SABER_NOISE_SEEDBYTES];
        for (int64_t i = 0; i < length; i++) {
            // generate secret
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            GenSecret(s[i], seed_s);

            memset(*b[i], 0, SABER_L * SABER_N * sizeof(uint16_t));
            RoundingMul(A, s[i], b[i], 0);
            // compress the vector b using Saber KEM 
            uint8_t compressed_b[SABER_POLYVECCOMPRESSEDBYTES];
            POLVECp2BS(compressed_b, b[i]);

            io->send_data(compressed_b, SABER_POLYVECCOMPRESSEDBYTES * sizeof(uint8_t));
        }
        // io->flush();

        for (int64_t i = 0; i < length; i++) {
            uint8_t compressed_bp [SABER_POLYVECCOMPRESSEDBYTES];
            io->recv_data(compressed_bp, SABER_POLYVECCOMPRESSEDBYTES * sizeof(uint8_t));
            //io->recv_data(*bp[i], SABER_L * SABER_N * sizeof(uint16_t));
            BS2POLVECp(compressed_bp, bp[i]);
        }
        io->flush();

        // compute v0 v1
        uint16_t *cm0 = new uint16_t[SABER_N];
        uint16_t *cm1 = new uint16_t[SABER_N];
        uint16_t *v0 = new uint16_t[SABER_N];
        uint16_t *v1 = new uint16_t[SABER_N];
        for (int64_t i = 0; i < length; i++) {
            memset(v0, 0, SABER_N * sizeof(uint16_t));
            memset(v1, 0, SABER_N * sizeof(uint16_t));
            memset(cm0, 0, SABER_N * sizeof(uint16_t));
            memset(cm1, 0, SABER_N * sizeof(uint16_t));
            for (int j = 0; j < SABER_L; j++) {
                for (int k = 0; k < SABER_N; k++) {
                    s[i][j][k] = Bits(s[i][j][k], SABER_EP, SABER_EP);
                }
            }
            InnerProd_plush1(bp[i], s[i], v0);
            for (int j = 0; j < SABER_L; j++) {
                for (int k = 0; k < SABER_N; k++) {
                    // bp[j][k] = bp[j][k] - b[j][k];
                    bp[i][j][k] = Bits(bp[i][j][k] - b[i][j][k], SABER_EP, SABER_EP);
                }
            }
            InnerProd_plush1(bp[i], s[i], v1);
            // compute cm0 cm1
            block m[2];
            for (int j = 0; j < SABER_N; j++) {
                cm0[j] = Bits(v0[j], SABER_EP - 1, SABER_ET);
                cm1[j] = Bits(v1[j], SABER_EP - 1, SABER_ET);
            }

            io->send_data(cm0, SABER_N * sizeof(uint16_t));
            // io->flush();
            io->send_data(cm1, SABER_N * sizeof(uint16_t));
            // io->flush();
            for (int j = 0; j < SABER_N; j++) {
                v0[j] = Bits(v0[j], SABER_EP, 1);
                v1[j] = Bits(v1[j], SABER_EP, 1);
            }
            m[0] = Hash::hash_for_block(v0, SABER_N * 2) ^ data0[i];
            m[1] = Hash::hash_for_block(v1, SABER_N * 2) ^ data1[i];
            
            io->send_data(m, 2 * sizeof(block));
            // io->flush();
        }
        for (int64_t i = 0; i < length; i++) {
            for (int j = 0; j < SABER_L; j++) {
                delete[] s[i][j];
            }
            delete[] s[i];
            delete[] *b[i];
            delete[] *bp[i];
            delete[] b[i];
            delete[] bp[i];
        }
        delete[] cm0;
        delete[] cm1;
        delete[] v0;
        delete[] v1;
        delete[] s;
        delete[] b;
        delete[] bp;
    }

    void recv(block* data, const bool* x, int64_t length) override {
        uint16_t A[SABER_L][SABER_L][SABER_N];
        GenMatrix(A, seed_A);

        uint16_t ***b = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            b[i] = new uint16_t*[SABER_L];
            *b[i] = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; j++) {
                b[i][j] = b[i][j - 1] + SABER_N;
            }
        }

        uint16_t ***bp = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            bp[i] = new uint16_t*[SABER_L];
            *bp[i] = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; j++) {
                bp[i][j] = bp[i][j - 1] + SABER_N;
            }
        }

        uint16_t ***sp = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            sp[i] = new uint16_t*[SABER_L];
            for (int j = 0; j < SABER_L; j++) {
                sp[i][j] = new uint16_t[SABER_N];
            }
        }
        
        uint16_t* vp = new uint16_t[SABER_N];
        uint16_t *cm0 = new uint16_t[SABER_N];
        uint16_t *cm1 = new uint16_t[SABER_N];

        uint8_t b_compressed[SABER_POLYVECCOMPRESSEDBYTES];

        for (int64_t i = 0; i < length; i++) {
            io->recv_data(b_compressed, SABER_POLYVECCOMPRESSEDBYTES * sizeof(uint8_t));
            BS2POLVECp(b_compressed, b[i]);
            // io->recv_data(*b[i], SABER_L * SABER_N * sizeof(uint16_t));
        }
        // io->flush();

        uint8_t seed_s[SABER_NOISE_SEEDBYTES];
        for (int64_t i = 0; i < length; i++) {
            // generate secret s'
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            GenSecret(sp[i], seed_s);

            memset(*bp[i], 0, SABER_L * SABER_N * sizeof(uint16_t));
            RoundingMul(A, sp[i], bp[i], 1);
            // keep constant time
            // if (x[i]) {
                for (int j = 0; j < SABER_L; j++) {
                    for (int k = 0; k < SABER_N; k++) {
                        // bp[j][k] = bp[j][k] + b[j][k];
                        bp[i][j][k] = Bits(bp[i][j][k] + x[i] * b[i][j][k], SABER_EP, SABER_EP);                    
                    }
                }
            // }
            uint8_t compressed_bp[SABER_POLYVECCOMPRESSEDBYTES];
            POLVECp2BS(compressed_bp, bp[i]);
            io->send_data(compressed_bp, SABER_POLYVECCOMPRESSEDBYTES * sizeof(uint8_t));
            //io->send_data(*bp[i], SABER_L * SABER_N * sizeof(uint16_t));
        }
        io->flush();

        for (int64_t i = 0; i < length; i++) {
            // compute vp
            memset(vp, 0, SABER_N * sizeof(uint16_t));
            for (int j = 0; j < SABER_L; ++j) {
                for (int k = 0; k < SABER_N; ++k) { 
                    sp[i][j][k] = Bits(sp[i][j][k], SABER_EP, SABER_EP);
                }
            }
            InnerProd_plush1(b[i], sp[i], vp);
            // recv cm0 cm1 e0 e1
            block m[2];
            memset(cm0, 0, SABER_N * sizeof(uint16_t));
            memset(cm1, 0, SABER_N * sizeof(uint16_t));
            io->recv_data(cm0, SABER_N * sizeof(uint16_t));
            io->recv_data(cm1, SABER_N * sizeof(uint16_t));
            // io->flush();
            io->recv_data(m, 2 * sizeof(block));
            // io->flush();
            if(x[i]) {
                for (int j = 0; j < SABER_N; ++j) {
                    cm1[j] = Bits(vp[j] - (cm1[j] << (SABER_EP - 1 - SABER_ET)) + h2, SABER_EP, 1);
                }
                data[i] = m[x[i]] ^ Hash::hash_for_block(cm1, SABER_N * 2);
            } else {
                for (int j = 0; j < SABER_N; ++j) {
                    cm0[j] = Bits(vp[j] - (cm0[j] << (SABER_EP - 1 - SABER_ET)) + h2, SABER_EP, 1);
                }
                data[i] = m[x[i]] ^ Hash::hash_for_block(cm0, SABER_N * 2);
            }
        }
        for (int64_t i = 0; i < length; i++) {
            for (int j = 0; j < SABER_L; j++) {
                delete[] sp[i][j];
            }
            delete[] sp[i];
            delete[] *b[i];
            delete[] *bp[i];
            delete[] b[i];
            delete[] bp[i];
        }
        delete[] cm0;
        delete[] cm1;
        delete[] vp;
        delete[] b;
        delete[] bp;
        delete[] sp;
    }
};
}

#endif