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


namespace emp {
template<typename IO>
class NPSaber2: public OT<IO> {
    public:
    IO* io;

    NPSaber2(IO* io, uint8_t* _seed_A = nullptr) {
        this->io = io;
    }

    ~NPSaber2() {

    }

    void send(const block* data0, const block* data1,int64_t length) override {
        uint8_t seed_A[SABER_SEEDBYTES];
        randombytes(seed_A, SABER_SEEDBYTES);
        shake128(seed_A, SABER_SEEDBYTES, seed_A, SABER_SEEDBYTES);
        // send seed_A
        io->send_data(seed_A, SABER_SEEDBYTES);
        // io->flush();

        // generate r and send
        uint16_t *r = new uint16_t[SABER_L * SABER_N];
        randombytes((reinterpret_cast<unsigned char *>(r)), SABER_L * SABER_N * sizeof(uint16_t));
        io->send_data(r, SABER_L * SABER_N * sizeof(uint16_t));
        io->flush();

        // generate A
        uint16_t A[SABER_L][SABER_L][SABER_N];
        GenMatrix(A, seed_A);    
        uint16_t ***s0 = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            s0[i] = new uint16_t*[SABER_L];
            for (int j = 0; j < SABER_L; j++) {
                s0[i][j] = new uint16_t[SABER_N];
            }
        }

        uint16_t ***b0p = new uint16_t**[length];
        uint16_t ***b1p = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            b0p[i] = new uint16_t*[SABER_L];
            b1p[i] = new uint16_t*[SABER_L];
            *b0p[i] = new uint16_t[SABER_L * SABER_N];
            *b1p[i] = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; j++) {
                b0p[i][j] = b0p[i][j - 1] + SABER_N;
                b1p[i][j] = b1p[i][j - 1] + SABER_N;
            }
        }

        uint16_t ***b0 = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            b0[i] = new uint16_t*[SABER_L];
            *b0[i] = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; j++) {
                b0[i][j] = b0[i][j - 1] + SABER_N;
            }
        }

        uint16_t *cm0 = new uint16_t[SABER_N];
        uint16_t *cm1 = new uint16_t[SABER_N];
        uint8_t seed_s[SABER_NOISE_SEEDBYTES];
        uint16_t *v0 = new uint16_t[SABER_N];
        uint16_t *v1 = new uint16_t[SABER_N];

        for (int64_t i = 0; i < length; ++i) {
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            GenSecret(s0[i], seed_s);

            io->recv_data(*b0p[i], SABER_L * SABER_N * sizeof(uint16_t));
            for (int j = 0; j < SABER_L; ++j) {
                for (int k = 0; k < SABER_N; ++k) {
                    b1p[i][j][k] = r[j * SABER_N + k] - b0p[i][j][k];
                }
            }

            memset(*b0[i], 0, SABER_L * SABER_N * sizeof(uint16_t));
            RoundingMul(A, s0[i], b0[i], 0);
            io->send_data(*b0[i], SABER_L * SABER_N * sizeof(uint16_t));
        }
        io->flush();
        
        block m[2];
        for (int64_t i = 0; i < length; ++i) {
            // compute ciphertexts and send
            memset(v0, 0, SABER_N * sizeof(uint16_t));
            memset(v1, 0, SABER_N * sizeof(uint16_t));
            memset(cm0, 0, SABER_N * sizeof(uint16_t));
            memset(cm1, 0, SABER_N * sizeof(uint16_t));
            for(int j = 0; j < SABER_L; ++j) {
                for(int k = 0; k < SABER_N; ++k) {
                    s0[i][j][k] = Bits(s0[i][j][k], SABER_EP, SABER_EP);
                }
            }
            InnerProd_plush1(b0p[i], s0[i], v0);
            InnerProd_plush1(b1p[i], s0[i], v1);

            for (int j = 0; j < SABER_N; ++j) {
                cm0[j] = Bits(v0[j], SABER_EP - 1, SABER_ET);
                cm1[j] = Bits(v1[j], SABER_EP - 1, SABER_ET);
            }

            io->send_data(cm0, SABER_N * sizeof(uint16_t));
            io->send_data(cm1, SABER_N * sizeof(uint16_t));

            for (int j = 0; j < SABER_N; ++j) {
                v0[j] = Bits(v0[j], SABER_EP, 1);
                v1[j] = Bits(v1[j], SABER_EP, 1);
            }
            m[0] = Hash::hash_for_block(v0, SABER_N * 2) ^ data0[i];
            m[1] = Hash::hash_for_block(v1, SABER_N * 2) ^ data1[i];
            io->send_data(m, 2 * sizeof(block));
        }

        for (int64_t i = 0; i < length; ++i) {
            for (int j = 0; j < SABER_L; ++j) {
                delete[] s0[i][j];
            }
            delete[] s0[i];
            delete[] *b0p[i];
            delete[] *b1p[i];
            delete[] *b0[i];
            delete[] b0p[i];
            delete[] b1p[i];
            delete[] b0[i];
        }
        delete[] s0;
        delete[] b0p;
        delete[] b1p;
        delete[] b0;
        delete[] v0;
        delete[] v1;
        delete[] cm0;
        delete[] cm1;
        delete[] r;
    }

    void recv(block* data, const bool* x, int64_t length) override {
        uint8_t seed_A[SABER_SEEDBYTES];
        io->recv_data(seed_A, SABER_SEEDBYTES);
        
        uint16_t *r = new uint16_t[SABER_L * SABER_N];
        io->recv_data(r, SABER_L * SABER_N * sizeof(uint16_t));
        io->flush();

        uint16_t A[SABER_L][SABER_L][SABER_N];
        GenMatrix(A, seed_A);

        uint16_t ***sp = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            sp[i] = new uint16_t*[SABER_L];
            for (int j = 0; j < SABER_L; j++) {
                sp[i][j] = new uint16_t[SABER_N];
            }
        }

        uint8_t seed_s[SABER_NOISE_SEEDBYTES];

        uint16_t ***b0p = new uint16_t**[length];
        uint16_t ***b1p = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            b0p[i] = new uint16_t*[SABER_L];
            b1p[i] = new uint16_t*[SABER_L];
            *b0p[i] = new uint16_t[SABER_L * SABER_N];
            *b1p[i] = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; j++) {
                b0p[i][j] = b0p[i][j - 1] + SABER_N;
                b1p[i][j] = b1p[i][j - 1] + SABER_N;
            }
        }
        
        uint16_t ***b0 = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            b0[i] = new uint16_t*[SABER_L];
            *b0[i] = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; j++) {
                b0[i][j] = b0[i][j - 1] + SABER_N;
            }
        }

        for (int64_t i = 0; i < length; ++i) {
            // memsset to 0!!!
            memset(*b0p[i], 0, SABER_L * SABER_N * sizeof(uint16_t));
            memset(*b1p[i], 0, SABER_L * SABER_N * sizeof(uint16_t));
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            GenSecret(sp[i], seed_s);
            if (x[i]) {
                RoundingMul(A, sp[i], b1p[i], 1);
                for (int j = 0; j < SABER_L; ++j) {
                    for (int k = 0; k < SABER_N; ++k) {
                        b0p[i][j][k] = r[j * SABER_N + k] - b1p[i][j][k];
                    }
                }
            } else {
                RoundingMul(A, sp[i], b0p[i], 1);
                for (int j = 0; j < SABER_L; ++j) {
                    for (int k = 0; k < SABER_N; ++k) {
                        b1p[i][j][k] = r[j * SABER_N + k] - b0p[i][j][k];
                    }
                }
            }

            io->send_data(*b0p[i], SABER_L * SABER_N * sizeof(uint16_t));

            io->recv_data(*b0[i], SABER_L * SABER_N * sizeof(uint16_t));
        }
        io->flush();

        block m[2];
        uint16_t *cm0 = new uint16_t[SABER_N];
        uint16_t *cm1 = new uint16_t[SABER_N];
        uint16_t* vp = new uint16_t[SABER_N];
        for (int64_t i = 0; i < length; ++i) {
            memset(cm0, 0, SABER_N * sizeof(uint16_t));
            memset(cm1, 0, SABER_N * sizeof(uint16_t));
            io->recv_data(cm0, SABER_N * sizeof(uint16_t));
            io->recv_data(cm1, SABER_N * sizeof(uint16_t));
            io->recv_data(m, 2 * sizeof(block));

            memset(vp, 0, SABER_N * sizeof(uint16_t));
            for (int j = 0; j < SABER_L; ++j) {
                for (int k = 0; k < SABER_N; ++k) { 
                    sp[i][j][k] = Bits(sp[i][j][k], SABER_EP, SABER_EP);
                }
            }
            if(x[i]) {
                InnerProd_plush1(b0[i], sp[i], vp);
                for (int j = 0; j < SABER_N; ++j) {
                    cm1[j] = Bits(vp[j] - (cm1[j] << (SABER_EP - 1 - SABER_ET)) + h2, SABER_EP, 1);
                }
                data[i] = m[x[i]] ^ Hash::hash_for_block(cm1, SABER_N * 2);
            }
            else {
                InnerProd_plush1(b0[i], sp[i], vp);
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
            delete[] *b0p[i];
            delete[] *b1p[i];
            delete[] *b0[i];
            delete[] b0p[i];
            delete[] b1p[i];
            delete[] b0[i];
        }
        delete[] sp;
        delete[] b0p;
        delete[] b1p;
        delete[] b0;
        delete[] cm0;
        delete[] cm1;
        delete[] vp;
        delete[] r;
    }
};

}



#endif