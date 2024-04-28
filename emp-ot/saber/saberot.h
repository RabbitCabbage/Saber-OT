#ifndef SaberOT_H__
#define SaberOT_H__
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
class SaberOT: public OT<IO> {
    public:
    IO* io;
    uint8_t* seed_A; // [SABER_SEEDBYTES];
    bool delete_seed_A = true;

    SaberOT(IO* io, uint8_t* _seed_A = nullptr) {
        this->io = io;
        if (_seed_A == nullptr) {
            seed_A = new uint8_t[SABER_SEEDBYTES];
            randombytes(seed_A, SABER_SEEDBYTES);
            shake128(seed_A, SABER_SEEDBYTES, seed_A, SABER_SEEDBYTES);
        } else {
            seed_A = _seed_A;
            delete_seed_A = false;
        }
    }

    ~SaberOT() {
        if (delete_seed_A) {
            delete[] seed_A;
        }
    }

    void send(const block* data0, const block* data1,int64_t length) override {
        // first get the shared vector r, send r
        // uint16_t r[SABER_L][SABER_N];
        // uint16_t **r = new uint16_t*[SABER_L];
        uint16_t *r = new uint16_t[SABER_L * SABER_N];
        // for (int i = 0; i < SABER_L; ++i) {
        //     r[i] = new uint16_t[SABER_N];
        // }
        randombytes((reinterpret_cast<unsigned char *>(r)), SABER_L * SABER_N * sizeof(uint16_t));
        // shake128((unsigned char *)r, SABER_L * SABER_N * sizeof(uint16_t), (unsigned char *)r, SABER_L * SABER_N * sizeof(uint16_t));
        io->send_data(r, SABER_L * SABER_N * sizeof(uint16_t));
        io->flush();

        // generate A
        uint16_t A[SABER_L][SABER_L][SABER_N];
        GenMatrix(A, seed_A);        

        for (int64_t i = 0; i < length; ++i) {
            // generate s_0, s_1 for the future v_0, v_1
            // initialize s_0, s_1
            uint16_t **s0 = new uint16_t*[SABER_L];
            uint16_t **s1 = new uint16_t*[SABER_L];
            for (int j = 0; j < SABER_L; ++j) {
                s0[j] = new uint16_t[SABER_N];
                s1[j] = new uint16_t[SABER_N];
            }
            uint8_t seed_s[SABER_NOISE_SEEDBYTES];
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            GenSecret(s0, seed_s);
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            GenSecret(s1, seed_s);

            // receive b'_0, compute b'_1
            uint16_t **b0p = new uint16_t*[SABER_L];
            uint16_t **b1p = new uint16_t*[SABER_L];
            *b0p = new uint16_t[SABER_L * SABER_N];
            *b1p = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; ++j) {
                b0p[j] = b0p[j - 1] + SABER_N;
                b1p[j] = b1p[j - 1] + SABER_N;
            }
            io->recv_data(*b0p, SABER_L * SABER_N * sizeof(uint16_t));
            io->flush();
            // b1p = r - b0p;
            // for (int j = 0; j < SABER_L; ++j) {
            //     for (int k = 0; k < SABER_N; ++k) {
            //         b1p[j][k] = r[j][k] - b0p[j][k];
            //     }
            // }
            for (int j = 0; j < SABER_L; ++j) {
                for (int k = 0; k < SABER_N; ++k) {
                    b1p[j][k] = r[j * SABER_N + k] - b0p[j][k];
                }
            }   

            // compute b_0, b_1, send b_0, b_1
            uint16_t **b0 = new uint16_t*[SABER_L];
            uint16_t **b1 = new uint16_t*[SABER_L];
            *b0 = new uint16_t[SABER_L * SABER_N];
            *b1 = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; ++j) {
                b0[j] = b0[j - 1] + SABER_N;
                b1[j] = b1[j - 1] + SABER_N;
            }
            RoundingMul(A, s0, b0, false);
            RoundingMul(A, s1, b1, false);
            io->send_data(*b0, SABER_L * SABER_N * sizeof(uint16_t));
            io->send_data(*b1, SABER_L * SABER_N * sizeof(uint16_t));
            io->flush();

            // compute ciphertexts and send
            uint16_t *cm0 = new uint16_t[SABER_N];
            uint16_t *cm1 = new uint16_t[SABER_N];
            uint16_t v0[SABER_N];
            uint16_t v1[SABER_N];
            // Bits(s0[j][k],SABER_EP,SABER_EP)
            for(int j = 0; j < SABER_L; ++j) {
                for(int k = 0; k < SABER_N; ++k) {
                    s0[j][k] = Bits(s0[j][k], SABER_EP, SABER_EP);
                    s1[j][k] = Bits(s1[j][k], SABER_EP, SABER_EP);
                }
            }
            InnerProd_plush1(b0p, s0, v0);
            InnerProd_plush1(b1p, s1, v1);

            block m[2];
            for (int j = 0; j < SABER_N; ++j) {
                cm0[j] = Bits(v0[j], SABER_EP - 1, SABER_ET);
                cm1[j] = Bits(v1[j], SABER_EP - 1, SABER_ET);
                io->send_data(cm0, SABER_N * sizeof(uint16_t));
                io->flush();
                io->send_data(cm1, SABER_N * sizeof(uint16_t));
                io->flush();
            }
            for (int j = 0; j < SABER_N; ++j) {
                v0[j] = Bits(v0[j], SABER_EP, 1);
                v1[j] = Bits(v1[j], SABER_EP, 1);
            }
            printf("send\n");
            printf("data0: %d\n", *reinterpret_cast<const int*>(&data0[i]));
            printf("data1: %d\n", *reinterpret_cast<const int*>(&data1[i]));
            block a0 = Hash::hash_for_block(v0, SABER_N * 2) ^ data0[i];
            block a1 = Hash::hash_for_block(v1, SABER_N * 2) ^ data1[i];
            printf("a0: %d\n", *reinterpret_cast<int*>(&a0));
            printf("a1: %d\n", *reinterpret_cast<int*>(&a1));
            m[0] = Hash::hash_for_block(v0, SABER_N * 2) ^ data0[i];
            m[1] = Hash::hash_for_block(v1, SABER_N * 2) ^ data1[i];
            printf("m0: %d\n", *reinterpret_cast<int*>(&m[0]));
            printf("m1: %d\n", *reinterpret_cast<int*>(&m[1]));
            io->send_data(m, 2 * sizeof(block));
            io->flush();
            for (int j = 0; j < SABER_L; ++j) {
                delete[] s0[j];
                delete[] s1[j];
            }
            delete[] s0;
            delete[] s1;
            delete[] *b0p;
            delete[] *b1p;
            delete[] *b0;
            delete[] *b1;
            delete[] b0p;
            delete[] b1p;
            delete[] b0;
            delete[] b1;
            delete[] cm0;
            delete[] cm1;
        }
        delete[] r;
    }

    void recv(block* data, const bool* x, int64_t length) override {
        // receive r
        // all the length cases use the same r, from NPOT
        uint16_t *r = new uint16_t[SABER_L * SABER_N];
        io->recv_data(r, SABER_L * SABER_N * sizeof(uint16_t));
        io->flush();

        // generate A
        uint16_t A[SABER_L][SABER_L][SABER_N];
        GenMatrix(A, seed_A);

        for (int64_t i = 0; i < length; ++i) {
            // generate two b' according to the bit x
            uint16_t **sp = new uint16_t*[SABER_L];
            uint8_t seed_s[SABER_NOISE_SEEDBYTES];
            uint16_t **b0p = new uint16_t*[SABER_L];
            uint16_t **b1p = new uint16_t*[SABER_L];
            *b0p = new uint16_t[SABER_L * SABER_N];
            *b1p = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; ++j) {
                b0p[j] = b0p[j - 1] + SABER_N;
                b1p[j] = b1p[j - 1] + SABER_N;
            }
            for (int j = 0; j < SABER_L; ++j) {
                sp[j] = new uint16_t[SABER_N];
            }
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            GenSecret(sp, seed_s);
            if (x[i]) {
                RoundingMul(A, sp, b1p, true);
                // b0p = r - b1p;
                for (int j = 0; j < SABER_L; ++j) {
                    for (int k = 0; k < SABER_N; ++k) {
                        b0p[j][k] = r[j * SABER_N + k] - b1p[j][k];
                    }
                }
            } else {
                RoundingMul(A, sp, b0p, true);
                // b1p = r - b0p;
                for (int j = 0; j < SABER_L; ++j) {
                    for (int k = 0; k < SABER_N; ++k) {
                        b1p[j][k] = r[j * SABER_N + k] - b0p[j][k];
                    }
                }
            }
            // send b'_0
            io->send_data(*b0p, SABER_L * SABER_N * sizeof(uint16_t));
            io->flush();

            // receive two b_0 and b_1
            uint16_t **b0 = new uint16_t*[SABER_L];
            uint16_t **b1 = new uint16_t*[SABER_L];
            *b0 = new uint16_t[SABER_L * SABER_N];
            *b1 = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; ++j) {
                b0[j] = b0[j - 1] + SABER_N;
                b1[j] = b1[j - 1] + SABER_N;
            }
            io->recv_data(*b0, SABER_L * SABER_N * sizeof(uint16_t));
            io->recv_data(*b1, SABER_L * SABER_N * sizeof(uint16_t));
            io->flush();

            // receive two message block
            block m[2];
            uint16_t *cm0 = new uint16_t[SABER_N];
            uint16_t *cm1 = new uint16_t[SABER_N];
            io->recv_data(cm0, SABER_N * sizeof(uint16_t));
            io->recv_data(cm1, SABER_N * sizeof(uint16_t));
            io->flush();
            io->recv_data(m, 2 * sizeof(block));
            io->flush();
            printf("recv\n");
            printf("m0: %d\n", *reinterpret_cast<int*>(&m[0]));
            printf("m1: %d\n", *reinterpret_cast<int*>(&m[1]));
            uint16_t* vp = new uint16_t[SABER_N];
            for (int j = 0; j < SABER_L; ++j) {
                for (int k = 0; k < SABER_N; ++k) { 
                    sp[j][k] = Bits(sp[j][k], SABER_EP, SABER_EP);
                }
            }
            if(x) {
                InnerProd_plush1(b1, sp, vp);
                for (int j = 0; j < SABER_N; ++j) {
                    cm1[j] = Bits(vp[j] - (cm1[j] << (SABER_EP - 1 - SABER_ET)) + h2, SABER_EP, 1);
                }
                block c0 = m[x[i]] ^ Hash::hash_for_block(cm1, SABER_N * 2);
                printf("c0: %d\n", *reinterpret_cast<int*>(&c0));
                data[i] = m[x[i]] ^ Hash::hash_for_block(cm1, SABER_N * 2);
            }
            else {
                InnerProd_plush1(b0, sp, vp);
                for (int j = 0; j < SABER_N; ++j) {
                    cm0[j] = Bits(vp[j] - (cm0[j] << (SABER_EP - 1 - SABER_ET)) + h2, SABER_EP, 1);
                }
                block c1 = m[x[i]] ^ Hash::hash_for_block(cm0, SABER_N * 2);
                printf("c1: %d\n", *reinterpret_cast<int*>(&c1));
                data[i] = m[x[i]] ^ Hash::hash_for_block(cm0, SABER_N * 2);
            }
            for (int j = 0; j < SABER_L; ++j) {
                delete[] sp[j];
            }
            delete[] sp;
            delete[] *b0p;
            delete[] *b1p;
            delete[] *b0;
            delete[] *b1;
            delete[] b0p;
            delete[] b1p;
            delete[] b0;
            delete[] b1;
            delete[] cm0;
            delete[] cm1;
            delete[] vp;
        }
        delete[] r;
    }
};

}



#endif