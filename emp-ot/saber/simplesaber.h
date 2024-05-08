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

    SimpleSaber(IO* io, uint8_t* _seed_A = nullptr){
        this->io = io;
    }

    ~SimpleSaber(){
    }

    void send(const block *data0, const block* data1, int64_t length) override {
        uint8_t seed_A[SABER_SEEDBYTES];
        randombytes(seed_A, SABER_SEEDBYTES);
        shake128(seed_A, SABER_SEEDBYTES, seed_A, SABER_SEEDBYTES);
        // send seed_A
        io->send_data(seed_A, SABER_SEEDBYTES);
        io->flush();

        // generate A
        uint16_t A[SABER_L][SABER_L][SABER_N];
        GenMatrix(A, seed_A);

        for (int64_t i = 0; i < length; i++) {
            // generate secret
            uint8_t seed_s[SABER_NOISE_SEEDBYTES];
            uint16_t **s = new uint16_t*[SABER_L];
            for (int j = 0; j < SABER_L; j++) {
                s[j] = new uint16_t[SABER_N];
            }
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            GenSecret(s, seed_s);

            // compute b
            uint16_t **b = new uint16_t*[SABER_L];
            *b = new uint16_t[SABER_L * SABER_N];
            memset(*b, 0, SABER_L * SABER_N * sizeof(uint16_t));
            for (int j = 1; j < SABER_L; j++) {
                b[j] = b[j - 1] + SABER_N;
            }
            RoundingMul(A, s, b, 0);
            io->send_data(*b, SABER_L * SABER_N * sizeof(uint16_t));
            io->flush();

            // receive b'
            uint16_t **bp = new uint16_t*[SABER_L];
            *bp = new uint16_t[SABER_L * SABER_N];
            memset(*bp, 0, SABER_L * SABER_N * sizeof(uint16_t));
            for (int j = 1; j < SABER_L; j++) {
                bp[j] = bp[j - 1] + SABER_N;
            }
            io->recv_data(*bp, SABER_L * SABER_N * sizeof(uint16_t));
            io->flush();

            // compute v0 v1
            uint16_t cm0[SABER_N];
            uint16_t cm1[SABER_N];
            uint16_t v0[SABER_N];
            uint16_t v1[SABER_N];
            memset(v0, 0, SABER_N * sizeof(uint16_t));
            memset(v1, 0, SABER_N * sizeof(uint16_t));
            memset(cm0, 0, SABER_N * sizeof(uint16_t));
            memset(cm1, 0, SABER_N * sizeof(uint16_t));
            for (int j = 0; j < SABER_L; j++) {
                for (int k = 0; k < SABER_N; k++) {
                    s[j][k] = Bits(s[j][k], SABER_EP, SABER_EP);
                }
            }
            InnerProd_plush1(bp, s, v0);
            for (int j = 0; j < SABER_L; j++) {
                for (int k = 0; k < SABER_N; k++) {
                    bp[j][k] = bp[j][k] - b[j][k];
                }
            }
            InnerProd_plush1(bp, s, v1);
            // compute cm0 cm1
            block m[2];
            for (int j = 0; j < SABER_N; j++) {
                cm0[j] = Bits(v0[j], SABER_EP - 1, SABER_ET);
                cm1[j] = Bits(v1[j], SABER_EP - 1, SABER_ET);
            }

            io->send_data(cm0, SABER_N * sizeof(uint16_t));
            io->flush();
            io->send_data(cm1, SABER_N * sizeof(uint16_t));
            io->flush();
            for (int j = 0; j < SABER_N; j++) {
                v0[j] = Bits(v0[j], SABER_EP, 1);
                v1[j] = Bits(v1[j], SABER_EP, 1);
            }
            m[0] = Hash::hash_for_block(v0, SABER_N * 2) ^ data0[i];
            m[1] = Hash::hash_for_block(v1, SABER_N * 2) ^ data1[i];
            
            io->send_data(m, 2 * sizeof(block));
            io->flush();
            for (int j = 0; j < SABER_L; ++j) {
                delete[] s[j];
            }
            delete[] s;
            delete[] *b;
            delete[] *bp;
            delete[] b;
            delete[] bp;
        }
    }

    void recv(block* data, const bool* x, int64_t length) override {
        uint8_t seed_A[SABER_SEEDBYTES];
        io->recv_data(seed_A, SABER_SEEDBYTES);
        io->flush();

        // generate A
        uint16_t A[SABER_L][SABER_L][SABER_N];
        GenMatrix(A, seed_A);

        for (int64_t i = 0; i < length; i++) {
            // receive b
            uint16_t **b = new uint16_t*[SABER_L];
            *b = new uint16_t[SABER_L * SABER_N];
            memset(*b, 0, SABER_L * SABER_N * sizeof(uint16_t));
            for (int j = 1; j < SABER_L; j++) {
                b[j] = b[j - 1] + SABER_N;
            }
            io->recv_data(*b, SABER_L * SABER_N * sizeof(uint16_t));
            io->flush();

            // generate secret s'
            uint8_t seed_s[SABER_NOISE_SEEDBYTES];
            uint16_t **sp = new uint16_t*[SABER_L];
            for (int j = 0; j < SABER_L; j++) {
                sp[j] = new uint16_t[SABER_N];
            }
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            GenSecret(sp, seed_s);

            // compute b'
            uint16_t **bp = new uint16_t*[SABER_L];
            *bp = new uint16_t[SABER_L * SABER_N];
            memset(*bp, 0, SABER_L * SABER_N * sizeof(uint16_t));
            for (int j = 1; j < SABER_L; j++) {
                bp[j] = bp[j - 1] + SABER_N;
            }
            RoundingMul(A, sp, bp, 1);
            if (x[i]) {
                for (int j = 0; j < SABER_L; j++) {
                    for (int k = 0; k < SABER_N; k++) {
                        bp[j][k] = bp[j][k] + b[j][k];
                    }
                }
            }
            io->send_data(*bp, SABER_L * SABER_N * sizeof(uint16_t));
            io->flush();

            // compute vp
            uint16_t* vp = new uint16_t[SABER_N];
            memset(vp, 0, SABER_N * sizeof(uint16_t));
            for (int j = 0; j < SABER_L; ++j) {
                for (int k = 0; k < SABER_N; ++k) { 
                    sp[j][k] = Bits(sp[j][k], SABER_EP, SABER_EP);
                }
            }
            InnerProd_plush1(b, sp, vp);
            // recv cm0 cm1 e0 e1
            block m[2];
            uint16_t *cm0 = new uint16_t[SABER_N];
            uint16_t *cm1 = new uint16_t[SABER_N];
            memset(cm0, 0, SABER_N * sizeof(uint16_t));
            memset(cm1, 0, SABER_N * sizeof(uint16_t));
            io->recv_data(cm0, SABER_N * sizeof(uint16_t));
            io->recv_data(cm1, SABER_N * sizeof(uint16_t));
            io->flush();
            io->recv_data(m, 2 * sizeof(block));
            io->flush();
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
            for (int j = 0; j < SABER_L; ++j) {
                delete[] sp[j];
            }
            delete[] sp;
            delete[] *b;
            delete[] *bp;
            delete[] b;
            delete[] bp;
            delete[] cm0;
            delete[] cm1;
            delete[] vp;
        }
    }
};
}

#endif