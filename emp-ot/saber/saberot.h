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

    SaberOT(IO* io, uint8_t* _seed_A = nullptr) {
        this->io = io;
    }

    ~SaberOT() {

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
        // print the matrix A, check the correctness
        // block a = Hash::hash_for_block(A, SABER_L * SABER_L * SABER_N * 2);
        // std::cout << "send hash of A: " << *(reinterpret_cast<const uint64_t *>(&a)) << std::endl;
        // generate s_0, s_1 for the future v_0, v_1
        // initialize s_0, s_1
        uint16_t ***s0 = new uint16_t**[length];
        uint16_t ***s1 = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            s0[i] = new uint16_t*[SABER_L];
            s1[i] = new uint16_t*[SABER_L];
            for (int j = 0; j < SABER_L; j++) {
                s0[i][j] = new uint16_t[SABER_N];
                s1[i][j] = new uint16_t[SABER_N];
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
        uint16_t ***b1 = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            b0[i] = new uint16_t*[SABER_L];
            b1[i] = new uint16_t*[SABER_L];
            *b0[i] = new uint16_t[SABER_L * SABER_N];
            *b1[i] = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; j++) {
                b0[i][j] = b0[i][j - 1] + SABER_N;
                b1[i][j] = b1[i][j - 1] + SABER_N;
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
            randombytes(seed_s, SABER_NOISE_SEEDBYTES);
            GenSecret(s1[i], seed_s);

            io->recv_data(*b0p[i], SABER_L * SABER_N * sizeof(uint16_t));
            // io->flush();
            // b1p = r - b0p;
            for (int j = 0; j < SABER_L; ++j) {
                for (int k = 0; k < SABER_N; ++k) {
                    b1p[i][j][k] = r[j * SABER_N + k] - b0p[i][j][k];
                }
            }

            // ================== debug ==================
            // block hash_b0p = Hash::hash_for_block(*b0p, SABER_L * SABER_N * 2);
            // block hash_b1p = Hash::hash_for_block(*b1p, SABER_L * SABER_N * 2);
            // std::cout << "sender b0p: " << *(reinterpret_cast<const uint64_t *>(&hash_b0p)) << std::endl;
            // std::cout << "sender b1p: " << *(reinterpret_cast<const uint64_t *>(&hash_b1p)) << std::endl;
            // ============================================

            // memsset to 0!!!
            memset(*b0[i], 0, SABER_L * SABER_N * sizeof(uint16_t));
            memset(*b1[i], 0, SABER_L * SABER_N * sizeof(uint16_t));
            RoundingMul(A, s0[i], b0[i], 0);
            RoundingMul(A, s1[i], b1[i], 0);
            io->send_data(*b0[i], SABER_L * SABER_N * sizeof(uint16_t));
            io->send_data(*b1[i], SABER_L * SABER_N * sizeof(uint16_t));
        }
        io->flush();

            // ================== debug ==================
            // block hash_b0 = Hash::hash_for_block(*b0, SABER_L * SABER_N * 2);
            // block hash_b1 = Hash::hash_for_block(*b1, SABER_L * SABER_N * 2);
            // std::cout << "sender b0: " << *(reinterpret_cast<const uint64_t *>(&hash_b0)) << std::endl;
            // std::cout << "sender b1: " << *(reinterpret_cast<const uint64_t *>(&hash_b1)) << std::endl;
            // ============================================
        
        block m[2];
        for (int64_t i = 0; i < length; ++i) {
            // compute ciphertexts and send
            // memsset to 0!!! v must be initialized to 0
            memset(v0, 0, SABER_N * sizeof(uint16_t));
            memset(v1, 0, SABER_N * sizeof(uint16_t));
            memset(cm0, 0, SABER_N * sizeof(uint16_t));
            memset(cm1, 0, SABER_N * sizeof(uint16_t));
            // Bits(s0[j][k],SABER_EP,SABER_EP)
            for(int j = 0; j < SABER_L; ++j) {
                for(int k = 0; k < SABER_N; ++k) {
                    s0[i][j][k] = Bits(s0[i][j][k], SABER_EP, SABER_EP);
                    s1[i][j][k] = Bits(s1[i][j][k], SABER_EP, SABER_EP);
                }
            }
            InnerProd_plush1(b0p[i], s0[i], v0);
            InnerProd_plush1(b1p[i], s1[i], v1);

            for (int j = 0; j < SABER_N; ++j) {
                cm0[j] = Bits(v0[j], SABER_EP - 1, SABER_ET);
                cm1[j] = Bits(v1[j], SABER_EP - 1, SABER_ET);
            }
            
            // ================== debug ==================
            // block a0 = Hash::hash_for_block(cm0, SABER_N * 2);
            // block a1 = Hash::hash_for_block(cm1, SABER_N * 2);
            // std::cout << "sender cm0: " << *(reinterpret_cast<const uint64_t *>(&a0)) << std::endl;
            // std::cout << "sender cm1: " << *(reinterpret_cast<const uint64_t *>(&a1)) << std::endl;
            // block hash_v0 = Hash::hash_for_block(v0, SABER_N * 2);
            // block hash_v1 = Hash::hash_for_block(v1, SABER_N * 2);
            // std::cout << "sender v0: " << *(reinterpret_cast<const uint64_t *>(&hash_v0)) << std::endl;
            // std::cout << "sender v1: " << *(reinterpret_cast<const uint64_t *>(&hash_v1)) << std::endl;
            // ============================================

            io->send_data(cm0, SABER_N * sizeof(uint16_t));
            io->send_data(cm1, SABER_N * sizeof(uint16_t));

            for (int j = 0; j < SABER_N; ++j) {
                v0[j] = Bits(v0[j], SABER_EP, 1);
                v1[j] = Bits(v1[j], SABER_EP, 1);
            }
            m[0] = Hash::hash_for_block(v0, SABER_N * 2) ^ data0[i];
            m[1] = Hash::hash_for_block(v1, SABER_N * 2) ^ data1[i];

            // ================== debug ==================
            // hash_v0 = Hash::hash_for_block(v0, SABER_N * 2);
            // hash_v1 = Hash::hash_for_block(v1, SABER_N * 2);
            // std::cout << "sender v0->2: " << *(reinterpret_cast<const uint64_t *>(&hash_v0)) << std::endl;
            // std::cout << "sender v1->2: " << *(reinterpret_cast<const uint64_t *>(&hash_v1)) << std::endl;
            // std::cout << "sender data0: " << *(reinterpret_cast<const uint64_t *>(&data0[i])) << std::endl;
            // std::cout << "sender data1: " << *(reinterpret_cast<const uint64_t *>(&data1[i])) << std::endl;
            // std::cout << "sender m[0]: " << *(reinterpret_cast<const uint64_t *>(&m[0])) << std::endl;
            // std::cout << "sender m[1]: " << *(reinterpret_cast<const uint64_t *>(&m[1])) << std::endl;
            // ============================================

            io->send_data(m, 2 * sizeof(block));
        }

        for (int64_t i = 0; i < length; ++i) {
            for (int j = 0; j < SABER_L; ++j) {
                delete[] s0[i][j];
                delete[] s1[i][j];
            }
            delete[] s0[i];
            delete[] s1[i];
            delete[] *b0p[i];
            delete[] *b1p[i];
            delete[] *b0[i];
            delete[] *b1[i];
            delete[] b0p[i];
            delete[] b1p[i];
            delete[] b0[i];
            delete[] b1[i];
        }
        delete[] s0;
        delete[] s1;
        delete[] b0p;
        delete[] b1p;
        delete[] b0;
        delete[] b1;
        delete[] v0;
        delete[] v1;
        delete[] cm0;
        delete[] cm1;
        delete[] r;
    }

    void recv(block* data, const bool* x, int64_t length) override {
        uint8_t seed_A[SABER_SEEDBYTES];
        // receive seed_A
        io->recv_data(seed_A, SABER_SEEDBYTES);
        // io->flush();
        
        // receive r
        // all the length cases use the same r, from NPOT
        uint16_t *r = new uint16_t[SABER_L * SABER_N];
        io->recv_data(r, SABER_L * SABER_N * sizeof(uint16_t));
        io->flush();

        // generate A
        uint16_t A[SABER_L][SABER_L][SABER_N];
        GenMatrix(A, seed_A);
        // block a = Hash::hash_for_block(A, SABER_L * SABER_L * SABER_N * 2);
        // std::cout << "recv hash of A: " << *(reinterpret_cast<const uint64_t *>(&a)) << std::endl;

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
        uint16_t ***b1 = new uint16_t**[length];
        for (int64_t i = 0; i < length; i++) {
            b0[i] = new uint16_t*[SABER_L];
            b1[i] = new uint16_t*[SABER_L];
            *b0[i] = new uint16_t[SABER_L * SABER_N];
            *b1[i] = new uint16_t[SABER_L * SABER_N];
            for (int j = 1; j < SABER_L; j++) {
                b0[i][j] = b0[i][j - 1] + SABER_N;
                b1[i][j] = b1[i][j - 1] + SABER_N;
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
                // b0p = r - b1p;
                for (int j = 0; j < SABER_L; ++j) {
                    for (int k = 0; k < SABER_N; ++k) {
                        b0p[i][j][k] = r[j * SABER_N + k] - b1p[i][j][k];
                    }
                }
            } else {
                RoundingMul(A, sp[i], b0p[i], 1);
                // b1p = r - b0p;
                for (int j = 0; j < SABER_L; ++j) {
                    for (int k = 0; k < SABER_N; ++k) {
                        b1p[i][j][k] = r[j * SABER_N + k] - b0p[i][j][k];
                    }
                }
            }

            // ================== debug ==================
            // block hash_b0p = Hash::hash_for_block(*b0p, SABER_L * SABER_N * 2);
            // block hash_b1p = Hash::hash_for_block(*b1p, SABER_L * SABER_N * 2);
            // std::cout << "recver b0p: " << *(reinterpret_cast<const uint64_t *>(&hash_b0p)) << std::endl;
            // std::cout << "recver b1p: " << *(reinterpret_cast<const uint64_t *>(&hash_b1p)) << std::endl;
            // ============================================

            // send b'_0
            io->send_data(*b0p[i], SABER_L * SABER_N * sizeof(uint16_t));
            // io->flush();

            io->recv_data(*b0[i], SABER_L * SABER_N * sizeof(uint16_t));
            io->recv_data(*b1[i], SABER_L * SABER_N * sizeof(uint16_t));
        }
        io->flush();

            // ================== debug ==================
            // block hash_b0 = Hash::hash_for_block(*b0, SABER_L * SABER_N * 2);
            // block hash_b1 = Hash::hash_for_block(*b1, SABER_L * SABER_N * 2);
            // std::cout << "recver b0: " << *(reinterpret_cast<const uint64_t *>(&hash_b0)) << std::endl;
            // std::cout << "recver b1: " << *(reinterpret_cast<const uint64_t *>(&hash_b1)) << std::endl;
            // ============================================

            // receive two message block
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

            // ================== debug ==================
            // block a0 = Hash::hash_for_block(cm0, SABER_N * 2);
            // block a1 = Hash::hash_for_block(cm1, SABER_N * 2);
            // std::cout << "recver cm0: " << *(reinterpret_cast<const uint64_t *>(&a0)) << std::endl;
            // std::cout << "recver cm1: " << *(reinterpret_cast<const uint64_t *>(&a1)) << std::endl;
            // std::cout << "recver m[0]: " << *(reinterpret_cast<const uint64_t *>(&m[0])) << std::endl;
            // std::cout << "recver m[1]: " << *(reinterpret_cast<const uint64_t *>(&m[1])) << std::endl;
            // ============================================

            memset(vp, 0, SABER_N * sizeof(uint16_t));
            for (int j = 0; j < SABER_L; ++j) {
                for (int k = 0; k < SABER_N; ++k) { 
                    sp[i][j][k] = Bits(sp[i][j][k], SABER_EP, SABER_EP);
                }
            }
            if(x[i]) {
                InnerProd_plush1(b1[i], sp[i], vp);
                for (int j = 0; j < SABER_N; ++j) {
                    cm1[j] = Bits(vp[j] - (cm1[j] << (SABER_EP - 1 - SABER_ET)) + h2, SABER_EP, 1);
                }
                data[i] = m[x[i]] ^ Hash::hash_for_block(cm1, SABER_N * 2);

                // ================== debug ==================
                // block hash_cm1 = Hash::hash_for_block(cm1, SABER_N * 2);
                // std::cout << "recver cm1->2: " << *(reinterpret_cast<const uint64_t *>(&hash_cm1)) << std::endl;
                // // block c1 = m[x[i]] ^ Hash::hash_for_block(cm1, SABER_N * 2);
                // // std::cout << "recv hash1: " << *(reinterpret_cast<const uint64_t *>(&c1)) << std::endl;
                // std::cout << "recver x = 1, data[i]: " << *(reinterpret_cast<const uint64_t *>(&data[i])) << std::endl;
                // ============================================
            }
            else {
                InnerProd_plush1(b0[i], sp[i], vp);
                for (int j = 0; j < SABER_N; ++j) {
                    cm0[j] = Bits(vp[j] - (cm0[j] << (SABER_EP - 1 - SABER_ET)) + h2, SABER_EP, 1);
                }
                data[i] = m[x[i]] ^ Hash::hash_for_block(cm0, SABER_N * 2);

                // ================== debug ==================
                // block hash_cm0 = Hash::hash_for_block(cm0, SABER_N * 2);
                // std::cout << "recver cm0->2: " << *(reinterpret_cast<const uint64_t *>(&hash_cm0)) << std::endl;
                // // block c0 = m[x[i]] ^ Hash::hash_for_block(cm0, SABER_N * 2);
                // // std::cout << "recv hash0: " << *(reinterpret_cast<const uint64_t *>(&c0)) << std::endl;
                // std::cout << "recver x = 0, data[i]: " << *(reinterpret_cast<const uint64_t *>(&data[i])) << std::endl;
                // ============================================
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
            delete[] *b1[i];
            delete[] b0p[i];
            delete[] b1p[i];
            delete[] b0[i];
            delete[] b1[i];
        }
        delete[] sp;
        delete[] b0p;
        delete[] b1p;
        delete[] b0;
        delete[] b1;
        delete[] cm0;
        delete[] cm1;
        delete[] vp;
        delete[] r;
    }
};

}



#endif