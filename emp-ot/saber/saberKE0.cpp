#include "poly.h"
#include "poly_mul.h"
#include "api.h"
#include "pack_unpack.h"
#include "cbd.h"
#include "fips202.h"
#include "SABER_params.h"
#include "rng.h"
#include <stdint.h>
#include <iostream>
#include <emp-tool/emp-tool.h>

int main(){
    uint8_t seed_A[SABER_SEEDBYTES];
    uint8_t seed_s[SABER_NOISE_SEEDBYTES];
    uint16_t A[SABER_L][SABER_L][SABER_N];
    uint16_t **s = new uint16_t*[SABER_L];
    uint16_t **b = new uint16_t*[SABER_L];
    uint16_t **sp = new uint16_t*[SABER_L];
    uint16_t **bp = new uint16_t*[SABER_L];
    for (int i = 0; i < SABER_L; i++) {
        s[i] = new uint16_t[SABER_N];
        b[i] = new uint16_t[SABER_N];
        sp[i] = new uint16_t[SABER_N];
        bp[i] = new uint16_t[SABER_N];
        // memset(b[i], 0, SABER_N * sizeof(uint16_t));
        // memset(bp[i], 0, SABER_N * sizeof(uint16_t));
    }
    // memset(b, 0, SABER_L * SABER_N * sizeof(uint16_t));
    // memset(bp, 0, SABER_L * SABER_N * sizeof(uint16_t));
    // print hash of b and bp
    emp::block hash_b = emp::Hash::hash_for_block(b, SABER_L * SABER_N * sizeof(uint16_t));
    emp::block hash_bp = emp::Hash::hash_for_block(bp, SABER_L * SABER_N * sizeof(uint16_t));
    std::cout << "hash_b: " << *reinterpret_cast<uint64_t *>(&hash_b) << std::endl;
    std::cout << "hash_bp: " << *reinterpret_cast<uint64_t *>(&hash_bp) << std::endl;

    uint16_t v[SABER_N];
    uint16_t vp[SABER_N];
    uint16_t cm[SABER_N];
    memset(v, 0, SABER_N * sizeof(uint16_t));
    memset(vp, 0, SABER_N * sizeof(uint16_t));
    // memset(cm, 0, SABER_N * sizeof(uint16_t));
    // the damn bug is from here, v and vp are not initialized

    // gen all vectors
    randombytes(seed_A, SABER_SEEDBYTES);
    shake128(seed_A, SABER_SEEDBYTES, seed_A, SABER_SEEDBYTES);
    GenMatrix(A, seed_A);
    randombytes(seed_s, SABER_NOISE_SEEDBYTES);
    GenSecret(s, seed_s);
    randombytes(seed_s, SABER_NOISE_SEEDBYTES);
    GenSecret(sp, seed_s);

    RoundingMul(A, s, b, 0);
    RoundingMul(A, sp, bp, 1);
    // print numbers of b and bp
    for (int i = 0; i < SABER_L; i++) {
        for (int j = 0; j < SABER_N; j++) {
            std::cout << b[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for (int i = 0; i < SABER_L; i++) {
        for (int j = 0; j < SABER_N; j++) {
            std::cout << bp[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    // print hash of b and bp
    hash_b = emp::Hash::hash_for_block(b, SABER_L * SABER_N * sizeof(uint16_t));
    hash_bp = emp::Hash::hash_for_block(bp, SABER_L * SABER_N * sizeof(uint16_t));
    std::cout << "hash_b: " << *reinterpret_cast<uint64_t *>(&hash_b) << std::endl;
    std::cout << "hash_bp: " << *reinterpret_cast<uint64_t *>(&hash_bp) << std::endl;


    for (int i = 0; i < SABER_L; i++) {
        for (int j = 0; j < SABER_N; j++) {
            s[i][j] = Bits(s[i][j], SABER_EP, SABER_EP);
            sp[i][j] = Bits(sp[i][j], SABER_EP, SABER_EP);
        }
    }
    InnerProd_plush1(b, sp, vp);
    InnerProd_plush1(bp, s, v);

    for (int i = 0; i < SABER_N; i++) {
        cm[i] = Bits(v[i], SABER_EP - 1, SABER_ET);
        v[i] = Bits(v[i], SABER_EP, 1);
        bool pp = Bits(vp[i] - (cm[i] << (SABER_EP - 1 - SABER_ET)) + h2, SABER_EP, 1);
        if (pp != v[i]) {
            std::cout << "Error!" << std::endl;
        }
    }
    return 0;
}