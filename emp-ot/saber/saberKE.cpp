#include "poly.h"
#include "poly_mul.h"
#include "api.h"
#include "pack_unpack.h"
#include "cbd.h"
#include "fips202.h"
#include "SABER_params.h"
#include "rng.h"
#include <emp-tool/emp-tool.h>
#include <iostream>

int main(){
    // simulate the sender
    uint8_t seed_A[SABER_SEEDBYTES];
    uint8_t seed_s[SABER_NOISE_SEEDBYTES];
    uint16_t A[SABER_L][SABER_L][SABER_N];
    // uint16_t s0[SABER_L][SABER_N];
    // uint16_t s1[SABER_L][SABER_N];
    // uint16_t b0[SABER_L][SABER_N];
    // uint16_t b1[SABER_L][SABER_N];
    uint16_t **s0 = new uint16_t*[SABER_L];
    uint16_t **s1 = new uint16_t*[SABER_L];
    uint16_t **b0 = new uint16_t*[SABER_L];
    uint16_t **b1 = new uint16_t*[SABER_L];
    for (int i = 0; i < SABER_L; i++) {
        s0[i] = new uint16_t[SABER_N];
        s1[i] = new uint16_t[SABER_N];
        b0[i] = new uint16_t[SABER_N];
        b1[i] = new uint16_t[SABER_N];
    }
    uint16_t r[SABER_L * SABER_N];

    // gen r
    randombytes(reinterpret_cast<unsigned char *>(r), SABER_N * SABER_L * sizeof(uint16_t));

    // Generate A
    randombytes(seed_A, SABER_SEEDBYTES);
    GenMatrix(A, seed_A);
    randombytes(seed_s, SABER_NOISE_SEEDBYTES);
    GenSecret(s0, seed_s);
    randombytes(seed_s, SABER_NOISE_SEEDBYTES);
    GenSecret(s1, seed_s);

    // simulate the receiver
    // use the same seed_A and matrix A
    // uint16_t sp[SABER_L][SABER_N];
    // uint16_t b1p[SABER_L][SABER_N];
    // uint16_t b0p[SABER_L][SABER_N];
    uint16_t **sp = new uint16_t*[SABER_L];
    uint16_t **b1p = new uint16_t*[SABER_L];
    uint16_t **b0p = new uint16_t*[SABER_L];
    for (int i = 0; i < SABER_L; i++) {
        sp[i] = new uint16_t[SABER_N];
        b1p[i] = new uint16_t[SABER_N];
        b0p[i] = new uint16_t[SABER_N];
    }
    randombytes(seed_s, SABER_NOISE_SEEDBYTES);
    GenSecret(sp, seed_s);
    bool x = 1;
    if (x) {
        RoundingMul(A, sp, b1p, 1);
        // b0p = r - b1p;
        for (int i = 0; i < SABER_L; i++) {
            for (int j = 0; j < SABER_N; j++) {
                b0p[i][j] = r[i * SABER_N + j] - b1p[i][j];
            }
        }
    } else {
        RoundingMul(A, sp, b0p, 1);
        // b1p = r - b0p;
        for (int i = 0; i < SABER_L; i++) {
            for (int j = 0; j < SABER_N; j++) {
                b1p[i][j] = r[i * SABER_N + j] - b0p[i][j];
            }
        }
    }

    // simulate the sender
    RoundingMul(A, s0, b0, 0);
    RoundingMul(A, s1, b1, 0);
    uint16_t cm0[SABER_N];
    uint16_t cm1[SABER_N];
    uint16_t v0[SABER_N];
    uint16_t v1[SABER_N];
    for(int j = 0; j < SABER_L; ++j) {
        for(int k = 0; k < SABER_N; ++k) {
            s0[j][k] = Bits(s0[j][k], SABER_EP, SABER_EP);
            s1[j][k] = Bits(s1[j][k], SABER_EP, SABER_EP);
        }
    }
    InnerProd_plush1(b0p, s0, v0);
    InnerProd_plush1(b1p, s1, v1);
    emp::block m[2];
    for (int j = 0; j < SABER_N; ++j) {
        cm0[j] = Bits(v0[j], SABER_EP - 1, SABER_ET);
        cm1[j] = Bits(v1[j], SABER_EP - 1, SABER_ET);
    }
    for (int j = 0; j < SABER_N; ++j) {
        v0[j] = Bits(v0[j], SABER_EP, 1);
        v1[j] = Bits(v1[j], SABER_EP, 1);
    }
    emp::block hash_v0 = emp::Hash::hash_for_block(v0, SABER_N * 2);
    emp::block hash_v1 = emp::Hash::hash_for_block(v1, SABER_N * 2);
    std::cout << "hash_v0: " << *reinterpret_cast<const uint64_t *>(v0) << std::endl;
    std::cout << "hash_v1: " << *reinterpret_cast<const uint64_t *>(v1) << std::endl;

    // simulate the receiver
    uint16_t vp[SABER_N];
    for (int j = 0; j < SABER_L; ++j) {
        for (int k = 0; k < SABER_N; ++k) { 
            sp[j][k] = Bits(sp[j][k], SABER_EP, SABER_EP);
        }
    }
    if (x) {
        InnerProd_plush1(b1, sp, vp);
        for (int j = 0; j < SABER_N; ++j) {
            cm1[j] = Bits(vp[j] - (cm1[j] << (SABER_EP - 1 - SABER_ET)) + h2, SABER_EP, 1);
        }
        emp::block hash_cm1 = emp::Hash::hash_for_block(cm1, SABER_N * 2);
        std::cout << "hash_cm1: " << *reinterpret_cast<const uint64_t *>(cm1) << std::endl;
        // check v1 and cm1
        for (int j = 0; j < SABER_N; ++j) {
            if (v1[j] != cm1[j]) {
                std::cout << j << std::endl;
            }
        }
    } else {
        InnerProd_plush1(b0, sp, vp);
        for (int j = 0; j < SABER_N; ++j) {
            cm0[j] = Bits(vp[j] - (cm0[j] << (SABER_EP - 1 - SABER_ET)) + h2, SABER_EP, 1);
        }
        emp::block hash_cm0 = emp::Hash::hash_for_block(cm0, SABER_N * 2);
        std::cout << "hash_cm0: " << *reinterpret_cast<const uint64_t *>(cm0) << std::endl;
        // check v0 and cm0
        for (int j = 0; j < SABER_N; ++j) {
            if (v0[j] != cm0[j]) {
                std::cout << j << std::endl;
            }
        }
    }
    return 0;
}