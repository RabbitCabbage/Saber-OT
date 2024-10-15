#include "test/test.h"
#include "emp-ot/emp-ot.h"
#include <iostream>
using namespace std;

const static int threads = 1;

int main(int argc, char** argv) {
    int length, port, party;

    if (argc <= 3)
        length = (1<<20) + 101;
    else
        length = (1<<atoi(argv[3])) + 101;

    parse_party_and_port(argv, &party, &port);
    NetIO * io = new NetIO(party==ALICE ? nullptr:"127.0.0.1", port);

    OTCO<NetIO> * co = new OTCO<NetIO>(io);
    double co_total_time = 0;
    for (int i = 0; i < 1000; ++i) {
        co_total_time += test_ot<OTCO<NetIO>>(co, io, party, 128);
    }
    cout << "128 COOTs (average over 1000 runs):\t" << (co_total_time / 1000) << " ms" << endl;

    OTNP<NetIO> * np = new OTNP<NetIO>(io);
    double np_total_time = 0;
    for (int i = 0; i < 1000; ++i) {
        np_total_time += test_ot<OTNP<NetIO>>(np, io, party, 128);
    }
    cout << "128 NPOTs (average over 1000 runs):\t" << (np_total_time / 1000) << " ms" << endl;

    MRKyber<NetIO> *mrkyber = new MRKyber<NetIO>(io);
    double mrkyber_total_time = 0;
    for (int i = 0; i < 1000; ++i) {
        mrkyber_total_time += test_ot<MRKyber<NetIO>>(mrkyber, io, party, 128);
    }
    cout << "128 MRKyber OTs (average over 1000 runs):\t" << (mrkyber_total_time / 1000) << " ms" << endl;
    delete mrkyber;
    delete np;
    delete co;
    delete io;
    return 0;
}