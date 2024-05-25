#include "test/test.h"
#include "emp-ot/emp-ot.h"
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
    cout <<"128 COOTs:\t"<<test_ot<OTCO<NetIO>>(co, io, party, 128)<<" ms"<<endl;

    OTNP<NetIO> * np = new OTNP<NetIO>(io); 
    cout <<"128 NPOTs:\t"<<test_ot<OTNP<NetIO>>(np, io, party, 128)<<" ms"<<endl;

    uint8_t* seed_A = new uint8_t[SABER_SEEDBYTES];
    uint16_t* r = new uint16_t[SABER_L * SABER_N];

    FILE *fp1 = fopen("random.txt", "rb");
    fread(r, 1, SABER_L * SABER_N * sizeof(uint16_t), fp1);
    fclose(fp1);
    FILE *fp2 = fopen("seedA.txt", "rb");
    fread(seed_A, 1, SABER_SEEDBYTES, fp2);
    fclose(fp2);

    // NPSaber1<NetIO> *npsaber1 = new NPSaber1<NetIO>(io, seed_A, r);
    // cout <<"128 NP Saber OTs(2 secrets):\t"<<test_ot<NPSaber1<NetIO>>(npsaber1, io, party, 128)<<" ms"<<endl;

    NPSaber2<NetIO> *npsaber2 = new NPSaber2<NetIO>(io, seed_A, r);
    cout <<"128 NP Saber OTs(1 secret):\t"<<test_ot<NPSaber2<NetIO>>(npsaber2, io, party, 128)<<" ms"<<endl;

    SimpleSaber<NetIO> *simplesaber = new SimpleSaber<NetIO>(io, seed_A);
    cout <<"128 Simplest Saber OTs:\t"<<test_ot<SimpleSaber<NetIO>>(simplesaber, io, party, 128)<<" ms"<<endl;

	// cout <<"Passive SABER OT\t"<<double(length)/test_ot<NPSaber1<NetIO>>(npsaber1, io, party, length)*1e6<<" OTps"<<endl;

    // ================== OTE ==================
	// cout <<"Passive SABER COT\t"<<double(length)/test_cot<NPSaber1<NetIO>>(npsaber1, io, party, length)*1e6<<" OTps"<<endl;
	// cout <<"Passive SABER ROT\t"<<double(length)/test_rot<NPSaber1<NetIO>>(npsaber1, io, party, length)*1e6<<" OTps"<<endl;
	// delete npsaber1;
	// npsaber1 = new NPSaber1<NetIO>(io);
	// cout <<"Active SABER OT\t"<<double(length)/test_ot<NPSaber1<NetIO>>(npsaber1, io, party, length)*1e6<<" OTps"<<endl;
	// // cout <<"Active SABER COT\t"<<double(length)/test_cot<NPSaber1<NetIO>>(npsaber1, io, party, length)*1e6<<" OTps"<<endl;
	// // cout <<"Active SABER ROT\t"<<double(length)/test_rot<NPSaber1<NetIO>>(npsaber1, io, party, length)*1e6<<" OTps"<<endl;
	// delete npsaber1;
    // ==========================================    

    delete np;
    // delete npsaber1;
    delete npsaber2;
    delete simplesaber;
    delete io;
    delete co;
    return 0;
}