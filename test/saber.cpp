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

    SaberOT<NetIO> *saberot = new SaberOT<NetIO>(io);
    cout <<"128 NP Saber OTs:\t"<<test_ot<SaberOT<NetIO>>(saberot, io, party, 128)<<" ms"<<endl;

    SimpleSaber<NetIO> *simplesaber = new SimpleSaber<NetIO>(io);
    cout <<"128 Simplest Saber OTs:\t"<<test_ot<SimpleSaber<NetIO>>(simplesaber, io, party, 128)<<" ms"<<endl;

	// cout <<"Passive SABER OT\t"<<double(length)/test_ot<SaberOT<NetIO>>(saberot, io, party, length)*1e6<<" OTps"<<endl;

    // ================== OTE ==================
	// cout <<"Passive SABER COT\t"<<double(length)/test_cot<SaberOT<NetIO>>(saberot, io, party, length)*1e6<<" OTps"<<endl;
	// cout <<"Passive SABER ROT\t"<<double(length)/test_rot<SaberOT<NetIO>>(saberot, io, party, length)*1e6<<" OTps"<<endl;
	// delete saberot;
	// saberot = new SaberOT<NetIO>(io);
	// cout <<"Active SABER OT\t"<<double(length)/test_ot<SaberOT<NetIO>>(saberot, io, party, length)*1e6<<" OTps"<<endl;
	// // cout <<"Active SABER COT\t"<<double(length)/test_cot<SaberOT<NetIO>>(saberot, io, party, length)*1e6<<" OTps"<<endl;
	// // cout <<"Active SABER ROT\t"<<double(length)/test_rot<SaberOT<NetIO>>(saberot, io, party, length)*1e6<<" OTps"<<endl;
	// delete saberot;
    // ==========================================    

    delete np;
    delete saberot;
    delete simplesaber;
    delete io;
    delete co;
    return 0;
}