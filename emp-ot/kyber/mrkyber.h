#ifndef MRKyber_h
#define MRKyber_h

extern "C"
{
#include <test/KyberOT/KyberOT.h>
}

namespace emp
{
template<typename IO>
class MRKyber: public OT<IO>{
    public:
    IO* io;

    MRKyber(IO* io) {
        this->io = io;
    }

    ~MRKyber(){
    }

    void recv (block* data, const bool* b, int64_t length) override {
        auto n = length;
		auto ot = std::vector<KyberOTRecver>{};
		auto pkBuff = std::vector<KyberOtRecvPKs>{};
		auto ctxts = std::vector<KyberOTCtxt>{};
		ot.resize(n);
		pkBuff.resize(n);
		ctxts.resize(n);

		for (uint64_t i = 0; i < n; ++i)
		{
			ot[i].b = b[i];

			//get receivers message and secret coins
			KyberReceiverMessage(&ot[i], &pkBuff[i]);
		}

        // send pkBuff as whole_msg to the sender
        auto whole_msg = new unsigned char [length * 2 * PKlength * sizeof(unsigned char)];
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < 2; ++j) {
                memcpy(whole_msg + i * 2 * PKlength + j * PKlength, pkBuff[i].keys[j], PKlength * sizeof(unsigned char));
            }
        }
        io->send_data(whole_msg, length * 2 * PKlength * sizeof(unsigned char));

        // receive ctxts from the sender
        auto whole_ctxt = new unsigned char[length * 2 * CTlength * sizeof(unsigned char)];
        io->recv_data(whole_ctxt, length * 2 * CTlength * sizeof(unsigned char));
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < 2; ++j) {
                memcpy(ctxts[i].sm[j], whole_ctxt + i * 2 * CTlength + j * CTlength, CTlength * sizeof(unsigned char));
            }
        }

        for (uint64_t i = 0; i < n; ++i)
		{
			KyberReceiverStrings(&ot[i], &ctxts[i]);
			memcpy(data + i, ot[i].rot, sizeof(block));
		}
    }

    void send (const block* data0, const block* data1, int64_t length) override {
        auto pkBuff = std::vector<KyberOtRecvPKs>{};
		auto ctxts = std::vector<KyberOTCtxt>{};
		auto ptxt = KyberOTPtxt{ };
		auto n = length;
		pkBuff.resize(n);
		ctxts.resize(n);

        // receive whole_msg from the receiver
        auto whole_msg = new unsigned char[length * 2 * PKlength * sizeof(unsigned char)];
        io->recv_data(whole_msg, length * 2 * PKlength * sizeof(unsigned char));
        // copy to pkBuff
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < 2; ++j) {
                memcpy(pkBuff[i].keys[j], whole_msg + i * 2 * PKlength + j * PKlength, PKlength * sizeof(unsigned char));
            }
        }

        for (uint64_t i = 0; i < n; ++i)
		{
			memcpy(ptxt.sot[0], &data0[i], sizeof(block));
			memset(ptxt.sot[0] + sizeof(block), 0, sizeof(ptxt.sot[0]) - sizeof(block));

			memcpy(ptxt.sot[1], &data1[i], sizeof(block));
			memset(ptxt.sot[1] + sizeof(block), 0, sizeof(ptxt.sot[1]) - sizeof(block));

			//get senders message, secret coins and ot strings
			KyberSenderMessage(&ctxts[i], &ptxt, &pkBuff[i]);
		}

        // send ctxt as whole_ctxt to the receiver
        auto whole_ctxt = new unsigned char[length * 2 * CTlength * sizeof(unsigned char)];
        for (int i = 0; i < length; ++i) {
            memcpy(whole_ctxt + i * 2 * CTlength, ctxts[i].sm[0], CTlength * sizeof(unsigned char));
            memcpy(whole_ctxt + i * 2 * CTlength + CTlength, ctxts[i].sm[1], CTlength * sizeof(unsigned char));
        }
        io->send_data(whole_ctxt, length * 2 * CTlength * sizeof(unsigned char));

        // Free allocated memory
        delete[] whole_msg;
        delete[] whole_ctxt;
    }
};
}

#endif// MRKyber_h