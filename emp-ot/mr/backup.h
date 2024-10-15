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

    void send (const block* data0, const block* data1, int64_t length) override {
        KyberOtRecvPKs* pkBuff = new KyberOtRecvPKs[length];
        KyberOTCtxt* ctxts = new KyberOTCtxt[length];
        KyberOTPtxt ptxt;
        unsigned char *whole_msg = new unsigned char[length * 2 * PKlength * sizeof(unsigned char)];
        // receive the KyberReceiverMessage from the receiver
        io->recv_data(whole_msg, length * 2 * PKlength * sizeof(unsigned char));

        // copy to pkBuff
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < 2; ++j) {
                memcpy(pkBuff[i].keys[j], whole_msg + i * 2 * PKlength + j * PKlength, PKlength * sizeof(unsigned char));
            }
        }
        
        // compute the KyberSenderMessage
        for (int i = 0; i < length; ++i)
		{
			memcpy(ptxt.sot[0], &data0[i], sizeof(block));
			memset(ptxt.sot[0] + sizeof(block), 0, sizeof(ptxt.sot[0]) - sizeof(block));

			memcpy(ptxt.sot[1], &data1[i], sizeof(block));
			memset(ptxt.sot[1] + sizeof(block), 0, sizeof(ptxt.sot[1]) - sizeof(block));

			//get senders message, secret coins and ot strings
			KyberSenderMessage(&ctxts[i], &ptxt, &pkBuff[i]);
		}
        // send the KyberSenderMessage to the receiver
        unsigned char *whole_ctxt = new unsigned char[length * 2 * CTlength * sizeof(unsigned char)];
        
        for (int i = 0; i < length; ++i) {
            memcpy(whole_ctxt + i * 2 * CTlength, ctxts[i].sm[0], CTlength * sizeof(unsigned char));
        }
        io->send_data(whole_ctxt, length * 2 * CTlength * sizeof(unsigned char));
        
        // Free allocated memory
        delete[] whole_msg;
        delete[] whole_ctxt;
        delete[] pkBuff;
        delete[] ctxts;
    }

    void recv (block* data, const bool* b, int64_t length) override {
        KyberOTRecver* ot = new KyberOTRecver[length];
        KyberOtRecvPKs* pkBuff = new KyberOtRecvPKs[length];
        KyberOTCtxt* ctxts = new KyberOTCtxt[length];

        // generate receiver message, KyberReceiverMessage, for all i in [0, length)
		for (int i = 0; i < length; ++i)
		{
			ot[i].b = b[i];
			//get receivers message and secret coins
			KyberReceiverMessage(&ot[i], &pkBuff[i]);
		}
        
        // send the KyberReceiverMessage to the sender
        unsigned char * whole_msg = new unsigned char [length * 2 * PKlength * sizeof(unsigned char)];
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < 2; ++j) {
                memcpy(whole_msg + i * 2 * PKlength + j * PKlength, pkBuff[i].keys[j], PKlength * sizeof(unsigned char));
            }
        }
        io->send_data(whole_msg, length * 2 * PKlength * sizeof(unsigned char));


        // receive the KyberSenderMessage from the sender
        unsigned char *whole_ctxt = new unsigned char[length * 2 * CTlength * sizeof(unsigned char)];
        io->recv_data(whole_ctxt, length * 2 * CTlength * sizeof(unsigned char));
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < 2; ++j) {
                memcpy(ctxts[i].sm[j], whole_ctxt + i * 2 * CTlength + j * CTlength, CTlength * sizeof(unsigned char));
            }
        }
        // for all i in [0, length), decrypt and get the message  (KyberReceiverString)
        for (int i = 0; i < length; ++i) {
            KyberReceiverStrings(&ot[i], &ctxts[i]);
            memcpy(data + i, ot[i].rot, sizeof(block));
        }
        // Free allocated memory
        delete[] whole_msg;
        delete[] whole_ctxt;
        delete[] ot;
        delete[] pkBuff;
        delete[] ctxts;
    }
};
}

#endif// MRKyber_h