#include "FileManager.hh"
static const int bftnseg=160;
static const int schnseg=64;
class EasirocManager: public FileManager{
	public:
		EasirocManager(){};
		TH2* DrawBftSlewing(int UorD, int seg){ int hn = 15000+1000*UorD+seg;return (TH2*)GetHistogram(hn); };
		TH2* DrawSchSlewing(int seg){ int hn = 23000+seg;return (TH2*)GetHistogram(hn); };
		TH1* DrawBftTDC(int UorD, int seg){ int hn = 11000+1000*UorD+seg;return (TH2*)GetHistogram(hn); };
		TH1* DrawSchTDC(int seg){ int hn = 21000+seg;return (TH2*)GetHistogram(hn); };

};
