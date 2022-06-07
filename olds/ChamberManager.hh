#include "Utils.hh"
#include "ParLimits.hh"
static const int BcNW = 64;
static const int Sdc1NW = 64;
static const int Sdc2XNW = 70;
static const int Sdc2YNW = 40;
static const int Sdc3NW = 128;
static const int Sdc4XNW = 96;
static const int Sdc4YNW = 64;

static const double BcDL = 1.5;
static const double Sdc1DL = 3;
static const double Sdc2DL = 5;
static const double Sdc3DL = 4.5;
static const double Sdc4DL = 10;
double Bcdt2 = 45;
double Sdc1dt2 = 80;
double Sdc2dt2 = 120;
double Sdc3dt2 = 110;
double Sdc4dt2 = 260;
int Bcid = 113;


enum dets{ BC3,BC4,Sdc1,Sdc2,Sdc3,Sdc4
};
TString Comments[8] = {"### SDC3-X1","### SDC3-X2", "### SDC3-Y1", "### SDC3-Y2","### SDC4-Y1","### SDC4-Y2","### SDC4-X1","###SDC4-X2"};
TString Sdc1Comments[10] = {"### SDC1-V1","### SDC1-V2", "### SDC1-X1", "### SDC1-X2","### SDC1-U1","### SDC1-U2","### SDC2-X1","###SDC2-X2","### SDC2-Y1","###SDC2-Y2"};
TString BcComments[12] = {"### BC3-X1","### BC3-X2", "### BC3-U1", "### BC3-U2","### BC3-V1","### BC3-V2","### BC4-U1","### BC4-U2","### BC4-V1","###BC4-V2","### BC4-X1","###BC4-X2"};
class ChamberManager{
	private:
		fstream pf;
		TFile* file;
	public:
		ChamberManager(){}
		void LoadFile(TString filename);
		void MakeParameterFile(string title);
		void WriteT0Parameter(int Id, int WireId, double p0, double p1);
		void WriteTParameter( double p0);
		void WriteComment(TString comment);
		void Layer( int layer);
		void Sdc1Layer( int layer);
		void BcLayer( int layer);
		void WriteDriftParameter(int Id, int WireId, int Type, int NParam, double* p);
		void LoadDriftParameter(string filename,int layer, double* p);
		TH1* GetHisto(int nl, int nw);
		TH1* GetHisto(int nl);
		TH1* GetDTHisto(int nl);
		TH2* GetDTDLHisto(int nl);
		TH2* GetDTDLLRHisto(int nl);
		TH1* GetDTTrackHisto(int nl);
		TGraphErrors* DLDTGraph(int layer, double MaxDL,double t1, double t2, bool track);
};
void ChamberManager::LoadFile(TString filename){
	file = new TFile(filename,"READ");
}
void ChamberManager::MakeParameterFile(string title){
	pf.open(title,fstream::out);
}
void ChamberManager::WriteComment(TString comment){
	pf<<comment<<endl;
}
void ChamberManager::Layer(int layer){
	WriteComment(Comments[layer]);
	WriteComment("###");
}
void ChamberManager::Sdc1Layer(int layer){
	WriteComment(Sdc1Comments[layer]);
	WriteComment("###");
}
void ChamberManager::BcLayer(int layer){
	WriteComment(BcComments[layer]);
	WriteComment("###");
}
void ChamberManager::LoadDriftParameter(string filename, int layer, double* p){
	fstream f;
	f.open(filename,fstream::in);
	double dats[10];
	for(int i=0;i<100;i++){
	ReadTSV(f,dats);
		if(dats[0]==layer){
			for(int j=0;j<6;j++){
				p[j]=dats[4+j];
			}
		continue;
		}
	}
};
void ChamberManager::WriteT0Parameter(int Id, int WireId, double p0, double p1){
	pf<<Id<<"\t";
	pf<<WireId<<"\t";
	pf<<p0<<"\t";
	pf<<p1<<endl;
}
void ChamberManager::WriteTParameter(double p0){
	pf<<"######"<<Form(" tmax = %f ",p0)<<"#####"<<endl;
}
void ChamberManager::WriteDriftParameter(int Id, int WireId, int Type, int NParam, double* p){
	pf<<Id<<"\t"	;
	pf<<WireId<<"\t"	;
	pf<<Type<<"\t"	;
	pf<<NParam	;
	for(int i=0;i<NParam;i++){
		pf<<"\t"<<p[i];
	}
	pf<<endl;
}
TH1* ChamberManager::GetHisto(int layer){
	int hn = layer*100+6;
	return (TH1*)file->Get(Form("h%d",hn));
}
TH1* ChamberManager::GetDTHisto(int layer){
	int hn = layer*100+3;
	return (TH1*)file->Get(Form("h%d",hn));
}
TH2* ChamberManager::GetDTDLHisto(int layer){
	int hn = layer*100+22;
	return (TH2*)file->Get(Form("h%d",hn));
}
TH2* ChamberManager::GetDTDLLRHisto(int layer){
	int hn = layer*100+19;
	return (TH2*)file->Get(Form("h%d",hn));
}
TH1* ChamberManager::GetDTTrackHisto(int layer){
	int hn = layer*100+12;
	return (TH1*)file->Get(Form("h%d",hn));
}

TH1* ChamberManager::GetHisto(int layer, int wire){
	int hn = layer*10000+wire;
	TH1* h = (TH1*)file->Get(Form("h%d",hn));
	int peak = h->GetMaximum();
	if(peak<200) h->Rebin(2);
	return h;
}
TGraphErrors* ChamberManager::DLDTGraph(int layer,double MaxDL, double t1,double t2, bool track){
	TH1* hist_dt;
	if(track) hist_dt = GetDTTrackHisto(layer);
	else	hist_dt = GetDTHisto(layer);
	const int nbin = hist_dt->GetNbinsX();
	double bc[nbin];
	double bcI[nbin];
	double cont[nbin];
	double countI[nbin];
	double bce[nbin];
	double cnte[nbin];
	double Integrated=0;
	double bw;
	int count=0;
	for(int i=1;i<nbin+1;++i){
		bc[i-1]=hist_dt->GetBinCenter(i);
		bw=hist_dt->GetBinWidth(i);
		cont[i-1]=hist_dt->GetBinContent(i);
		if(bc[i-1]<t2&&bc[i-1]>t1){
			Integrated+=cont[i-1];
			bcI[count]=bc[i-1];
			bce[count]=bw/sqrt(12);
			countI[count]=Integrated;
			count++;
		}
	}
	cout<<count<<endl;
	countI[count]=(1+1/count);
	for(int i=0;i<count;++i){
		countI[i]=MaxDL*countI[i]/Integrated;
		cnte[i]=(MaxDL*countI[i+1]/Integrated-countI[i])/sqrt(12);
	}
	TString ht = Form("DLDT_Hist%d",layer);
	TGraphErrors* graph = new TGraphErrors(count,bcI, countI,bce,cnte);
	return graph;
}
