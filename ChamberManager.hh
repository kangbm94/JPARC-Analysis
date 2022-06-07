#include "Utils.hh"
#include "ParLimits.hh"
#include "FileManager.hh"
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
double Ztof[2] = {1565.3,1604.15};
double W0tof=12.5;
double Widthtof=80;
double dWdXtof=75;

TString xTof(int i){
	return Form("x0[0]+u0[0]*%f",Ztof[i]);
}
TString yTof(int i){
	return Form("y0[0]+v0[0]*%f",Ztof[i]);
}
enum dets{ BC3,BC4,Sdc1,Sdc2,Sdc3,Sdc4
};
TString Comments[8] = {"### SDC3-X1","### SDC3-X2", "### SDC3-Y1", "### SDC3-Y2","### SDC4-Y1","### SDC4-Y2","### SDC4-X1","###SDC4-X2"};
TString Sdc1Comments[10] = {"### SDC1-V1","### SDC1-V2", "### SDC1-X1", "### SDC1-X2","### SDC1-U1","### SDC1-U2","### SDC2-X1","###SDC2-X2","### SDC2-Y1","###SDC2-Y2"};
TString BcComments[12] = {"### BC3-X1","### BC3-X2", "### BC3-U1", "### BC3-U2","### BC3-V1","### BC3-V2","### BC4-U1","### BC4-U2","### BC4-V1","###BC4-V2","### BC4-X1","###BC4-X2"};


class ChamberManager: public FileManager{
	public:
		ChamberManager(){}
	
		void LoadBcOut(){LoadChain("bcout");}
		void LoadSdcIn(){LoadChain("sdcin");}
		void LoadSdcOut(){LoadChain("sdcout");}


		void LoadDriftParameter(string filename,int layer, double* p);
		TH1* GetTDCHisto(int layer){
			int hn = layer*100+6;
			return (TH1*)GetHistogram(hn);}
		TH1* GetTDCHisto(int layer,int wire);
		TH2* GetDTDLHisto(int layer){
			int hn = layer*100+22;
			return (TH2*)GetHistogram(hn);}
		TH2* GetLeftDTDLHisto(int layer){
			int hn = layer*100+23;
			return (TH2*)GetHistogram(hn);}
		TH2* GetRightDTDLHisto(int layer){
			int hn = layer*100+24;
			return (TH2*)GetHistogram(hn);}
		TH2* GetDTDLLRHisto(int layer){
			int hn = layer*100+19;
			return (TH2*)GetHistogram(hn);}
		TH1* GetDTTrackHisto(int layer){
			int hn = layer*100+12;
			return (TH1*)GetHistogram(hn);}
		TH1* GetDTHisto(int layer){
			int hn = layer*100+30;
			return (TH1*)GetHistogram(hn);}
		TGraphErrors* DLDTGraph(int layer, double MaxDL,double t1, double t2, bool track);
		void FillToFHitDistribution(TString ht, int seg);
		void WriteDriftParameter(int DetID, int WireID, int Type, int Npar, vector<double> par);
};
void ChamberManager:: WriteDriftParameter(int DetID, int WireID, int Type, int Npar, vector<double> par){
	vector<int> ID = {DetID,WireID,Type,Npar};
	WriteParameter(ID,par);
}

void ChamberManager::FillToFHitDistribution(TString ht, int seg){
			TH2D* H = new TH2D(ht,ht,500,-1000,1000,500,-1000,1000);
			int UorD = (seg)%2;
			TString arg = yTof(UorD)+":"+xTof(UorD)+">>"+ht;
			TCut SegCut = Form("TofSeg[0]==%d&&nhTof==1",seg);
	//		TCut SegCut = Form("TofSeg==%d",seg);
			if(seg==0){
				SegCut = Form("TofSeg[0]!=%d&&nhTof==2",seg);
			}
			TCut cut=SegCut;
			double center = (seg-W0tof)*dWdXtof;
			TLine* left = new TLine(center-Widthtof/2,-1000,center-Widthtof/2,1000);
			TLine* right = new TLine(center+Widthtof/2,-1000,center+Widthtof/2,1000);
			left->SetLineColor(kRed);
			left->SetLineWidth(2);
			right->SetLineColor(kRed);
			right->SetLineWidth(2);
			DataChain->Draw(arg,cut,"col");
			left->Draw("same");
			right->Draw("same");
}
void ChamberManager::LoadDriftParameter(string filename, int layer, double* p){
	fstream ParamFile;
	ParamFile.open(filename,fstream::in);
	double dats[10];
	for(int i=0;i<100;i++){
	ReadTSV(ParamFile,dats);
		if(dats[0]==layer){
			for(int j=0;j<6;j++){
				p[j]=dats[4+j];
			}
		continue;
		}
	}
};

TH1* ChamberManager::GetTDCHisto(int layer, int wire){
	int hn = layer*10000+wire;	
	TH1* h = (TH1*)GetHistogram(hn);
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
