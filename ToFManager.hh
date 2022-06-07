#include "Utils.hh"
#include "FileManager.hh"
#include "PhysicalConstants.hh"
static const int tofnseg = 24;
TF1* f_gaus = new TF1("f_gaus","gaus");
TF1* f_landau = new TF1("f_landau","landau");
TF1* slew_func= new TF1("slew_func",SlewFunc,0,5,3);
TCut VertexCut(vector <double> Origin, vector <double> VertexSize){
	TCut CutX = Form("abs(vtx-%f)<%f",Origin[0],VertexSize[0]);
	TCut CutY = Form("abs(vty-%f)<%f",Origin[1],VertexSize[1]);
	TCut CutZ = Form("abs(vtx-%f)<%f",Origin[2],VertexSize[2]);
	return CutX&&CutY&&CutZ;
}
TF1* f_gauss[6];



class ToFManager:public FileManager{
	private:
		TFile* hist_file;
	public:
		ToFManager(){};
		void SaveHisto(TString filename);
		void LoadKHodo(){ LoadChain("khodo");}
		virtual void LoadKK(){ LoadChain("kk");}
		TH1* GetHisto(int seg, int UD);
		TH1* GetADCHisto(int seg, int UD, int particle){int hn = 33001+100*seg+particle+3*UD;return (TH1*)GetHistogram(hn);}
		TH1* GetPEDHisto(int seg, int UD){int hn = 30001+100*seg+UD;return (TH1*)GetHistogram(hn);}
		TH1* GetStofHisto(int seg, int particle){
			int hn = 10000+100*seg+particle;
			return (TH1*)GetHistogram(hn);}
		TH2D* GetPHCHisto(int seg, int UD,int particle);
		TH2D* GetPHCedHisto(int seg, int UD,int particle);
		void FitTimewalk(int seg, int UD,int particle,bool mod);
		void LoadOldVTOF(TString filename,double* oldVTOF);
		void LoadOldCableOffset(TString filename,double* oldOffset);
};
void ToFManager::LoadOldVTOF(TString filename,double* oldVTOF){
	LoadParameterFile(filename);
	double data[7];
	for(int i=0;i<tofnseg;++i){
		ReadTSV(OldParameterFile,data);
		oldVTOF[i]=data[5];
	}
}
void ToFManager::LoadOldCableOffset(TString filename,double* oldOffset){
	LoadParameterFile(filename);
	double data[7];
	for(int i=0;i<8;++i){
		ReadTSV(OldParameterFile,data);
		oldOffset[i]=data[5];
	}
}
void ToFManager::SaveHisto(TString filename){
	hist_file = new TFile(filename,"recreate");
}
TH1* ToFManager::GetHisto(int seg, int UD){// 0 -> Up, 1 -> Down
	int hn = 6*10000+100*seg + 3+UD;
	TH1* h =(TH1*)GetHistogram(hn);
	int peak = h->GetBinContent(h->GetMaximumBin());
	int peak_ref = 20;
	int nbin = 200/peak + 1;
	if(nbin> 1)	h->Rebin(nbin);
	return h;
}
TH2D* ToFManager::GetPHCHisto(int seg, int UD,int particle){
	int hn = 4*10000+100*seg +10*UD+ particle+1;
	TH2D* h =(TH2D*)GetHistogram(hn);
	return h;
};
TH2D* ToFManager::GetPHCedHisto(int seg, int UD,int particle){
	int hn = 5*10000+100*seg +10*UD+ particle+1;
	TH2D* h =(TH2D*)GetHistogram(hn);
	return h;
};

