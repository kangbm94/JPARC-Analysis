#include "HodoscopeManager.hh"
const int bh1segs=11,bh2segs=8;
int mod=1;//1 for bh1, 2 for bh2,4 for htof
void PHC(){
	cout<<"PHC(int run)"<<endl;
}
void PHC(int run){
	TChain* chain = new TChain("tree");
	TString filename=Filename(run);
	chain->Add(filename);
	HodoscopeManager h;
	int segs[3]={0,bh1segs,bh2segs};
	mod = 1;//BH1
     
	h.MakePHCFile(Form("/home/had/kangbm/HodoManager/PHCFiles/HodoPHCParam_run0%d",run),mod,run);
	h.LoadPHCFile("/home/had/kangbm/HodoManager/PHC_Format");
	for(int i=0;i<segs[mod];i++){
		h.FitTimewalk(chain,0,mod,i+1,0);
	}
	h.BH1_U2D();
	for(int i=0;i<segs[mod];i++){
		h.FitTimewalk(chain,0,mod,i+1,1);
	}
	h.BH1_2_BH2();
	mod=2;//BH2
	for(int i=0;i<segs[mod];i++){
		h.FitTimewalk(chain,0,mod,i+1,0);
	}
	h.BH2_U2D();
	for(int i=0;i<segs[mod];i++){
		h.FitTimewalk(chain,0,mod,i+1,1);
	}
	h.OtherLine();
}
