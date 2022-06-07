#include "DstChain.hh"
#include "/Users/MIN/ROOTSharedLibs/LinesInSpace.hh"
void LoadGeometry(){
}

void DstAnal(){
	cout<<"DstAnal(int run1, int run2)"<<endl;
	cout<<"DstAnal(int run)"<<endl;
	LoadGeometry();
}
int runlist[100]={0};
void DstAnal(int run1, int run2){
	TString filename = "DstKKAna0";
	TString root = ".root";
	int run=0;
	for(int i=0;i<run2-run1+1;i++){
		runlist[i]=run1+i;//
		run=runlist[i];
		kkchain->Add(filename+to_string(run)+root);
	}
	SetBranchAddress(kkchain);
	cout<<"BranchAddressSet"<<endl;
	int ent = kkchain->GetEntries();
	cout<<ent<<endl;
	double dist[10][10];
	double thta[10][10];
	double sim[10][10];
	int SdcInTrack,SdcOutTrack;
	double SdcInx0,SdcIny0,SdcInu0,SdcInv0,SdcInChi2;
	double SdcOutx0,SdcOuty0,SdcOutu0,SdcOutv0,SdcOutChi2;
	double distance,theta,similarity;
	double SdcInChiCut=5,SdcOutChiCut=20;
	//	int nev;
	TFile* file = new TFile("Distance0"+to_string(run)+root,"RECREATE");
	TTree* tree = new TTree("tree","tree");
	tree->Branch("evnum",&evnum,"evnum/I");
	tree->Branch("SdcInTrack",&SdcInTrack,"SdcInTrack/I");
	tree->Branch("SdcInx0",&SdcInx0,"SdcInx0/D");	
	tree->Branch("SdcIny0",&SdcIny0,"SdcIny0/D");	
	tree->Branch("SdcInu0",&SdcInu0,"SdcInu0/D");	
	tree->Branch("SdcInv0",&SdcInv0,"SdcInv0/D");	
	tree->Branch("SdcInChi2",&SdcInChi2,"SdcInChi2/D");	
	tree->Branch("SdcOutTrack",&SdcOutTrack,"SdcOutTrack/I");
	tree->Branch("SdcOutx0",&SdcOutx0,"SdcOutx0/D");	
	tree->Branch("SdcOuty0",&SdcOuty0,"SdcOuty0/D");	
	tree->Branch("SdcOutu0",&SdcOutu0,"SdcOutu0/D");	
	tree->Branch("SdcOutv0",&SdcOutv0,"SdcOutv0/D");	
	tree->Branch("distance",&distance,"distancce/D");	
	tree->Branch("SdcOutChi2",&SdcOutChi2,"SdcOutChi2/D");	
	tree->Branch("theta",&theta,"theta/D");	
	tree->Branch("similarity",&similarity,"similarity/D");	
	cout<<"Output File Loaded"<<endl;
	cout<<"ChainEntryGot"<<endl;
	//	ent = 100;
	double zin = -500,zout=1000;
	for(int i=0;i<ent;i++){
		Indicator(i,ent);
		kkchain->GetEntry(i);
		for(int j=0;j<10;j++){
			for(int k=0;k<10;k++){
				//				dist[j][k]=9e6;
				thta[j][k]=4;
			}
		}

		for(int j=0;j<ntSdcIn;j++){
			Line TrIn=Line(x0SdcIn[j],y0SdcIn[j],u0SdcIn[j],v0SdcIn[j]);
			for(int k=0;k<ntSdcOut;k++){
				Line TrOut=Line(x0SdcOut[k],y0SdcOut[k],u0SdcOut[k],v0SdcOut[k]);
				dist[j][k] = TrOut.Distance(TrIn);
				thta[j][k] = TrOut.Angle(TrIn);
				sim[j][k] = TrOut.Similarity(TrIn,zin,zout);
			}
		}
		similarity = sim[0][0];
		for(int j=0;j<ntSdcIn;j++){
			for(int k=0;k<ntSdcOut;k++){
				similarity = Min(similarity,sim[j][k]);
			}
		}
		if(similarity<0){
			SdcInx0=-9999;
			SdcIny0=-9999;
			SdcInu0=-9999;
			SdcInv0=-9999;
			SdcOutx0=-9999;
			SdcOuty0=-9999;
			SdcOutu0=-9999;
			SdcOutv0=-9999;
			SdcInTrack=-9999;
			SdcOutTrack=-9999;
			distance=-9999;
			SdcInChi2=-9999;
			SdcOutChi2=-9999;
			theta=-9999;
			similarity=-9999;
		}
		else{
			for(int j=0;j<ntSdcIn;j++){
				for(int k=0;k<ntSdcOut;k++){
					if(similarity==sim[j][k]){
						SdcInTrack = j;SdcOutTrack=k;
						SdcInx0=x0SdcIn[j];
						SdcIny0=y0SdcIn[j];
						SdcInu0=u0SdcIn[j];
						SdcInv0=v0SdcIn[j];
						SdcInChi2=chisqrSdcIn[j];
						SdcOutx0=x0SdcOut[k];
						SdcOuty0=y0SdcOut[k];
						SdcOutu0=u0SdcOut[k];
						SdcOutv0=v0SdcOut[k];
						SdcOutChi2=chisqrSdcOut[j];
						distance = dist[j][k];
						theta = thta[j][k];
						similarity = sim[j][k];
						continue;
					}
				}
			}
		}
		tree->Fill();
	}
	file->Write();
}
void DstAnal(int run){
	DstAnal(run,run);
}
