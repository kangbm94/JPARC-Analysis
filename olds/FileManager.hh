class FileManager{
	private:
		TFile* DataFile;

	protected:
		fstream ParameterFile;
		fstream ParameterTemplate;
		TChain* DataChain;

	public:
		FileManager(){}
		void LoadFile(TString FileName){DataFile = new TFile(FileName,"READ");}
		void LoadChain(TString ChainName){DataChain = (TChain*)DataFile->Get(ChainName);}
		TObject* GetHistogram(int HistNumber);
		void Draw(TString Argument, TCut Cut, TString Options = "");
		void MakeParameterFile(TString FileName){ParameterFile.open(FileName,fstream::out);}
		void LoadParameterFile(TString FileName){ParameterTemplate.open(FileName,fstream::in);}
		void WriteParameter(vector<int> ID,vector<double> Param);
		void WriteComment(TString Comment){ParameterFile<<Comment<<endl;}
};

TObject* FileManager::GetHistogram(int HistNumber){
	TString HistName = Form("h%d",HistNumber);
	return (TObject*)DataFile->Get(HistName);
}

void FileManager::Draw(TString Argument, TCut Cut, TString Options = ""){
	DataChain->Draw(Argument, Cut, Options);
}

void FileManager::WriteParameter(vector<int> ID,vector<double> Param){
	for(int i=0;i<ID.size();++i){
		ParameterFile<<ID[i]<<"\t";
	}
	for(int i=0;i<Param.size()-1;++i){
		ParameterFile<<Param[i]<<"\t";
	}
	ParameterFile<<Param[Param.size()-1]<<endl;
}
