class UBFitter{
	private:
		static vector<double>datapoint;
		static UBFitter* currentInstance;
	public:
		UBFitter(){
			//currentInstance = this;
		};
		UBFitter(vector<double> dat){
			SetDatapoint(dat);
		}
		void SetDatapoint(const vector<double>& dat){
//			datapoint.clear();
			datapoint = dat;
		}
		static void DoFit(){
			TMinuit Min(1);
			double arglist[10];
			Min.SetFCN(FCN);
			int ierflg = 0;
			Min.mnparm(0,"p1",0.1,0.01,-1,1,ierflg);
			arglist[0]=500;
			arglist[1]=1;
			Min.mnexcm("MIGRAD",arglist,2,ierflg);
			double param;
			double parE;
			Min.GetParameter(0,param,parE);
			}
		double GetPar(){
			return 1;
		}
};
