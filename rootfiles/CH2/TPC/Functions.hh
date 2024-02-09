
double dLInner = 3;
double dLOuter = 20;
double dLBoundary = 15;
int MaxLd = 15;
int LdBin(double len){
	if(len < dLBoundary){
		return len / dLInner;
	}
	else{
		return 6+ (len-30) / dLOuter;
	}
}
double LdLen(int bin){
	int BoundaryBin = LdBin(dLBoundary+1e-9);
	if(bin <= BoundaryBin){
		return dLInner * bin;
	}
	else{
		return (bin - BoundaryBin)* dLOuter + dLBoundary;
	}
}
double PPiMax = 0.23;
double PPiMin = 0.02;
int MaxBinPPi = 21;//MeV
int PPiBin(double mom){
	if(mom < PPiMin) return 0;
	if(mom > PPiMax) return MaxBinPPi+1;
	double dPPi = (PPiMax - PPiMin)/MaxBinPPi;
	int bin = (mom - PPiMin) / dPPi + 1;
	return bin;
}
double PtPiMom(int bin){
	double dPPi = (PPiMax - PPiMin)/MaxBinPPi;
	if(bin == 0) return PPiMin - dPPi/2;
	if(bin > MaxBinPPi) return PPiMax+ dPPi/2;
	return dPPi * bin + PPiMin - dPPi/2;
}
double PPMax = 1.;
double PPMin = 0.2;
int MaxBinPP = 16;//MeV
int PPBin(double mom){
	if(mom < PPMin) return 0;
	if(mom > PPMax) return MaxBinPP+1;
	double dPP = (PPMax - PPMin)/MaxBinPP;
	int bin = (mom - PPMin) / dPP + 1;
	return bin;
}
double PtPMom(int bin){
	double dPP = (PPMax - PPMin)/MaxBinPP;
	if(bin == 0) return PPMin - dPP/2;
	if(bin > MaxBinPP) return PPMax+ dPP/2;
	return dPP * bin + PPMin - dPP/2;
}
double PxPiMax = 0.18;
double PxPiMin = 0.;
int MaxBinPxPi = 18; 
int PxPiBin(double mom){
	double dPxPi = (PxPiMax-PxPiMin)/MaxBinPxPi;
	if(mom < PxPiMin) return 0;
	if(mom > PxPiMax) return MaxBinPxPi+1;
	int bin = (mom - PxPiMin) / dPxPi + 1;
	return bin;
}
double PxPiMom(int bin){
	double dPxPi = (PxPiMax-PxPiMin)/MaxBinPxPi;
	if(bin == 0) return PxPiMin - dPxPi/2;
	if(bin > MaxBinPxPi) return PxPiMax+ dPxPi/2;
	return dPxPi * bin + PxPiMin - dPxPi/2;
}
double ResP(double* res,double* par){
/*
 For a Multivariable function F of Ms, F(M), 
 the Covariance of F, V(F) should satisfy
 V(F) = J^T V(M) J, where J = dFdM is a Jacobian of F respect to M.
	Here, We set F = P = abs(Pt*1/sin(th));
	and obtain the Variance of P from Pt,Th,Ph.
 */
	double VPt = res[0]*res[0];
	double VTh = res[1]*res[1];
	double VPh = res[2]*res[2];

	double Pt = par[0];
	double Th = par[1];
	double Ph = par[2];
	double k = cos(Th);
	double dPdPt = 1./abs(sin(Th));
	double dPdTh = Pt*cos(Th)/sin(Th)/sin(Th);
	double dPdPh = 0;
	/*	
	double dPt2 = j*j/(1-j*j); //Square of Derivatives of Py respect to Pt;
	double dTh2 = Pt*Pt*pow((1-2*j*j),2)/pow(1-j*j,3)*pow(cos(Th)*sin(Ph),2);//Square of Derivatives of Py respect to Th;
	double dPh2 = Pt*Pt*pow((1-2*j*j),2)/pow(1-j*j,3)*pow(sin(Th)*cos(Th),2);//Square of Derivatives of Py respect to Th;
*/
//	double VPy = VPt * dPt2 + VTh * dTh2 + VPh * dPh2;
//	double VP = VPy + VPt;
	double dPt2  = dPdPt*dPdPt;
	double dTh2  = dPdTh*dPdTh;
	double dPh2  = dPdPh*dPdPh;
	double VP = VPt * dPt2 + VTh * dTh2 + VPh * dPh2;
	return sqrt(VP);
}
void GetPResolution(double* input,int* inputI,vector<vector<double>>Rarray,double* scale,double* output){
	double px = input[0];
	double py = input[1];
	double pz = input[2];
	TVector3 HV(-px,pz,py);
	double HTh = HV.Theta();
	double HPh = HV.Phi();
	int nh = inputI[0];
	int MinHit = inputI[1];
	int MaxHit = inputI[2];
	double pt = hypot(px,pz);
	double ps = scale[2];
	double ths = scale[3];
	double phs = scale[4];
	auto RPt = Rarray.at(0); 
	auto RTh = Rarray.at(1); 
	auto RPh = Rarray.at(2); 
	if(nh > MaxHit)nh=MaxHit;
	if(nh < MinHit)nh=MinHit;

	int pbin = PPBin(pt);
	if( pbin > MaxBinPP) pbin = MaxBinPP+1;
	double ResPt = RPt.at(nh-MinHit);
	double ResTh = RTh.at(pbin)*ths;
	double ResPh = RPh.at(pbin)*phs;
	double resP[3] = {pt*ResPt,ResTh,ResPh};
	double parP[3] = {pt,HTh,HPh};
	double ResPP = ResP(resP,parP)*ps;
	output[0] = ResPP;
	output[1] = ResTh;
	output[2] = ResPh;
}
void GetPiResolution(double* input,int* inputI,vector<vector<double>>Rarray,double* scale,double* output){
	double px = input[0];
	double py = input[1];
	double pz = input[2];
	TVector3 HV(-px,pz,py);
	double HTh = HV.Theta();
	double HPh = HV.Phi();
	int nh = inputI[0];
	int MinHit = inputI[1];
	int MaxHit = inputI[2];
	double pt = hypot(px,pz);
	double ps = scale[5];
	double ths = scale[6];
	double phs = scale[7];
	auto RPt = Rarray.at(0); 
	auto RTh = Rarray.at(1); 
	auto RPh = Rarray.at(2); 
	if(nh > MaxHit)nh=MaxHit;
	if(nh < MinHit)nh=MinHit;

	int pbin = PPiBin(pt);
	if( pbin > MaxBinPPi) pbin = MaxBinPPi+1;
	double ResPt = RPt.at(nh-MinHit);
	double ResTh = RTh.at(pbin)*ths;
	double ResPh = RPh.at(pbin)*phs;
	double resP[3] = {pt*ResPt,ResTh,ResPh};
	double parP[3] = {pt,HTh,HPh};
	double ResPP = ResP(resP,parP)*ps;
	output[0] = ResPP;
	output[1] = ResTh;
	output[2] = ResPh;
}

TMatrixD MergedMatrix(TMatrixD A, TMatrixD B){
	int rA = A.GetNrows();
	int cA = A.GetNcols();
	int rB = B.GetNrows();
	int cB = B.GetNcols();
	int rC = rA + rB;
	int cC = cA + cB;
	TMatrixD C(rC,cC);
	TMatrixDSub(C, 0, rA -1, 0, cA -1) = A;
	TMatrixDSub(C, rA, rC -1, cA, cC -1) = B;
	return C;
}
TMatrixD MergedOffDiagonalMatrix(TMatrixD A, TMatrixD B){
	auto C = MergedMatrix(A,B);
	for(int i=0;i<C.GetNcols();++i){
		C(i,i)=0;
	}
	return C;
}
