#ifndef Math_hh
#define Math_hh 1
double XiParam[2] = {-0.401,-2.1};
double AXi = XiParam[0];
double BXi = sqrt(1 - AXi*AXi) * sin(XiParam[1]*acos(-1)/180);
double CXi = sqrt(1 - AXi*AXi) * cos(XiParam[1]*acos(-1)/180);

double LdParam[2] = {0.732,-6.5};
double ALd = LdParam[0];
double BLd = sqrt(1 - ALd*ALd) * cos(LdParam[1]*acos(-1)/180);
double CLd = sqrt(1 - ALd*ALd) * sin(LdParam[1]*acos(-1)/180);
double mP = .938272;

TVector3 NormalCross(TVector3 A, TVector3 B){
	auto C = A.Cross(B);
	C = C*(1./C.Mag());
	return C;
}
double NormalDot(TVector3 A, TVector3 B){
	auto nA = A.Mag();
	auto nB = B.Mag();
	return A*B / nA / nB;
}

double MeV = 1000;
TLorentzVector ToHelix(TLorentzVector B){
		double x = -B.X();
		double y = B.Z();
		double z = B.Y();
		double t = B.T();
		return TLorentzVector(x,y,z,t);
}
TVector3 ToHelix(TVector3 B){
		double x = -B.X();
		double y = B.Z();
		double z = B.Y();
		return TVector3(x,y,z);
}
TVector3 ToGlobal(TVector3 B){
		double x = -B.X();
		double y = B.Z();
		double z = B.Y();
		return TVector3(x,y,z);
}
TMatrixD MinorMatrix(TMatrixD A,int row, int col){
	int Arow = A.GetNrows();
	int Acol = A.GetNcols();
	if(row < 0 or row > Arow -1 or col > Acol -1 or col < 0){
		cout<<"Failed to get Minor:";
		A.Print();
		return A; 
	}
	TMatrixD B(Arow-1,Acol-1);
	int Brow = 0;
	for(int ir = 0;ir < Arow; ++ir){
		if(ir == row) continue;
		int Bcol = 0;
		for(int ic = 0;ic < Acol; ++ic){
			if(ic == col) continue;
			B(Brow,Bcol) = A(ir,ic);
			Bcol++;
		}
		Brow++;
	}
	return B;
}
void CovThPh(double* res, double* par, double* cov){
	//* Theta and phi in global coordinate
	// will be respresented in terms of helix coordinate:
	// (x,y,z) -> (z,-x,y) where
	// x = r sin(Th)cos(Ph) = -r sin(th)sin(ph)
	// y = r sin(Th)sin(Ph) = r cos(th)
	// z = r cos(Th) = rsin(th)cos(ph)
	// where th = acos( y/r ) = acos[sin(Th)sin(Ph)]
	// and ph = atan2( -x/z ) = atan2[ sin(Th)cos(Ph) / cos(Th) ] = atan2[ tan(Th)cos(Ph)];
	// Then, Th = acos( z/r ) = acos[ sin(th)cos(ph)]
	// and Ph = atan2( y/x) = atan2[cos(th)/sin(th)/sin(ph)] = atan2[cot(th)/sin(ph)]
	// Then the covariance will be
	// V(Th,Ph) = J^T V(th,ph) J  where J = (dTh/dth , dTh/dph)
	// 																			(dPh/dth , dPh/dph)
	//and we have dTh/dth = - cos(th)sin(ph)/sqrt(1-sin^2(th)sin^2(ph));
	//						dTh/dph = - sin(th)cos(ph)/sqrt(1-sin^2(th)sin^2(ph));
	//						dPh/dth = cos(ph)/(sin^2 ph + cot^2 th);
	//						dPh/dph = cot(th)/(sin^2 ph + cot^2 th);
	double sth=res[0];
	double sph=res[1];
	
	double th=par[0];
	double ph=par[1];
	double dTdt = - cos(th)*sin(ph) / sqrt(1 - pow(sin(th)*sin(ph),2));
	double dTdp = - sin(th)*cos(ph) / sqrt(1 - pow(sin(th)*sin(ph),2));
	double dPdt = cos(ph) / (sin(ph)*sin(ph)+1./tan(th)/tan(th));
	double dPdp = 1./ (sin(ph)*sin(ph)+1./tan(ph)/tan(ph)) / tan(th)/tan(th);

	double VJ11 = sth*dTdt;
	double VJ12 = sth*dTdp;
	double VJ21 = sph*dPdt;
	double VJ22 = sph*dPdp;

	double VTVJ11 = dTdt*VJ11 + dPdt*VJ21;
	double VTVJ12 = dTdt*VJ12 + dPdt*VJ22;
	double VTVJ21 = dTdp*VJ11 + dPdp*VJ21;
	double VTVJ22 = dTdp*VJ12 + dPdp*VJ22;

	cov[0] = VTVJ11;
	cov[1] = VTVJ12;
	cov[2] = VTVJ21;
	cov[3] = VTVJ22;
}
void WeightedFill(TH1D* h,double dat, double w){
	int bin = h->FindBin(dat);
	h->AddBinContent(bin,w);
}
int GetBin(int nbin, double lb,double hb, double cont){
	double bw = (hb - lb)/nbin;
	double range = cont - lb;
	int bin = floor(range / bw);
	if (bin<0) return -1;
	if (bin >= nbin) return nbin;
	else return bin;
}
double GetBinCenter(int nbin, double lb, double hb, int bin){
	double bw = (hb - lb)/nbin;
	double offset = lb + bw/2; 
	if(bin < 0) return lb - bw/2;
	if(bin >= nbin) return hb + bw/2;
	return offset + bw * bin;	
}




#endif
