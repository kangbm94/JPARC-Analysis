#include "DstBranches.hh"
TChain* kkchain = new TChain("kk");
void SetBranchAddress(TChain* chain){	
	chain->SetBranchAddress("evnum",&evnum);
	chain->SetBranchAddress("nlBcOut",&nlBcOut);
	chain->SetBranchAddress("ntBcOut",&ntBcOut);
	chain->SetBranchAddress("nhBcOut",&nhBcOut);
	chain->SetBranchAddress("chisqrBcOut",&chisqrBcOut);
	chain->SetBranchAddress("x0BcOut",x0BcOut);
	chain->SetBranchAddress("y0BcOut",y0BcOut);
	chain->SetBranchAddress("u0BcOut",u0BcOut);
	chain->SetBranchAddress("v0BcOut",v0BcOut);
	chain->SetBranchAddress("xtgtBcOut",xtgtBcOut);
	chain->SetBranchAddress("ytgtBcOut",ytgtBcOut);
//	chain->SetBranchAddress("xbh2gtBcOut",xbh2gtBcOut);
//	chain->SetBranchAddress("ybh2gtBcOut",ybh2gtBcOut);

	chain->SetBranchAddress("nlSdcIn",&nlSdcIn);
	chain->SetBranchAddress("ntSdcIn",&ntSdcIn);
	chain->SetBranchAddress("nhSdcIn",&nhSdcIn);
	chain->SetBranchAddress("chisqrSdcIn",&chisqrSdcIn);
	chain->SetBranchAddress("x0SdcIn",x0SdcIn);
	chain->SetBranchAddress("y0SdcIn",y0SdcIn);
	chain->SetBranchAddress("u0SdcIn",u0SdcIn);
	chain->SetBranchAddress("v0SdcIn",v0SdcIn);

	chain->SetBranchAddress("nlSdcOut",&nlSdcOut);
	chain->SetBranchAddress("ntSdcOut",&ntSdcOut);
	chain->SetBranchAddress("nhSdcOut",&nhSdcOut);
	chain->SetBranchAddress("chisqrSdcOut",&chisqrSdcOut);
	chain->SetBranchAddress("x0SdcOut",x0SdcOut);
	chain->SetBranchAddress("y0SdcOut",y0SdcOut);
	chain->SetBranchAddress("u0SdcOut",u0SdcOut);
	chain->SetBranchAddress("v0SdcOut",v0SdcOut);
}
