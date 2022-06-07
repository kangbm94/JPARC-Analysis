//Parameters.
//Id	Name	Xk	Yk	Zk	TA	RA1	RA2	L	Res	W0	dXdW	Ofs
double sdc1_x=0,sdc1_y=0;
double sdc1_z=-932.7;

double sdc2_x=0,sdc2_y=0;
double sdc2_z=-480.66;

double sdc3_x=0,sdc3_y=0;
double sdc3_z=894.7;

double sdc4_x=0,sdc4_y=0;
double sdc4_z=1120.6;

//Common//
int trigtpat[32];
int evnum;
//Beam//
double CTime0,CBtof0;
//Sch
const int SchNL= 100;//AddressArrayLength
int nhSch;
int csSch[SchNL];double tSch[SchNL];double wSch[SchNL];double SchPos[SchNL];double SchSeg[SchNL];
//BC//
const int BcOutNT=10;
int nlBcOut,ntBcOut;
int nhBcOut[BcOutNT];
double chisqrBcOut[BcOutNT];
double x0BcOut[BcOutNT];double y0BcOut[BcOutNT];double u0BcOut[BcOutNT];double v0BcOut[BcOutNT];
double xtgtBcOut[BcOutNT];double ytgtBcOut[BcOutNT];double xbh2BcOut[BcOutNT];double ybh2BcOut[BcOutNT];
//SdcIn//
const int SdcInNT=10;
int nlSdcIn,ntSdcIn;
int nhSdcIn[SdcInNT];
double chisqrSdcIn[SdcInNT];
double x0SdcIn[SdcInNT];double y0SdcIn[SdcInNT];double u0SdcIn[SdcInNT];double v0SdcIn[SdcInNT];
double xtgtSdcIn[SdcInNT];double ytgtSdcIn[SdcInNT];
//SdcOut//
const int SdcOutNT=10;
int nlSdcOut,ntSdcOut;
int nhSdcOut[SdcOutNT];
double chisqrSdcOut[SdcOutNT];
double x0SdcOut[SdcOutNT];double y0SdcOut[SdcOutNT];double u0SdcOut[SdcOutNT];double v0SdcOut[SdcOutNT];
double xtgtSdcOut[SdcOutNT];double ytgtSdcOut[SdcOutNT];

//KuramaHodo//
