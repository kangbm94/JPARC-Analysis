double m_pi = PionMass()/1000;//MeV to GeV
double m_pi2= m_pi*m_pi;
double tof_to_m2(double stof, double p, double path){
	Double_t beta = path/stof/LightSpeed();//LightSpeed: mm / ns
	return p*p*(1.-beta*beta)/beta/beta;//stof: ns,  p: GeV/c, path: mm
}
double m2_to_tof(double m2, double p, double path){
	return path*sqrt(m2+p*p)/p/LightSpeed();
}

double stof_before[24]={
19.8423	,
18.8674	,
18.031	,
16.9369	,
16.3764	,
14.666	,
15.4117	,
13.9669	,
13.6375	,
13.9333	,
14.2121	,
13.5435	,
13.6645	,
13.7424	,
4.4016	,
14.4078	,
15.3851	,
16.255	,
16.8695	,
16.8441	,
17.4235	,
17.6772	,
19.4738	,
19.9883	,
};
