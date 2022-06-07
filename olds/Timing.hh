double dgaus(double* x, double* p){
	double g1 = Gaussian(x[0],p[0],p[1],p[2]);//mean,sigma,amplitude
	double g2 = Gaussian(x[0],p[0]+p[3],p[4],p[5]);
	return g1+g2;	
}


double Time(){
	return 1;
}
