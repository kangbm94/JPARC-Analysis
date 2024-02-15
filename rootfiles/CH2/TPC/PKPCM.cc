void PKPCM(){
}
void PKPCM(double p){
	double mk = 0.493;
	double mp = 0.938;
	TLorentzVector B(0,0,p,hypot(p,mk));
	TLorentzVector T(0,0,0,mp);
	auto C = B+T;
	cout<<C.Mag()<<endl;
}
