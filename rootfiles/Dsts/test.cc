#include "/Users/MIN/ROOTSharedLibs/LinesInSpace.hh"
void test(){
	Line tr1 = Line(0.1,0,0,0.1);
	Line tr2 = Line(0.095,0,0,0.095);
	double dist = tr1.Angle(tr2);
	cout<<dist<<endl;
}
