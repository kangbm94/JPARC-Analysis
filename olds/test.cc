#include "Utils.hh"
void test(){
	double data[10]={0};
	fstream file;
	file.open("ToFParam_default");
	cout<<"First"<<endl;
	ReadTSV(file,data);
	for(int i=0;i<7;i++){
		cout<<data[i]<<endl;
	}
	cout<<"second"<<endl;
	ReadTSV(file,data);
	for(int i=0;i<7;i++){
		cout<<data[i]<<endl;
	}
	cout<<"Third"<<endl;
	ReadTSV(file,data);
	for(int i=0;i<7;i++){
		cout<<data[i]<<endl;
	}
}
