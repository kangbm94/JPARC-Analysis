void LogReader(){
	fstream f;
	int nl=1342;
	//	nl=370;
	f.open("E42RunSummary2021May.csv");
	string cont[100];
	ofstream fo;
//	fo.open("HRuns.txt");
	fo.open("CH2Runs.txt");
	int cnt=0;
	for(int i=0;i<nl;i++){
//	while	(ReadingCSV(f,cont)){
	ReadCSV(f,cont);
		if(i==0){
			continue;
		}
		if(cont[3].at(0)=='C'&&cont[3].at(1)=='H'){
			fo<<cont[2]<<","<<stoi(cont[7])<<","<<stoi(cont[8])<<","<<stoi(cont[9])<<endl;
			//		cout<<"H Production: "<<cont[2]<<endl;
		}
	}
fo.close();
}
