#include "KKManager.hh"
TCut VertexCut(vector <double> Origin, vector <double> Size,int n){
	TCut CutX = Form("abs(vtx[%d]-%f)<%f",n,Origin[0],Size[0]);
	TCut CutY = Form("abs(vty[%d]-%f)<%f",n,Origin[1],Size[1]);
	TCut CutZ = Form("abs(vtx[%d]-%f)<%f",n,Origin[2],Size[2]);
	return CutX&&CutY&&CutZ;
}
