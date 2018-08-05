#include <cmath>
#include "LocalTimeStep.h"
#include "../InputOutput/Mesh.h"
#include "../InputOutput/SteadyReader.h"
#include "../Inputs.h"
using namespace std;

LocalTimeStep::LocalTimeStep(Mesh* mp,SteadyReader* srp,double* dt) {
	int celltot=mp->GetCelltot();
	double* minDis=mp->GetMinDis();
	double** qs=srp->GetPrimitiveDomain();
	double u,v,V,a,S;
	for(int i=0;i<celltot;i++){
		u=qs[i][1];	v=qs[i][2];
		V=sqrt(u*u+v*v);
		a=sqrt(gamma*qs[i][3]/qs[i][0]);
		S=a+V;
		dt[i]=minDis[i]/S;
	}
}
