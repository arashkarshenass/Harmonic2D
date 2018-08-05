#include "SteadyJacobian.h"

#include "../InputOutput/Mesh.h"
#include "../InputOutput/SteadyReader.h"
#include "../Utility/Utility.h"
#include "../Inputs.h"
#include <cmath>

using namespace std;

SteadyJacobian::SteadyJacobian(Mesh* mpT,SteadyReader* srp,Utility* up) {

	mp=mpT;
	int cellTotalDomain=mp->GetCelltot();			//number of domain cells
	int cellTotalWall=mp->GetCellTotalWall();				//number of wall ghost cells
	int cellTotalInlet=mp->GetCellTotalInlet();			//number of inlet ghost cells
	int cellTotalOutlet=mp->GetCellTotalOutlet();			//number of outlet ghost cells
	int edgeSize=mp->GetEdgeNumber();				//number of faces of each domian cells
	connectivity=mp->GetConnectivity();				//map pf interior cells
	mapWall=mp->GetMapWall();						//map of wall ghost cells
	mapInlet=mp->GetMapInlet();						//map of inlet ghost cells
	mapOutlet=mp->GetMapOutlet();					//map of outlet ghost cells
	qsDom=srp->GetPrimitiveDomain();				//primitive variable of domain cells (in xy)
	unMatDom=srp->GetUnDom();						//normal velocity of each face of each cell of domain
	utMatDom=srp->GetUtDom();						//normal velocity of each face of each cell of domain
	nxMat=mp->GetNx();								//nx of each face of each cell of domain
	nyMat=mp->GetNy();								//ny of each face of each cell of domain
	lam=new double[4];								//eigen value lamda
	g1=gamma-1;										//gamma-1


	//allacating memory for jacobians of domain from current cell side
	jacobCurrentSide=new double***[cellTotalDomain];			//for each cell
	for (int i=0;i<cellTotalDomain;i++){
		jacobCurrentSide[i]=new double**[edgeSize];			//for each face
		for (int j=0;j<4;j++){
			jacobCurrentSide[i][j]=new double*[4];			//rows
			for (int k=0;k<4;k++)
				jacobCurrentSide[i][j][k]=new double[4];	//column
		}
	}

	//allacating memory for jacobians of domain from neighbor cell side
	jacobNeighborSide=new double***[cellTotalDomain];			//for each cell
	for (int i=0;i<cellTotalDomain;i++){
		jacobNeighborSide[i]=new double**[edgeSize];			//for each face
		for (int j=0;j<4;j++){
			jacobNeighborSide[i][j]=new double*[4];			//rows
			for (int k=0;k<4;k++)
				jacobNeighborSide[i][j][k]=new double[4];	//column
			}
		}

	//primitive variables of wall ghost cell (in NT)
	qsWall=new double*[cellTotalWall];
	for (int i=0;i<cellTotalWall;i++)
		qsWall[i]=new double[4];

	//primitive variables of inlet ghost cell (in NT)
	qsInlet=new double*[cellTotalInlet];
	for (int i=0;i<cellTotalInlet;i++)
		qsInlet[i]=new double[4];

	//primitive variables of outlet ghost cell (in NT)
	qsOutlet=new double*[cellTotalOutlet];
	for (int i=0;i<cellTotalOutlet;i++)
		qsOutlet[i]=new double[4];

	//wall ghost cell primitive variables using wall BC
	for (int wallNum=0;wallNum<cellTotalWall;wallNum++)
		SteadyWallBC(wallNum);

	//inlet ghost cell primitive variables using inlet BC
	for (int inletNum=0;inletNum<cellTotalInlet;inletNum++)
		SteadyInletBC(inletNum);

	//outlet ghost cell primitive variables using inlet BC
	for (int outletNum=0;outletNum<cellTotalOutlet;outletNum++)
		SteadyOutletBC(outletNum);

	//Steady Jacobian Coefficient
	for (int cellNum=0;cellNum<cellTotalDomain;cellNum++)
		for (int faceNum=0;faceNum<edgeSize;faceNum++)
			CalculateCoefficient(cellNum,faceNum);

	delete[] lam;
}

void SteadyJacobian::CalculateCoefficient(int cn,int fn){
	double un,ut,a,s2;

	//-------jacobian from current side
	a=sqrt(gamma*qsDom[cn][3]/qsDom[cn][0]);		s2=qsDom[cn][1]*qsDom[cn][1]+qsDom[cn][2]*qsDom[cn][2];
	un=unMatDom[cn][fn];		ut=utMatDom[cn][fn];
	lam[0]=0.5*(un+abs(un));
	lam[1]=lam[0];
	lam[2]=0.5*((un+a)+abs(un+a));
	lam[3]=0.5*((un-a)+abs(un-a));
	JacobianMatrix(cn,fn,un,ut,a,s2,jacobCurrentSide);


	//-------jacobian from neighbor side
	//find neighbour cell
	int neighCell=connectivity[cn][fn];
	int wallNum,inletNum,outletNum,neighFace;
	switch(neighCell){
	case -1:		//neighcell is wall ghost cell
			//which wall ghost cell?
			wallNum=0;
			while(mapWall[wallNum][0]!=cn)
				wallNum=wallNum+1;
			//primituve variables of that wall ghost cell (in n-t)
			un=qsWall[wallNum][1];
			ut=qsWall[wallNum][2];
			a=sqrt(gamma*qsWall[wallNum][3]/qsWall[wallNum][0]);
			break;
	case -2:		//neighcell is inlet ghost cell
			//which inlet ghost cell?
			inletNum=0;
			while(mapInlet[inletNum][0]!=cn)
				inletNum=inletNum+1;
			//primituve variables of that inlet ghost cell (in n-t)
			un=qsInlet[inletNum][1];
			ut=qsInlet[inletNum][2];
			a=sqrt(gamma*qsInlet[inletNum][3]/qsInlet[inletNum][0]);
			break;
	case -3:		//neighcell is outlet ghost cell
			//which outlet ghost cell?
			outletNum=0;
			while(mapOutlet[outletNum][0]!=cn)
				outletNum=outletNum+1;
			//primituve variables of that outlet ghost cell (in n-t)
			un=qsOutlet[outletNum][1];
			ut=qsOutlet[outletNum][2];
			a=sqrt(gamma*qsOutlet[outletNum][3]/qsOutlet[outletNum][0]);
			break;
	default:		//neighbour is interior cell with number neigh
			//which face of the neighbour cell?
			neighFace=mp->GetNeighborFaceNum(cn,fn,neighCell);
			//primituve variables of that neighbor cell (in n-t)
			un=-unMatDom[neighCell][neighFace];
			ut=-utMatDom[neighCell][neighFace];
			a=sqrt(gamma*qsDom[neighCell][3]/qsDom[neighCell][0]);
			break;
	}//end of switch
	s2=un*un+ut*ut;

	lam[0]=0.5*(un-abs(un));
	lam[1]=lam[0];
	lam[2]=0.5*((un+a)-abs(un+a));
	lam[3]=0.5*((un-a)-abs(un-a));
	JacobianMatrix(cn,fn,un,ut,a,s2,jacobNeighborSide);
}

void SteadyJacobian::SteadyWallBC(int wallNum){

	//cell number to which wall is connected
	int cellNumDomain=mapWall[wallNum][0];
	//face number of this cell
	int faceNumDomain=mapWall[wallNum][1];

	//primitive variables of domain cell
	double rhoDom=qsDom[cellNumDomain][0];			double pDom=qsDom[cellNumDomain][3];
	double unDom=unMatDom[cellNumDomain][faceNumDomain];	double utDom=utMatDom[cellNumDomain][faceNumDomain];

	//primitive variables of the wall ghost cell
	double rhoGst=rhoDom;					double pGst=pDom;
	double unGst=-unDom;					double utGst=utDom;

	qsWall[wallNum][0]=rhoGst;	qsWall[wallNum][3]=pGst;
	qsWall[wallNum][1]=unGst;	qsWall[wallNum][2]=utGst;
}

void SteadyJacobian::SteadyInletBC(int inletNum){
	//Note that inlet boundaru condition relations are based on x-y values but NOT n-t.

	//cell number to which inlet is connected
	int cellNumDomain=mapInlet[inletNum][0];
	//face number of this cell
	int faceNumDomain=mapInlet[inletNum][1];

	//variables at domain (in x-y)
	double rhoDom=qsDom[cellNumDomain][0];		double pDom=qsDom[cellNumDomain][3];
	double uDom=qsDom[cellNumDomain][1];		double vDom=qsDom[cellNumDomain][2];
	double aDom=sqrt(gamma*pDom/rhoDom);

	//blazek method
	double unDom=unMatDom[cellNumDomain][faceNumDomain];
	double sDom=sqrt(uDom*uDom+vDom/vDom);
	double rMinus=unDom-(2*aDom/g1);
	double cosThetaDom=-unDom/sDom;
	double aStag2=aDom*aDom+0.5*g1*sDom*sDom;
	double aBoundary=-rMinus*g1/(g1*cosThetaDom*cosThetaDom+2)*(1+cosThetaDom*sqrt(((g1*cosThetaDom*cosThetaDom+2)*aStag2)/(g1*rMinus*rMinus)-g1/2));
	double aGhost=2*aBoundary-aDom;
	double tGhost=tStagIn*(aGhost*aGhost/aStag2);
	double pGhost=pStagIn*pow(tGhost/tStagIn,gamma/g1);
	double rhoGhost=pGhost/R/tGhost;
	double sGhost=sqrt(2*cp*(tStagIn-tGhost));
	double uGhost=sGhost*cos(thetaIn);
	double vGhost=sGhost*sin(thetaIn);

	//variables at ghost (in n-t)
	double unGhost=uGhost*nxMat[cellNumDomain][faceNumDomain]+vGhost*nyMat[cellNumDomain][faceNumDomain];
	double utGhost=-uGhost*nyMat[cellNumDomain][faceNumDomain]+vGhost*nxMat[cellNumDomain][faceNumDomain];

	qsInlet[inletNum][0]=rhoGhost;		qsInlet[inletNum][1]=unGhost;	qsInlet[inletNum][2]=utGhost;	qsInlet[inletNum][3]=pGhost;
}

void SteadyJacobian::SteadyOutletBC(int outletNum){
	//Note that inlet boundaru condition relations are based on x-y values but NOT n-t.

	//cell number to which inlet is connected
	int cellNumDomain=mapOutlet[outletNum][0];
	//face number of this cell
	int faceNumDomain=mapOutlet[outletNum][1];

	//variables at domain (in x-y)
	double rhoDom=qsDom[cellNumDomain][0];		double pDom=qsDom[cellNumDomain][3];
	double uDom=qsDom[cellNumDomain][1];		double vDom=qsDom[cellNumDomain][2];
	double aDom=sqrt(gamma*pDom/rhoDom);

	//variables at ghost (in x-y)
	//predictor
	double pGhost=stdPrsOut;
	double rhoGhost=rhoDom+(pGhost-pDom)/aDom/aDom;
	double uGhost=abs(uDom-(pGhost-pDom)/rhoDom/aDom);
	double aGhost=sqrt(gamma*pGhost/rhoGhost);

	//corrector
	double roAvg=(rhoGhost+rhoDom)/2;
	double aAvg=(aGhost+aDom)/2;
	rhoGhost=rhoDom+(pGhost-pDom)/aAvg/aAvg;
	uGhost=abs(uDom-(pGhost-pDom)/roAvg/aAvg);
	double vGhost=0;

	//variables at ghost (in n-t)
	double unGhost=uGhost*nxMat[cellNumDomain][faceNumDomain]+vGhost*nyMat[cellNumDomain][faceNumDomain];
	double utGhost=-uGhost*nyMat[cellNumDomain][faceNumDomain]+vGhost*nxMat[cellNumDomain][faceNumDomain];

	qsOutlet[outletNum][0]=rhoGhost;	qsOutlet[outletNum][1]=unGhost;		qsOutlet[outletNum][2]=utGhost;		qsOutlet[outletNum][3]=pGhost;
}

void SteadyJacobian::JacobianMatrix(int cn,int fn,double un,double ut,double a,double s2,double****jacob){
	jacob[cn][fn][0][0]=lam[0]* (1-g1*s2/2/a/a) + lam[2]* ((g1*s2-2*un*a)/4/a/a) + lam[3]* ((g1*s2+2*un*a)/4/a/a);
	jacob[cn][fn][0][1]=lam[0]* (g1*un/a/a) + lam[2]* ((a-g1*un)/2/a/a) + lam[3]* (-(a+g1*un)/2/a/a);
	jacob[cn][fn][0][2]=g1*ut/2/a/a* (2*lam[0] - lam[2] - lam[3]);
	jacob[cn][fn][0][3]=-g1/2/a/a* (2*lam[0] - lam[2] - lam[3]);
	jacob[cn][fn][1][0]=lam[0]* ((2*a*a-g1*s2)*un/2/a/a) + lam[2]* ((un+a)*(g1*s2-2*a*un)/4/a/a) + lam[3]* ((un-a)*(g1*s2+2*a*un)/4/a/a);
	jacob[cn][fn][1][1]=lam[0]* (g1*un*un/a/a) + lam[2]* ((a+un)*(a-g1*un)/2/a/a) + lam[3]* ((a-un)*(a+g1*un)/2/a/a);
	jacob[cn][fn][1][2]=lam[0]* (g1*un*ut/a/a) + lam[2]* (-g1*(un+a)*ut/2/a/a) + lam[3]* (-g1*(un-a)*ut/2/a/a);
	jacob[cn][fn][1][3]=lam[0]* (-g1*un/a/a)+ lam[2]* (g1*(un+a)/2/a/a) + lam[3]* (g1*(un-a)/2/a/a);
	jacob[cn][fn][2][0]=lam[0]* (-g1*ut*s2/2/a/a) + lam[2]* ((g1*s2-2*un*a)*ut/4/a/a) + lam[3]*  ((g1*s2+2*un*a)*ut/4/a/a);
	jacob[cn][fn][2][1]=lam[0]* (g1*un*ut/a/a) + lam[2]* ((a-g1*un)*ut/2/a/a) + lam[3]* (-(a+g1*un)*ut/2/a/a) ;
	jacob[cn][fn][2][2]=g1*ut*ut/2/a/a* (2*lam[0] - lam[2] - lam[3]) + lam[0];
	jacob[cn][fn][2][3]=-g1*ut/2/a/a* (2*lam[0] - lam[2] - lam[3]);
	jacob[cn][fn][3][0]=lam[0]* ((-g1*s2*s2+2*a*a*(un*un-ut*ut))/4/a/a) + lam[2]* ((g1*s2-2*a*un)*(2*a*a/g1+s2+2*a*un)/8/a/a) + lam[3]* ((g1*s2+2*a*un)*(2*a*a/g1+s2-2*a*un)/8/a/a);
	jacob[cn][fn][3][1]=lam[0]* (g1*un*s2/2/a/a) + lam[2]* ((2*a*a/g1+s2+2*a*un)*(a-g1*un)/4/a/a) + lam[3]* (-(2*a*a/g1+s2-2*a*un)*(a+g1*un)/4/a/a);
	jacob[cn][fn][3][2]=lam[0]* ((g1*s2+2*a*a)*ut/2/a/a) + lam[2]* (-ut*g1*(2*a*a/g1+s2+2*un*a)/4/a/a) + lam[3]* (-ut*g1*(2*a*a/g1+s2-2*un*a)/4/a/a);
	jacob[cn][fn][3][3]=lam[0]* (-g1*s2/2/a/a) + lam[2]* (g1*(2*a*a/g1+s2+2*un*a)/4/a/a) + lam[3]* (g1*(2*a*a/g1+s2-2*un*a)/4/a/a);
}

double** SteadyJacobian::GetQsWall(){
	return qsWall;
}

double** SteadyJacobian::GetQsInlet(){
	return qsInlet;
}

double** SteadyJacobian::GetQsOutlet(){
	return qsOutlet;
}

double**** SteadyJacobian::GetJacobianCurrentSide(){
	return jacobCurrentSide;
}

double**** SteadyJacobian::GetJacobianNeighborSide(){
	return jacobNeighborSide;
}





