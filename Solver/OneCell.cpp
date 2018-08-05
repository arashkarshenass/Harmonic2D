/*=================================================
                Flux Class source
=================================================*/


#include "OneCell.h"
#include "../InputOutput/Mesh.h"
#include "../InputOutput/SteadyReader.h"
#include "../BoundaryCondition/Inlet.h"
#include "../BoundaryCondition/Outlet.h"
#include "../BoundaryCondition/Wall.h"
#include "../FaceFlux/StegerWarming.h"
#include "../Inputs.h"
#include "SteadyJacobian.h"

Inlet inlet;						//create inlet object
Outlet outlet;						//create outlet object
Wall wall;							//create wall object
StegerWarming sw;					//create upwind object

OneCell::OneCell(double*cfr,double*cfi,double** wRTemp,double** wImTemp ,Mesh* mpT,SteadyReader* srp,SteadyJacobian* sjp){
	mp=mpT;
	cellFluxReal=cfr;
	cellFluxImag=cfi;
	wrDomain=wRTemp;									//real conservative harmonic variables of domain's all cells
	wiDomain=wImTemp;								//imag conservative harmonic variables of domain's all cells
	wrCellCurr=new double[4];							//real conservative harmonic variables of current cell
	wiCellCurr=new double[4];							//imag conservative harmonic variables of current cell
	wrFaceCurr=new double[4];							//real conservative harmonic variables of current cell
	wiFaceCurr=new double[4];							//imag conservative harmonic variables of current cell
	wrCellNeigh=new double[4];							//real conservative harmonic variables of neighbor cell
	wiCellNeigh=new double[4];							//imag conservative harmonic variables of neighbor cell
	wrFaceNeigh=new double[4];							//real conservative harmonic variables of neighbor cell
	wiFaceNeigh=new double[4];							//imag conservative harmonic variables of neighbor cell
	faceFluxReal=new double[4];					//real fluxes passing a face
	faceFluxImag=new double[4];					//imag flux passing a face
	jacobPlus=new double*[4];					//the face jacobian from current side
	jacobMinus=new double*[4];					//the face jacobian from neighbor side
	for(int i=0;i<4;i++){
		jacobPlus[i]=new double[4];
		jacobMinus[i]=new double[4];
	}
	jacobCurrentSideAllDomain=sjp->GetJacobianCurrentSide();		//face jaconbian of domain's all cells from current side
	jacobNeighborSideAllDomain=sjp->GetJacobianNeighborSide();		//face jacobians of domain's all cells from neighbor side
	qs=srp->GetPrimitiveDomain();									//steady primitive variables of domain's all cells (in xy)
	connectivity=mp->GetConnectivity();					//get connectuvity of cells from mesh obj
	edgeNumber=mp->GetEdgeNumber();						//get edge number of cells from mesh obj
	nxMat=mp->GetNx();
	nyMat=mp->GetNy();
	dlMat=mp->GetDl();


	wall.Initializer(wrFaceCurr,wiFaceCurr,wrFaceNeigh,wiFaceNeigh);					//initializing wall object
	inlet.Initializer(wrCellCurr,wiCellCurr,wrCellNeigh,wiCellNeigh,mp,srp,sjp);				//initializing inlet object
	outlet.Initializer(wrCellCurr,wiCellCurr,wrCellNeigh,wiCellNeigh,mp,srp,sjp);				//initializing outletaksel object
	sw.Initializer(wrFaceCurr,wiFaceCurr,wrFaceNeigh,wiFaceNeigh,faceFluxReal,faceFluxImag,jacobPlus,jacobMinus);
}

void OneCell::ClaculateCellFlux(int cn){
	cellNum=cn;
	wrCellCurr[0]=wrDomain[cellNum][0];		wrCellCurr[1]=wrDomain[cellNum][1];					//real conservative variables current cell (in xy)
	wrCellCurr[2]=wrDomain[cellNum][2];		wrCellCurr[3]=wrDomain[cellNum][3];
	wiCellCurr[0]=wiDomain[cellNum][0];		wiCellCurr[1]=wiDomain[cellNum][1];					//imag conservative variables current cell (in xy)
	wiCellCurr[2]=wiDomain[cellNum][2];		wiCellCurr[3]=wiDomain[cellNum][3];

	for(faceNum=0;faceNum<edgeNumber;faceNum++){			//loop for faces of a cell

		XYtoNT(wrCellCurr,wiCellCurr,wrFaceCurr,wiFaceCurr);									//calculate real and imag conservative variables current cell (in nt)

		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				jacobPlus[i][j]=jacobCurrentSideAllDomain[cellNum][faceNum][i][j];				//jacobian plus of face
				jacobMinus[i][j]=jacobNeighborSideAllDomain[cellNum][faceNum][i][j];			//jacobian minus forface
			}
		}

		int neigh=connectivity[cellNum][faceNum];					//neighbour cell number
		switch(neigh){
		case -1:
				wall.CalculateGhost();															//calculates harmonic variable of wall ghost cell (in nt)
				break;
		case -2:
				inlet.CalculateGhost(cellNum,faceNum);											//calculates harmonic variable of inlet ghost cell (in xy)
				XYtoNT(wrCellNeigh,wiCellNeigh,wrFaceNeigh,wiFaceNeigh);						//transfer into nt
				break;
		case -3:
				outlet.CalculateGhost(cellNum,faceNum);											//calculates harmonic variable of outlet ghost cell (in xy)
				XYtoNT(wrCellNeigh,wiCellNeigh,wrFaceNeigh,wiFaceNeigh);						//transfer into nt
				break;
		default:
				wrCellNeigh[0]=wrDomain[neigh][0];	wrCellNeigh[1]=wrDomain[neigh][1];			//neighbor cell (in xy)
				wrCellNeigh[2]=wrDomain[neigh][2];	wrCellNeigh[3]=wrDomain[neigh][3];
				wiCellNeigh[0]=wiDomain[neigh][0];	wiCellNeigh[1]=wiDomain[neigh][1];
				wiCellNeigh[2]=wiDomain[neigh][2];	wiCellNeigh[3]=wiDomain[neigh][3];
				XYtoNT(wrCellNeigh,wiCellNeigh,wrFaceNeigh,wiFaceNeigh);						//transfer into nt
				break;
		}//end of switch


		sw.CalculateFaceFlux(dlMat[cellNum][faceNum]);			//calculating face flux in n-t coordinate
		NTtoXY();												//transfer face flux from n-t tot x-y

		cellFluxReal[0]+=faceFluxReal[0];				//add face real fluxes into cell real fluxes
		cellFluxReal[1]+=faceFluxReal[1];
		cellFluxReal[2]+=faceFluxReal[2];
		cellFluxReal[3]+=faceFluxReal[3];

		cellFluxImag[0]+=faceFluxImag[0];				//add face imag fluxes into cell imag fluxes
		cellFluxImag[1]+=faceFluxImag[1];
		cellFluxImag[2]+=faceFluxImag[2];
		cellFluxImag[3]+=faceFluxImag[3];
		}//end of face loop
}//end of function

void OneCell::XYtoNT(double*wrCell,double*wiCell,double*wrFace,double*wiFace){
	wrFace[0]=wrCell[0];
	wrFace[1]=wrCell[1]*nxMat[cellNum][faceNum]+wrCell[2]*nyMat[cellNum][faceNum];
	wrFace[2]=-wrCell[1]*nyMat[cellNum][faceNum]+wrCell[2]*nxMat[cellNum][faceNum];
	wrFace[3]=wrCell[3];

	wiFace[0]=wiCell[0];
	wiFace[1]=wiCell[1]*nxMat[cellNum][faceNum]+wiCell[2]*nyMat[cellNum][faceNum];
	wiFace[2]=-wiCell[1]*nyMat[cellNum][faceNum]+wiCell[2]*nxMat[cellNum][faceNum];
	wiFace[3]=wiCell[3];
}

void OneCell::NTtoXY(){
	double faceFluxRealN=faceFluxReal[1];
	double faceFluxRealT=faceFluxReal[2];
	faceFluxReal[1]=faceFluxRealN*nxMat[cellNum][faceNum]-faceFluxRealT*nyMat[cellNum][faceNum];
	faceFluxReal[2]=faceFluxRealN*nyMat[cellNum][faceNum]+faceFluxRealT*nxMat[cellNum][faceNum];

	double faceFluxImagN=faceFluxImag[1];
	double faceFluxImagT=faceFluxImag[2];
	faceFluxImag[1]=faceFluxImagN*nxMat[cellNum][faceNum]-faceFluxImagT*nyMat[cellNum][faceNum];
	faceFluxImag[2]=faceFluxImagN*nyMat[cellNum][faceNum]+faceFluxImagT*nxMat[cellNum][faceNum];
}
