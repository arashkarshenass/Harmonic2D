/*=================================================
                Wall Boundary Class source
=================================================*/

//here u is Un (normal velocity) and v is Ut (tangential velocity)

#include "Wall.h"
#include "../InputOutput/Mesh.h"
#include "../InputOutput/SteadyReader.h"
#include "../Inputs.h"
#include "../Solver/SteadyJacobian.h"

void Wall::Initializer(double*wrFaceCurrTemp,double*wiFaceCurrTemp,double*wrFaceNeighTemp,double*wiFaceNeighTemp)
{
	wrFaceDomain=wrFaceCurrTemp;
	wiFaceDomain=wiFaceCurrTemp;
	wrFaceGhost=wrFaceNeighTemp;
	wiFaceGhost=wiFaceNeighTemp;
}

void Wall::CalculateGhost()
{
	wrFaceGhost[0]=wrFaceDomain[0];
	wrFaceGhost[1]=-wrFaceDomain[1];
	wrFaceGhost[2]=wrFaceDomain[2];
	wrFaceGhost[3]=wrFaceDomain[3];

	wiFaceGhost[0]=wiFaceDomain[0];
	wiFaceGhost[1]=-wiFaceDomain[1];
	wiFaceGhost[2]=wiFaceDomain[2];
	wiFaceGhost[3]=wiFaceDomain[3];

}
