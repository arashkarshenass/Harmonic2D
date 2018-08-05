/*=================================================
               Cell Flux Class header
=================================================*/

#ifndef ONECELL_H_
#define ONECELL_H_

class Mesh;
class SteadyReader;
class SteadyJacobian;
class Utility;

class OneCell
{
public:
	OneCell(double*,double*,double**,double**,Mesh*,SteadyReader*,SteadyJacobian*);
	void ClaculateCellFlux(int);
	void XYtoNT(double*,double*,double*,double*);
	void NTtoXY();

private:
	Mesh* mp=nullptr;
	double* cellFluxReal=nullptr;
	double* cellFluxImag=nullptr;
	double** wrDomain=nullptr;
	double** wiDomain=nullptr;
	double* wrCellCurr=nullptr;
	double* wiCellCurr=nullptr;
	double* wrCellNeigh=nullptr;
	double* wiCellNeigh=nullptr;
	double* wrFaceCurr=nullptr;
	double* wiFaceCurr=nullptr;
	double* wrFaceNeigh=nullptr;
	double* wiFaceNeigh=nullptr;
	double* faceFluxReal=nullptr;
	double* faceFluxImag=nullptr;
	double** jacobPlus=nullptr;
	double** jacobMinus=nullptr;
	double**** jacobCurrentSideAllDomain=nullptr;
	double**** jacobNeighborSideAllDomain=nullptr;
	double** qs=nullptr;
	int** connectivity=nullptr;
	int edgeNumber=0;
	double** nxMat=nullptr;
	double** nyMat=nullptr;
	double** dlMat=nullptr;
	int cellNum=0;
	int faceNum=0;





};

#endif
