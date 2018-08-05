#ifndef SOLVER_STEADYJACOBIAN_H_
#define SOLVER_STEADYJACOBIAN_H_

class Mesh;
class SteadyReader;
class Utility;

class SteadyJacobian {
public:
	SteadyJacobian(Mesh*,SteadyReader*,Utility*);
	void CalculateCoefficient(int,int);
	void JacobianMatrix(int,int,double,double,double,double,double****);
	void SteadyWallBC(int);
	void SteadyInletBC(int);
	void SteadyOutletBC(int);
	double** GetQsWall();
	double** GetQsInlet();
	double** GetQsOutlet();
	double**** GetJacobianCurrentSide();
	double**** GetJacobianNeighborSide();


private:
	Mesh* mp=nullptr;
	int** connectivity=nullptr;
	int** mapWall=nullptr;
	int** mapInlet=nullptr;
	int** mapOutlet=nullptr;
	double** qsDom=nullptr;
	double** unMatDom=nullptr;
	double** utMatDom=nullptr;
	double** nxMat=nullptr;
	double** nyMat=nullptr;
	double** qsWall=nullptr;
	double** qsInlet=nullptr;
	double** qsOutlet=nullptr;
	double**** jacobCurrentSide=nullptr;
	double**** jacobNeighborSide=nullptr;
	double* lam=nullptr;
	double g1=0;
};

#endif
