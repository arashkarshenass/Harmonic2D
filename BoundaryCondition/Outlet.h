/*=================================================
                Outlet Boundary Class header
=================================================*/
#ifndef OUTLET_H_
#define OUTLET_H_
class Mesh;
class SteadyReader;
class SteadyJacobian;
class Outlet
{
public:
	void Initializer (double*,double*,double*,double*,Mesh*,SteadyReader*,SteadyJacobian*);
	void CalculateGhost(int,int);

private:
	double* wRealDomain=nullptr;
	double* wImagDomain=nullptr;
	double* wRealGhost=nullptr;
	double* wImagGhost=nullptr;
	double** qSteadyDomain=nullptr;
	double** qSteadyOutlet=nullptr;
	int**    mapOutlet=nullptr;
	double** nxMat=nullptr;
	double** nyMat=nullptr;
	double  pRealGhost=0;
	double  pImagGhost=0;
	double  g1=0;
};

#endif
