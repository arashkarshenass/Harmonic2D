/*=================================================
                Inlet Boundary Class header
=================================================*/
#ifndef INLET_H_
#define INLET_H_

class Mesh;
class SteadyReader;
class SteadyJacobian;

class Inlet
{
public:
	void Initializer (double*,double*,double*,double*,Mesh*,SteadyReader*,SteadyJacobian*);
	void CalculateGhost(int,int);
private:

	double*  wRealDomain=nullptr;
	double*  wImagDomain=nullptr;
	double*	 wRealGhost=nullptr;
	double*  wImagGhost=nullptr;
	double** qSteadyDomain=nullptr;
	double** qSteadyInlet=nullptr;
	int** 	 mapInlet=nullptr;
	double** nxMat=nullptr;
	double** nyMat=nullptr;
	double 	 g1=0;
};

#endif
