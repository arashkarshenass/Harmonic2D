/*=================================================
                Outlet Boundary Class source
=================================================*/

//here u is Un (normal velocity) and v is Ut (tangential velocity)

#include "Outlet.h"
#include "../InputOutput/Mesh.h"
#include "../InputOutput/SteadyReader.h"
#include "../Solver/SteadyJacobian.h"
#include "../Inputs.h"

#include <cmath>
#include <iostream>
using namespace std;

void Outlet::Initializer(double*wrCellCurr,double*wiCellCurr,double*wrCellNeigh,double*wiCellNeigh,Mesh*mp,SteadyReader*srp,SteadyJacobian*sjp)
{

	wRealDomain=wrCellCurr;
	wImagDomain=wiCellCurr;
	wRealGhost=wrCellNeigh;
	wImagGhost=wiCellNeigh;
	qSteadyDomain=srp->GetPrimitiveDomain();
	qSteadyOutlet=sjp->GetQsOutlet();
	mapOutlet=mp->GetMapOutlet();
	nxMat=mp->GetNx();
	nyMat=mp->GetNy();
	g1=gamma-1;
	pRealGhost=dP*stdPrsOut;														//real part of physical boundary condition
	pImagGhost=0;																	//imaginary part of physical boundary condition
}

void Outlet::CalculateGhost(int cellNumDomain,int faceNumDomain){

	//Numeric boundary condition
	wRealGhost[0]=wRealDomain[0];															//real numeric boundary condition 1
	wImagGhost[0]=wImagDomain[0];															//imag numeric boundary condition 1
	wRealGhost[1]=wRealDomain[1];															//real numeric boundary condition 2
	wImagGhost[1]=wImagDomain[1];															//imag numeric boundary condition 2
	wRealGhost[2]=wRealDomain[2];															//real numeric boundary condition 3
	wImagGhost[2]=wImagDomain[2];															//imag numeric boundary condition 3

	//find associated outlet ghost cell
	int outletNum=0;
	while (mapOutlet[outletNum][0]!=cellNumDomain){
		outletNum=outletNum+1;
	}

	//steady variables of ghost
	double rhoSteadyGhost=qSteadyOutlet[outletNum][0];																			//steady rho ghost
	double unSteadyGhost=qSteadyOutlet[outletNum][1];																			//steady Un ghost
	double utSteadyGhost=qSteadyOutlet[outletNum][2];																			//steady Ut ghost
	double uSteadyGhost=unSteadyGhost*nxMat[cellNumDomain][faceNumDomain]-utSteadyGhost*nyMat[cellNumDomain][faceNumDomain];	//steady u ghost
	double vSteadyGhost=unSteadyGhost*nyMat[cellNumDomain][faceNumDomain]+utSteadyGhost*nxMat[cellNumDomain][faceNumDomain];	//steady v ghost

	//calculate conservative harmonic variables at ghost
	double rhoRealGhost=wRealGhost[0];																							//rho real ghost
	double rhoImagGhost=wImagGhost[0];																							//rho imag ghost
	double urg=(wRealGhost[1]-rhoRealGhost*uSteadyGhost)/rhoSteadyGhost;														//u real ghost
	double uig=(wImagGhost[1]-rhoImagGhost*uSteadyGhost)/rhoSteadyGhost;														//u imag ghost
	double vrg=(wRealGhost[2]-rhoRealGhost*vSteadyGhost)/rhoSteadyGhost;														//v real ghost
	double vig=(wImagGhost[2]-rhoImagGhost*vSteadyGhost)/rhoSteadyGhost;														//v imagy ghost
	wRealGhost[3]=pRealGhost/g1 + rhoRealGhost/2*(uSteadyGhost*uSteadyGhost+vSteadyGhost*vSteadyGhost) + rhoSteadyGhost*(uSteadyGhost*urg+vSteadyGhost*vrg);	//real conservative harmonic 3 at ghost
	wImagGhost[3]=pImagGhost/g1 + rhoImagGhost/2*(uSteadyGhost*uSteadyGhost+vSteadyGhost*vSteadyGhost) + rhoSteadyGhost*(uSteadyGhost*uig+vSteadyGhost*vig);	//imag conservative harmonic 3 at ghost
}

