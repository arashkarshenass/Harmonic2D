/*=================================================
                Inlet Boundary Class source
=================================================*/

#include "Inlet.h"
#include "../InputOutput/Mesh.h"
#include "../InputOutput/SteadyReader.h"
#include "../Solver/SteadyJacobian.h"
#include "../Inputs.h"

#include <cmath>
#include <iostream>
using namespace std;

void Inlet::Initializer(double*wrCellCurr,double*wiCellCurr,double*wrCellNeigh,double*wiCellNeigh,Mesh*mp,SteadyReader*srp,SteadyJacobian*sjp)
{
	wRealDomain=wrCellCurr;
	wImagDomain=wiCellCurr;
	wRealGhost=wrCellNeigh;
	wImagGhost=wiCellNeigh;
	qSteadyDomain=srp->GetPrimitiveDomain();
	qSteadyInlet=sjp->GetQsInlet();
	mapInlet=mp->GetMapInlet();
	nxMat=mp->GetNx();
	nyMat=mp->GetNy();



	g1=gamma-1;
}


void Inlet::CalculateGhost(int cellNumDomain, int faceNumDomain)
{

	//steady variables of Domian
	double rhoSteadyDomain=qSteadyDomain[cellNumDomain][0];															//steady rho domain
	double uSteadyDomain=qSteadyDomain[cellNumDomain][1];																//steady u domain
	double vSteadyDomain=qSteadyDomain[cellNumDomain][2];																//steady v domain
	double pSteadyDomain=qSteadyDomain[cellNumDomain][3];																//steady P domain
	double aSteadyDomain=sqrt(gamma*pSteadyDomain/rhoSteadyDomain);												//steadu a domain
	double sSteadyDomain=sqrt(uSteadyDomain*uSteadyDomain+vSteadyDomain*vSteadyDomain);							//steady s domain

	//harmonic primitive variables domain
	double rhoRealDomain=wRealDomain[0];																			//rho real domain
	double rhoImagDomain=wImagDomain[0];																			//rho imag domian
	double uRealDomain=(wRealDomain[1]-rhoRealDomain*uSteadyDomain)/rhoSteadyDomain;								//u real domain
	double uImagDomain=(wImagDomain[1]-rhoImagDomain*uSteadyDomain)/rhoSteadyDomain;								//u imag domain
	double vRealDomain=(wRealDomain[2]-rhoRealDomain*vSteadyDomain)/rhoSteadyDomain;								//v real domian
	double vImagDomain=(wImagDomain[2]-rhoImagDomain*vSteadyDomain)/rhoSteadyDomain;								//v imag domian
	double pRealDomain=g1*(wRealDomain[3]-rhoSteadyDomain*(uSteadyDomain*uRealDomain+vSteadyDomain*vRealDomain)-0.5*rhoRealDomain*(uSteadyDomain*uSteadyDomain+vSteadyDomain*vSteadyDomain));		//p real domian
	double pImagDomain=g1*(wImagDomain[3]-rhoSteadyDomain*(uSteadyDomain*uImagDomain+vSteadyDomain*vImagDomain)-0.5*rhoImagDomain*(uSteadyDomain*uSteadyDomain+vSteadyDomain*vSteadyDomain));		//p imag domian
	double aRealDomain=gamma/2/aSteadyDomain*(pRealDomain/rhoSteadyDomain-pSteadyDomain/rhoSteadyDomain/rhoSteadyDomain*rhoRealDomain);					//a real domain
	double aImagDomain=gamma/2/aSteadyDomain*(pImagDomain/rhoSteadyDomain-pSteadyDomain/rhoSteadyDomain/rhoSteadyDomain*rhoImagDomain);					//a imag domian
	double sRealDomain=(uSteadyDomain*uRealDomain+vSteadyDomain*vRealDomain)/sSteadyDomain;						//s real domain
	double sImagDomain=(uSteadyDomain*uImagDomain+vSteadyDomain*vImagDomain)/sSteadyDomain;						//s imag domain
	double mRealDomain=sRealDomain/aSteadyDomain-sSteadyDomain*aRealDomain/aSteadyDomain/aSteadyDomain;			//M real domain
	double mImagDomain=sImagDomain/aSteadyDomain-sSteadyDomain*aImagDomain/aSteadyDomain/aSteadyDomain;			//M imag domain

	//numeric boundary condition
	double mRealGhost=mRealDomain;																				//zero order extrapolation
	double mImagGhost=mImagDomain;																				//zero order extrapolation

	//find associated inlet ghost cell
	int inletNum=0;
	while (mapInlet[inletNum][0]!=cellNumDomain){
		inletNum=inletNum+1;
	}

	//steady variables of ghost
	double rhoSteadyGhost=qSteadyInlet[inletNum][0];												//steady rho ghost
	double unSteadyGhost=qSteadyInlet[inletNum][1];													//steady Un ghost
	double utSteadyGhost=qSteadyInlet[inletNum][2];													//steady Ut ghost
	double uSteadyGhost=unSteadyGhost*nxMat[cellNumDomain][faceNumDomain]-utSteadyGhost*nyMat[cellNumDomain][faceNumDomain];	//steady u ghost
	double vSteadyGhost=unSteadyGhost*nyMat[cellNumDomain][faceNumDomain]+utSteadyGhost*nxMat[cellNumDomain][faceNumDomain];	//steady v ghost
	double pSteadyGhost=qSteadyInlet[inletNum][3];													//steady P ghost
	double aSteadyGhost=sqrt(gamma*pSteadyGhost/rhoSteadyGhost);									//steady a ghost
	double mSteadyGhost=sqrt(uSteadyGhost*uSteadyGhost+vSteadyGhost*vSteadyGhost)/aSteadyGhost;		//steady M ghost




	//physical boundary condition 1
	double pRealGhost=-gamma*pStagIn*mSteadyGhost*mRealGhost/pow(1+0.5*g1*mSteadyGhost*mSteadyGhost,gamma/g1+1);				//static real P at ghost
	double pImagGhost=-gamma*pStagIn*mSteadyGhost*mImagGhost/pow(1+0.5*g1*mSteadyGhost*mSteadyGhost,gamma/g1+1);				//static imaginary P at ghost

	//physical boundary condition 2
	double roStagIn=pStagIn/tStagIn/R;																							//stagnation rho at ghost (working with rho is easier than temp)
	double rhoRealGhost=-roStagIn*mSteadyGhost*mRealGhost/pow(1+0.5*g1*mSteadyGhost*mSteadyGhost,gamma/g1);						//static real rho at ghost
	double rhoImagGhost=-roStagIn*mSteadyGhost*mImagGhost/pow(1+0.5*g1*mSteadyGhost*mSteadyGhost,gamma/g1);						//static imaginary rho at ghost

	//other harmonic primitive variables at ghost cell
	double aRealGhost=gamma/2/aSteadyGhost*(pRealGhost/rhoSteadyGhost-pSteadyGhost/rhoSteadyGhost/rhoSteadyGhost*rhoRealGhost);	//real sound speed at ghost
	double aImagGhost=gamma/2/aSteadyGhost*(pImagGhost/rhoSteadyGhost-pSteadyGhost/rhoSteadyGhost/rhoSteadyGhost*rhoImagGhost);	//imaginary sound speed at ghost
	double sRealGhost=mSteadyGhost*aRealGhost+mRealGhost*aSteadyGhost;															//real absolute velocity at ghost
	double sImagGhost=mSteadyGhost*aImagGhost+mImagGhost*aSteadyGhost;															//imaginary absolute velocity at ghost
	double uRealGhost=sRealGhost*cos(thetaIn);																							//real u at ghots
	double uImagGhost=sImagGhost*cos(thetaIn);																							//imaginary u at ghost
	double vRealGhost=sRealGhost*sin(thetaIn);																							//real v at ghost
	double vImagGhost=sImagGhost*sin(thetaIn);																							//imaginary v at ghost

	//harmonic conservative variables at ghsot cell
	wRealGhost[0]=rhoRealGhost;
	wImagGhost[0]=rhoImagGhost;

	wRealGhost[1]=rhoSteadyGhost*uRealGhost+rhoRealGhost*uSteadyGhost;
	wImagGhost[1]=rhoSteadyGhost*uImagGhost+rhoImagGhost*uSteadyGhost;

	wRealGhost[2]=rhoSteadyGhost*vRealGhost+rhoRealGhost*vSteadyGhost;
	wImagGhost[2]=rhoSteadyGhost*vImagGhost+rhoImagGhost*vSteadyGhost;

	wRealGhost[3]=pRealGhost/g1 + rhoRealGhost/2*(uSteadyGhost*uSteadyGhost+vSteadyGhost*vSteadyGhost) + rhoSteadyGhost*(uSteadyGhost*uRealGhost+vSteadyGhost*vRealGhost);
	wImagGhost[3]=pImagGhost/g1 + rhoImagGhost/2*(uSteadyGhost*uSteadyGhost+vSteadyGhost*vSteadyGhost) + rhoSteadyGhost*(uSteadyGhost*uImagGhost+vSteadyGhost*vImagGhost);
}


