#include "SolutionCalculator.h"

#include "Mesh.h"
#include "SteadyReader.h"
#include "../Solver/HarmonicSolver.h"
#include "../Inputs.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <valarray>
using namespace std;


SolutionCalculator::SolutionCalculator(Mesh*mp,SteadyReader*srp,HarmonicSolver*hsp) {
	double g1=gamma-1;
	int cellTotal=mp->GetCelltot();
	double** qSteady=srp->GetPrimitiveDomain();
	double** wHatReal=hsp->GetConservativeReal();
	double** wHatImag=hsp->GetConservativeImag();

	//solution steady
	double rhoSteady[cellTotal]={0};
	double uSteady[cellTotal]={0};
	double vSteady[cellTotal]={0};
	double pSteady[cellTotal]={0};
	double sSteady[cellTotal]={0};
	double aSteady[cellTotal]={0};
	double mSteady[cellTotal]={0};

	for(int i=0;i<cellTotal;i++){
		rhoSteady[i]=qSteady[i][0];
		uSteady[i]=qSteady[i][1];
		vSteady[i]=qSteady[i][2];
		pSteady[i]=qSteady[i][3];
		sSteady[i]=sqrt(uSteady[i]*uSteady[i]+vSteady[i]*vSteady[i]);
		aSteady[i]=sqrt(gamma*pSteady[i]/rhoSteady[i]);
		mSteady[i]=sSteady[i]/aSteady[i];
	}

	//solution int frequency domain in real and imgag format
	double rhoHatReal[cellTotal]={0};
	double rhoHatImag[cellTotal]={0};
	double uHatReal[cellTotal]={0};
	double uHatImag[cellTotal]={0};
	double vHatReal[cellTotal]={0};
	double vHatImag[cellTotal]={0};
	double pHatReal[cellTotal]={0};
	double pHatImag[cellTotal]={0};
	double sHatReal[cellTotal]={0};
	double sHatImag[cellTotal]={0};
	double aHatReal[cellTotal]={0};
	double aHatImag[cellTotal]={0};
	double mHatReal[cellTotal]={0};
	double mHatImag[cellTotal]={0};

	for(int i=0;i<cellTotal;i++){
		rhoHatReal[i]=wHatReal[i][0];
		rhoHatImag[i]=wHatImag[i][0];
		uHatReal[i]=(wHatReal[i][1]-rhoHatReal[i]*uSteady[i])/rhoSteady[i];
		uHatImag[i]=(wHatImag[i][1]-rhoHatImag[i]*uSteady[i])/rhoSteady[i];
		vHatReal[i]=(wHatReal[i][2]-rhoHatReal[i]*vSteady[i])/rhoSteady[i];
		vHatImag[i]=(wHatImag[i][2]-rhoHatImag[i]*vSteady[i])/rhoSteady[i];
		pHatReal[i]=g1*(wHatReal[i][3]-0.5*(2*(uSteady[i]*wHatReal[i][1]+vSteady[i]*wHatReal[i][2])-wHatReal[i][0]*(uSteady[i]*uSteady[i]+vSteady[i]*vSteady[i])));
		pHatImag[i]=g1*(wHatImag[i][3]-0.5*(2*(uSteady[i]*wHatImag[i][1]+vSteady[i]*wHatImag[i][2])-wHatImag[i][0]*(uSteady[i]*uSteady[i]+vSteady[i]*vSteady[i])));
		sHatReal[i]=(uSteady[i]*uHatReal[i]+vSteady[i]*vHatReal[i])/sSteady[i];
		sHatImag[i]=(uSteady[i]*uHatImag[i]+vSteady[i]*vHatImag[i])/sSteady[i];
		aHatReal[i]=gamma/2/aSteady[i]*(pHatReal[i]/rhoSteady[i]-pSteady[i]*rhoHatReal[i]/rhoSteady[i]/rhoSteady[i]);
		aHatReal[i]=gamma/2/aSteady[i]*(pHatImag[i]/rhoSteady[i]-pSteady[i]*rhoHatImag[i]/rhoSteady[i]/rhoSteady[i]);
		mHatReal[i]=sHatReal[i]/aSteady[i]-sSteady[i]*aHatReal[i]/aSteady[i]/aSteady[i];
		mHatImag[i]=sHatImag[i]/aSteady[i]-sSteady[i]*aHatImag[i]/aSteady[i]/aSteady[i];
	}

	//solution in frequeny domain in maginutude and phase angle
	amplitude.w1=new double[cellTotal];			//structure to store amplitudes
	amplitude.w2=new double[cellTotal];
	amplitude.w3=new double[cellTotal];
	amplitude.w4=new double[cellTotal];
	amplitude.rho=new double[cellTotal];
	amplitude.u=new double[cellTotal];
	amplitude.v=new double[cellTotal];
	amplitude.P=new double[cellTotal];
	amplitude.M=new double[cellTotal];

	for (int cn=0;cn<cellTotal;cn++){
		amplitude.w1[cn]=sqrt(wHatReal[cn][0]*wHatReal[cn][0]+wHatImag[cn][0]*wHatImag[cn][0]);
		amplitude.w2[cn]=sqrt(wHatReal[cn][1]*wHatReal[cn][1]+wHatImag[cn][1]*wHatImag[cn][1]);
		amplitude.w3[cn]=sqrt(wHatReal[cn][2]*wHatReal[cn][2]+wHatImag[cn][2]*wHatImag[cn][2]);
		amplitude.w4[cn]=sqrt(wHatReal[cn][3]*wHatReal[cn][3]+wHatImag[cn][3]*wHatImag[cn][3]);
		amplitude.rho[cn]=sqrt(rhoHatReal[cn]*rhoHatReal[cn]+rhoHatImag[cn]*rhoHatImag[cn]);
		amplitude.u[cn]=sqrt(uHatReal[cn]*uHatReal[cn]+uHatImag[cn]*uHatImag[cn]);
		amplitude.v[cn]=sqrt(vHatReal[cn]*vHatReal[cn]+vHatImag[cn]*vHatImag[cn]);
		amplitude.P[cn]=sqrt(pHatReal[cn]*pHatReal[cn]+pHatImag[cn]*pHatImag[cn]);
		amplitude.M[cn]=sqrt(mHatReal[cn]*mHatReal[cn]+mHatImag[cn]*mHatImag[cn]);
	}

	phi.w1=new double[cellTotal];			//structure to store phase angle
	phi.w2=new double[cellTotal];
	phi.w3=new double[cellTotal];
	phi.w4=new double[cellTotal];
	phi.rho=new double[cellTotal];
	phi.u=new double[cellTotal];
	phi.v=new double[cellTotal];
	phi.P=new double[cellTotal];
	phi.M=new double[cellTotal];

	for (int cn=0;cn<cellTotal;cn++){
		phi.w1[cn]=PhiCalculator(wHatReal[cn][0],wHatImag[cn][0]);
		phi.w2[cn]=PhiCalculator(wHatReal[cn][1],wHatImag[cn][1]);
		phi.w3[cn]=PhiCalculator(wHatReal[cn][2],wHatImag[cn][2]);
		phi.w4[cn]=PhiCalculator(wHatReal[cn][3],wHatImag[cn][3]);
		phi.rho[cn]=PhiCalculator(rhoHatReal[cn],rhoHatImag[cn]);
		phi.u[cn]=PhiCalculator(uHatReal[cn],uHatImag[cn]);
		phi.v[cn]=PhiCalculator(vHatReal[cn],vHatImag[cn]);
		phi.P[cn]=PhiCalculator(pHatReal[cn],pHatImag[cn]);
		phi.M[cn]=PhiCalculator(mHatReal[cn],mHatImag[cn]);
	}

	cout<<"\n"<<endl;
	cout<<"|wHat1|  w1Phi  |wHat2|  w2Phi  |wHat3|  w3Phi  |wHat4|  w4Phi  |rhoHat|  rhoPhi   |uHat|   uPhi   |vHat|   vPhi   |pHat|   pPhi   |mHat|   mPhi"<<endl;
	cout<<"------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
	for (int cn=0;cn<cellTotal;cn++){
		cout<<setw(6)<<fixed<<left<<cn<<
		setw(15)<<setprecision(8)<<fixed<<left<<amplitude.w1[cn]<<setw(15)<<setprecision(8)<<fixed<<left<<phi.w1[cn]<<
		setw(15)<<setprecision(8)<<fixed<<left<<amplitude.w2[cn]<<setw(15)<<setprecision(8)<<fixed<<left<<phi.w2[cn]<<
		setw(15)<<setprecision(8)<<fixed<<left<<amplitude.w3[cn]<<setw(15)<<setprecision(8)<<fixed<<left<<phi.w3[cn]<<
		setw(15)<<setprecision(8)<<fixed<<left<<amplitude.w4[cn]<<setw(15)<<setprecision(8)<<fixed<<left<<phi.w4[cn]<<
		setw(15)<<setprecision(8)<<fixed<<left<<amplitude.rho[cn]<<setw(15)<<setprecision(8)<<fixed<<left<<phi.rho[cn]<<
		setw(15)<<setprecision(8)<<fixed<<left<<amplitude.u[cn]<<setw(15)<<setprecision(8)<<fixed<<left<<phi.u[cn]<<
		setw(15)<<setprecision(8)<<fixed<<left<<amplitude.v[cn]<<setw(15)<<setprecision(8)<<fixed<<left<<phi.v[cn]<<
		setw(15)<<setprecision(8)<<fixed<<left<<amplitude.P[cn]<<setw(15)<<setprecision(8)<<fixed<<left<<phi.P[cn]<<
		setw(15)<<setprecision(8)<<fixed<<left<<amplitude.M[cn]<<setw(15)<<setprecision(8)<<fixed<<left<<phi.M[cn]<<endl;
	}
	cout<<"\n"<<endl;


	//convert solution into time domain
	double periodTime=1/redF;						//period of oscilation in seconds
	instantTotal=periodDevision+1;				//number of time instants to calculate solution at
	double radialFreq=2*pi/periodTime;				//radial frequency
	dt=periodTime/periodDevision;			//it includes time=0, time=inst1, time=inst2, ... and time =T
	double time[instantTotal]={0};					//instants in which solution is required 	0 <= time <= T
	for (int i=0;i<instantTotal;i++)
		time[i]=i*dt;


	harmonicSinus=new TimeSolution[instantTotal];		//array structure to store sinus harmonic solution at different times
	for (int t=0;t<instantTotal;t++){
		harmonicSinus[t].w1=new double[cellTotal];
		harmonicSinus[t].w2=new double[cellTotal];
		harmonicSinus[t].w3=new double[cellTotal];
		harmonicSinus[t].w4=new double[cellTotal];
		harmonicSinus[t].rho=new double[cellTotal];
		harmonicSinus[t].u=new double[cellTotal];
		harmonicSinus[t].v=new double[cellTotal];
		harmonicSinus[t].P=new double[cellTotal];
		harmonicSinus[t].M=new double[cellTotal];
	}

	harmonicCosinus=new TimeSolution[instantTotal];		//array structure to store cosinus harmonic solution at different times
	for (int t=0;t<instantTotal;t++){
		harmonicCosinus[t].w1=new double[cellTotal];
		harmonicCosinus[t].w2=new double[cellTotal];
		harmonicCosinus[t].w3=new double[cellTotal];
		harmonicCosinus[t].w4=new double[cellTotal];
		harmonicCosinus[t].rho=new double[cellTotal];
		harmonicCosinus[t].u=new double[cellTotal];
		harmonicCosinus[t].v=new double[cellTotal];
		harmonicCosinus[t].P=new double[cellTotal];
		harmonicCosinus[t].M=new double[cellTotal];
	}


	for (int t=0;t<instantTotal;t++){
		harmonicSinus[t].time=time[t];
		harmonicCosinus[t].time=time[t];
		for(int cn=0;cn<cellTotal;cn++){
			harmonicSinus[t].w1[cn]=amplitude.w1[cn]*sin(radialFreq*time[t]+phi.w1[cn]);
			harmonicSinus[t].w2[cn]=amplitude.w2[cn]*sin(radialFreq*time[t]+phi.w2[cn]);
			harmonicSinus[t].w3[cn]=amplitude.w3[cn]*sin(radialFreq*time[t]+phi.w3[cn]);
			harmonicSinus[t].w4[cn]=amplitude.w4[cn]*sin(radialFreq*time[t]+phi.w4[cn]);
			harmonicSinus[t].rho[cn]=amplitude.rho[cn]*sin(radialFreq*time[t]+phi.rho[cn]);
			harmonicSinus[t].u[cn]=amplitude.u[cn]*sin(radialFreq*time[t]+phi.u[cn]);
			harmonicSinus[t].v[cn]=amplitude.v[cn]*sin(radialFreq*time[t]+phi.v[cn]);
			harmonicSinus[t].P[cn]=amplitude.P[cn]*sin(radialFreq*time[t]+phi.P[cn]);
			harmonicSinus[t].M[cn]=amplitude.M[cn]*sin(radialFreq*time[t]+phi.M[cn]);
			harmonicCosinus[t].w1[cn]=amplitude.w1[cn]*cos(radialFreq*time[t]+phi.w1[cn]);
			harmonicCosinus[t].w2[cn]=amplitude.w2[cn]*cos(radialFreq*time[t]+phi.w2[cn]);
			harmonicCosinus[t].w3[cn]=amplitude.w3[cn]*cos(radialFreq*time[t]+phi.w3[cn]);
			harmonicCosinus[t].w4[cn]=amplitude.w4[cn]*cos(radialFreq*time[t]+phi.w4[cn]);
			harmonicCosinus[t].rho[cn]=amplitude.rho[cn]*cos(radialFreq*time[t]+phi.rho[cn]);
			harmonicCosinus[t].u[cn]=amplitude.u[cn]*cos(radialFreq*time[t]+phi.u[cn]);
			harmonicCosinus[t].v[cn]=amplitude.v[cn]*cos(radialFreq*time[t]+phi.v[cn]);
			harmonicCosinus[t].P[cn]=amplitude.P[cn]*cos(radialFreq*time[t]+phi.P[cn]);
			harmonicCosinus[t].M[cn]=amplitude.M[cn]*cos(radialFreq*time[t]+phi.M[cn]);
		}
	}


	//combine steady and harmonic to get unsteady
	unsteadySinus=new TimeSolution[instantTotal];		//array structure to store unsteady sinus solution at different times
	for (int t=0;t<instantTotal;t++){
		unsteadySinus[t].w1=new double[cellTotal];
		unsteadySinus[t].w2=new double[cellTotal];
		unsteadySinus[t].w3=new double[cellTotal];
		unsteadySinus[t].w4=new double[cellTotal];
		unsteadySinus[t].rho=new double[cellTotal];
		unsteadySinus[t].u=new double[cellTotal];
		unsteadySinus[t].v=new double[cellTotal];
		unsteadySinus[t].P=new double[cellTotal];
		unsteadySinus[t].M=new double[cellTotal];
	}

	unsteadyCosinus=new TimeSolution[instantTotal];		//array structure to store unsteady cosinus solution at different times
	for (int t=0;t<instantTotal;t++){
		unsteadyCosinus[t].w1=new double[cellTotal];
		unsteadyCosinus[t].w2=new double[cellTotal];
		unsteadyCosinus[t].w3=new double[cellTotal];
		unsteadyCosinus[t].w4=new double[cellTotal];
		unsteadyCosinus[t].rho=new double[cellTotal];
		unsteadyCosinus[t].u=new double[cellTotal];
		unsteadyCosinus[t].v=new double[cellTotal];
		unsteadyCosinus[t].P=new double[cellTotal];
		unsteadyCosinus[t].M=new double[cellTotal];
	}


	for (int t=0;t<instantTotal;t++){
		unsteadySinus[t].time=time[t];
		unsteadyCosinus[t].time=time[t];
		for(int cn=0;cn<cellTotal;cn++){
			unsteadySinus[t].w1[cn]=rhoSteady[cn] + harmonicSinus[t].w1[cn];
			unsteadySinus[t].w2[cn]=rhoSteady[cn]*uSteady[cn] + harmonicSinus[t].w2[cn];
			unsteadySinus[t].w3[cn]=rhoSteady[cn]*vSteady[cn] + harmonicSinus[t].w3[cn];
			unsteadySinus[t].w4[cn]=(pSteady[cn]/(gamma-1) + rhoSteady[cn]*sSteady[cn]*sSteady[cn]) + harmonicSinus[t].w4[cn];
			unsteadySinus[t].rho[cn]=rhoSteady[cn] + harmonicSinus[t].rho[cn];
			unsteadySinus[t].u[cn]=uSteady[cn] + harmonicSinus[t].u[cn];
			unsteadySinus[t].v[cn]=vSteady[cn] + harmonicSinus[t].v[cn];
			unsteadySinus[t].P[cn]=pSteady[cn] + harmonicSinus[t].P[cn];
			unsteadySinus[t].M[cn]=mSteady[cn] + harmonicSinus[t].M[cn];
			unsteadyCosinus[t].w1[cn]=rhoSteady[cn] + harmonicCosinus[t].w1[cn];
			unsteadyCosinus[t].w2[cn]=rhoSteady[cn]*uSteady[cn] + harmonicCosinus[t].w2[cn];
			unsteadyCosinus[t].w3[cn]=rhoSteady[cn]*vSteady[cn] + harmonicCosinus[t].w3[cn];
			unsteadyCosinus[t].w4[cn]=(pSteady[cn]/(gamma-1) + rhoSteady[cn]*sSteady[cn]*sSteady[cn]) + harmonicCosinus[t].w4[cn];
			unsteadyCosinus[t].rho[cn]=rhoSteady[cn] + harmonicCosinus[t].rho[cn];
			unsteadyCosinus[t].u[cn]=uSteady[cn] + harmonicCosinus[t].u[cn];
			unsteadyCosinus[t].v[cn]=vSteady[cn] + harmonicCosinus[t].v[cn];
			unsteadyCosinus[t].P[cn]=pSteady[cn] + harmonicCosinus[t].P[cn];
			unsteadyCosinus[t].M[cn]=mSteady[cn] + harmonicCosinus[t].M[cn];
		}
	}


}

double SolutionCalculator::PhiCalculator(double real,double imag){
	double phi;

	if (real==0){
			if(imag>0){
				phi=pi/2;
				return phi;
			}
			else if(imag<0){
				phi=3/2*pi;
				return phi;
			}
			else{		//imag=0
				phi=0;
				cout<<"Warning: both real and imag values got zero !"<<endl;
				return phi;
			}
	}

	if (imag==0){
			if (real>0){
				phi=0;
				return phi;
			}
			else{		//real<0
				phi=pi;
				return phi;
			}
	}

	if(real>0)
		phi=atan(imag/real);
	else
		phi=pi+atan(imag/real);

	return phi;
}



Frequency* SolutionCalculator::GetAmplitude(){
	return &amplitude;
}

Frequency* SolutionCalculator::GetPhi(){
	return &phi;
}

TimeSolution* SolutionCalculator::GetHarmonicSinus(){
	return harmonicSinus;
}

TimeSolution* SolutionCalculator::GetHarmonicCosinus(){
	return harmonicCosinus;
}

TimeSolution* SolutionCalculator::GetUnsteadySinus(){
	return unsteadySinus;
}

TimeSolution* SolutionCalculator::GetUnsteadyCosinus(){
	return unsteadyCosinus;
}

int SolutionCalculator::GetInstantNumber(){
	return instantTotal;
}

double SolutionCalculator::GetDt(){
	return dt;
}
