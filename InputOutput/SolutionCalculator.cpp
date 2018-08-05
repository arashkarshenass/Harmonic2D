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


	harmonic=new TimeSolution[instantTotal];		//array structure to store harmonic solution at different times
	for (int t=0;t<instantTotal;t++){
		harmonic[t].w1=new double[cellTotal];
		harmonic[t].w2=new double[cellTotal];
		harmonic[t].w3=new double[cellTotal];
		harmonic[t].w4=new double[cellTotal];
		harmonic[t].rho=new double[cellTotal];
		harmonic[t].u=new double[cellTotal];
		harmonic[t].v=new double[cellTotal];
		harmonic[t].P=new double[cellTotal];
		harmonic[t].M=new double[cellTotal];
	}

	for (int t=0;t<instantTotal;t++){
		harmonic[t].time=time[t];
		if(harmonicType==1)
			for(int cn=0;cn<cellTotal;cn++){
				harmonic[t].w1[cn]=amplitude.w1[cn]*sin(radialFreq*time[t]+phi.w1[cn]);
				harmonic[t].w2[cn]=amplitude.w2[cn]*sin(radialFreq*time[t]+phi.w2[cn]);
				harmonic[t].w3[cn]=amplitude.w3[cn]*sin(radialFreq*time[t]+phi.w3[cn]);
				harmonic[t].w4[cn]=amplitude.w4[cn]*sin(radialFreq*time[t]+phi.w4[cn]);
				harmonic[t].rho[cn]=amplitude.rho[cn]*sin(radialFreq*time[t]+phi.rho[cn]);
				harmonic[t].u[cn]=amplitude.u[cn]*sin(radialFreq*time[t]+phi.u[cn]);
				harmonic[t].v[cn]=amplitude.v[cn]*sin(radialFreq*time[t]+phi.v[cn]);
				harmonic[t].P[cn]=amplitude.P[cn]*sin(radialFreq*time[t]+phi.P[cn]);
				harmonic[t].M[cn]=amplitude.M[cn]*sin(radialFreq*time[t]+phi.M[cn]);
			}
		else if(harmonicType==2)
			for(int cn=0;cn<cellTotal;cn++){
				harmonic[t].w1[cn]=amplitude.w1[cn]*cos(radialFreq*time[t]+phi.w1[cn]);
				harmonic[t].w2[cn]=amplitude.w2[cn]*cos(radialFreq*time[t]+phi.w2[cn]);
				harmonic[t].w3[cn]=amplitude.w3[cn]*cos(radialFreq*time[t]+phi.w3[cn]);
				harmonic[t].w4[cn]=amplitude.w4[cn]*cos(radialFreq*time[t]+phi.w4[cn]);
				harmonic[t].rho[cn]=amplitude.rho[cn]*cos(radialFreq*time[t]+phi.rho[cn]);
				harmonic[t].u[cn]=amplitude.u[cn]*cos(radialFreq*time[t]+phi.u[cn]);
				harmonic[t].v[cn]=amplitude.v[cn]*cos(radialFreq*time[t]+phi.v[cn]);
				harmonic[t].P[cn]=amplitude.P[cn]*cos(radialFreq*time[t]+phi.P[cn]);
				harmonic[t].M[cn]=amplitude.M[cn]*cos(radialFreq*time[t]+phi.M[cn]);
			}
	}


	//combine steady and harmonic to get unsteady
	unsteady=new TimeSolution[instantTotal];		//array structure to store unsteady solution at different times
	for (int t=0;t<instantTotal;t++){
		unsteady[t].w1=new double[cellTotal];
		unsteady[t].w2=new double[cellTotal];
		unsteady[t].w3=new double[cellTotal];
		unsteady[t].w4=new double[cellTotal];
		unsteady[t].rho=new double[cellTotal];
		unsteady[t].u=new double[cellTotal];
		unsteady[t].v=new double[cellTotal];
		unsteady[t].P=new double[cellTotal];
		unsteady[t].M=new double[cellTotal];
	}

	for (int t=0;t<instantTotal;t++){
		unsteady[t].time=time[t];
		for(int cn=0;cn<cellTotal;cn++){
			unsteady[t].w1[cn]=rhoSteady[cn] + harmonic[t].w1[cn];
			unsteady[t].w2[cn]=rhoSteady[cn]*uSteady[cn] + harmonic[t].w2[cn];
			unsteady[t].w3[cn]=rhoSteady[cn]*vSteady[cn] + harmonic[t].w3[cn];
			unsteady[t].w4[cn]=(pSteady[cn]/(gamma-1) + rhoSteady[cn]*sSteady[cn]*sSteady[cn]) + harmonic[t].w4[cn];
			unsteady[t].rho[cn]=rhoSteady[cn] + harmonic[t].rho[cn];
			unsteady[t].u[cn]=uSteady[cn] + harmonic[t].u[cn];
			unsteady[t].v[cn]=vSteady[cn] + harmonic[t].v[cn];
			unsteady[t].P[cn]=pSteady[cn] + harmonic[t].P[cn];
			unsteady[t].M[cn]=mSteady[cn] + harmonic[t].M[cn];
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

TimeSolution* SolutionCalculator::GetHarmonic(){
	return harmonic;
}

TimeSolution* SolutionCalculator::GetUnsteady(){
	return unsteady;
}

int SolutionCalculator::GetInstantNumber(){
	return instantTotal;
}

double SolutionCalculator::GetDt(){
	return dt;
}
