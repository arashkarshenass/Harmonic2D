			/*=======================================
 	   Harmonic Solver Class Source
========================================*/

#include "HarmonicSolver.h"
#include "../InputOutput/Mesh.h"
#include "../InputOutput/SteadyReader.h"
#include "SteadyJacobian.h"
#include "../Accelerator/LocalTimeStep.h"
#include "OneCell.h"
#include "Residue.h"
#include "../Utility/Utility.h"
#include "../Inputs.h"

#include <iostream>
#include <iomanip>
#include <cmath>


using namespace std;

HarmonicSolver::HarmonicSolver(Mesh* mp,SteadyReader* srp,SteadyJacobian* sjp,Utility* upTemp){

	up=upTemp;										//utility pointer
	celltot=mp->GetCelltot();						//total cell number
	vol=mp->GetArea();								//cell volume
	dt=new double[celltot];											//allocate dt vector (time steps from local time step class)
	wr=new double*[celltot];										//allocate real conservative perturbed variables inside cells
	wi=new double*[celltot];										//allocate img conservative perturbed variables inside cells
	wrOld=new double*[celltot];									//allocate old real conservative perturbed variables inside cells (from previous iteration)
	wiOld=new double*[celltot];								//allocate old img conservative perturbed variables inside cells (from previous iteration)
	dwr=new double*[celltot];										//allocate change of real conservative variables inside cells
	dwi=new double*[celltot];									//allocate change of img conservative variables inside cells
	cellFluxReal=new double[4];										//allocate total real fluxes passing a cell
	cellFluxImag=new double[4];										//allocate total img fluxes passing a cell
	resR=new double[4];												//allocate maximum real residues of cells during an iteration
	resIm=new double[4];											//allocate maximum img residues of cells during an iteration

	for(int i=0;i<celltot;i++){
		wr[i]=new double[4];
		wi[i]=new double[4];
		wrOld[i]=new double[4];
		wiOld[i]=new double[4];
		dwr[i]=new double[4];
		dwi[i]=new double[4];
	}

	up->MatrixInitializer(wrOld,celltot,4,0);					//set vectors and matrices to zero
	up->MatrixInitializer(wiOld,celltot,4,0);
	up->MatrixInitializer(dwr,celltot,4,0);
	up->MatrixInitializer(dwi,celltot,4,0);
	up->VectorInitializer(cellFluxReal,4,0);
	up->VectorInitializer(cellFluxImag,4,0);
	up->VectorInitializer(resR,4,0);
	up->VectorInitializer(resIm,4,0);

	for(int i=0;i<celltot;i++){										//initializing real and img conservative variables
		wr[i][0]=wr1;			wi[i][0]=wi1;
		wr[i][1]=wr2;			wi[i][1]=wi2;
		wr[i][2]=wr3;			wi[i][2]=wi3;
		wr[i][3]=wr4;			wi[i][3]=wi4;
	}

	cout<<"\nSolution is initialized as below for all cells"<<endl;
	cout<<"W1Real    W1Imag    W2Real   W2Imag     W3Real    W3Imag    W4Real     W4Imag"<<endl;
	cout<<"--------------------------------------------------------------------------------------------"<<endl;
	cout<<wr1<<"     "<<wi1<<"     "<<wr2<<"     "<<wi2<<"     "<<wr3<<"     "<<wi3<<"     "<<wr4<<"     "<<wi4<<endl;

	omega=2*pi*F;
	cout<<"\nReduced frequency is "<<F<<" Hz"<<endl;
	cout<<"\nRedial frequency is "<<omega<<" rad/sec"<<endl;



	LocalTimeStep lts(mp,srp,dt);															//create local time step object and calculate time step

	OneCell oneCell(cellFluxReal,cellFluxImag,wr,wi,mp,srp,sjp);					//create cellFlux object with inputs to its structor
	ocp=&oneCell;																				//address of cellFlux object

	Residue residue(dwr,dwi,resR,resIm,up,celltot);								//create residue object and send addresses to it
	rp=&residue;																					//address of residue object

	Itterator();

	cout<<"\nIteration is finished with "<<itr<<" cycles."<<endl;;
	cout<<"\nCell#"<<"   W1_R    "<<"    W1_I    "<<"   W2_R    "<<"   W2_I    "<<"   W3_R    "<<"   W3_I    "<<"   W4_R    "<<"   W4_I    "<<endl;
	cout<<"--------------------------------------------------------------------------------------------------------------------------------"<<endl;
	for(int i=0;i<celltot;i++){
		cout<<setw(7)<<fixed<<left<<i<<
		setw(18)<<setprecision(8)<<fixed<<left<<wr[i][0]<<setw(18)<<setprecision(8)<<fixed<<left<<wi[i][0]<<
		setw(18)<<setprecision(8)<<fixed<<left<<wr[i][1]<<setw(18)<<setprecision(8)<<fixed<<left<<wi[i][1]<<
		setw(18)<<setprecision(8)<<fixed<<left<<wr[i][2]<<setw(18)<<setprecision(8)<<fixed<<left<<wi[i][2]<<
		setw(18)<<setprecision(8)<<fixed<<left<<wr[i][3]<<setw(18)<<setprecision(8)<<fixed<<left<<wi[i][3]<<endl;
	}
}

void HarmonicSolver::Itterator(){
	cout<<"\nIteration is started."<<endl;
	cout<<"\n-------------------------------------Resuduals---------------------------------"<<endl;
	cout<<"#------W1Real------W1Imag------W2Real------W2Imag------W3Real------W3Imag------W4Real------W4Imag"<<endl;
	double alpha[]={0.666,0.666,1};
	double maxResidue=100;
	itr=0;

	while (maxResidue>resRef){									//psudo-time loop
		 itr+=1;
		 up->MatrixEqualizer(wrOld,wr,celltot,4);				//set w real old to w real at the begining of a time iteration
		 up->MatrixEqualizer(wiOld,wi,celltot,4);				//set w imag old to w imag at the begining of a time iteration

		 for(int rk=0;rk<3;rk++){								//runge-kutta loop
			up->MatrixInitializer(dwr,celltot,4,0);			//set dw real to zero at the begining of each r-k iteration
			up->MatrixInitializer(dwi,celltot,4,0);				//set dw imag to zero at the begining of each r-k iteration

			for(int cn=0;cn<celltot;cn++){						//cell loop
				up->VectorInitializer(cellFluxReal,4,0);			//set cellFlux real to zero before calculating real fluxes of a cell
				up->VectorInitializer(cellFluxImag,4,0);			//set cellFlux imag to zero before calculating imag fluxes of a cell
				ocp->ClaculateCellFlux(cn);						//calculating real and imag fluxes of a cell

				for (int j=0;j<4;j++){								//updating conservative variables
					dwr[cn][j]=(omega*wi[cn][j]-cellFluxReal[j]/vol[cn])*dt[cn];			//variation in wreal[i] of cell cn
					dwi[cn][j]=(-omega*wr[cn][j]-cellFluxImag[j]/vol[cn])*dt[cn];			//variation in wimag[i] of cell cn
					wr[cn][j]=wrOld[cn][j]+alpha[rk]*dwr[cn][j];			//updating wreal cell cn
					wi[cn][j]=wiOld[cn][j]+alpha[rk]*dwi[cn][j];			//updatong wimag cell cn
				}	//end up variable update loop
			}	//end of cell loop
		 }	//end of rk loop


		 switch(itr){												//prepare scaling of residues at the first iteration and print them
		 case 1:
		 	  	switch(scaleRes){
		 	  	case 1:
		 	  	  		rp->FirstIteration();
		 	  	  		break;
		 	  	}
		 	  	break;
		 }

		 if(itr%nfres==0 && itr!=1){
			rp->CalculateResidue(itr);			//calculate scaled residues and printing them
			maxResidue=rp->findMax();			//find maximum residue for time loop control
		 }
	}	//end of time while loop
} //end of function

double** HarmonicSolver::GetConservativeReal() const{
	return wr;
}

double** HarmonicSolver::GetConservativeImag() const{
	return wi;
}
