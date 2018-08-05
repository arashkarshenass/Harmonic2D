#include "TecPlotEr.h"

#include "Mesh.h"
#include "SolutionCalculator.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

TecPlotEr::TecPlotEr(Mesh* mp,SolutionCalculator* scp){
	int celltot=mp->GetCelltot();
	int nodetot=mp->GetNodetot();
	double* x=mp->GetX();
	double* y=mp->GetY();
	int** cellmap=mp->GetCellmap();
	Frequency* amplitude=scp->GetAmplitude();
	Frequency* phi=scp->GetPhi();
	TimeSolution* harmonic=scp->GetHarmonic();
	TimeSolution* unsteady=scp->GetUnsteady();
	int instantNumber=scp->GetInstantNumber();
	double dt=scp->GetDt();


	puts("\nSaving solution files with TecPlot Format ...");

	//--------------------------------------------------writing frequency solution
	ofstream FrequencySolutionFile;
	FrequencySolutionFile.open("FrequencySolution.dat");
	FrequencySolutionFile<<"TITLE=\"FrequencySolution\"\nFILETYPE=FULL"<<endl;
	FrequencySolutionFile<<"VARIABLES=\"X\" \"Y\" \"Amplitude w1\" \"Amplitude w2\" \"Amplitude w3\" \"Amplitude w4\" \"Amplitude rho\" \"Amplitude u\" \"Amplitude v\" \"Amplitude P\" \"Amplitude M\" \"Phi w1\" \"Phi w2\" \"Phi w3\" \"Phi w4\" \"Phi rho\" \"Phi u\" \"Phi v\" \"Phi P\" \"Phi M\""<<endl;
	FrequencySolutionFile<<"ZONE\nT=Interior\nNODES="<<nodetot<<endl;
	FrequencySolutionFile<<"ELEMENTS="<<celltot<<endl;
	FrequencySolutionFile<<"DATAPACKING=BLOCK\nZONETYPE=FEQUADRILATERAL"<<endl;
	FrequencySolutionFile<<"VARLOCATION=([3-20]=CELLCENTERED)"<<endl;
	FrequencySolutionFile<<setprecision(6)<<fixed<<scientific;

	FrequencySolutionFile<<"#x-coord of vertices"<<endl;
	for(int nn=0;nn<nodetot;nn++)
		FrequencySolutionFile<<x[nn]<<endl;

	FrequencySolutionFile<<"#y-coord of vertices"<<endl;
	for(int nn=0;nn<nodetot;nn++)
		FrequencySolutionFile<<y[nn]<<endl;

	FrequencySolutionFile<<"#w1 Amplitude of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<amplitude->w1[cn]<<endl;

	FrequencySolutionFile<<"#w2 Amplitude of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<amplitude->w2[cn]<<endl;

	FrequencySolutionFile<<"#w3 Amplitude of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<amplitude->w3[cn]<<endl;

	FrequencySolutionFile<<"#w4 Amplitude of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<amplitude->w4[cn]<<endl;

	FrequencySolutionFile<<"#Density Amplitude of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<amplitude->rho[cn]<<endl;

	FrequencySolutionFile<<"#X-Velocity Amplitude of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<amplitude->u[cn]<<endl;

	FrequencySolutionFile<<"#Y-Velocity Amplitude of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<amplitude->v[cn]<<endl;

	FrequencySolutionFile<<"#Static Pressure Amplitude of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<amplitude->P[cn]<<endl;

	FrequencySolutionFile<<"#Mach Number Amplitude of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<amplitude->M[cn]<<endl;

	FrequencySolutionFile<<"#w1 Phase Angle of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<phi->w1[cn]<<endl;

	FrequencySolutionFile<<"#w2 Phase Angle of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<phi->w2[cn]<<endl;

	FrequencySolutionFile<<"#w3 Phase Angle of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<phi->w3[cn]<<endl;

	FrequencySolutionFile<<"#w4 Phase Angle of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<phi->w4[cn]<<endl;

	FrequencySolutionFile<<"#Density Phase Angle of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<phi->rho[cn]<<endl;

	FrequencySolutionFile<<"#X-Velocity Phase Angle of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<phi->u[cn]<<endl;

	FrequencySolutionFile<<"#Y-Velocity Phase Angle of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<phi->v[cn]<<endl;

	FrequencySolutionFile<<"#Static Pressure Phase Angle of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<phi->P[cn]<<endl;

	FrequencySolutionFile<<"#Mach Number Phase Angle of cell centers"<<endl;
	for(int cn=0;cn<celltot;cn++)
		FrequencySolutionFile<<phi->M[cn]<<endl;

	FrequencySolutionFile<<"#Cellmap of cells"<<endl;
	for(int i=0;i<celltot;i++)
		FrequencySolutionFile<<cellmap[i][0]+1<<" "<<cellmap[i][1]+1<<" "<<cellmap[i][2]+1<<" "<<cellmap[i][3]+1<<endl;

	FrequencySolutionFile.close();
	cout<<"Harmonic Solution in frequency domain is saved successfully (1 file)"<<endl;

	//----------------------------------------------writing harmonic solution in time
	char fileName[200];
	ofstream TimeSolutionFile;
	for (int tn=0;tn<instantNumber;tn++){
		snprintf(fileName,200,"Harmonic_Solution_at_t=%4.2f.dat",dt*tn);
		TimeSolutionFile.open(fileName);
		TimeSolutionFile<<"TITLE=\"Harmonic Solution at time= "<<dt*tn<<" seconds\""<<endl;
		TimeSolutionFile<<"FILETYPE=FULL"<<endl;
		TimeSolutionFile<<"VARIABLES=\"X\" \"Y\" \"w1\" \"w2\" \"w3\" \"w4\" \"rho\" \"u\" \"v\" \"P\" \"M\""<<endl;
		TimeSolutionFile<<"ZONE\nT=Interior\nNODES="<<nodetot<<endl;
		TimeSolutionFile<<"ELEMENTS="<<celltot<<endl;
		TimeSolutionFile<<"DATAPACKING=BLOCK\nZONETYPE=FEQUADRILATERAL"<<endl;
		TimeSolutionFile<<"VARLOCATION=([3-11]=CELLCENTERED)"<<endl;
		TimeSolutionFile<<setprecision(6)<<fixed<<scientific;

		TimeSolutionFile<<"#x-coord of vertices"<<endl;
		for(int nn=0;nn<nodetot;nn++)
			TimeSolutionFile<<x[nn]<<endl;

		TimeSolutionFile<<"#y-coord of vertices"<<endl;
		for(int nn=0;nn<nodetot;nn++)
			TimeSolutionFile<<y[nn]<<endl;

		TimeSolutionFile<<"#w1 Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<harmonic[tn].w1[cn]<<endl;

		TimeSolutionFile<<"#w2 Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<harmonic[tn].w2[cn]<<endl;

		TimeSolutionFile<<"#w3 Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<harmonic[tn].w3[cn]<<endl;

		TimeSolutionFile<<"#w4 Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<harmonic[tn].w4[cn]<<endl;

		TimeSolutionFile<<"#Density Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<harmonic[tn].rho[cn]<<endl;

		TimeSolutionFile<<"#X-Velocity Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<harmonic[tn].u[cn]<<endl;

		TimeSolutionFile<<"#Y-Velocity Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<harmonic[tn].v[cn]<<endl;

		TimeSolutionFile<<"#Static Pressure Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<harmonic[tn].P[cn]<<endl;

		TimeSolutionFile<<"#Mach Number Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<harmonic[tn].M[cn]<<endl;

		TimeSolutionFile<<"#Cellmap of cells"<<endl;
		for(int i=0;i<celltot;i++)
			TimeSolutionFile<<cellmap[i][0]+1<<" "<<cellmap[i][1]+1<<" "<<cellmap[i][2]+1<<" "<<cellmap[i][3]+1<<endl;

		TimeSolutionFile.close();
	}	//enf of for loop

	cout<<"Harmonic solution in time domain is saved successfully ("<<instantNumber<<" files)"<<endl;


	//------------------------------------------writing unsteady solution in time
	for (int tn=0;tn<instantNumber;tn++){
		snprintf(fileName,200,"Unsteady_Solution_at_t=%4.2f.dat",dt*tn);
		TimeSolutionFile.open(fileName);
		TimeSolutionFile<<"TITLE=\"Unsteady Solution at time= "<<dt*tn<<" seconds\""<<endl;
		TimeSolutionFile<<"FILETYPE=FULL"<<endl;
		TimeSolutionFile<<"VARIABLES=\"X\" \"Y\" \"w1\" \"w2\" \"w3\" \"w4\" \"rho\" \"u\" \"v\" \"P\" \"M\""<<endl;
		TimeSolutionFile<<"ZONE\nT=Interior\nNODES="<<nodetot<<endl;
		TimeSolutionFile<<"ELEMENTS="<<celltot<<endl;
		TimeSolutionFile<<"DATAPACKING=BLOCK\nZONETYPE=FEQUADRILATERAL"<<endl;
		TimeSolutionFile<<"VARLOCATION=([3-11]=CELLCENTERED)"<<endl;
		TimeSolutionFile<<setprecision(6)<<fixed<<scientific;

		TimeSolutionFile<<"#x-coord of vertices"<<endl;
		for(int nn=0;nn<nodetot;nn++)
			TimeSolutionFile<<x[nn]<<endl;

		TimeSolutionFile<<"#y-coord of vertices"<<endl;
		for(int nn=0;nn<nodetot;nn++)
			TimeSolutionFile<<y[nn]<<endl;

		TimeSolutionFile<<"#w1 Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<unsteady[tn].w1[cn]<<endl;

		TimeSolutionFile<<"#w2 Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<unsteady[tn].w2[cn]<<endl;

		TimeSolutionFile<<"#w3 Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<unsteady[tn].w3[cn]<<endl;

		TimeSolutionFile<<"#w4 Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<unsteady[tn].w4[cn]<<endl;

		TimeSolutionFile<<"#Density Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<unsteady[tn].rho[cn]<<endl;

		TimeSolutionFile<<"#X-Velocity Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<unsteady[tn].u[cn]<<endl;

		TimeSolutionFile<<"#Y-Velocity Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<unsteady[tn].v[cn]<<endl;

		TimeSolutionFile<<"#Static Pressure Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<unsteady[tn].P[cn]<<endl;

		TimeSolutionFile<<"#Mach Number Harmonic of cell centers"<<endl;
		for(int cn=0;cn<celltot;cn++)
			TimeSolutionFile<<unsteady[tn].M[cn]<<endl;

		TimeSolutionFile<<"#Cellmap of cells"<<endl;
		for(int i=0;i<celltot;i++)
			TimeSolutionFile<<cellmap[i][0]+1<<" "<<cellmap[i][1]+1<<" "<<cellmap[i][2]+1<<" "<<cellmap[i][3]+1<<endl;

		TimeSolutionFile.close();
		}	//enf of for loop
		cout<<"Unsteady solution in time domain is saved successfully ("<<instantNumber<<" files)"<<endl;

}
