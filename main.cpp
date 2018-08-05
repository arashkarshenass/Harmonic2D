/*
 ==================================================
          Linearized 2D Euler Harmonic NT
 ==================================================
   2D Euler Linearized Harmonic Flow Solver
   	   	   	   based on unstructured Grids
                  developed by
               Arash Karshenass
       Middle East Technical University
                 Ankara- TURKEY
 ================================================== */


#include "InputOutput/Mesh.h"
#include "InputOutput/SteadyReader.h"
#include "Solver/SteadyJacobian.h"
#include "Solver/HarmonicSolver.h"
#include "InputOutput/SolutionCalculator.h"
#include "InputOutput/TecPlotEr.h"
#include "Utility/Utility.h"
#include <cstdio>
#include <iostream>



using namespace std;

//--------Inputs----------
//through input header file


int main()
{
	cout<<"------Program is started-----"<<endl;
	char meshName[]="NASA_2D_C++.su2";     		//name of mesh file
	char steadyName[]="Code_NT_AUSM.dat";			//name of steady solutin file

	Utility utility;

	Mesh mesh(meshName);									//create object mesh with constructor to read mesh file
	mesh.Geometry(&utility);										//calculates mesh geometry details
	mesh.Connectivity();											//calculate connectivity
	mesh.GhostCellStructure();
	SteadyReader readObj(steadyName,&mesh);				//create object std with constructor to read steady results
	SteadyJacobian jacObj(&mesh,&readObj,&utility);
	HarmonicSolver hrmslvObj(&mesh,&readObj,&jacObj,&utility);		//create object hrmnSlv with constructor to prepare harmonic solver
	SolutionCalculator finalSolution(&mesh,&readObj,&hrmslvObj);
	TecPlotEr tecplotSave(&mesh,&finalSolution);
	return 0;
}
