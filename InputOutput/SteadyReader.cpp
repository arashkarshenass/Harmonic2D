/*======================================
 	   Steady Solution Class Source
=======================================*/

#include "SteadyReader.h"

#include "Mesh.h"
#include "../Inputs.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <string>
using namespace std;

SteadyReader::SteadyReader(char pub_name_in_func[],Mesh* mp){

	int celltot=mp->GetCelltot();		//get total cell number from mesh
	int nodetot=mp->GetNodetot();		//get total node number from mesh
	int edgeNum=mp->GetEdgeNumber();	//get edgenumber from mesh
	double** nxMat=mp->GetNx();			//get Nx from mesh
	double** nyMat=mp->GetNy();			//get Ny from mesh

	qsDom=new double*[celltot];			//steady primitive variables of domain cells
	for (int i=0;i<celltot;i++)
		qsDom[i]=new double[4];

	UnMatDom=new double*[celltot];			//Un of faces
	UtMatDom=new double*[celltot];			//Ut of faces
	for (int i=0;i<celltot;i++){
		UnMatDom[i]=new double[edgeNum];
		UtMatDom[i]=new double[edgeNum];
	}

//----------------------------------------------Openning file-----------------------------------
	char fileName[100];
	for(int i=0;i<100;i++)											//get steady solution file name
			fileName[i]=pub_name_in_func[i];

	FILE *myFile;														//Reading steady solution file
	myFile=fopen(fileName,"r");

	while(myFile==NULL){										//file is found or NOT?{
			printf("Steady Solution file cannot be found!!\n");
			printf("Enter file name: ");
			cin>>fileName;   										//or you can use cin>>filename;
			myFile=fopen(fileName,"r");
	}

	//----------------------------------------------Reading file-----------------------------------
	for(int i=0;i<10;i++)								//skip initial headers
		fscanf(myFile,"%*[^\n]\n",nullptr);

	for (int i=0;i<nodetot+1;i++)						//skip x-coordinates
		fscanf(myFile,"%*[^\n]\n",nullptr);

	for (int i=0;i<nodetot+1;i++)						//skip y-coordinates
		fscanf(myFile,"%*[^\n]\n",nullptr);

	fscanf(myFile,"%*[^\n]\n",nullptr);					//skip rho header

	for(int i=0;i<celltot;i++)							//reading density
		fscanf(myFile,"%lf\n",*(qsDom+i)+0);

	fscanf(myFile,"%*[^\n]\n",nullptr);					//skip u header

	for(int i=0;i<celltot;i++)							//reading u
		fscanf(myFile,"%lf\n",*(qsDom+i)+1);

	fscanf(myFile,"%*[^\n]\n",nullptr);					//skip v header

	for(int i=0;i<celltot;i++)							//reading v
		fscanf(myFile,"%lf\n",*(qsDom+i)+2);

	fscanf(myFile,"%*[^\n]\n",nullptr);					//skip p header

	for(int i=0;i<celltot;i++)							//reading p
		fscanf(myFile,"%lf\n",*(qsDom+i)+3);

	fclose(myFile);
	cout<<"\nSteady Solution read successfully"<<endl;

	//convert u and v into Un and Ut and store them
	for (int cn=0;cn<celltot;cn++)
		for(int fn=0;fn<edgeNum;fn++){
			UnMatDom[cn][fn]=qsDom[cn][1]*nxMat[cn][fn]+qsDom[cn][2]*nyMat[cn][fn];
			UtMatDom[cn][fn]=-qsDom[cn][1]*nyMat[cn][fn]+qsDom[cn][2]*nxMat[cn][fn];
	}

	if(prntstd==1){													//printing steady solution (primitives and conservatives)
		cout<<"\n   Steady Solution results (Primitive values"<<endl;
		cout<<"#       ro        u         v         P"<<endl;
		cout<<"-------------------------------------------------------------"<<endl;
		for(int i=0;i<celltot;i++){
			cout<<left<<setw(7)<<i<<left<<setw(11)<<setprecision(5)<<fixed<<qsDom[i][0]
			<<left<<setw(11)<<setprecision(4)<<fixed<<qsDom[i][1]<<left<<setw(11)<<
			setprecision(4)<<fixed<<qsDom[i][2]<<left<<setw(10)<<setprecision(1)<<
			fixed<<qsDom[i][3]<<endl;
		}
	}
}

double** SteadyReader::GetPrimitiveDomain() const{
	return qsDom;
}

double** SteadyReader::GetUnDom() const{
	return UnMatDom;
}

double** SteadyReader::GetUtDom() const{
	return UtMatDom;
}





