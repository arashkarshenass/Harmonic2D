/*=================================================
                Residue Class source
=================================================*/
#include "Residue.h"
#include "../Utility/Utility.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
using namespace std;

Residue::Residue(double** dwRTemp,double**dwImTemp,double* resRTemp,double* resImTemp,Utility* upTemp,int celltotTemp){
	celltot=celltotTemp;
	dwR=dwRTemp;
	dwIm=dwImTemp;
	resR=resRTemp;
	resIm=resImTemp;
	dwRFirst=new double[4];
	dwImFirst=new double[4];
	up=upTemp;

	dwRFirst[0]=1;
	dwRFirst[1]=1;
	dwRFirst[2]=1;
	dwRFirst[3]=1;
	dwImFirst[0]=1;
	dwImFirst[1]=1;
	dwImFirst[2]=1;
	dwImFirst[3]=1;
}

void Residue::FirstIteration(){
	up->VectorInitializer(resR,4,0);					//set residuesR to 0
	up->VectorInitializer(resIm,4,0);				//set residuesIm to 0

	double totR[4]={0};		double totIm[4]={0};

	//average resuduesR of the 1st iteration to be used for scaling of next iterations
	for (int j=0;j<celltot;j++){
		totR[0]=totR[0]+dwR[j][0]*dwR[j][0];
		totR[1]=totR[1]+dwR[j][1]*dwR[j][1];
		totR[2]=totR[2]+dwR[j][2]*dwR[j][2];
		totR[3]=totR[3]+dwR[j][3]*dwR[j][3];
	}

	dwRFirst[0]=sqrt(totR[0])/celltot;
	dwRFirst[1]=sqrt(totR[1])/celltot;
	dwRFirst[2]=sqrt(totR[2])/celltot;
	dwRFirst[3]=sqrt(totR[3])/celltot;

	resR[0]=log10(dwRFirst[0]);
	resR[1]=log10(dwRFirst[1]);
	resR[2]=log10(dwRFirst[2]);
	resR[3]=log10(dwRFirst[3]);

	//average resuduesIm of the 1st iteration to be used for scaling of next iterations
	for (int j=0;j<celltot;j++){
		totIm[0]=totIm[0]+dwIm[j][0]*dwIm[j][0];
		totIm[1]=totIm[1]+dwIm[j][1]*dwIm[j][1];
		totIm[2]=totIm[2]+dwIm[j][2]*dwIm[j][2];
		totIm[3]=totIm[3]+dwIm[j][3]*dwIm[j][3];
	}

	dwImFirst[0]=sqrt(totIm[0])/celltot;
	dwImFirst[1]=sqrt(totIm[1])/celltot;
	dwImFirst[2]=sqrt(totIm[2])/celltot;
	dwImFirst[3]=sqrt(totIm[3])/celltot;

	resIm[0]=log10(dwImFirst[0]);
	resIm[1]=log10(dwImFirst[1]);
	resIm[2]=log10(dwImFirst[2]);
	resIm[3]=log10(dwImFirst[3]);

	//print first iteration residue
	cout<<setw(7)<<fixed<<left<<1
	<<setw(12)<<setprecision(3)<<fixed<<left<<resR[0]<<setw(12)<<setprecision(3)<<fixed<<left<<resIm[0]
	<<setw(12)<<setprecision(3)<<fixed<<left<<resR[1]<<setw(12)<<setprecision(3)<<fixed<<left<<resIm[1]
	<<setw(12)<<setprecision(3)<<fixed<<left<<resR[2]<<setw(12)<<setprecision(3)<<fixed<<left<<resIm[2]
	<<setw(12)<<setprecision(3)<<fixed<<left<<resR[3]<<setw(12)<<setprecision(3)<<fixed<<left<<resIm[3]<<endl;

}

void Residue::CalculateResidue(int itr){

	up->VectorInitializer(resR,4,0);			//set residues to 0
	up->VectorInitializer(resIm,4,0);			//set residues to 0

	double scaledDwR[celltot][4]={0};
	double scaledDwI[celltot][4]={0};

	//convert residues to scaled residue
	for (int n=0;n<celltot;n++){
		scaledDwR[n][0]=abs(dwR[n][0]/dwRFirst[0]);
		scaledDwR[n][1]=abs(dwR[n][1]/dwRFirst[1]);
		scaledDwR[n][2]=abs(dwR[n][2]/dwRFirst[2]);
		scaledDwR[n][3]=abs(dwR[n][3]/dwRFirst[3]);
		scaledDwI[n][0]=abs(dwIm[n][0]/dwImFirst[0]);
		scaledDwI[n][1]=abs(dwIm[n][1]/dwImFirst[1]);
		scaledDwI[n][2]=abs(dwIm[n][2]/dwImFirst[2]);
		scaledDwI[n][3]=abs(dwIm[n][3]/dwImFirst[3]);
	}

	//find the maximum scaled rasidue among all cells
	for(int n=0;n<celltot;n++){
		if(abs(scaledDwR[n][0])>resR[0])
			resR[0]=scaledDwR[n][0];
		if(abs(scaledDwR[n][1])>resR[1])
			resR[1]=scaledDwR[n][1];
		if(abs(scaledDwR[n][2])>resR[2])
			resR[2]=scaledDwR[n][2];
		if(abs(scaledDwR[n][3])>resR[3])
			resR[3]=scaledDwR[n][3];

		if(abs(scaledDwI[n][0])>resIm[0])
			resIm[0]=scaledDwI[n][0];
		if(abs(scaledDwI[n][1])>resIm[1])
			resIm[1]=scaledDwI[n][1];
		if(abs(scaledDwI[n][2])>resIm[2])
			resIm[2]=scaledDwI[n][2];
		if(abs(scaledDwI[n][3])>resIm[3])
			resIm[3]=scaledDwI[n][3];
	}

	//take logarithm of maximum scaled residue
	resR[0]=log10(resR[0]);
	resR[1]=log10(resR[1]);
	resR[2]=log10(resR[2]);
	resR[3]=log10(resR[3]);
	resIm[0]=log10(resIm[0]);
	resIm[1]=log10(resIm[1]);
	resIm[2]=log10(resIm[2]);
	resIm[3]=log10(resIm[3]);

	cout<<setw(7)<<fixed<<left<<itr
	<<setw(12)<<setprecision(3)<<fixed<<left<<resR[0]<<setw(12)<<setprecision(3)<<fixed<<left<<resIm[0]
	<<setw(12)<<setprecision(3)<<fixed<<left<<resR[1]<<setw(12)<<setprecision(3)<<fixed<<left<<resIm[1]
	<<setw(12)<<setprecision(3)<<fixed<<left<<resR[2]<<setw(12)<<setprecision(3)<<fixed<<left<<resIm[2]
	<<setw(12)<<setprecision(3)<<fixed<<left<<resR[3]<<setw(12)<<setprecision(3)<<fixed<<left<<resIm[3]<<endl;
}

double Residue::findMax(){
	double maxR=-100;
	for(int i=0;i<4;i++)
		if (resR[i]>maxR)
			maxR=resR[i];

	double maxIm=-100;
	for(int i=0;i<4;i++)
		if (resIm[i]>maxIm)
			maxIm=resIm[i];

	double max=maxR;
	if(maxIm>maxR)
		max=maxIm;

	return max;
}
