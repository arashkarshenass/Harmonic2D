/*=================================================
                Steger Class source
=================================================*/
#include "StegerWarming.h"
#include "../InputOutput/Mesh.h"
#include "../Inputs.h"

#include<cmath>
using namespace std;

void StegerWarming::Initializer(double*wrFaceCurrTemp,double*wiFaceCurrTemp,double*wrFaceNeighTemp,double*wiFaceNeighTemp,double*faceFluxRealTemp,double*faceFluxImagTemp,double**jacobPlusTemp,double**jacobMinusTemp)
{
	wrFaceCurr=wrFaceCurrTemp;
	wiFaceCurr=wiFaceCurrTemp;
	wrFaceNeigh=wrFaceNeighTemp;
	wiFaceNeigh=wiFaceNeighTemp;
	faceFluxReal=faceFluxRealTemp;
	faceFluxImag=faceFluxImagTemp;
	jacobPlus=jacobPlusTemp;
	jacobMinus=jacobMinusTemp;

}

void StegerWarming::CalculateFaceFlux(double dl){

	//F real +
	double fRealPlus[4]={0};
	fRealPlus[0]=jacobPlus[0][0]*wrFaceCurr[0] + jacobPlus[0][1]*wrFaceCurr[1] + jacobPlus[0][2]*wrFaceCurr[2] + jacobPlus[0][3]*wrFaceCurr[3];
	fRealPlus[1]=jacobPlus[1][0]*wrFaceCurr[0] + jacobPlus[1][1]*wrFaceCurr[1] + jacobPlus[1][2]*wrFaceCurr[2] + jacobPlus[1][3]*wrFaceCurr[3];
	fRealPlus[2]=jacobPlus[2][0]*wrFaceCurr[0] + jacobPlus[2][1]*wrFaceCurr[1] + jacobPlus[2][2]*wrFaceCurr[2] + jacobPlus[2][3]*wrFaceCurr[3];
	fRealPlus[3]=jacobPlus[3][0]*wrFaceCurr[0] + jacobPlus[3][1]*wrFaceCurr[1] + jacobPlus[3][2]*wrFaceCurr[2] + jacobPlus[3][3]*wrFaceCurr[3];

	//F real -
	double fRealMinus[4]={0};
	fRealMinus[0]=jacobMinus[0][0]*wrFaceNeigh[0] + jacobMinus[0][1]*wrFaceNeigh[1] + jacobMinus[0][2]*wrFaceNeigh[2] + jacobMinus[0][3]*wrFaceNeigh[3];
	fRealMinus[1]=jacobMinus[1][0]*wrFaceNeigh[0] + jacobMinus[1][1]*wrFaceNeigh[1] + jacobMinus[1][2]*wrFaceNeigh[2] + jacobMinus[1][3]*wrFaceNeigh[3];
	fRealMinus[2]=jacobMinus[2][0]*wrFaceNeigh[0] + jacobMinus[2][1]*wrFaceNeigh[1] + jacobMinus[2][2]*wrFaceNeigh[2] + jacobMinus[2][3]*wrFaceNeigh[3];
	fRealMinus[3]=jacobMinus[3][0]*wrFaceNeigh[0] + jacobMinus[3][1]*wrFaceNeigh[1] + jacobMinus[3][2]*wrFaceNeigh[2] + jacobMinus[3][3]*wrFaceNeigh[3];

	//F real
	faceFluxReal[0]=(fRealPlus[0]+fRealMinus[0])*dl;
	faceFluxReal[1]=(fRealPlus[1]+fRealMinus[1])*dl;
	faceFluxReal[2]=(fRealPlus[2]+fRealMinus[2])*dl;
	faceFluxReal[3]=(fRealPlus[3]+fRealMinus[3])*dl;

	//F imag +
	double fImagPlus[4]={0};
	fImagPlus[0]=jacobPlus[0][0]*wiFaceCurr[0] + jacobPlus[0][1]*wiFaceCurr[1] + jacobPlus[0][2]*wiFaceCurr[2] + jacobPlus[0][3]*wiFaceCurr[3];
	fImagPlus[1]=jacobPlus[1][0]*wiFaceCurr[0] + jacobPlus[1][1]*wiFaceCurr[1] + jacobPlus[1][2]*wiFaceCurr[2] + jacobPlus[1][3]*wiFaceCurr[3];
	fImagPlus[2]=jacobPlus[2][0]*wiFaceCurr[0] + jacobPlus[2][1]*wiFaceCurr[1] + jacobPlus[2][2]*wiFaceCurr[2] + jacobPlus[2][3]*wiFaceCurr[3];
	fImagPlus[3]=jacobPlus[3][0]*wiFaceCurr[0] + jacobPlus[3][1]*wiFaceCurr[1] + jacobPlus[3][2]*wiFaceCurr[2] + jacobPlus[3][3]*wiFaceCurr[3];

	//F imag -
	double fImagMinus[4]={0};
	fImagMinus[0]=jacobMinus[0][0]*wiFaceNeigh[0] + jacobMinus[0][1]*wiFaceNeigh[1] + jacobMinus[0][2]*wiFaceNeigh[2] + jacobMinus[0][3]*wiFaceNeigh[3];
	fImagMinus[1]=jacobMinus[1][0]*wiFaceNeigh[0] + jacobMinus[1][1]*wiFaceNeigh[1] + jacobMinus[1][2]*wiFaceNeigh[2] + jacobMinus[1][3]*wiFaceNeigh[3];
	fImagMinus[2]=jacobMinus[2][0]*wiFaceNeigh[0] + jacobMinus[2][1]*wiFaceNeigh[1] + jacobMinus[2][2]*wiFaceNeigh[2] + jacobMinus[2][3]*wiFaceNeigh[3];
	fImagMinus[3]=jacobMinus[3][0]*wiFaceNeigh[0] + jacobMinus[3][1]*wiFaceNeigh[1] + jacobMinus[3][2]*wiFaceNeigh[2] + jacobMinus[3][3]*wiFaceNeigh[3];

	//F imag
	faceFluxImag[0]=(fImagPlus[0]+fImagMinus[0])*dl;
	faceFluxImag[1]=(fImagPlus[1]+fImagMinus[1])*dl;
	faceFluxImag[2]=(fImagPlus[2]+fImagMinus[2])*dl;
	faceFluxImag[3]=(fImagPlus[3]+fImagMinus[3])*dl;




}
