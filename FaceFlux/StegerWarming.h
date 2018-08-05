/*=================================================
                Steger Class header
=================================================*/
#ifndef FACEFLUX_STEGERWARMING_H_
#define FACEFLUX_STEGERWARMING_H_

class Mesh;
class StegerWarming {
public:
	void Initializer(double*,double*,double*,double*,double*,double*,double**,double**);
	void CalculateFaceFlux(double);


private:
	double* wrFaceCurr=nullptr;
	double* wiFaceCurr=nullptr;
	double* wrFaceNeigh=nullptr;
	double* wiFaceNeigh=nullptr;
	double* faceFluxReal=nullptr;
	double* faceFluxImag=nullptr;
	double** jacobPlus=nullptr;
	double** jacobMinus=nullptr;
};

#endif
