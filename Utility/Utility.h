

#ifndef UTILITY_UTILITY_H_
#define UTILITY_UTILITY_H_

class Utility {
public:
	void VectorInitializer(double*,int,double);
	void MatrixInitializer(double**,int,int,double);
	void VectorEqualizer(double*,double*,int);
	void MatrixEqualizer(double**,double**,int,int);
	void Monitoring2D(int,int,double**);
	double VectorProduct2D(double,double,double,double);
	double* MatrixXvector(double**,double*,int,int);
};

#endif
