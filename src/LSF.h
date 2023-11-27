#pragma once

#include <iostream>
#include "vectortemplate.h"
#include "matrixtemplate.h"
#include "matrixtoolbox.h"
#include <exception>

class LSF {
public:
	LSF();

	void fit(int terms, 
		const CVector<double>& X, 
		const CVector<double>& y, 
		CVector<double>& coeff, 
		double& residuals, 
		double& r2);
	
private:
	Matrix<double> transform_X(const CVector<double>& X, int terms);
	void calc_residuals(const CVector<double>& X, const CVector<double>& y, CVector<double>& coeff, double& residuals, double& r2);

};
LSF::LSF(){
}

Matrix<double> LSF::transform_X(const CVector<double>& X, int terms) {
	Matrix<double> transformed(terms, X.GetSize(), 1);
	transformed(1) = CVector<double>(X.GetSize(), 1);
	transformed(2) = CVector<double>(X);
	for (int i = 3; i <= terms; i++) {
		CVector<double> temp_vec(X);
		MATRIXTOOLBOX::apply(temp_vec, pow, static_cast<double>(i-1));
		transformed(i) = temp_vec;
	}
	Matrix<double> transformed_transposed(X.GetSize(), terms);
	MATRIXTOOLBOX::Transpose(transformed, transformed_transposed);
	return transformed_transposed;
}

void LSF::calc_residuals(const CVector<double>& X, const CVector<double>& y, CVector<double>& coeff, double& residuals, double& r2) {
	
	try{
		CVector<double > y_hat = transform_X(X, coeff.GetSize()) * coeff;
		CVector<double> temp = y - y_hat;
		MATRIXTOOLBOX::apply(temp, pow, 2.0);;
		double SSreg = MATRIXTOOLBOX::sum(temp);
		double mean = MATRIXTOOLBOX::sum(y) / y.GetSize();
		temp = y - mean;
		MATRIXTOOLBOX::apply(temp, pow, 2.0);
		double SStot = MATRIXTOOLBOX::sum(temp);
		residuals = SSreg;
		r2 = 1 - SSreg / SStot;
	}
	catch (std::exception) {
		throw;
	}

}

void LSF::fit(int terms,
	const CVector<double>& X,
	const CVector<double>& y,
	CVector<double>& coeff,
	double& residuals,
	double& r2) {
	try {
		if (X.GetSize() != y.GetSize()) {
			throw std::invalid_argument("Arrays need to have the same length");
		}
		int length = X.GetSize();
		Matrix<double> transformed = transform_X(X, terms);
		Matrix<double> transformed_transposed(terms, X.GetSize());
		MATRIXTOOLBOX::Transpose(transformed, transformed_transposed);

		CVector<double> XTy = transformed_transposed * y;
		Matrix<double> XTX = transformed_transposed * transformed;


		MATRIXTOOLBOX::AxEqb(XTX, coeff, XTy);
		calc_residuals(X, y, coeff, residuals, r2);
	}
	catch(std::exception){
		throw;
	}
}
	