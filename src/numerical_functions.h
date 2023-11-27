#pragma once
#define EPSILON 0.000001
#include <functional>
#include <limits>
#include "print.h"
#include <cmath>
#include <iomanip>


enum Method { None, ForwardDiff, BackwardDiff, CentralDiff, FourthOrderCentralDiff};

double sign(double x) {
	if (x >= 0.0) {
		return 1.0;
	}
	else {
		return -1.0;
	}	
}


double ForwardMethod(int FUNCNO, double central_point, double h, double(*userfunc)(int FUNCNO, double central_point)) {
	double fxph = userfunc(FUNCNO, central_point + h);
	double fx = userfunc(FUNCNO, central_point);
	double dDeriv = (fxph - fx) / h;
	return dDeriv;
}

double BackwardMethod(int FUNCNO, double central_point, double h, double(*userfunc)(int FUNCNO, double central_point)) {
	double fx = userfunc(FUNCNO, central_point);
	double fxmh = userfunc(FUNCNO, central_point - h);
	double dDeriv = (fx - fxmh) / h;
	return dDeriv;
}

double CentralMethod(int FUNCNO, double central_point, double h, double(*userfunc)(int FUNCNO, double central_point)) {
	double fxph = userfunc(FUNCNO, central_point + h);
	double fxmh = userfunc(FUNCNO, central_point - h);
	double dDeriv = (fxph - fxmh) / (2 * h);
	return dDeriv;
}

double FourthOrderCentralMethod(int FUNCNO, double central_point, double h, double(*userfunc)(int FUNCNO, double central_point)) {
	double fxp2h = userfunc(FUNCNO, central_point + 2.0*h);
	double fxm2h = userfunc(FUNCNO, central_point - 2.0 * h);
	double fxph = userfunc(FUNCNO, central_point + h);
	double fxmh = userfunc(FUNCNO, central_point - h);
	double dDeriv = (-1 * fxp2h + 8 * fxph - 8 * fxmh + fxm2h) / (12 * h);
	return dDeriv;
}

double Trapezodial(int FUNCNO, int n, double upper_limit, double lower_limit,
	double(*userfunc)(int FUNCNO, double x)) {
	double interval = (upper_limit - lower_limit) / static_cast<double>(n);
	double area_under_curve = 0;
	for (double i = lower_limit; i <= upper_limit; i += interval) {
		if ((i == lower_limit) || (i == upper_limit)) {
			area_under_curve += userfunc(FUNCNO, i);
		}
		else {
			area_under_curve += 2*userfunc(FUNCNO, i);
		}
		
	}
	return area_under_curve * 1 / 2 * interval;
	
}

double Simpson(int FUNCNO, int n, double upper_limit, double lower_limit,
	double(*userfunc)(int FUNCNO, double x)) {
	double interval = (upper_limit - lower_limit) / static_cast<double>(n);
	double area_under_curve = 0;
	for (double i = lower_limit; i <= upper_limit; i += interval) {
		if ((i == lower_limit) || (i == upper_limit)) {
			area_under_curve += userfunc(FUNCNO, i);
		}
		else if (static_cast<int>((i - lower_limit)/interval)%2 == 1) {
			area_under_curve += 4 * userfunc(FUNCNO, i);
		}
		else {
			area_under_curve += 2 * userfunc(FUNCNO, i);
		}

	}
	return area_under_curve * 1 / 3 * interval;

}

bool inRange(double low, double high, double x) {
	return ((x - high) * (x - low) <= 0);
}


double deriv_func(int FUNCNO, double x) {

	switch (FUNCNO) {

	case 1:
		return pow(x, 3.0) - 2 * pow(x, 2.0) + 10 * x - 5;

	case 2:
		return pow(x, 4.0) - 2 * pow(x, 3.0) - 6 * x * x + x + 5;

	case 3:
		return (2 * pow(x, 3.0) - 4 * x + 10) / ((x - 1) * (x - 2));

	default:
		return 0;
	}

}
void test_deriv() {

	using namespace std;

	const int W = 20, S = 7, P = 15;
	enum Method { None, ForwardDiff, BackwardDiff, CentralDiff, FourthOrderCentralDiff};
	string method[] = { "None", "ForwardMethod", "BackwardMethod", "CentralMethod", "FourthOrderCentralMethod"};
	const int FUNC_NO = 2;
	double h_trial[] = { 1e-15, 1.0e-10, 1.0e-8, 1.0e-5, 1.0e-2, 1.0e-1 };
	double central_point = 2;
	double exact = -15;

	print("\n");

	print(setw(W), "   h      ");
	print(setw(S), "Method");
	print(setw(W), "  f'(x)   ");
	print(setw(W), "Rel. Error", "\n");
	print(string(3 * W + S, '-'), '\n');

	double best_error = std::numeric_limits<double>::max();
	double best_deriv, best_h;
	int elements = sizeof(h_trial) / sizeof(h_trial[0]);
	Method BestMethod = None;

	for (int i = 0; i < elements; i++) {
		double h = h_trial[i];
		double deriv_fd = ForwardMethod(FUNC_NO, central_point, h, deriv_func);
		double rel_error_fd = fabs((exact - deriv_fd) / exact);
		if (rel_error_fd < best_error) {
			best_error = rel_error_fd;
			best_deriv = deriv_fd;
			best_h = h;
			BestMethod = ForwardDiff;
		}
		double deriv_bd = BackwardMethod(FUNC_NO, central_point, h, deriv_func);
		double rel_error_bd = fabs((exact - deriv_bd) / exact);
		if (rel_error_bd < best_error) {
			best_error = rel_error_bd;
			best_deriv = deriv_bd;
			best_h = h;
			BestMethod = BackwardDiff;
		}
		double deriv_cd = CentralMethod(FUNC_NO, central_point, h, deriv_func);
		double rel_error_cd = fabs((exact - deriv_cd) / exact);
		if (rel_error_cd < best_error) {
			best_error = rel_error_cd;
			best_deriv = deriv_cd;
			best_h = h;
			BestMethod = CentralDiff;
		}
		double deriv_focd = FourthOrderCentralMethod(FUNC_NO, central_point, h, deriv_func);
		double rel_error_focd = fabs((exact - deriv_focd) / exact);
		if (rel_error_focd < best_error) {
			best_error = rel_error_focd;
			best_deriv = deriv_focd;
			best_h = h;
			BestMethod = FourthOrderCentralDiff;
		}

		std::cout.setf(ios::fixed); cout.precision(P);
		print(setw(W), h);
		print(setw(S), "FD");
		print(setw(W), deriv_fd);
		print(setw(W), rel_error_fd, "\n", "");
		print(setw(W), h);
		print(setw(S), "BD");
		print(setw(W), deriv_bd);
		print(setw(W), rel_error_bd, "\n", "");
		print(setw(W), h);
		print(setw(S), "CD");
		print(setw(W), deriv_cd);
		print(setw(W), rel_error_cd, "\n", "");
		print(setw(W), h);
		print(setw(S), "FOCD");
		print(setw(W), deriv_focd);
		print(setw(W), rel_error_focd, "\n", "");


	}

	print('\n', '\n');
	print("Best derivative estimate = ", best_deriv, "\n");
	print("          Coresponding h = ", best_h, "\n");
	print("          Method Used = ", method[BestMethod], "\n");
}

double brent_solver_mine(double(*func)(const int FUNC_NO, double x), int FUNC_NO, double a, double b, double c, int count = 0) {

	const int W = 20, S = 7, P = 15;



	using namespace std;

	if (count < 1) {
		print("\n");

		print(setw(W), "   a      ");
		print(setw(W), "   b      ");
		print(setw(W), "   c      ");
		print(setw(W), "  f(b)  ");
		print(setw(W), "  m  ");
		print("\n");
		print(string(5 * W, '-'), '\n');
	}



	double fa = func(FUNC_NO, a);
	double fb = func(FUNC_NO, b);
	double m = (a + b) / 2;

	cout.setf(ios::fixed); cout.precision(P);
	print(setw(W), a);
	print(setw(W), b);
	print(setw(W), c);
	print(setw(W), fb);
	print(setw(W), m, "\n", "");


	if (fabs(fb) < EPSILON || fabs(m) < EPSILON) {
		return b;
	}
	count++;
	if (count > 100) {
		print("Could not converge\n");
		return 0;
	}

	
	double s = b - fb * (b-a) / (fb - fa);

	if (inRange(m, b, s)) {
		c = b;
		b = s;
		if (fa * func(FUNC_NO, b) > 0) {
			a = c;
		}
	}
	else {
		c = b;
		b = m;
		if (func(FUNC_NO, b) * fa > 0) {
			a = c;
		}
		else {
			b = m;
		}
	}

	return brent_solver_mine(func, FUNC_NO, a, b, c, count);
	

}

double brent_solver(double(*func)(const int FUNC_NO, double x), int FUNC_NO, double a, double b, double c, int count = 0) {

	const int W = 20, S = 7, P = 15;

	using namespace std;

	if (count < 1) {
		print("\n");
		print(setw(W), "   a      ");
		print(setw(W), "  b  ");
		print(setw(W), "  m  ");
		print("\n");
		print(string(3 * W, '-'), '\n');
	}

	double m = (c - b) / 2;
	double bp;
	double bpp;
	/*
	
	double s = b - func(FUNC_NO, b) * (b - a) / (func(FUNC_NO, b) - func(FUNC_NO, a));
	double c;
	*/
	cout.setf(ios::fixed); cout.precision(P);
	print(setw(W), a);
	print(setw(W), b);
	print(setw(W), m, "\n", "");

	double fa = func(FUNC_NO, a);
	double fb = func(FUNC_NO, b);

	if ((fa * fb) > 0) {
		print("Not Bracketed", '\n');
		return 0;
	}
	
	if (fabs(fa) < EPSILON || fabs(m) < EPSILON) {
		return b;
	}

	double s = b - fb * (b - a) / (fb - fa);
	if (inRange(b, b + m, s)) {
		bp = s;
		
	}
	else {
		bpp = b + m;
	}
	if (fabs(b - bpp) > EPSILON) {
		bp = bpp;
	}
	else {
		bp = b + EPSILON * sign(m);
	}
	
	double bnew = bp;
	a = b;
	double fbnew = func(FUNC_NO, bnew);

	if (fbnew*fb < 0){
		c = b;
	}
	b = bnew;
	count++;
	if (count > 100) {
		print("Could not converge\n");
		return 0;
	}
	if (a > b) {
		return brent_solver(func, FUNC_NO, b, a, count);
	}

	else {
		return brent_solver(func, FUNC_NO, a, b, count);
	}

}


double defaultInit[2] = { 1,0 };
double quadrature_eval(double points[], double weights[], int n, int FUNC_NO, 
	double(*integral_func)(int FUNC_NO, double x), double mapper[2] = defaultInit) {

	double eval = 0.0;
	for (int i = 0; i < n; i++) {
		double mapped_x = mapper[0] * points[i] + mapper[1];
		eval += weights[i] * integral_func(FUNC_NO, mapped_x) * mapper[0];
	}
	return eval;

}

double quadrature(int FUNC_NO, int n, double upper_limit, double lower_limit,
	double(*integral_func)(int FUNC_NO, double x)) {
	double mapper[2];
	if ((lower_limit != -1) || (upper_limit != 1)) {
		double m = (upper_limit - lower_limit) / 2;
		mapper[0] = m;
		mapper[1] = upper_limit - m;
		upper_limit = 1;
		lower_limit = -1;
	}
	switch (n) {

	case 1: {
		double points[1] = { 0.0 }, weights[1] = { 2.0 };
		return quadrature_eval(points, weights, n, 1, integral_func, mapper);
	}

	case 2: {
		double points[2] = { 1.0 / sqrt(3), -1.0 / sqrt(3) }, weights[2] = { 1.0, 1.0 };
		return quadrature_eval(points, weights, n, 1, integral_func, mapper);
	}


	case 3: {
		double points[3] = { 0.0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0) },
			weights[3] = { 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 };
		return quadrature_eval(points, weights, n, 1, integral_func, mapper);
	}

	case 4: {
		double points[4] = { sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)) },
			weights[4] = { (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 };
		return quadrature_eval(points, weights, n, 1, integral_func, mapper);
	}
	case 5: {

		double points[5] = { 0.0, 1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)), -1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)), 1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)), -1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) },
			weights[5] = { 128.0 / 225.0, (322 + 13 * sqrt(70)) / 900 , (322 + 13 * sqrt(70)) / 900 , (322 - 13 * sqrt(70)) / 900 , (322 - 13 * sqrt(70)) / 900 };
		return quadrature_eval(points, weights, n, 1, integral_func, mapper);
	}


	default:
		return  std::numeric_limits<double>::quiet_NaN();;
	}

}



void bracket_finder(double(*func)(const int FUNCNO, double x), int FUNCNO, double& a, double& b) {
	a = -100;
	b = -101;


	while (func(FUNCNO, a) * func(FUNCNO, b) >=0) {
		a += 2;
		b -=2;
	}
	print("out");
	while (func(FUNCNO, a) * func(FUNCNO, b) <= 0) {
		print("fa = ", func(FUNCNO, a));
		print(" fb = ", func(FUNCNO, b), "\n");
		b++;
	}
	b--;
	a++;
	

}
