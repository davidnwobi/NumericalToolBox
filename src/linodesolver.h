#pragma once
#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <cfenv>
#include "vectortemplate.h"
#include "matrixtemplate.h"
#include "matrixtoolbox.h"


/**
 * @namespace LINODESOLVER
 * @brief The LINODESOLVER namespace provides numerical methods to solve ordinary differential equations (ODEs).
 *
 * The LINODESOLVER namespace contains functions and classes that implement numerical methods to solve systems of ordinary differential equations (ODEs).
 * It includes the ode4 method, a fourth-order Runge-Kutta algorithm, which is widely used for solving initial value problems (IVPs) of ODEs.
 * The namespace also contains utility functions for handling floating-point exceptions and error conditions.
 *
 * The key components in the LINODESOLVER namespace are:
 * - FunctionContainer: A class that holds a collection of functions to be solved in a system of ODEs.
 * - FunctionPtr: A function pointer type representing a single variable function of type double.
 * - ode4: The main function for solving a system of ODEs using the ode4 method.
 * - ode4wrapperLinear: A helper function to solve a linear ODE using the ode4 method with specified boundary conditions.
 *
 * Example usage of the LINODESOLVER namespace:
 * @code{.cpp}
 * #include "linodesolver.h"
 *
 * // Define the ODE functions for a second-order linear ODE
 * double y1_prime(CVector<double> points) {
 *     // y' = y2 (the second derivative of y is represented as y2)
 *     return points(3);
 * }
 *
 * double y2_prime(CVector<double> points) {
 *     // y'' = -y' - y (the first derivative of y is represented as y1)
 *     return -points(3) - points(2);
 * }
 *
 * int main() {
 *     try {
 *         LINODESOLVER::FunctionContainer functions;
 *         functions.addFunction(&y1_prime);
 *         functions.addFunction(&y2_prime);
 *
 *         // Set up the initial conditions: y(0) = 0 and y'(0) = 1
 *         Matrix<double> solution;
 *         solution = LINODESOLVER::ode4wrapperLinear(functions, 0.0, 0.2, 0.001, { 0.0, 1.0, 1.0 });
 *
 *         // Print the solution matrix
 *         std::cout << solution;
 *     } catch (std::exception e) {
 *         std::cout << e.what();
 *     }
 *
 *     return 0;
 * }
 * @endcode
 *
 * Note: The above example demonstrates solving a second-order linear ODE using ode4 with specified boundary conditions at x = 0.
 *       For different ODEs or systems of ODEs, you can define appropriate functions and initial conditions accordingly.
 */

namespace LINODESOLVER {

	/**
	 * @brief Function pointer type for representing a single variable function of type double.
	 */
	using FunctionPtr = double (*)(CVector<double>);

	/**
	 * @brief A class that holds a collection of functions to be solved.
	 */
	class FunctionContainer {
	public:
		/**
		 * @brief A vector of std::function objects representing the functions to be solved.
		 */
		std::vector<std::function<double(CVector<double>)>> functions;

		/**
		 * @brief Adds a function to the container.
		 * @param func A function pointer to be added to the container.
		 */
		void addFunction(FunctionPtr func) {
			functions.emplace_back(func);
		}
	};

	/**
	 * @brief Solve a system of ordinary differential equations using the Runge-Kutta-Fehlberg numerical method.
	 * @param funcs The FunctionContainer containing the functions to be solved.
	 * @param x_points The vector of x points at which to evaluate the solution.
	 * @param init_points The initial values for the dependent variables at x_points(1).
	 * @return A Matrix<double> containing the solution of the ODE system at the given x points.
	 */
	Matrix<double> ode45(FunctionContainer funcs,
		double start_point, double end_point,
		CVector<double> init_points,
		double TOL = 10e-2);

	/**
	 * @brief Solve a linear ordinary differential equation using the ode45 numerical method.
	 * @param functions The FunctionContainer containing the linear differential equation function.
	 * @param start_point The starting x value.
	 * @param end_point The ending x value.
	 * @param boundary_conditions An initializer list containing the boundary conditions (initial values) for the dependent variable.
	 * @return A Matrix<double> containing the solution of the linear ODE at the given x points.
	 * @throw std::exception If the number of boundary conditions does not match the number of equations or if the step size is invalid.
	 */
	Matrix<double> ode45wrapperLinear(FunctionContainer functions, double start_point, double end_point, std::initializer_list<double> boundary_conditions, double tolerance);

	/**
	 * @brief This namespace contains internal helper functions used in the LINODESOLVER namespace.
	 *
	 * The `detail` namespace contains utility functions that are used internally in the LINODESOLVER namespace.
	 * These functions are not intended to be used directly by the user but provide crucial support for the main
	 * ODE solver functions.
	 */
	namespace detail {

		/**
		 * @brief Enumeration of error types.
		 */
		enum class ERROR { BOUNDARYCONDITIONANDEQUATIONMISMATCH, FLOATEXCEPTION, INVALIDHVALUE };

		/**
		 * @brief Handle different types of errors and throw appropriate exceptions.
		 * @param err The error type.
		 * @throw std::exception Depending on the error type.
		 */
		void ErrorHandler(ERROR err);

		/**
		 * @brief Check if floating-point exception occurred due to overflow.
		 * @param x The floating-point number to check.
		 * @return std::runtime_error
		 */
		void FloatException(double x);

		/**
		 * @brief Compute the k1 coefficients for the ode45 method.
		 * @param funcs The FunctionContainer containing the functions to be solved.
		 * @param init_points The initial values for the dependent variables.
		 * @return A CVector<double> containing the k1 coefficients.
		 */
		CVector<double> compute_k1_ODE45(FunctionContainer funcs, CVector<double> init_points, double step_size);

		/**
			* @brief Compute the k2 coefficients for the ode45 method.
			* @param funcs The FunctionContainer containing the functions to be solved.
			* @param init_points The initial values for the dependent variables.
			* @param k_coeff The matrix of k coefficents.
			* @param step_size The step size for the numerical integration.
			* @param rkf45_tableau The RKF45 Butcher's Tableau
			* @return A CVector<double> containing the k2 coefficients.
			*/
		CVector<double> compute_k2_ODE45(FunctionContainer funcs, CVector<double> init_points,
			const Matrix<double>& k_coeff, double step_size, const Matrix<double>& rkf45_tableau);

		/**
		 * @brief Compute the k3 coefficients for the ode45 method.
		 * @param funcs The FunctionContainer containing the functions to be solved.
		 * @param init_points The initial values for the dependent variables.
		 * @param k_coeff The matrix of k coefficents.
		 * @param step_size The step size for the numerical integration.
		 * @param rkf45_tableau The RKF45 Butcher's Tableau
		 * @return A CVector<double> containing the k3 coefficients.
		 */
		CVector<double> compute_k3_ODE45(FunctionContainer funcs, CVector<double> init_points, 
			const Matrix<double>& k_coeff, double step_size, const Matrix<double>& rkf45_tableau);

		/**
		 * @brief Compute the k4 coefficients for the ode45 method.
		 * @param funcs The FunctionContainer containing the functions to be solved.
		 * @param init_points The initial values for the dependent variables.
		 * @param k_coeff The matrix of k coefficents.
		 * @param step_size The step size for the numerical integration.
		 * @param rkf45_tableau The RKF45 Butcher's Tableau
		 * @return A CVector<double> containing the k3 coefficients.
		 */
		CVector<double> compute_k4_ODE45(FunctionContainer funcs, CVector<double> init_points,
			const Matrix<double>& k_coeff, double step_size, const Matrix<double>& rkf45_tableau);

		/**
		 * @brief Compute the k5 coefficients for the ode45 method.
		 * @param funcs The FunctionContainer containing the functions to be solved.
		 * @param init_points The initial values for the dependent variables.
		 * @param k_coeff The matrix of k coefficents.
		 * @param step_size The step size for the numerical integration.
		 * @param rkf45_tableau The RKF45 Butcher's Tableau
		 * @return A CVector<double> containing the k3 coefficients.
		 */
		CVector<double> compute_k5_ODE45(FunctionContainer funcs, CVector<double> init_points,
			const Matrix<double>& k_coeff, double step_size, const Matrix<double>& rkf45_tableau);

		/**
		 * @brief Compute the k6 coefficients for the ode45 method.
		 * @param funcs The FunctionContainer containing the functions to be solved.
		 * @param init_points The initial values for the dependent variables.
		 * @param k_coeff The matrix of k coefficents.
		 * @param step_size The step size for the numerical integration.
		 * @param rkf45_tableau The RKF45 Butcher's Tableau
		 * @return A CVector<double> containing the k3 coefficients.
		 */
		CVector<double> compute_k6_ODE45(FunctionContainer funcs, CVector<double> init_points,
			const Matrix<double>& k_coeff, double step_size, const Matrix<double>& rkf45_tableau);

		/**
		 * @brief Compute the coefficients k1, k2, k3, k4, k5, k6 for the ode45 method.
		 * @param funcs The FunctionContainer containing the functions to be solved.
		 * @param init_points The initial values for the dependent variables.
		 * @param step_size The step size for the numerical integration.
		 * @return A Matrix<double> containing the k coefficients.
		 */
		Matrix<double> compute_k_ODE45(FunctionContainer funcs, CVector<double> init_points, 
			double step_size, const Matrix<double>& rkf45_tableau);


		/**
		 * @brief Check for overflow exception and throw an exception if detected.
		 */
		void overflow_detected();

		/**
		 * @brief Check for underflow exception and throw an exception if detected.
		 */
		void underflow_detected();
	}

}

// Compute the k1 coefficients for the ode4 method.
CVector<double> LINODESOLVER::detail::compute_k1_ODE45(FunctionContainer funcs, CVector<double> init_points, double step_size) {
	// Note the step size is used within this function as per the algorithim used
	// It could just as easily be applied when the calculations are done as with ODE4 
	// It is kept here for consistency
	CVector<double> k1(funcs.functions.size());
	size_t i = 1;

	// Evaluate each function in the FunctionContainer at the initial points and store the results in k1.
	for (const auto& func : funcs.functions) {
		k1(i) = func(init_points) * step_size;

		 FloatException(k1(i));

		i++;
	}
	detail::overflow_detected();
	detail::underflow_detected();
	return k1;
}

CVector<double> LINODESOLVER::detail::compute_k2_ODE45(FunctionContainer funcs, CVector<double> init_points,
	const Matrix<double>& kcoeff, double step_size, const Matrix<double>& rkf45_tableau) {

	CVector<double> k2(funcs.functions.size());
	
	CVector<double> points_for_calc(init_points);

	// Compute the intermediate points for k2 using the current points and k1.
	points_for_calc(1) += step_size * rkf45_tableau(1, 1);

	 FloatException(points_for_calc(1));

	for (int i = 2; i <= points_for_calc.GetSize(); i++) {

		points_for_calc(i) += (kcoeff(1, i - 1) * rkf45_tableau(1, 2));

		 FloatException(points_for_calc(i));
	}

	size_t i = 1;
	// Evaluate each function in the FunctionContainer at the intermediate points and store the results in k2.
	for (const auto& func : funcs.functions) {
		k2(i) = func(points_for_calc) * step_size;

		 FloatException(k2(i));

		i++;
	}
	detail::overflow_detected();
	detail::underflow_detected();
	return k2;
}

CVector<double> LINODESOLVER::detail::compute_k3_ODE45(FunctionContainer funcs, CVector<double> init_points,
	const Matrix<double>& kcoeff, double step_size, const Matrix<double>& rkf45_tableau) {

	CVector<double> k3(funcs.functions.size());
	size_t i = 1;
	CVector<double> points_for_calc(init_points);

	// Compute the intermediate points for k3 using the current points k1, k2.
	points_for_calc(1) += step_size * rkf45_tableau(2, 1);

	 FloatException(points_for_calc(1));

	for (int i = 2; i <= points_for_calc.GetSize(); i++) {
		points_for_calc(i) += (kcoeff(1, i - 1) * rkf45_tableau(2, 2) + kcoeff(2, i - 1) * rkf45_tableau(2, 3));

		 FloatException(points_for_calc(i));
	}

	// Evaluate each function in the FunctionContainer at the intermediate points and store the results in k2.
	for (const auto& func : funcs.functions) {
		k3(i) = func(points_for_calc) * step_size;

		 FloatException(k3(i));

		i++;
	}
	detail::overflow_detected();
	detail::underflow_detected();
	return k3;
}

CVector<double> LINODESOLVER::detail::compute_k4_ODE45(FunctionContainer funcs, CVector<double> init_points,
	const Matrix<double>& kcoeff, double step_size, const Matrix<double>& rkf45_tableau) {

	CVector<double> k4(funcs.functions.size());
	size_t i = 1;
	CVector<double> points_for_calc(init_points);

	// Compute the intermediate points for k4 using the current points k1, k2 and k3.
	points_for_calc(1) += step_size * rkf45_tableau(3, 1);

	 FloatException(points_for_calc(1));

	for (int i = 2; i <= points_for_calc.GetSize(); i++) {
		points_for_calc(i) += (kcoeff(1, i - 1) * rkf45_tableau(3, 2) +
			kcoeff(2, i - 1) * rkf45_tableau(3, 3) + kcoeff(3, i - 1) * rkf45_tableau(3, 4));

		 FloatException(points_for_calc(i));
	}

	// Evaluate each function in the FunctionContainer at the intermediate points and store the results in k2.
	for (const auto& func : funcs.functions) {
		k4(i) = func(points_for_calc) * step_size;

		 FloatException(k4(i));

		i++;
	}
	detail::overflow_detected();
	detail::underflow_detected();
	return k4;
}

CVector<double> LINODESOLVER::detail::compute_k5_ODE45(FunctionContainer funcs, CVector<double> init_points,
	const Matrix<double>& kcoeff, double step_size, const Matrix<double>& rkf45_tableau) {

	CVector<double> k5(funcs.functions.size());
	size_t i = 1;
	CVector<double> points_for_calc(init_points);

	// Compute the intermediate points for k5 using the current points k1, k2, k3 and k4.
	points_for_calc(1) += step_size * rkf45_tableau(4, 1);

	 FloatException(points_for_calc(1));

	for (int i = 2; i <= points_for_calc.GetSize(); i++) {
		points_for_calc(i) += (kcoeff(1, i - 1) * rkf45_tableau(4, 2) +
			kcoeff(2, i - 1) * rkf45_tableau(4, 3) + kcoeff(3, i - 1) * rkf45_tableau(4, 4) +
			kcoeff(4, i - 1) * rkf45_tableau(4, 5));

		 FloatException(points_for_calc(i));
	}

	// Evaluate each function in the FunctionContainer at the intermediate points and store the results in k2.
	for (const auto& func : funcs.functions) {
		k5(i) = func(points_for_calc) * step_size;

		 FloatException(k5(i));

		i++;
	}
	detail::overflow_detected();
	detail::underflow_detected();
	return k5;
}

CVector<double> LINODESOLVER::detail::compute_k6_ODE45(FunctionContainer funcs, CVector<double> init_points,
	const Matrix<double>& kcoeff, double step_size, const Matrix<double>& rkf45_tableau) {

	CVector<double> k6(funcs.functions.size());
	size_t i = 1;
	CVector<double> points_for_calc(init_points);

	// Compute the intermediate points for k5 using the current points k1, k2, k3 and k4.
	points_for_calc(1) += step_size * rkf45_tableau(5, 1);

	 FloatException(points_for_calc(1));

	for (int i = 2; i <= points_for_calc.GetSize(); i++) {
		points_for_calc(i) += (kcoeff(1, i - 1) * rkf45_tableau(5, 2) +
			kcoeff(2, i - 1) * rkf45_tableau(5, 3) + kcoeff(3, i - 1) * rkf45_tableau(5, 4) +
			kcoeff(4, i - 1) * rkf45_tableau(5, 5) + kcoeff(5, i - 1) * rkf45_tableau(5, 6));

		 FloatException(points_for_calc(i));
	}

	// Evaluate each function in the FunctionContainer at the intermediate points and store the results in k2.
	for (const auto& func : funcs.functions) {
		k6(i) = func(points_for_calc) * step_size;

		 FloatException(k6(i));

		i++;
	}
	detail::overflow_detected();
	detail::underflow_detected();
	return k6;
}

// Compute the coefficients k1, k2, k3, k4, k5, k6 for the ode45 method.
Matrix<double>  LINODESOLVER::detail::compute_k_ODE45(FunctionContainer funcs, CVector<double> init_points, 
	double step_size, const Matrix<double>& rkf45_tableau) {

	Matrix<double>  k_coeff(6, funcs.functions.size());
	k_coeff(1) = compute_k1_ODE45(funcs, init_points, step_size);
	k_coeff(2) = compute_k2_ODE45(funcs, init_points, k_coeff, step_size, rkf45_tableau);
	k_coeff(3) = compute_k3_ODE45(funcs, init_points, k_coeff, step_size, rkf45_tableau);
	k_coeff(4) = compute_k4_ODE45(funcs, init_points, k_coeff, step_size, rkf45_tableau);
	k_coeff(5) = compute_k5_ODE45(funcs, init_points, k_coeff, step_size, rkf45_tableau);
	k_coeff(6) = compute_k6_ODE45(funcs, init_points, k_coeff, step_size, rkf45_tableau);

	return k_coeff;
}

Matrix<double> LINODESOLVER::ode45(FunctionContainer funcs, double start_point, double end_point,
	CVector<double> init_points, double TOL) {
	std::vector<CVector<double>> points;
	points.push_back(init_points);
	Matrix<double> rkf45_tableau(5, 6);
	rkf45_tableau(1) = { 1 / 4, 1 / 4, 0, 0, 0, 0 };
	rkf45_tableau(2) = { 3 / 8, 3 / 32, 9 / 32, 0, 0, 0 };
	rkf45_tableau(3) = { 12 / 13, 1932 / 2197, -7200 / 2197, 7296 / 2197, 0, 0 };
	rkf45_tableau(4) = { 1, 439 / 216, -8, 3680 / 513, -845 / 4104, 0 };
	rkf45_tableau(5) = { 1 / 2, -8 / 27, 2, -3544 / 2565, 1859 / 4104, -11 / 40 };
	double step_size = end_point - start_point;
	while (points.back()(1) < end_point){

		CVector<double> prev_points;
		prev_points = points.back();
		CVector<double> curr_points(prev_points.GetSize());
		Matrix<double> k_coeffs;

		double rk4 = 0;
		double rkf45 = 0;
		
		while (true){
			// Compute the coefficients k1, k2, k3, and k4 for the current step.
			k_coeffs = detail::compute_k_ODE45(funcs, prev_points, step_size, rkf45_tableau);
			// Update the points for the current step using the k coefficients and the step size.
			CVector<double> errors(curr_points.GetSize()-1);
			for (int j = 2; j <= funcs.functions.size() + 1; j++) {

				{
					curr_points(1) = points.back()(1) + step_size;
					rk4 = points.back()(j) + (k_coeffs(1, j - 1) * 25 / 216 + k_coeffs(3, j - 1) * 1408 / 2565 +
						k_coeffs(4, j - 1) * 2197 / 4101 - k_coeffs(5, j - 1) * 1 / 5);
					rkf45 = points.back()(j) + (k_coeffs(1, j - 1) * 16 / 135 + k_coeffs(3, j - 1) * 6656 / 12825 +
						k_coeffs(4, j - 1) * 28561 / 56430 - k_coeffs(5, j - 1) * 9 / 50 + k_coeffs(6, j - 1) * 2 / 55);

					errors(j - 1) = rkf45 - rk4;
				}

				curr_points(j) = rk4;
				detail::FloatException(curr_points(j));

			}
			double error = MATRIXTOOLBOX::TwoNorm(errors);
			if (error > (TOL)) {

				double s = pow(TOL / error, 1.0 / 4.0) * 0.84;
				step_size *= std::min(2.0, std::max(0.1, s));
				k_coeffs = detail::compute_k_ODE45(funcs, prev_points, step_size, rkf45_tableau);
			}
			else if (error == 0) {
				step_size *= 2;
			}
			else {
				break;
			}
		}
		
		
		detail::overflow_detected();
		detail::underflow_detected();
		points.push_back(curr_points);
	}
	Matrix<double> points_matrix(points.size(), points.back().GetSize());
	for (int i = 1; i <= points_matrix.get_rows(); i++) {
		points_matrix(i) = std::move(points.at(i - 1));
	}
	return points_matrix;
}

// Wrapper function to solve a linear ODE using ode4 with specified boundary conditions.
Matrix<double> LINODESOLVER::ode45wrapperLinear(FunctionContainer functions, double start_point, double end_point, std::initializer_list<double> boundary_conditions, double tolerance = 0.01) {
	try {
		if (boundary_conditions.size() - 1 != functions.functions.size()) {
			detail::ErrorHandler(detail::ERROR::BOUNDARYCONDITIONANDEQUATIONMISMATCH);
		}
		CVector<double> init_points;
		init_points = boundary_conditions;
		Matrix<double> solution;
		//solution = LINODESOLVER::ode4(functions, x_points, init_points, step_size);
		Matrix<double>solution2(LINODESOLVER::ode45(functions, start_point, end_point, init_points, tolerance));
		return solution2;
	}
	catch (std::exception e) {
		throw;
	}
}

// Function to check for overflow.
void LINODESOLVER::detail::overflow_detected() {
	if ((bool)std::fetestexcept(FE_OVERFLOW)) {
		std::cerr << "Warning: Underflow has occurred" << std::endl;
	}
}

// Function to check for underflow.
void LINODESOLVER::detail::underflow_detected() {
	if ((bool)std::fetestexcept(FE_UNDERFLOW)) {
		std::cerr << "Warning: Overflow has occurred" << std::endl;
	}
}

// Function to handle errors.
void LINODESOLVER::detail::ErrorHandler(ERROR err) {
	switch (err) {
	case LINODESOLVER::detail::ERROR::BOUNDARYCONDITIONANDEQUATIONMISMATCH:
		throw std::invalid_argument("The number of equations provided is not compatible with the number of boundary conditions");
	case LINODESOLVER::detail::ERROR::FLOATEXCEPTION:
		throw std::invalid_argument("Overflow or underflow condition during program execution.");
	case LINODESOLVER::detail::ERROR::INVALIDHVALUE:
		throw std::invalid_argument("Step size must be greater than zero");
	default:
		throw std::invalid_argument("Unknown error");
		break;
	}
}

void LINODESOLVER::detail::FloatException(double x)
{
	if (x == 0.0) {
		// Treat zero as a normal value (no exception thrown for zero).
		return;
	}

	switch (std::fpclassify(x))
	{
	case FP_INFINITE:  // Infinite, e.g. 1.0/0.0, 2*DBL_MAX
		throw std::runtime_error("Floating-point exception: Infinite value encountered. Value: " + std::to_string(x));
	case FP_NAN:       // Not-a-Number, e.g. 0.0/0.0, sqrt(-1.0)
		throw std::runtime_error("Floating-point exception: Not-a-Number (NaN) encountered. Value: " + std::to_string(x));
	case FP_SUBNORMAL: // Subnormal, e.g. DBL_MIN/2
		throw std::runtime_error("Floating-point exception: Subnormal value encountered. Value: " + std::to_string(x));
	case FP_NORMAL:    // normal value
		// Do nothing, no exception thrown for normal values.
		break;
	default:
		throw std::runtime_error("Unknown floating-point exception. Value: " + std::to_string(x));
	}
}

