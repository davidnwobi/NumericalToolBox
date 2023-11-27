#pragma once
#include "vectortemplate.h"
#include "matrixtemplate.h"
#include <cmath>
#include <exception>
#include <stdexcept>
#include<cstdlib>
#include <variant>

# define M_PI           3.14159265358979323846 
/**
 * @namespace MATRIXTOOLBOX
 * @brief Namespace containing various matrix and vector manipulation functions.
 */
namespace MATRIXTOOLBOX {


    /**
     * @brief Enum defining possible errors in the matrix toolbox.
     */
    enum Error {
        INCOMPATIBLEVECTORS,
        INCOMPATIBLEMATRICES,
        INCOMPATIBLEVECTORSCROSSPRODUCT,
        SQUAREMATRIXREQUIRED,
        DEPENDENTROWS
    };

    /**
    * @enum Mode
    * @brief Enumeration representing the mode of operation for solving linear equations.
    */
    enum Mode {
        DEFAULT,  /**< Default mode */
        LU,       /**< LU factorization mode */
        LDLT      /**< LDLT factorization mode */
    };

    enum RandMode {
        NORMAL,
        SYMMETRIC,
        POSITIVEDEF
    };

    /**
     * @brief Generates a random vector of type `CVector<double>`.
     * @param rows The number of rows in the vector.
     * @param start The starting value for random number generation (inclusive).
     * @param end The ending value for random number generation (exclusive).
     * @param seed The seed value for the random number generator.
     * @return A random vector of type `CVector<double>`.
     */
    CVector<double> random_vector(int rows, int start = 1, int end = 10, int seed = 1);

    /**
     * @brief Generates a random matrix of type `Matrix<double>`.
     * @param rows The number of rows in the matrix.
     * @param cols The number of columns in the matrix.
     * @param start The starting value for random number generation (inclusive).
     * @param end The ending value for random number generation (exclusive).
     * @param seed The seed value for the random number generator.
     * @return A random matrix of type `Matrix<double>`.
     */
    Matrix<double> random_matrix(int rows, int cols, int start = 1, int end = 10, RandMode randmode= NORMAL, int seed = 1);

    Matrix<double> identity(int n);

    template <class T, class... Args>
    void apply(CVector<T>& A, T (*func)(T, Args...), Args... args);

    template <class T>
    T sum(const CVector<T>& A);

    /**
     * @brief Calculates the dot product of two vectors.
     * @tparam T The data type of the vectors.
     * @param A The first vector.
     * @param B The second vector.
     * @return The dot product of the two vectors.
     */
    template <class T>
    T DotProduct(const CVector<T>& A, const CVector<T>& B);

    /**
     * @brief Calculates the 2-norm (Euclidean norm) of a vector.
     * @tparam T The data type of the vector.
     * @param A The vector.
     * @return The 2-norm of the vector.
     */
    template <class T>
    double TwoNorm(const CVector<T>& A);

    /**
     * @brief Normalizes a vector by dividing each element by its 2-norm.
     * @tparam T The data type of the vector.
     * @param A The vector to normalize.
     */
    template <class T>
    void Normalize(CVector<T>& A);

    /**
     * @brief Calculates the maximum value in a vector.
     * @tparam T The data type of the vector.
     * @param A The vector.
     * @return The maximum value in the vector.
     */
    template <class T>
    T MaxValue(const CVector<T>& A);

    /**
     * @brief Calculates the minimum value in a vector.
     * @tparam T The data type of the vector.
     * @param A The vector.
     * @return The minimum value in the vector.
     */
    template <class T>
    T MinValue(const CVector<T>& A);

    /**
     * @brief Calculates the maximum norm (infinity norm) of a vector.
     * @tparam T The data type of the vector.
     * @param A The vector.
     * @return The maximum norm of the vector.
     */
    template <class T>
    T MaxNorm(const CVector<T>& A);

    /**
     * @brief Calculates the cross product of two vectors in 3D space.
     * @tparam T The data type of the vectors.
     * @param A The first vector.
     * @param B The second vector.
     * @param C The resulting cross product vector.
     */
    template <class T>
    void CrossProduct(const CVector<T>& A, const CVector<T>& B, CVector<T>& C);

    /**
     * @brief Performs row reduction on a matrix and a vector.
     * @tparam T The data type of the matrix and vector.
     * @param A The matrix to be row reduced.
     * @param B The vector to be row reduced.
     * @param TOL The tolerance value for numerical comparisons.
     */
    template <class T>
    void RowReducer(Matrix<T>& A, CVector<T>& B, boost::optional<std::vector<std::tuple<int, int>>&> swaps = boost::none, T TOL = 1e-6);

    /**
     * @brief Performs row reduction on a matrix.
     * @tparam T The data type of the matrix.
     * @param A The matrix to be row reduced.
     * @param TOL The tolerance value for numerical comparisons.
     */
    template <class T>
    void RowReducer(Matrix<T>& A, T TOL = 1e-6);

    /**
     * @brief Calculates the determinant of a matrix.
     * @tparam T The data type of the matrix.
     * @param A The matrix.
     * @param c The determinant value (output parameter).
     */
    template <class T>
    void Determinant(const Matrix<T>& A, T& c);

    /**
     * @brief Performs forward substitution to solve a system of linear equations.
     * @tparam T The data type of the matrix and vectors.
     * @param A The matrix of coefficients.
     * @param x The solution vector.
     * @param b The constant vector.
     * @param mode The mode of substitution (DEFAULT, LU, LDLT).
     * @param TOL The tolerance value for numerical comparisons.
     */
    template <class T>
    void forward_substitution(const Matrix<T>& A, CVector<T>& x, const CVector<T>& b, Mode mode, T TOL = 1e-6);

    /**
     * @brief Performs backward substitution to solve a system of linear equations.
     * @tparam T The data type of the matrix and vectors.
     * @param A The matrix of coefficients.
     * @param x The solution vector.
     * @param b The constant vector.
     * @param mode The mode of substitution (DEFAULT, LU, LDLT).
     * @param TOL The tolerance value for numerical comparisons.
     */
    template <class T>
    void backward_substitution(const Matrix<T>& A, CVector<T>& x, const CVector<T>& b, Mode mode = Mode::DEFAULT, T TOL = 1e-6);

    /**
     * @brief Calculates the maximum norm (infinity norm) of a matrix.
     * @tparam T The data type of the matrix.
     * @param A The matrix.
     * @return The maximum norm of the matrix.
     */

    template <class T>
    Matrix<double> inverse_of_matrix(const Matrix<T>& P, double TOL= 10e-6);

    template <class T>
    T MaxNorm(const Matrix<T>& A);

    /**
     * @brief Transposes a matrix.
     * @tparam T The data type of the matrix.
     * @param A The matrix to transpose.
     * @param B The transposed matrix (output parameter).
     */
    template <class T>
    void Transpose(const Matrix<T>& A, Matrix<T>& B);

    /**
     * @brief Transposes a matrix.
     * @tparam T The data type of the matrix.
     * @param A The matrix to transpose.
     * @return The transposed matrix (output parameter).
     */
    template <class T>
    Matrix<T> Transpose(const Matrix<T>& A);

    /**
     * @brief Performs LU factorization on a matrix.
     * @tparam T The data type of the matrix.
     * @param A The matrix to be factorized.
     * @param TOL The tolerance value for numerical comparisons.
     */
    template <class T>
    void LUFactorization(Matrix<T>& A, T TOL = 1e-6);

    /**
     * @brief Solves a system of linear equations using LU factorization.
     * @tparam T The data type of the matrix and vectors.
     * @param A The matrix of coefficients.
     * @param x The solution vector.
     * @param b The constant vector.
     */
    template <class T>
    void LUSolve(const Matrix<T>& A, CVector<T>& x, const CVector<T>& b, double TOL=10e-6);

    /**
     * @brief Solves a system of linear equations using matrix equation Ax = b.
     * @tparam T The data type of the matrix and vectors.
     * @param A The matrix of coefficients.
     * @param x The solution vector.
     * @param b The constant vector.
     * @param TOL The tolerance value for numerical comparisons.
     */
    template <class T>
    void AxEqb(Matrix<T>& A, CVector<T>& x, CVector<T>& b, boost::optional<std::vector<std::tuple<int, int>>&> swaps = boost::none, T TOL = 1e-6);

    /**
     * @brief Performs LDLT factorization on a matrix.
     * @tparam T The data type of the matrix.
     * @param A The matrix to be factorized.
     * @param TOL The tolerance value for numerical comparisons.
     */
    template <class T>
    void LDLTFactorization(Matrix<T>& A, T TOL = 1e-6);

    /**
     * @brief Solves a system of linear equations using LDLT factorization.
     * @tparam T The data type of the matrix and vectors.
     * @param A The matrix of coefficients.
     * @param x The solution vector.
     * @param b The constant vector.
     * @param TOL The tolerance value for numerical comparisons.
     */
    template <class T>
    void LDLTSolve(const Matrix<T>& A, CVector<T>& x, const CVector<T>& b, T TOL = 1e-6);

    /**
     * @brief Handles errors in the matrix toolbox by throwing exceptions.
     * @param err The error code.
     */
    void ErrorHandler(Error err);

    namespace EIGEN {

        void vector_iter_inverse(const Matrix<double>& A, CVector<double>& x, std::variant<Matrix<double>, double> estimate, double& lambda, int max_iter, double TOL = 1e-6);

        void eigen_Jacobi(Matrix<double>& A, Matrix<double>& P, int max_iter, bool verbose=false, double TOL = 1e-6);

        void generalized_eigen_Jacobi(Matrix<double>& K,  Matrix<double>& X, Matrix<double>& LAMDBA, 
            Matrix<double>& M, int max_iter, bool verbose=false, boost::optional<Matrix<double>&>residuals=boost::none, double TOL = 1e-6);

        namespace detail {
            Matrix<double> jrot(Matrix<double>& A, int p, int q);
            Matrix<double> gjrot(Matrix<double>& K, Matrix<double>& M, int p, int q, double TOL=1e-6);
            void find_max_abs_val_loc(Matrix<double>& A, int& i, int& j);
            bool diagonalized(Matrix<double>& A, double TOL);
        }
    }

}

CVector<double> MATRIXTOOLBOX::random_vector(int rows, int start, int end, int seed) {
    CVector<double> A(rows, 0.0);
    try {
        if (seed == 1) {
            seed = (unsigned)time(NULL);
        }
        std::srand(seed);

        

        for (int i = 1; i <= A.GetSize(); i++) {
            A(i) = start + (static_cast<double>(rand()) / RAND_MAX) * (end - start);
        }
        return A;
    }
    catch (...) {
        throw;
    }
    return A;
}

Matrix<double> MATRIXTOOLBOX::random_matrix(int rows, int cols, int start, int end, RandMode randmode, int seed ) {
    Matrix<double> A(rows, cols);
    try {
        if (seed == 1) {
            seed = (unsigned)time(NULL);
        }
        std::srand(seed);

        

        for (int i = 1; i <= A.get_rows(); i++) {
            for (int j = 1; j <= A.get_columns(); j++) {
                A(i,j) = start + (static_cast<double>(rand()) / RAND_MAX) * (end - start);
            }
            
        }

        switch (randmode)
        {
        case MATRIXTOOLBOX::NORMAL:
            return A;
            break;
        case MATRIXTOOLBOX::SYMMETRIC:
            if (rows != cols) {
                throw std::invalid_argument("Must be a sqaure matrix");
            }
            return (A + MATRIXTOOLBOX::Transpose(A)) * 0.5;
            break;
        case MATRIXTOOLBOX::POSITIVEDEF:
            if (rows != cols) {

                throw std::invalid_argument("Must be a sqaure matrix");
            }
            A += MATRIXTOOLBOX::Transpose(A);
            A += 0.5;
            A += identity(rows) * rows;
            return A;
            break;
        default:
            break;
        }


    }
    catch(std::exception){
        throw;
    }
    return A;

}

template <class T, class... Args>
void MATRIXTOOLBOX::apply(CVector<T>& A, T (*func)(T, Args...), Args... args) {
    int size = A.GetSize();
    for (int i = 1; i <= size; i++) {
        A(i) = func(A(i), args...);
    }
}

template <class T>
T MATRIXTOOLBOX::sum(const CVector<T>& A) {
    T sum = 0;
    int size = A.GetSize();
    for (int i = 1; i <= size; i++) {
        sum += A(i);
    }
    return sum;
}


template <class T>
T MATRIXTOOLBOX::DotProduct(const CVector<T>& A, const CVector<T>& B) {
    if (A.GetSize() != B.GetSize()) {
        ErrorHandler(Error::INCOMPATIBLEVECTORS);
        return 0;
    }
    int sum = 0;
    for (int i = 1; i <= A.GetSize(); i++) {
        sum += (A(i) * B(i));
    } 
    return sum;
}

template <class T>
double MATRIXTOOLBOX::TwoNorm(const CVector<T>& A) {
    T sum = 0;
    for (int i = 1; i <= A.GetSize(); i++) {
        sum += (A(i) * A(i));
    }
    return sqrt(sum);

}

template <class T>
void MATRIXTOOLBOX::Normalize(CVector<T>& A) {
    double norm = MATRIXTOOLBOX::TwoNorm(A);
    A /= norm;
}

template <class T>
T MATRIXTOOLBOX::MaxValue(const CVector<T>& A) {
    T max = A(1);
    for (int i = 2; i <= A.GetSize(); i++) {
        if (A(i) > max) {
            max = A(i);
        }
    }
    return max;
}

template <class T>
T MATRIXTOOLBOX::MinValue(const CVector<T>& A) {
    T min = A(1);
    for (int i = 2; i <= A.GetSize(); i++) {
        if (A(i) < min) {
            min = A(i);
        }
    }
    return min;
}


template <class T>
T MATRIXTOOLBOX::MaxNorm(const CVector<T>& A) {
    T max = fabs(A(1));
    for (int i = 2; i <= A.GetSize(); i++) {
        if (fabs(A(i)) > max) {
            max = fabs(A(i));
        }
    }
    return max;
}


template <class T>
void MATRIXTOOLBOX::CrossProduct(const CVector<T>& A, const CVector<T>& B, CVector<T>& C)
{
    if (A.GetSize() != 3 || B.GetSize() != 3 || C.GetSize() != 3)
    {
        // Handle error: Vectors A, B, and C must have size 3 for cross product calculation
        ErrorHandler(MATRIXTOOLBOX::Error::INCOMPATIBLEVECTORSCROSSPRODUCT);
    }

    C.SetSize(3);
    C(1) = A(2) * B(3) - A(3) * B(2);
    C(2) = A(3) * B(1) - A(1) * B(3);
    C(3) = A(1) * B(2) - A(2) * B(1);
}

template <class T>
int max_val_in_col(Matrix<T>& A, int row, T TOL) {
    double max_val = -INFINITY;
    int max_j = row;
    for (int j = row; j <= A.get_rows(); j++) {
        if (fabs(A(row, j)) > max_val && A(row, j) < TOL) {
            max_val = fabs(A(row, j));
            max_j = j;
        }
    }
    return max_j;
}
template <class T>
void MATRIXTOOLBOX::RowReducer( Matrix<T>& A, CVector<T>& B, boost::optional<std::vector<std::tuple<int, int>>&> swaps, T TOL) {
    if (A.get_rows() != A.get_columns()) {
        ErrorHandler(Error::SQUAREMATRIXREQUIRED);
    }
    if (A.get_rows() != B.GetSize()) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
    }
    int N = A.get_rows();
    for (int k = 1; k <= N - 1; k++) {// Pivot row
        //check if this has the max value in col
        int max_j = max_val_in_col(A, k, TOL);
        if (max_j != k) {
            A.SwapRows(max_j, k);
            
            T temp = B(max_j);
            B(max_j) = B(k);
            B(k) = temp;
            if (swaps.has_value()) {
                std::tuple temp = { max_j, k };
                swaps.get().emplace_back(temp);
            }
        }
        if (A(k,k) != 0.0) {
            for (int i = k + 1; i <= N; i++) { // Current row
                if (A(i, k) != 0) {
                    double ratio = A(i, k) / A(k, k);

                    for (int j = k; j <= N; j++) { // iterate through cols
                        A(i, j) -= (ratio * A(k, j));
                        if (fabs(A(i, j)) < TOL) {
                            A(i, j) = 0;
                        }
                    }
                    B(i) -= (ratio * B(k));
                    if (fabs(B(i)) < TOL) {
                        B(i) = 0;
                    }
                }

            }
        }
    }
}

template <class T>
void MATRIXTOOLBOX::RowReducer(Matrix<T>& A, T TOL) {
    if (A.get_rows() != A.get_columns()) {
        ErrorHandler(Error::SQUAREMATRIXREQUIRED);
    }
    int N = A.get_rows();
    for (int k = 1; k <= N - 1; k++) {// Pivot row
        for (int i = k + 1; i <= N; i++) { // Current row
            if (A(i, k) != 0) {
                double ratio = A(i, k) / A(k, k);
                for (int j = k; j <= N; j++) { // iterate through cols
                    A(i, j) -= (ratio * A(k, j));
                    if (fabs(A(i, j)) < TOL) {
                        A(i, j) = 0;
                    }
                }
            }

        }
    }
}

template <class T>
void MATRIXTOOLBOX::Determinant(const Matrix<T>& A, T& c) {
    try {
        Matrix<T> tempA(A);
        RowReducer(tempA);
        c = 1;
        for (int i = 1; i <= A.get_rows(); i++) {
            c *= tempA(i, i);
        }
    }
    catch (...) {
        throw;
    }
}

template <class T>
void MATRIXTOOLBOX::forward_substitution(const Matrix<T>& A, CVector<T>& x, const CVector<T>& B, Mode mode, T TOL) {
    if (A.get_rows() != A.get_columns()) {
        ErrorHandler(Error::SQUAREMATRIXREQUIRED);
    }
    if (A.get_rows() != B.GetSize()) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
    }

    x.SetSize(B.GetSize());
    x = B;

    int N = A.get_rows();
    switch (mode)
    {
    case MATRIXTOOLBOX::DEFAULT:
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= i - 1; j++) {
                x(i) -= A(i, j) * x(j);
                if (A(i, i) == 0) {
                    ErrorHandler(Error::DEPENDENTROWS);
                }
            }
            x(i) /= A(i, i);
            if (fabs(x(i)) < TOL) {
                x(i) = 0;
            }
        }
        break;
    case MATRIXTOOLBOX::LU:
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= i - 1; j++) {
                x(i) -= A(i, j) * x(j);
                if (A(i, i) == 0) {
                    ErrorHandler(Error::DEPENDENTROWS);
                }
            }
            if (fabs(x(i)) < TOL) {
                x(i) = 0;
            }
        }
        break;
    case MATRIXTOOLBOX::LDLT:
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= i - 1; j++) {
                x(i) -= A(i, j) * x(j);
                if (A(i, i) == 0) {
                    ErrorHandler(Error::DEPENDENTROWS);
                }
            }
            if (fabs(x(i)) < TOL) {
                x(i) = 0;
            }
        }
        break;
    default:
        break;
    }



}

template <class T>
void MATRIXTOOLBOX::backward_substitution(const Matrix<T>& A, CVector<T>& x, const CVector<T>& B, Mode mode, T TOL) {
    try{ 
            if (A.get_rows() != A.get_columns()) {
                ErrorHandler(Error::SQUAREMATRIXREQUIRED); 
            }
            if (A.get_rows() != B.GetSize()) {
                ErrorHandler(Error::INCOMPATIBLEMATRICES);
            }

            x.SetSize(B.GetSize());
            x = B;

            int N = A.get_rows();


        switch (mode){
            case MATRIXTOOLBOX::DEFAULT:

                for (int i = N; i >= 1; i--) {
                    for (int j = i+1; j <= N; j++) {
                        x(i) -= A(i, j) * x(j);
                    }
                    if (A(i, i) == 0) {
                        
                        ErrorHandler(Error::DEPENDENTROWS);
                        
                    }
                    
                    x(i) /= A(i, i);
                    if (fabs(x(i)) < TOL) {
                        x(i) = 0;
                    }
                }
                break;
            default:
                break;
            }
        
    }
        
    catch(std::exception){
        throw;
    }
}

template <class T>
Matrix<double> MATRIXTOOLBOX::inverse_of_matrix(const Matrix<T>& P, double TOL) {

    Matrix<double> inverse_mat(P.get_rows(), P.get_rows());
    size_t row = 1;
    Matrix<T> large_mat(P.get_rows() * P.get_rows(), P.get_rows() * P.get_rows(), 0);
    for (int i = 1; i <= P.get_rows(); i++) {
        for (int j = 1; j <= P.get_rows(); j++) {
            for (int k = 1; k <= P.get_rows(); k++) {
                large_mat(row, ((k - 1) * P.get_columns() + j)) = P(i, k);

            }
            row++;
        }
    }
    CVector<double> inverse_vec;
    CVector<double> identity_matrix_vector(P.get_rows() * P.get_rows(), 0);
    for (int i = 1; i <= P.get_rows(); i++) {
        identity_matrix_vector((i - 1) * P.get_columns() + i) = 1;
    }
    std::vector<std::tuple<int, int>> swaps;
    AxEqb(large_mat, inverse_vec, identity_matrix_vector, swaps, TOL);

    for (int i = 1; i <= P.get_rows(); i++) {
        for (int j = 1; j <= P.get_rows(); j++) {
            inverse_mat(i, j) = inverse_vec((i - 1) * P.get_columns() + j);
            if (fabs(inverse_mat(i, j)) < TOL){
                inverse_mat(i, j) = 0;
            }
        }
    }
    return inverse_mat;
}

template <class T>
void MATRIXTOOLBOX::AxEqb(Matrix<T>& A, CVector<T>& x, CVector<T>& b, 
    boost::optional<std::vector<std::tuple<int, int>>&> swaps, T TOL) {
    try {

        if (swaps.has_value()) {
            MATRIXTOOLBOX::RowReducer(A, b, swaps.get(), TOL);
        }
        else {
            MATRIXTOOLBOX::RowReducer(A, b, boost::none, TOL);
        }
    
        MATRIXTOOLBOX::backward_substitution(A, x, b, Mode::DEFAULT, TOL);


    }
    catch (std::exception) {
        throw;
    }
}

//double MatrixMaxNorm(const std::vector<std::vector<T>>& matrix) {
//    double maxNorm = 0.0;
//
//    for (const auto& row : matrix) {
//        double rowSum = 0.0;
//        for (const auto& element : row) {
//            rowSum += std::abs(element);
//        }
//        maxNorm = std::max(maxNorm, rowSum);
//    }
//
//    return maxNorm;
//}

template <class T>
T MATRIXTOOLBOX::MaxNorm(const Matrix<T>& A) {
    try{
        T maxNorm = 0.0;
        for (int i = 1; i <= A.get_rows(); i++) {
            T temp_norm = MaxNorm(A(i));
            if (temp_norm > maxNorm) {
                maxNorm = temp_norm;
            }
        }
        return maxNorm;
    }
    catch (...) {
        throw;
    }
}

template <class T>
void MATRIXTOOLBOX::Transpose(const Matrix<T>& A, Matrix<T>& B) {
    B.SetSize(A.get_columns(), A.get_rows());
    for (int i = 1; i <= A.get_rows(); i++) {
        for (int j = 1; j <= A.get_columns(); j++) {

               B(j, i) = A(i, j);

        }
    }
}

template <class T>
Matrix<T> MATRIXTOOLBOX::Transpose(const Matrix<T>& A) {
    Matrix<T> B;
    B.SetSize(A.get_columns(), A.get_rows());
    for (int i = 1; i <= A.get_rows(); i++) {
        for (int j = 1; j <= A.get_columns(); j++) {

            B(j, i) = A(i, j);

        }
    }
    return B;
}


template <class T>
void MATRIXTOOLBOX::LUFactorization(Matrix<T>& A, T TOL) {
    try
    {
        if (A.get_rows() != A.get_columns()) {
            ErrorHandler(Error::SQUAREMATRIXREQUIRED);
        }

        for (int j = 1; j <= A.get_rows(); j++) {
            for (int i = 1; i <= j; i++) {

                for (int k = 1; k <= i - 1; k++) {
                    if (i == k) {
                        A(i, j) -= 1 * A(k, j);
                    }
                    else {
                        A(i, j) -= A(i, k) * A(k, j);
                    }

                }
                if (fabs(A(j, j)) < TOL) {
                    std::cout << j << "\n";
                    std::cerr << "The equations are not linearly independent" << "\n";
                }
            }
            for (int i = j + 1; i <= A.get_rows(); i++) {
                for (int k = 1; k <= j - 1; k++) {
                    A(i, j) -= A(i, k) * A(k, j);
                }
                A(i, j) /= A(j, j);
                
            }

        }
    }
    catch (...) {
        throw;
    }
}

template <class T>
void MATRIXTOOLBOX::LUSolve(const Matrix<T>& A, CVector<T>& x, const CVector<T>& b, double TOL) {

    try
    {
        Matrix<T> LU(A);
        MATRIXTOOLBOX::LUFactorization(LU, TOL);
        
        CVector<T> y;
        forward_substitution(LU, y, b, Mode::LU, TOL);
        backward_substitution(LU, x, y, Mode::DEFAULT, TOL);
    }
    catch (...) {
        throw;
    }


}

template <class T>
void MATRIXTOOLBOX::LDLTFactorization(Matrix<T>& A, T TOL) {
    try {
        if (A.get_rows() != A.get_columns()) {
            ErrorHandler(Error::SQUAREMATRIXREQUIRED);
        }
        Matrix<T> L(A.get_rows(), A.get_columns(), 0);
        Matrix<T> D(A.get_rows(), A.get_columns(), 0);
        for (int i = 1; i<= A.get_rows(); i++) {
            D(i, i) = A(i, i);
            for (int j = 1; j <= i - 1; j++) {
                if (i == j) {
                    A(i, i) -= A(j, j);
                }
                else {
                    A(i, i) -= A(i, j) * A(i, j) * A(j, j);
                }
            }
            for (int j = i + 1; j <= A.get_rows(); j++) {
                L(j, i) = A(j, i);
                for (int k = 1; k <= i - 1; k++) {
                    A(j, i) -= A(j, k) * A(k, k) * A(i, k);
                }
                A(j, i) /= A(i, i);

            }
        }
    }
    catch (...) {
        throw;
    }

}

template <class T>
void MATRIXTOOLBOX::LDLTSolve(const Matrix<T>& A, CVector<T>& x, const CVector<T>& b, T TOL) {
    try
    {
        if (A.get_rows() != A.get_columns()) {
            ErrorHandler(Error::SQUAREMATRIXREQUIRED);
        }

        Matrix<T> LD(A);
        LDLTFactorization(LD, TOL);

        CVector<T> y;
        forward_substitution(LD, y, b, Mode::LDLT);

        Matrix<T> LT(LD.get_rows(), LD.get_columns());
        Transpose(LD, LT);

        for (int i = 1; i < A.get_rows(); i++) {
            for (int j = i+1; j <= A.get_columns(); j++) {
                LT(i, j) = LT(i, i) * LT(i, j);

            }
        }

        backward_substitution(LD, x, y, Mode::DEFAULT);

        std::cout << x;
    }
    catch (...)
    {
        throw;
    }
}


void MATRIXTOOLBOX::ErrorHandler(MATRIXTOOLBOX::Error err) {

    switch (err)
    {
    case MATRIXTOOLBOX::INCOMPATIBLEVECTORS:
        throw std::runtime_error("Vector operation cannot take place due to incompatible vectors.");
        break;
    case MATRIXTOOLBOX::INCOMPATIBLEMATRICES:
        throw std::runtime_error("Matrix operation cannot take place due to incompatible matrices.");
        break;
    case MATRIXTOOLBOX::INCOMPATIBLEVECTORSCROSSPRODUCT:
        throw std::runtime_error("Invalid vector size for cross product calculation.");
        break;
    case MATRIXTOOLBOX::SQUAREMATRIXREQUIRED:
        throw std::runtime_error("A square matrix is required for this operation");
        break;
    case MATRIXTOOLBOX::DEPENDENTROWS:
        throw std::runtime_error("The equations are not linearly independent");
        break;
    default:
        throw std::runtime_error("You've fucked something up. Write good code");
        break;
    }

}

Matrix<double> MATRIXTOOLBOX::identity(int n) {
    Matrix<double> ident(n,n,0);
    for (int i = 1; i <= ident.get_rows(); i++) {
        ident(i, i) = 1;
    }
    return ident;
}


void MATRIXTOOLBOX::EIGEN::vector_iter_inverse(const Matrix<double>& A, CVector<double>& x, std::variant<Matrix<double>, double> estimate, double& lambda, int max_iter, double TOL) {
    if (A.get_rows() != A.get_columns()) {
        ErrorHandler(Error::SQUAREMATRIXREQUIRED);
    }
    CVector<double> x_1(x);
    size_t k = 1;
    Matrix<double> M (identity(A.get_rows()));
    Matrix<double> K(A - identity(A.get_rows()) * lambda);
    double lambda_1 = lambda;
    
    double est_dob;
    if (std::holds_alternative<Matrix<double>>(estimate)) {
        do
        {
            K = A - (identity(A.get_rows()) * std::get<Matrix<double>>(estimate));
            AxEqb(K, x_1, x);
            Normalize(x_1);
            Matrix<double>x_mat(x_1);
            lambda = (Transpose(x_mat) * (A * x_1))(1);
            x = x_1;

            k++;
        } while ((fabs((lambda_1 - lambda) / lambda_1) > TOL) && (k < max_iter));
        x = x * (1 / x(x.GetSize()));
    }
    else {
        do
        {
            K = A - (identity(A.get_rows()) * std::get<double>(estimate));
            AxEqb(K, x_1, x);
            Normalize(x_1);
            Matrix<double>x_mat(x_1);
            lambda = (Transpose(x_mat) * (A * x_1))(1);
            x = x_1;

            k++;
        } while ((fabs((lambda_1 - lambda) / lambda_1) > TOL) && (k < max_iter));
        x = x * (1 / x(x.GetSize()));
    }
    
    
}

void MATRIXTOOLBOX::EIGEN::eigen_Jacobi(Matrix<double>& A, Matrix<double>& P,
    int max_iter, bool verbose, double TOL) {
    if (A.get_rows() != A.get_columns()) {
        ErrorHandler(Error::SQUAREMATRIXREQUIRED);
        return;
    }
    size_t k = 0;
    int n = A.get_rows();
    Matrix<double> R(n, n);
    int max_i = 0;
    int max_j = 0;
    P = identity(n);
    bool done;
    do {
        detail::find_max_abs_val_loc(A, max_i, max_j);
        R = detail::jrot(A, max_i, max_j);
        P = P * R;
        // verify inverse

        A = (Transpose(R) * A) * R;

        k++;
        done = detail::diagonalized(A, TOL);
        
    } while ((k < max_iter) && !done);
    if (verbose) {
        P = Transpose(P);
        for (int i = 1; i <= n; i++) {
            std::cout << "{lamdba_" << i << " = " << A(i, i) << " (v_" << i << " = " << P(i) * (1 / P(i)(P(i).GetSize())) << ")}" << "\n";
        }

        P = Transpose(P);
    }
    
    

}

void MATRIXTOOLBOX::EIGEN::generalized_eigen_Jacobi(Matrix<double>& K, Matrix<double>& X, Matrix<double>& LAMDBA, 
    Matrix<double>& M, int max_iter, bool verbose, boost::optional<Matrix<double>&>residuals, double TOL) {
    if (K.get_rows() != K.get_columns()) {
        ErrorHandler(Error::SQUAREMATRIXREQUIRED);
        return;
    }
    size_t k = 0;
    int n = K.get_rows();
    if (n == 1)
        return;
    Matrix<double> R(n, n);
    int max_i = 1;
    int max_j = 1;
    bool done;
    detail::find_max_abs_val_loc(K, max_i, max_j);
    R = detail::gjrot(K, M, max_i, max_j);
    X = identity(n);
    std::cout << X;
    do
    {   
        if (((K(max_i, max_j) * K(max_i, max_j)) / (K(max_i, max_i) * K(max_j, max_j))) < TOL) {
            detail::find_max_abs_val_loc(K, max_i, max_j);
            //std::cout << "max_i: " << max_i << " max_j: " << max_j << "\n";
            R = detail::gjrot(K, M, max_i, max_j);
        }
       
        X = X* R;
        
        Matrix<double> K_hat ((Transpose(R) * K) * R);
        M = (Transpose(R) * M) * R;
        Matrix<double> M_inverse_square(M);
        Matrix<double> M_inverse(n, n, 0);
        for (int i = 1; i <= n; i++) {
            for (int j=1; j <= n; j++) {
                M_inverse_square(i, j) = 1 / sqrt(M(i, i));
                
            } 
            M_inverse(i, i) = 1 / M(i, i);
        }

        //std::cout << X;
        LAMDBA = M_inverse * K_hat;
        done = true;
        for (int i = 1; i <= n; i++) {
            if ((K_hat(i, i) - K(i, i)) / K_hat(i, i)>TOL) {
                done = false;
                break;
            }
        }

        K = std::move(K_hat);
        //std::cout << "M: " << M << "M_inverse_square: " << M_inverse_square << "\n";
        k++;
    } while ((k < max_iter) && !done);
    if (residuals.has_value()) {
        residuals.get() = K * X - LAMDBA * M * X;
    }
    if (verbose) {
        X = Transpose(X);
        for (int i = 1; i <= n; i++) {
            std::cout << "{lamdba_" << i << " = " << LAMDBA(i, i) << " (v_" << i << " = " << X(i) * (1 / X(i)(X(i).GetSize())) << ")}" << "\n";
        }
        X = Transpose(X);
        std::cout << "Residuals: " << K * X - LAMDBA * M * X;
    }
    
}

bool MATRIXTOOLBOX::EIGEN::detail::diagonalized(Matrix<double>& A, double TOL) {
    int n = A.get_rows();
    bool done = true;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if ((i != j) && (fabs(A(i, j)) > TOL)) {
                done = false;
                break;
            }
        }
    }
    return done;
}

void MATRIXTOOLBOX::EIGEN::detail::find_max_abs_val_loc(Matrix<double>& A, int& i, int& j) {
    double max_val = 0;
    int max_i = 1;
    int max_j = 1;
    int n = A.get_rows();
    for (int r = 1; r <= n; r++) {
        for (int c = 1; c <= n; c++) {
            if ((r != c) && (fabs(A(r, c)) > fabs(max_val))) {
                max_val = A(r, c);
                max_i = r;
                max_j = c;
            }
        }
    }
    i = max_i;
    j = max_j;

}

Matrix<double> MATRIXTOOLBOX::EIGEN::detail::jrot(Matrix<double>& A, int p, int q) {
    double theta = 0;
    if (A(p, p) == A(q, q)) {
        theta = M_PI / 4;
    }
    else {
        theta = 0.5 * atan(2 * A(p, q) / (A(p, p) - A(q, q)));
    }
    Matrix<double>R (identity(A.get_rows()));
    R(p, p) = cos(theta);
    R(q, q) = cos(theta);
    R(p, q) = -sin(theta);
    R(q, p) = sin(theta);

    return R;
}

Matrix<double> MATRIXTOOLBOX::EIGEN::detail::gjrot(Matrix<double>& K, Matrix<double>& M, int p, int q, double TOL) {
    double a = K(p, p) * M(p, q) - M(p, p) * K(p, q);
    double b = K(q, q) * M(p, q) - M(q, q) * K(p, q);
    double c = K(p, p) * M(q, q) - M(p, p) * K(q, q);
    double alpha, beta;
    if (fabs(a)< TOL ) {
        alpha = K(p, q) / K(q, q);
        beta = 0.0;
    }
    else if (fabs(b) < TOL) {
        beta = -K(p, q) / K(q, q);
        alpha = 0.0;
    }
    else {
        alpha = (-0.5 * c + (c / fabs(c)) * sqrt(0.25 * c * c + a * b)) / a;
        beta = -a * alpha / b;
    }
    Matrix<double>R(MATRIXTOOLBOX::identity(K.get_rows()));

    R(p, q) = alpha;
    R(q, p) = beta;
    
    //std::cout << std::setprecision(15) << "alpha: " << alpha << " beta: " << beta << "\n";
    return R;
}
