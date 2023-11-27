#include <gtest/gtest.h>
#include "vectortemplate.h"
#include "matrixtemplate.h"


class MatrixTest : public ::testing::Test {
protected:

    void SetUp() override {
        // Optional setup before each test

    }

    void TearDown() override {
        // Optional teardown after each test
    }
};

TEST_F(MatrixTest, DefaultConstructor) {
    Matrix<int> matrix;
    EXPECT_EQ(matrix.get_rows(), 0);
    EXPECT_EQ(matrix.get_columns(), 0);
}

// Test constructor: Matrix(nRows, nCols)
TEST_F(MatrixTest, ConstructorWithSize) {
    Matrix<int> matrix(3, 4);
    EXPECT_EQ(matrix.get_rows(), 3);
    EXPECT_EQ(matrix.get_columns(), 4);

}



// Test constructor: Matrix(nRows, nCols, init_val)
TEST_F(MatrixTest, ConstructorWithSizeAndValue) {
    Matrix<int> matrix(2, 3, 5);
    EXPECT_EQ(matrix.get_rows(), 2);
    EXPECT_EQ(matrix.get_columns(), 3);

    // Verify that all elements are initialized to 5
    for (int i = 1; i <= matrix.get_rows(); i++) {
        for (int j = 1; j <= matrix.get_columns(); j++) {
            EXPECT_EQ(matrix(i, j), 5);
        }
    }
}

// Test constructor: Matrix(const CVector<T>& vector)
TEST_F(MatrixTest, ConstructorFromVector) {
    CVector<int> vector(3, 1);
    vector(1) = 10;
    vector(2) = 20;
    vector(3) = 30;

    Matrix<int> matrix(vector);
    EXPECT_EQ(matrix.get_rows(), 3);
    EXPECT_EQ(matrix.get_columns(), 1);

    // Verify that matrix is created from the vector correctly
    for (int i = 1; i <= matrix.get_rows(); i++) {
        EXPECT_EQ(matrix(i, 1), vector(i));
    }
}

// Test matrix equality: operator==
TEST_F(MatrixTest, MatrixEquality) {
    Matrix<int> matrix1(2, 3, 5);
    Matrix<int> matrix2(2, 3, 5);
    Matrix<int> matrix3(2, 3, 10);

    // Verify that matrix1 is equal to matrix2, but not equal to matrix3
    EXPECT_EQ(matrix1, matrix2);
    EXPECT_NE(matrix1, matrix3);
}

// Test matrix inequality: operator!=
TEST_F(MatrixTest, MatrixInequality) {
    Matrix<int> matrix1(2, 3, 5);
    Matrix<int> matrix2(2, 3, 5);
    Matrix<int> matrix3(2, 3, 10);

    // Verify that matrix1 is not equal to matrix3
    EXPECT_NE(matrix1, matrix3);
}

TEST_F(MatrixTest, AssignmentOperator) {
    Matrix<int> matrix1(2, 3, 5);
    Matrix<int> matrix2;
    matrix2 = matrix1;

    // Verify that matrix2 is a deep copy of matrix1
    EXPECT_EQ(matrix1, matrix2);
}

// Test copy constructor
TEST_F(MatrixTest, CopyConstructor) {
    Matrix<int> matrix1(2, 3, 5);
    Matrix<int> matrix2(matrix1);

    // Verify that matrix2 is a deep copy of matrix1
    EXPECT_EQ(matrix1, matrix2);
}

// Test assignment operator


// Test element-wise addition: operator+
TEST_F(MatrixTest, ElementWiseAddition) {
    Matrix<int> matrix1(2, 3, 5);
    Matrix<int> matrix2(2, 3, 10);
    Matrix<int> result = matrix1 + matrix2;

    // Verify that each element of the result is the sum of the corresponding elements in matrix1 and matrix2
    for (int i = 1; i <= result.get_rows(); i++) {
        for (int j = 1; j <= result.get_columns(); j++) {
            EXPECT_EQ(result(i, j), matrix1(i, j) + matrix2(i, j));
        }
    }
}

// Test element-wise addition with a scalar: operator+
TEST_F(MatrixTest, ElementWiseAdditionWithScalar) {
    Matrix<int> matrix(2, 3, 5);
    int scalar = 10;
    Matrix<int> result = matrix + scalar;

    // Verify that each element of the result is the sum of the corresponding element in the matrix and the scalar
    for (int i = 1; i <= result.get_rows(); i++) {
        for (int j = 1; j <= result.get_columns(); j++) {
            EXPECT_EQ(result(i, j), matrix(i, j) + scalar);
        }
    }
}

// Test compound addition: operator+=
TEST_F(MatrixTest, CompoundAddition) {
    Matrix<int> matrix1(2, 3, 5);
    Matrix<int> matrix2(2, 3, 10);
    matrix1 += matrix2;

    // Verify that matrix1 is modified correctly after compound addition with matrix2
    for (int i = 1; i <= matrix1.get_rows(); i++) {
        for (int j = 1; j <= matrix1.get_columns(); j++) {
            EXPECT_EQ(matrix1(i, j), 15);
        }
    }
}

// Test compound addition with a scalar: operator+=
TEST_F(MatrixTest, CompoundAdditionWithScalar) {
    Matrix<int> matrix(2, 3, 5);
    int scalar = 10;
    matrix += scalar;

    // Verify that each element of the matrix is modified correctly after compound addition with the scalar
    for (int i = 1; i <= matrix.get_rows(); i++) {
        for (int j = 1; j <= matrix.get_columns(); j++) {
            EXPECT_EQ(matrix(i, j), 15);
        }
    }
}

// Test element-wise subtraction: operator-
TEST_F(MatrixTest, ElementWiseSubtraction) {
    Matrix<int> matrix1(2, 3, 10);
    Matrix<int> matrix2(2, 3, 5);
    Matrix<int> result = matrix1 - matrix2;

    // Verify that each element of the result is the difference of the corresponding elements in matrix1 and matrix2
    for (int i = 1; i <= result.get_rows(); i++) {
        for (int j = 1; j <= result.get_columns(); j++) {
            EXPECT_EQ(result(i, j), matrix1(i, j) - matrix2(i, j));
        }
    }
}

// Test element-wise subtraction with a scalar: operator-
TEST_F(MatrixTest, ElementWiseSubtractionWithScalar) {
    Matrix<int> matrix(2, 3, 10);
    int scalar = 5;
    Matrix<int> result = matrix - scalar;

    // Verify that each element of the result is the difference of the corresponding element in the matrix and the scalar
    for (int i = 1; i <= result.get_rows(); i++) {
        for (int j = 1; j <= result.get_columns(); j++) {
            EXPECT_EQ(result(i, j), matrix(i, j) - scalar);
        }
    }
}

// Test compound subtraction: operator-=
TEST_F(MatrixTest, CompoundSubtraction) {
    Matrix<int> matrix1(2, 3, 10);
    Matrix<int> matrix2(2, 3, 5);
    matrix1 -= matrix2;

    // Verify that matrix1 is modified correctly after compound subtraction with matrix2
    for (int i = 1; i <= matrix1.get_rows(); i++) {
        for (int j = 1; j <= matrix1.get_columns(); j++) {
            EXPECT_EQ(matrix1(i, j), 5);
        }
    }
}

// Test compound subtraction with a scalar: operator-=
TEST_F(MatrixTest, CompoundSubtractionWithScalar) {
    Matrix<int> matrix(2, 3, 10);
    int scalar = 5;
    matrix -= scalar;

    // Verify that each element of the matrix is modified correctly after compound subtraction with the scalar
    for (int i = 1; i <= matrix.get_rows(); i++) {
        for (int j = 1; j <= matrix.get_columns(); j++) {
            EXPECT_EQ(matrix(i, j), 5);
        }
    }
}

// Test matrix multiplication with another matrix: operator*
TEST_F(MatrixTest, MatrixMultiplication) {
    Matrix<int> matrix1(2, 3, 2);
    Matrix<int> matrix2(3, 2, 3);
    Matrix<int> result = matrix1 * matrix2;

    // Verify the result of matrix multiplication
    EXPECT_EQ(result.get_rows(), 2);
    EXPECT_EQ(result.get_columns(), 2);
    EXPECT_EQ(result(1, 1), 18);
    EXPECT_EQ(result(1, 2), 18);
    EXPECT_EQ(result(2, 1), 18);
    EXPECT_EQ(result(2, 2), 18);
}

// Test matrix multiplication with a vector: operator*
TEST_F(MatrixTest, VectorMultiplication) {
    Matrix<int> matrix(2, 3, 2);
    CVector<int> vector(3, 1);
    vector(1) = 2;
    vector(2) = 3;
    vector(3) = 4;

    CVector<int> result = matrix * vector;

    // Verify the result of matrix multiplication with a vector
    EXPECT_EQ(result.GetSize(), 2);
    EXPECT_EQ(result(1), 18);
    EXPECT_EQ(result(2), 18);
}

// Test scalar multiplication: operator*
TEST_F(MatrixTest, ScalarMultiplication) {
    Matrix<int> mat(2, 3, 2);
    int scalar = 5;
    Matrix<int> result = mat * scalar;

    // Verify that each element of the result is the product of the corresponding element in the matrix and the scalar
    for (int i = 1; i <= result.get_rows(); i++) {
        for (int j = 1; j <= result.get_columns(); j++) {
            EXPECT_EQ(result(i, j), mat(i, j) * scalar);
        }
    }
}

// Test compound scalar multiplication: operator*=
TEST_F(MatrixTest, CompoundScalarMultiplication) {
    Matrix<int> mat(2, 3, 2);
    int scalar = 5;
    mat *= scalar;

    // Verify that each element of the matrix is modified correctly after compound scalar multiplication
    for (int i = 1; i <= mat.get_rows(); i++) {
        for (int j = 1; j <= mat.get_columns(); j++) {
            EXPECT_EQ(mat(i, j), 10);
        }
    }
}

// Test scalar division: operator/
TEST_F(MatrixTest, ScalarDivision) {
    Matrix<int> matrix(2, 3, 10);
    int scalar = 2;
    Matrix<int> result = matrix / scalar;

    // Verify that each element of the result is the division of the corresponding element in the matrix by the scalar
    for (int i = 1; i <= result.get_rows(); i++) {
        for (int j = 1; j <= result.get_columns(); j++) {
            EXPECT_EQ(result(i, j), matrix(i, j) / scalar);
        }
    }
}

// Test compound scalar division: operator/=
TEST_F(MatrixTest, CompoundScalarDivision) {
    Matrix<int> matrix(2, 3, 10);
    int scalar = 2;
    matrix /= scalar;

    // Verify that each element of the matrix is modified correctly after compound scalar division
    for (int i = 1; i <= matrix.get_rows(); i++) {
        for (int j = 1; j <= matrix.get_columns(); j++) {
            EXPECT_EQ(matrix(i, j), 5);
        }
    }
}

// Test get_rows() and get_columns() functions
TEST_F(MatrixTest, GetRowsAndColumns) {
    Matrix<int> matrix(4, 6);
    EXPECT_EQ(matrix.get_rows(), 4);
    EXPECT_EQ(matrix.get_columns(), 6);
}

// Test matrix element access with invalid indices
TEST_F(MatrixTest, InvalidElementAccess) {
    Matrix<int> matrix(3, 3, 5);

    // Attempt to access elements with invalid indices and verify that the appropriate exception is thrown
    EXPECT_THROW(matrix(0, 2), std::exception);
    EXPECT_THROW(matrix(4, 2), std::exception);
    EXPECT_THROW(matrix(2, 0), std::exception);
    EXPECT_THROW(matrix(2, 4), std::exception);
}

// Test display function
// Test << operator for Matrix class
TEST_F(MatrixTest, StreamOperator) {
    // Create a 2x2 matrix
    Matrix<int> matrix(2, 2);
    matrix(1, 1) = 1;
    matrix(1, 2) = 2;
    matrix(2, 1) = 3;
    matrix(2, 2) = 4;

    // Redirect cout to a stringstream to capture the output
    std::stringstream ss;
    std::streambuf* cout_buf = std::cout.rdbuf();
    std::cout.rdbuf(ss.rdbuf());

    // Use the << operator to output the matrix to the stringstream
    std::cout << matrix;

    // Restore cout buffer
    std::cout.rdbuf(cout_buf);

    // Verify that the output matches the expected format
    std::string expected_output = "2x2 Matrix:\n[ \n  [ 1 2 ]\n  [ 3 4 ]\n ]\n";
    EXPECT_EQ(ss.str(), expected_output);
}

TEST_F(MatrixTest, InvalidSizeError) {
    // Attempt to create a matrix with invalid size (negative rows and columns)
    ASSERT_THROW(Matrix<int> matrix(-2, 3), std::exception);

    // Attempt to create a matrix with invalid size (zero rows and columns)
    ASSERT_THROW(Matrix<int> matrix(0, 0), std::exception);
}

// Test error handling for invalid index access
TEST_F(MatrixTest, InvalidIndexError) {
    // Create a 2x2 matrix
    Matrix<int> matrix(2, 2);

    // Attempt to access an element with an invalid index (out of range)
    ASSERT_THROW(matrix(3, 1), std::exception);
    ASSERT_THROW(matrix(1, 3), std::exception);
}

// Test error handling for incompatible matrices
TEST_F(MatrixTest, IncompatibleMatricesError) {
    // Create two matrices of different sizes
    Matrix<int> matrix1(2, 2);
    Matrix<int> matrix2(3, 3);

    // Attempt to perform element-wise addition with incompatible matrices
    ASSERT_THROW(matrix1 + matrix2, std::exception);

    // Attempt to perform matrix multiplication with incompatible matrices
    ASSERT_THROW(matrix1 * matrix2, std::exception);
}

// Test error handling for division by zero
TEST_F(MatrixTest, DivisionByZeroError) {
    // Create a 2x2 matrix
    Matrix<int> matrix(2, 2);

    // Attempt to perform scalar division by zero
    ASSERT_THROW(matrix / 0, std::exception);
}


