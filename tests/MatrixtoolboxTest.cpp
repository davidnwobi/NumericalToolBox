#include <gtest/gtest.h>
#include "vectortemplate.h"
#include "matrixtemplate.h"
#include "matrixtoolbox.h"

class MatrixToolboxTest : public ::testing::Test {

    
protected:
    Matrix<double> matrix;
    CVector<double> vector;

    void SetUp() override {
        matrix.SetSize(3, 3);
        matrix(1, 1) = 1.0;
        matrix(1, 2) = 2.0;
        matrix(1, 3) = 3.0;
        matrix(2, 1) = 2.0;
        matrix(2, 2) = 1.0;
        matrix(2, 3) = 3.0;
        matrix(3, 1) = 3.0;
        matrix(3, 2) = 2.0;
        matrix(3, 3) = 1.0;
        // Refactor
        vector.SetSize(3);
        vector(1) = 14;
        vector(2) = 13;
        vector(3) = 10;
    }

    
};

// Test random_vector function
TEST_F(MatrixToolboxTest, RandomVectorTest) {
    // Test generating a random vector of size 5 with default range and seed
    CVector<double> vec = MATRIXTOOLBOX::random_vector(5);
    ASSERT_EQ(vec.GetSize(), 5);
    // Test generating a random vector with custom range and seed
    CVector<double> vec_custom = MATRIXTOOLBOX::random_vector(3, 0, 100, 42);
    ASSERT_EQ(vec_custom.GetSize(), 3);
}

// Test random_matrix function
TEST_F(MatrixToolboxTest, RandomMatrixTest) {
    // Test generating a random matrix of size 3x3 with default range and seed
    Matrix<double> mat = MATRIXTOOLBOX::random_matrix(3, 3);
    ASSERT_EQ(mat.get_rows(), 3);
    ASSERT_EQ(mat.get_columns(), 3);
    // Test generating a random matrix with custom range and seed
    Matrix<double> mat_custom = MATRIXTOOLBOX::random_matrix(2, 2, -50, 50);
    ASSERT_EQ(mat_custom.get_rows(), 2);
    ASSERT_EQ(mat_custom.get_columns(), 2);
}

int incrementByTwo(int x, int y) { return x + 2; };
// Test apply function
TEST_F(MatrixToolboxTest, ApplyTest) {
    CVector<int> vec(5);
    for (int i = 1; i <= 5; ++i) {
        vec(i) = i;
    }

    // Define a function to apply (in this case, increment by 2)
    
    MATRIXTOOLBOX::apply(vec, incrementByTwo, 0);

    ASSERT_EQ(vec(1), 3);
    ASSERT_EQ(vec(2), 4);
    ASSERT_EQ(vec(3), 5);
    ASSERT_EQ(vec(4), 6);
    ASSERT_EQ(vec(5), 7);
}

// Test the sum function
TEST_F(MatrixToolboxTest, SumTest) {
    // Create a CVector with some values
    CVector<int> vec(5);
    vec(1) = 1;
    vec(2) = 2;
    vec(3) = 3;
    vec(4) = 4;
    vec(5) = 5;

    // Call the sum function
    int result = MATRIXTOOLBOX::sum(vec);

    // Check the expected result
    EXPECT_EQ(result, 15); // The sum of 1 + 2 + 3 + 4 + 5 should be 15
}

TEST_F(MatrixToolboxTest, DotProductTest) {
    // Create two CVectors with some values
    CVector<int> vecA(3);
    vecA(1) = 1;
    vecA(2) = 2;
    vecA(3) = 3;

    CVector<int> vecB(3);
    vecB(1) = 2;
    vecB(2) = 3;
    vecB(3) = 4;

    // Call the DotProduct function
    int result = MATRIXTOOLBOX::DotProduct(vecA, vecB);

    // Check the expected result
    EXPECT_EQ(result, 20); // The dot product of {1, 2, 3} and {2, 3, 4} should be 20
}

TEST_F(MatrixToolboxTest, TwoNormTest) {
    CVector<int> vec(3);
    vec(1) = 1;
    vec(2) = 2;
    vec(3) = -2;

    double result = MATRIXTOOLBOX::TwoNorm(vec);

    EXPECT_DOUBLE_EQ(result, 3.0); // sqrt(1^2 + 2^2 + 2^2) = 3.0
}

// Test the Normalize function
TEST_F(MatrixToolboxTest, NormalizeTest) {
    CVector<double> vec(3);
    vec(1) = 1.0;
    vec(2) = 2.0;
    vec(3) = 2.0;

    MATRIXTOOLBOX::Normalize(vec);
    // Check that the vector is normalized (its 2-norm is now 1.0)
    EXPECT_DOUBLE_EQ(MATRIXTOOLBOX::TwoNorm(vec), 1.0);
}

// Test the MaxValue function
TEST_F(MatrixToolboxTest, MaxValueTest) {
    CVector<int> vec(5);
    vec(1) = 3;
    vec(2) = 7;
    vec(3) = 2;
    vec(4) = 9;
    vec(5) = 5;

    int result = MATRIXTOOLBOX::MaxValue(vec);

    EXPECT_EQ(result, 9); // The maximum value in the vector is 9
}

// Test the MinValue function
TEST_F(MatrixToolboxTest, MinValueTest) {
    CVector<int> vec(5);
    vec(1) = 3;
    vec(2) = 7;
    vec(3) = 2;
    vec(4) = 9;
    vec(5) = 5;

    int result = MATRIXTOOLBOX::MinValue(vec);

    EXPECT_EQ(result, 2); // The minimum value in the vector is 2
}

// Test the MaxNorm function
TEST_F(MatrixToolboxTest, MaxNormTest) {
    CVector<int> vec(5);
    vec(1) = -3;
    vec(2) = 7;
    vec(3) = -2;
    vec(4) = -9;
    vec(5) = 5;

    int result = MATRIXTOOLBOX::MaxNorm(vec);

    EXPECT_EQ(result, 9); // The maximum absolute value in the vector is 9
}

// Test the CrossProduct function
TEST_F(MatrixToolboxTest, CrossProductTest) {
    CVector<int> vecA(3);
    vecA(1) = 1;
    vecA(2) = 2;
    vecA(3) = 3;

    CVector<int> vecB(3);
    vecB(1) = 2;
    vecB(2) = 3;
    vecB(3) = 4;

    CVector<int> vecC(3);

    MATRIXTOOLBOX::CrossProduct(vecA, vecB, vecC);

    // Check the expected result
    EXPECT_EQ(vecC(1), -1); // 2*4 - 3*3 = -1
    EXPECT_EQ(vecC(2), 2); // 3*2 - 1*4 = 2
    EXPECT_EQ(vecC(3), -1); // 1*3 - 2*2 = -1
}

TEST_F(MatrixToolboxTest, RowReducerTest) {
    Matrix<double> row_reduced_sol(3, 3);

    row_reduced_sol(1, 1) = 1.0;
    row_reduced_sol(1, 2) = 2.0;
    row_reduced_sol(1, 3) = 3.0;
    row_reduced_sol(2, 1) = 0;
    row_reduced_sol(2, 2) = -3;
    row_reduced_sol(2, 3) = -3;
    row_reduced_sol(3, 1) = 0.0;
    row_reduced_sol(3, 2) = 0.0;
    row_reduced_sol(3, 3) = -4;


    CVector<double> row_reduced_vector(3);
    row_reduced_vector(1) = 14.0;
    row_reduced_vector(2) = -15.0;
    row_reduced_vector(3) = -12.0;

    MATRIXTOOLBOX::RowReducer(matrix, vector);

    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            EXPECT_FLOAT_EQ(row_reduced_sol(i, j), matrix(i, j));
        }
        EXPECT_FLOAT_EQ(row_reduced_vector(i), vector(i));
    }

    
}

TEST_F(MatrixToolboxTest, RowReducerTestMatrixOnly) {
    Matrix<double> row_reduced_sol(3, 3);

    row_reduced_sol(1, 1) = 1.0;
    row_reduced_sol(1, 2) = 2.0;
    row_reduced_sol(1, 3) = 3.0;
    row_reduced_sol(2, 1) = 0;
    row_reduced_sol(2, 2) = -3;
    row_reduced_sol(2, 3) = -3;
    row_reduced_sol(3, 1) = 0.0;
    row_reduced_sol(3, 2) = 0.0;
    row_reduced_sol(3, 3) = -4;


    MATRIXTOOLBOX::RowReducer(matrix);

    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            EXPECT_FLOAT_EQ(row_reduced_sol(i, j), matrix(i, j));
        }
    }


}

TEST_F(MatrixToolboxTest, Determinant) {
    Matrix<double> row_reduced_sol(3, 3);

    double sol;


    MATRIXTOOLBOX::Determinant(matrix, sol);

    EXPECT_EQ(sol, 12);


}

TEST_F(MatrixToolboxTest, ForwardSubstitution) {
    Matrix<double> A(3, 3);

    A(1, 1) = 1;
    A(1, 2) = 0;
    A(1, 3) = 0;
    A(2, 1) = 2;
    A(2, 2) = 4;
    A(2, 3) = 0;
    A(3, 1) = 3;
    A(3, 2) = 2;
    A(3, 3) = 5;

    CVector<double> b(3);
    b(1) = 1;
    b(2) = 6;
    b(3) = 10;

    CVector<double> x(3);
    
    MATRIXTOOLBOX::forward_substitution(A, x, b, MATRIXTOOLBOX::Mode::DEFAULT);

    for (int i = 1; i <= 3; i++) {
        EXPECT_FLOAT_EQ(x(i), 1);
    }


}

TEST_F(MatrixToolboxTest, BackwardSubstitution) {
    MATRIXTOOLBOX::RowReducer(matrix, vector);

    CVector<double> x(3);
    
    MATRIXTOOLBOX::backward_substitution(matrix, x, vector, MATRIXTOOLBOX::Mode::DEFAULT);
    for (int i = 1; i <= 3; i++) {
        EXPECT_FLOAT_EQ(x(i), i);
    }


}

