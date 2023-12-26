
#include <gtest/gtest.h>
#include "vectortemplate.h"


class CVectorTest : public ::testing::Test {
protected:
    CVector<double> vector;

    void SetUp() override {
        // Optional setup before each test
        vector.SetSize(5);
        vector.Set(0); // Initialize the vector with zeros
    }

    void TearDown() override {
        // Optional teardown after each test
    }
};

// Test default constructor
TEST(VectorTest, DefaultConstructor) {
    CVector<int> vector;
    EXPECT_EQ(vector.GetSize(), 0);
}

// Test constructor with size
TEST(VectorTest, ConstructorWithSize) {
    CVector<double> vector(10);
    EXPECT_EQ(vector.GetSize(), 10);
}

// Test constructor with size and initial value
TEST(VectorTest, ConstructorWithSizeAndInitialValue) {
    CVector<int> vector(3, 100);
    for (int i = 1; i <= vector.GetSize(); i++) {
        EXPECT_EQ(vector(i), 100);
    }
}

// Test copy constructor
TEST(VectorTest, CopyConstructor) {
    CVector<int> original(5, 50);
    CVector<int> copy(original);
    for (int i = 1; i <= copy.GetSize(); i++) {
        EXPECT_EQ(copy(i), 50);
    }
}

// Test setting and getting vector size
TEST(VectorTest, SetGetSize) {
    CVector<double> vector;
    vector.SetSize(8);
    EXPECT_EQ(vector.GetSize(), 8);
}

// Test vector element access
TEST(VectorTest, ElementAccess) {
    CVector<int> vector(5, 10);
    EXPECT_EQ(vector(3), 10);
    vector(3) = 20;
    EXPECT_EQ(vector(3), 20);
}

// Test vector addition
TEST_F(CVectorTest, VectorAddition) {
    CVector<double> other(5, 5.0);
    CVector<double> result = vector + other;
    for (int i = 1; i <= result.GetSize(); i++) {
        EXPECT_EQ(result(i), 5.0);
    }
}

TEST_F(CVectorTest, VectorSubtraction) {
    CVector<double> other(5, 2.0);
    CVector<double> result = vector - other;
    for (int i = 1; i <= result.GetSize(); i++) {
        EXPECT_EQ(result(i), -2.0);
    }
}

// Test scalar subtraction
TEST_F(CVectorTest, ScalarSubtraction) {
    double scalar = 3.5;
    CVector<double> result = vector - scalar;
    for (int i = 1; i <= result.GetSize(); i++) {
        EXPECT_EQ(result(i), -3.5);
    }
}

// Test vector multiplication by scalar
TEST_F(CVectorTest, VectorMultiplicationByScalar) {
    double scalar = 2.5;
    vector.Set(2.0);
    CVector<double> result = vector * scalar;
    for (int i = 1; i <= result.GetSize(); i++) {
        EXPECT_EQ(result(i), 2.5 * 2.0);
    }
}

TEST_F(CVectorTest, VectorMultiplicationByVectorElementWise) {
    CVector<double> vec2(vector.GetSize(), 4.0);
    vector.Set(2.0);
    CVector<double> result = vector * vec2;
    for (int i = 1; i <= result.GetSize(); i++) {
        EXPECT_EQ(result(i), 8.0);
    }
}


// Test vector division by scalar
TEST_F(CVectorTest, VectorDivisionByScalar) {
    double scalar = 0.5;
    CVector<double> result = vector / scalar;
    for (int i = 1; i <= result.GetSize(); i++) {
        EXPECT_EQ(result(i), 0.0);
    }
}

TEST_F(CVectorTest, VectorDivisionByVectorElementWise) {
    CVector<double> vec2(vector.GetSize(), 4.0);
    vector.Set(2.0);
    CVector<double> result;
    result = vector / vec2;
    for (int i = 1; i <= result.GetSize(); i++) {
        EXPECT_EQ(result(i), 0.5);
    }
}

// Test assignment operator
TEST_F(CVectorTest, AssignmentOperator) {
    CVector<double> other(5, 100);
    vector = other;
    for (int i = 1; i <= vector.GetSize(); i++) {
        EXPECT_EQ(vector(i), 100);
    }
}

// Test initializer list
TEST_F(CVectorTest, InitializerList) {
    std::initializer_list<double> v = { 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0 };
    vector = v;
    EXPECT_EQ(vector.GetSize(), 10);
    for (int i = 1; i <= vector.GetSize(); i++) {
        EXPECT_EQ(vector(i), i);
    }
    std::cout << vector;
}

// Test ostream operator (<<) overload
TEST_F(CVectorTest, OstreamOperatorOverload) {
    std::stringstream ss;
    ss << vector;
    EXPECT_EQ(ss.str(), "[ 0 0 0 0 0 ]");
}

// Test scalar addition
TEST_F(CVectorTest, ScalarAddition) {
    double scalar = 3.5;
    CVector<double> result = vector + scalar;
    for (int i = 1; i <= result.GetSize(); i++) {
        EXPECT_EQ(result(i), 3.5);
    }
}

TEST_F(CVectorTest, CompoundAddition) {
    CVector<double> other(5, 5);
    vector += other;
    for (int i = 1; i <= vector.GetSize(); i++) {
        EXPECT_EQ(vector(i), 5);
    }
}

// Test compound subtraction
TEST_F(CVectorTest, CompoundSubtraction) {
    CVector<double> other(5, 2);
    vector -= other;
    for (int i = 1; i <= vector.GetSize(); i++) {
        EXPECT_EQ(vector(i), -2);
    }
}

TEST_F(CVectorTest, CompoundMultiplication) {
    vector.Set(2);
    vector *= 5;
    for (int i = 1; i <= vector.GetSize(); i++) {
        EXPECT_EQ(vector(i), 10.0);
    }
}

// Test compound subtraction
TEST_F(CVectorTest, CompoundDivision) {
    vector.Set(2);
    vector /= 1.25;
    for (int i = 1; i <= vector.GetSize(); i++) {
        EXPECT_EQ(vector(i), 2.0 / 1.25);
    }
}

TEST_F(CVectorTest, CompoundVectorMultiplicationByVectorElementWise) {
    CVector<double> vec2(vector.GetSize(), 4.0);
    vector.Set(2.0);
    vector *= vec2;
    for (int i = 1; i <= vector.GetSize(); i++) {
        EXPECT_EQ(vector(i), 8.0);
    }
}

TEST_F(CVectorTest, CompoundVectorDivisionByVectorElementWise) {
    CVector<double> vec2(vector.GetSize(), 4.0);
    vector.Set(2.0);
    vector /= vec2;
    for (int i = 1; i <= vector.GetSize(); i++) {
        EXPECT_EQ(vector(i), 0.5);
    }
}

// Test incompatible vectors for addition
TEST_F(CVectorTest, IncompatibleVectorsAddition) {
    CVector<double> other(3, 5); // Create a vector of different size
    ASSERT_THROW(vector + other, std::exception);
}

// Test incompatible vectors for subtraction
TEST_F(CVectorTest, IncompatibleVectorsSubtraction) {
    CVector<double> other(3, 2); // Create a vector of different size
    ASSERT_THROW(vector - other, std::exception);
}

// Test for divide by Zero
TEST_F(CVectorTest, DivideByZero) {

    ASSERT_THROW(vector / 0, std::exception);
}

// Test for divide by Zero
TEST_F(CVectorTest, CompoundDivideByZero) {
    ASSERT_THROW(vector /= 0, std::exception);
}

// Test for bad vector creation
TEST_F(CVectorTest, Set) {
    ASSERT_THROW(CVector<double>(-1), std::exception);
}