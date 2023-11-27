# Numerical Toolbox

Welcome to the Numerical Toolbox, a collection of tools designed to provide numerical solutions to various mathematical problems. Whether you're working on linear algebra, solving systems of equations, finding eigenvalues, fitting curves, or tackling differential equations, this toolbox has you covered.

## Basic Vector and Matrix Operations

### CVector: Basic Vector Operations

```cpp
CVector<double> v1{ 1, 2, 3, 4, 5 };
std::cout << v1 << "\n";
CVector<double> v2{ 1, 2, 3, 4, 5 };
std::cout << v1 + v2 << "\n";
std::cout << v1 - v2 << "\n";
std::cout << v1 * v2 << "\n";
```

### Matrix: Basic Matrix Operations
```cpp
Matrix<double> m1 = { 3,3, 1, 2, 3, 4, 5, 6, 7, 8, 9 }; // First two numbers are the dimensions of the matrix
std::cout << m1 << "\n";
Matrix<double> m2 = { 3,3, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
std::cout << m1 + m2 << "\n";
std::cout << m1 - m2 << "\n";
std::cout << m1 * m2 << "\n";

CVector<double> v3 = { 1, 2, 3 };
std::cout << m1 * v3 << "\n";
```

#### You can assign a vector to a row of a matrix. Must be the same size as the row
```cpp
m1(2) = v3;
std::cout << m1 << "\n";
```

## Linear Algebra

### System of Linear Equations

#### x + y + z = 3
#### 2x + 3y + 4z = 20
#### x + y - z = 0

```cpp
Matrix<double> A { 3, 3, 1, 1, 1, 
                        2, 3, 4, 
                        1, 1, -1 };
CVector<double> b { 3, 20, 0 };
CVector<double> x;
MATRIXTOOLBOX::AxEqb(A, x, b);
std::cout << x << "\n";

x = {0, 0, 0};
MATRIXTOOLBOX::LUSolve(A, x, b);
std::cout << x << "\n";

x = {0, 0, 0};
MATRIXTOOLBOX::LDLTSolve(A, x, b);
std::cout << x << "\n";
```

// Random Matrix
// Generate a random matrix 
Matrix<double> rand_mat = MATRIXTOOLBOX::random_matrix(4, 4, -500, 500, MATRIXTOOLBOX::RandMode::POSITIVEDEF, 10);
std::cout << rand_mat << "\n";

// Identity Matrix
// Generate an identity matrix
Matrix<double> identity_mat = MATRIXTOOLBOX::identity(4);

// Transpose
// Transpose a matrix
Matrix<double> test = {3,3,3,2,1,
                            4,2,2,
                            5,6,2};
std::cout << MATRIXTOOLBOX::Transpose(test) << "\n";

//Determinant
// Calculate the determinant of a matrix
double det{};
MATRIXTOOLBOX::Determinant(test, det) ;
//-6
std::cout << det<< "\n";

// Use Inverse Iteration to find the nearest eigenvalue to a of a matrix
A = { 3, 3, 2, 1, 1,
            1, 3, 1, 
            1, 1, 4};
CVector<double> eigenvector = { 1, 1, 1 };
eigenvector/=sqrt(3);
double eigenvalue{5.0};
MATRIXTOOLBOX::EIGEN::vector_iter_inverse(A, eigenvector, 5.0, eigenvalue, 100);
// 5.214319743184
std::cout << "Eigenvalue: " << eigenvalue << "\n";
std::cout << "Eigenvector: " << eigenvector << "\n";


// Solve Eigenvalue problem of the form BX = LAMDBA*X	
Matrix <double> B = {2, 2 , 3,-1,
                            -1, 3};
Matrix <double> LAM;
// add true as the last argument to print the eigenvalues and eigenvectors
MATRIXTOOLBOX::EIGEN::eigen_Jacobi(A, LAM, 100, true);
//Diagonals are the eigenvalues
std::cout << LAM << "\n";


// Solve Eigenvalue problem of the form KX = LAMDBA*M*X
Matrix<double> K = {3,3,3,2,1,2,2,1,1,1,1};
Matrix<double> M(MATRIXTOOLBOX::identity(K.get_rows()));
Matrix<double> X;
Matrix<double> LAMDBA;
std::cout << "Generalized Eigenvalue Problem\n";
MATRIXTOOLBOX::EIGEN::generalized_eigen_Jacobi(K, X, LAMDBA, M, 100, true);
// std::cout << K << LAMDBA * M;


// Least Squares Fit
std::cout << "Least Squares Fit\n";
std::cout << "Function is x^2\n";
LSF lsf;
int terms = 3;
x = CVector<double> {1, 2, 3, 4, 5, 6, 7};
CVector<double> y {1, 4, 9, 16, 25, 36, 49};

CVector<double> coeff(terms);
double residuals;
double r2;
lsf.fit(terms, x, y, coeff, residuals, r2);
std::cout << coeff << "\n";
std::cout << "residuals = " << residuals << "\n";
std::cout << "r^2 = " << r2 << "\n";

//Modify it a bit
std::cout << "Add a bit of noise\n";
x = CVector<double> {1.01, 1.99, 3.01, 3.99, 5.01, 5.99, 7.01};
y = CVector<double> {1.01, 3.98, 9.02, 15.99, 25.01, 35.99, 49.01};

lsf.fit(terms, x, y, coeff, residuals, r2);
std::cout << coeff << "\n";
std::cout << "residuals = " << residuals << "\n";
std::cout << "r^2 = " << r2 << "\n";

//ODE Solver
std::cout << "Solve ODEs\n";
std::cout << "Solve Y''' + Y'' + Y' + Y = 0; Initial conditions: x = 0; Y(0) = 0, Y'(0) = 0, Y''(0) = 0, Y''(0) = 10" << "\n";
// // Decouple into 3 first order ODEs

// // Y1' = Y2
// // Y2' = Y3
// // Y3' = -Y1 - Y2 - Y3

LINODESOLVER::FunctionContainer functions;
functions.addFunction([](CVector<double> points) {return points(3); });
functions.addFunction([](CVector<double> points) {return points(4); });
functions.addFunction([](CVector<double> points) {return -points(2) - points(3) - points(4); });
Matrix<double> solution;
solution = LINODESOLVER::ode45wrapperLinear(functions, 0, 10, { 0, 0, 0, 10 }, 0.1);
std::cout << std::setprecision(15) << solution(solution.get_rows()) << "\n";
std::cout << "Number of points: " << solution.get_rows() << "\n";

//Send to file. Can be plotted in MATLAB or gnuplot if you have it
//solution.to_file("e^2.csv");

// Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");
// solution.plot(true);

// //MATRIXTOOLBOX::Transpose(solution).plot(gp);

// Integaration

auto integral_func = [](int FUNCNO, double x) {return x*x*x; };
auto solution_func = [](double x) {return x*x*x*x/4; };

double test_limits[] = { 1, 1024 * 1024 };
const int FUNC_NO = 1;
const double upper_limit = 79;
const double lower_limit = 10;

print(std::setprecision(15), "Solution: ", solution_func(upper_limit) - solution_func(lower_limit), "\n");
print("Trapezodial: ", Trapezodial(FUNC_NO, 16356, upper_limit, lower_limit, integral_func), "\n");
print("Simpson: ", Simpson(FUNC_NO, 16356, upper_limit, lower_limit, integral_func), "\n");
print("Gauss-Legendre: ", quadrature(FUNC_NO, 5, upper_limit, lower_limit, integral_func), "\n");

//Differentiation

auto diff_func = [](int FUNCNO, double x) {return 1000*atan(x); };
auto diff_solution_func = [](double x) {return 1000/(1+x*x); };

double central_point = 100;
double h = 0.01;
print(std::setprecision(15), "Solution: ", diff_solution_func(central_point), "\n");
print("Forward Difference: ", ForwardMethod(FUNC_NO, central_point, h, diff_func), "\n");
print("Backward Difference: ", BackwardMethod(FUNC_NO, central_point, h, diff_func), "\n");
print("Central Difference: ", CentralMethod(FUNC_NO, central_point, h, diff_func), "\n");
print("Fourth Order Central Difference: ", FourthOrderCentralMethod(FUNC_NO, central_point, h, diff_func), "\n");