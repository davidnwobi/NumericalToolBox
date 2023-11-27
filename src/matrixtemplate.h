#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <iomanip>
#include <tuple>
#include <boost/optional.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <stdexcept>
#include <string>
#include <vector>
#include "vectortemplate.h"
#include "gnuplot-iostream.h"

template <class T>
class Matrix;

template <class T>
std::ostream& operator<<(std::ostream&, const Matrix<T>&);

/**
 * @brief Represents a mathematical matrix.
 *
 * @tparam T The type of elements in the matrix.
 */
template <class T>
class Matrix {
public:
    enum class Error {
        INVALIDSIZE,
        MEMORYALLOCATION,
        INVALIDINDEX,
        INCOMPATIBLEMATRICES,
        DIVISIONBYZERO,
        NODIMMATRIX,
        BADFILEOPENING,
        BADFILEREAD,
        BADFILEWRITE
    };

    Matrix();
    Matrix(int nRows, int nCols);
    Matrix(int nRows, int nCols, T init_val);
    Matrix(const CVector<T>& vector);
    Matrix(const Matrix& mat);
    Matrix(Matrix&& mat) noexcept;
    Matrix(std::initializer_list<T> elements);
    ~Matrix();

    void operator=(const Matrix& mat);
    void operator=(Matrix&& mat) noexcept;
    void operator=(const T& val);
    void operator=(std::initializer_list<T> elements);

    T& operator()(int i, int j);
    T operator()(int i, int j) const;
    CVector<T>& operator()(int i);
    CVector<T> operator()(int i) const;

    Matrix operator+(const Matrix& mat) const;
    Matrix operator+(Matrix&& mat) const noexcept;
    Matrix operator+(const T& val) const;
    void operator+=(const Matrix& mat);
    void operator+=(const T& val);

    Matrix operator-(const Matrix& mat) const;
    Matrix operator-(Matrix&& mat) const noexcept;
    Matrix operator-(const T& val) const;
    void operator-=(const Matrix& mat);
    void operator-=(const T& val);

    CVector<T> operator*(const CVector<T>& vec) const;
    Matrix operator*(const Matrix& mat) const ;
    Matrix operator*(const T& val) const;
    void operator*=(const CVector<T>& vec);
    void operator*=(const Matrix& mat);
    void operator*=(const T& val);

    Matrix operator/(const T& val) const;
    void operator/=(const T& val);


    bool operator==(const Matrix& mat) const;
    bool operator!=(const Matrix& mat) const;

    int get_rows() const;
    int get_columns() const;



    friend std::ostream& operator<< <T>(std::ostream& os, const Matrix<T>& matrix);
    void SetSize(int nRows, int nCols);
    void to_file(const std::string& loc);
    inline Matrix from_file(const std::string& loc);
    void plot(bool first_is_x = false, boost::optional<int> col = boost::none);
    std::vector<std::vector<T>> convert_to_vector();
    void SwapRows(int from,int to);

private:
    CVector<CVector<T>> m_pData;
    int m_nRows;
    int m_nCols;
    void InitializeData(int nRows, int nCols, boost::optional<T> init_val = boost::none);
    void CopyData(const Matrix& mat);
    void ErrorHandler(Error err) const;
    
};


template <class T>
Matrix<T> Matrix<T>::from_file(const std::string& loc) {
    std::ifstream infile(loc);
    if (!infile.is_open()) {
        ErrorHandler(Error::BADFILEOPENING);
    }
    std::vector<CVector<double>> temp_mat;
    while (!infile.eof() && !infile.fail())
     {
        std::vector<double> data;
        std::string line;
        std::getline(infile, line);
        if (!line.empty() && line.find_first_not_of(" \t\n\v\f\r") != std::string::npos) {
            std::vector<std::string> tokens;
            boost::split(tokens, line, boost::is_any_of(","));

            for (const auto& token : tokens) {
                try {
                    double value = boost::lexical_cast<double>(token);
                    data.push_back(value);
                }
                catch (const boost::bad_lexical_cast& e) {
                    throw std::runtime_error("Error: Invalid data format in the file.");
                }
            }
            CVector<double> temp(data.size());
            size_t i = 1;
            for (auto& x : data) {
                temp(i) = data.at(i);
            }
            temp_mat.push_back(std::move(temp));
        }
        
    }
    Matrix<T> mat(temp_mat.size(), temp_mat.front().GetSize());
    size_t i = 1;
    for (auto& x : temp_mat) {
        mat(i) = std::move(temp_mat.at(i));
    }
    
    infile.close();
    return mat;

}

template <class T>
void Matrix<T>::to_file(const std::string& loc) {
    std::ofstream Outfile(loc, std::ios::out | std::ios::trunc);

    if (!Outfile.is_open()) {
        ErrorHandler(Error::BADFILEOPENING); // Handle file opening error
    }

    for (int i = 1; i <= m_nRows; i++) {
        Outfile << std::setiosflags(std::ios::left);
        for (int j = 1; j <= m_nCols; j++) {
            Outfile << m_pData(i)(j) << ",";
        }

        // Check for write errors after writing each row
        if (Outfile.fail()) {
            ErrorHandler(Error::BADFILEWRITE); // Handle file write error
        }

        Outfile << "\n";

        // Check for write errors after writing the newline character
        if (Outfile.fail()) {
            ErrorHandler(Error::BADFILEWRITE); // Handle file write error
        }
    }

    Outfile.close();
}

template <class T>
void Matrix<T>::InitializeData(int nRows, int nCols, boost::optional<T> init_val) {
    try {
        if ((nRows < 1) || (nCols < 1)) {
            ErrorHandler(Error::INVALIDSIZE);
        }
        else {
            m_nRows = nRows;
            m_nCols = nCols;
            if (init_val.has_value()) {
                m_pData = CVector<CVector<T>>(nRows);
                for (int i = 1; i <= m_nRows; i++) {
                    CVector<T> temp(nCols, init_val.get());
                    m_pData(i) = temp;
                }
            }
            else {
                m_pData = CVector<CVector<T>>(nRows);
                for (int i = 1; i <= m_nRows; i++) {
                    CVector<T> temp(nCols);
                    m_pData(i) = temp;
                }
            }
        }
    }
    catch (std::bad_alloc) {
        ErrorHandler(Error::MEMORYALLOCATION);
    }
}
template <class T>
void Matrix<T>::CopyData(const Matrix& mat) {
    try {
        m_nRows = mat.get_rows();
        m_nCols = mat.get_columns();
        for (int i = 1; i <= m_nRows; i++) {
            m_pData = mat;
        }

    }
    catch (std::bad_alloc) {
        ErrorHandler(Error::MEMORYALLOCATION);
    }
}
template <class T>
void Matrix<T>::SetSize(int nRows, int nCols) {
    InitializeData(nRows, nCols);
}

template <class T>
int Matrix<T>::get_rows() const {
    return m_nRows;
}

template <class T>
int Matrix<T>::get_columns() const {
    return m_nCols;
}

template <class T>
Matrix<T>::Matrix() : m_pData(CVector<CVector<T>>()), m_nRows(0), m_nCols(0) {}

template <class T>
Matrix<T>::Matrix(int nRows, int nCols) {
    InitializeData(nRows, nCols);
}

template <class T>
Matrix<T>::Matrix(int nRows, int nCols, T init_val) {
    InitializeData(nRows, nCols, init_val);
}

template <class T>
Matrix<T>::Matrix(const CVector<T>& vector) {
    if (vector.GetSize() < 1) {
        ErrorHandler(Error::INVALIDSIZE);
    }
    m_nRows = vector.GetSize();
    m_nCols = 1;
    m_pData = CVector<CVector<T>>(m_nRows);
    for (int i = 1; i <= m_nRows; i++) {
        m_pData(i) = CVector<T>(1, vector(i));
    }
}


template <class T>
Matrix<T>::Matrix(const Matrix& mat) {

    m_nRows = mat.m_nRows;
    m_nCols = mat.m_nCols;
    m_pData = mat.m_pData;
}

template <class T>
Matrix<T>::Matrix(Matrix&& mat) noexcept {
    m_nRows = mat.m_nRows;
    m_nCols = mat.m_nCols;
    m_pData = std::move(mat.m_pData);
    mat.m_pData = CVector<CVector<T>> ();
    mat.m_nRows = 0;
    mat.m_nCols = 0;
}
template <class T>
Matrix<T>::Matrix(std::initializer_list<T> elements) {
    std::vector<T> vec_temp = elements;
    // Check if the first two elements are the dimensions of the matrix
    if (static_cast<int>(vec_temp[0]) * static_cast<int>(vec_temp[1]) + 2 != elements.size()) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
    }
    Matrix<T> temp(static_cast<int>(vec_temp[0]), static_cast<int>(vec_temp[1])); // check if this is possible
    for (size_t i = 1; i <= temp.get_rows(); i++) {
        for (size_t j = 1; j <= temp.get_columns(); j++) {
            temp(i, j) = vec_temp[2 + (i - 1) * temp.get_columns() + (j - 1)];
        }
    }
    *this = std::move(temp);

}


template <class T>
Matrix<T>::~Matrix() {
   
}

template <class T>
void Matrix<T>::operator=(const Matrix& mat) {
    
    m_nRows = mat.m_nRows;
    m_nCols = mat.m_nCols;
    m_pData = mat.m_pData;
}

template <class T>
void Matrix<T>::operator=(Matrix&& mat) noexcept{
    
    m_nRows = mat.m_nRows;
    m_nCols = mat.m_nCols;
    m_pData = std::move(mat.m_pData);
    mat.m_pData = CVector<CVector<T>>();
    mat.m_nRows = 0;
    mat.m_nCols = 0;
}

template <class T>
void Matrix<T>::operator=(const T& val) {
    if ((m_nRows < 1) || (m_nCols < 1)) {
        ErrorHandler(Error::NODIMMATRIX);
    }
    for (int i = 1; i <= m_nRows * m_nCols; i++) {
        m_pData->operator()(i) = val;
    }
}
template <class T>
void Matrix<T>::operator=(std::initializer_list<T> elements) {
    std::vector<T> vec_temp = elements;
    if (static_cast<int>(vec_temp[0]) * static_cast<int>(vec_temp[1])+2 != elements.size()) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
    }
    Matrix<T> temp(static_cast<int>(vec_temp[0]), static_cast<int>(vec_temp[1])); // check if this is possible
    for (size_t i = 1; i <= temp.get_rows(); i++) {
        for (size_t j = 1; j <= temp.get_columns(); j++) {
            temp(i, j) = vec_temp[2 + (i - 1) * temp.get_columns() + (j-1)];
        }
    }
    *this = std::move(temp);

}

template <class T>
T& Matrix<T>::operator()(int i, int j) {
    return m_pData(i)(j);
}

template <class T>
T Matrix<T>::operator()(int i, int j) const {
    return m_pData(i)(j);
}

template <class T>
CVector<T> Matrix<T>::operator()(int i) const {
    return m_pData(i);
}

template <class T>
CVector<T>& Matrix<T>::operator()(int i) {
    return m_pData(i);
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
    os << matrix.get_rows() << "x" << matrix.get_columns() << " Matrix:" << std::endl;
    os << "[ " << std::endl;
    int rows = matrix.get_rows();
    int cols = matrix.get_columns();
    for (int i = 1; i <= rows; i++) {
        os << std::setprecision(4) << std::setw(4) << matrix.m_pData(i) << "\n";
    }
    os << " ]" << std::endl;
    return os;
}

template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix& mat) const {
    if ((m_nRows != mat.m_nRows) || (m_nCols != mat.m_nCols)) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
        return Matrix<T>();
    }
    Matrix<T> temp(*this);
    for (int i = 1; i <= m_nRows; i++) {
        temp(i) += mat(i);
    }
    return temp;
}

template <class T>
Matrix<T> Matrix<T>::operator+(Matrix&& mat) const noexcept {
    if ((m_nRows != mat.m_nRows) || (m_nCols != mat.m_nCols)) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
        return Matrix<T>();
    }
    Matrix<T> temp(std::move(mat));
    for (int i = 1; i <= m_nRows; i++) {
        temp(i) += this->operator()(i);
    }
    return temp;
}

template <class T>
Matrix<T> Matrix<T>::operator+(const T& val) const {
    Matrix<T> temp(*this);
    for (int i = 1; i <= m_nRows; i++) {
        for (int j = 1; j <= m_nCols; j++)
            temp(i, j) += val;
    }
    return temp;
}

template <class T>
void Matrix<T>::operator+=(const Matrix& mat) {
    if ((m_nRows != mat.m_nRows) || (m_nCols != mat.m_nCols)) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
        return;
    }
    for (int i = 1; i <= m_nRows; i++) {
        this->operator()(i) += mat(i);
    }

}
template <class T>
void Matrix<T>::operator+=(const T& val) {
    for (int i = 1; i <= m_nRows; i++) {
        for (int j = 1; j <= m_nCols; j++) {
            this->operator()(i, j) += val;
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix& mat) const{
    if ((m_nRows != mat.m_nRows) || (m_nCols != mat.m_nCols)) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
        return Matrix<T>();
    }
    Matrix<T> temp(*this);
    for (int i = 1; i <= m_nRows; i++) {
        temp(i) -=  mat(i);
    }
    return temp;
}

template <class T>
Matrix<T> Matrix<T>::operator-(Matrix&& mat) const noexcept {
    if ((m_nRows != mat.m_nRows) || (m_nCols != mat.m_nCols)) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
        return Matrix<T>();
    }
    Matrix<T> temp(std::move(mat));
    for (int i = 1; i <= m_nRows; i++) {
        temp(i) -= this->operator()(i);
    }
    return temp;
}

template <class T>
Matrix<T> Matrix<T>::operator-(const T& val) const {
    Matrix<T> temp(*this);
    for (int i = 1; i <= m_nRows; i++) {
        for (int j = 1; j <= m_nCols; j++)
            temp(i, j) -= val;
    }
    return temp;
}

template <class T>
void Matrix<T>::operator-=(const Matrix& mat) {
    if ((m_nRows != mat.m_nRows) || (m_nCols != mat.m_nCols)) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
        return;
    }
    for (int i = 1; i <= m_nRows; i++) {
        this->operator()(i) -= mat(i);
    }

}
template <class T>
void Matrix<T>::operator-=(const T& val) {
    for (int i = 1; i <= m_nRows; i++) {
        for (int j = 1; j <= m_nCols; j++) {
            this->operator()(i, j) -= val;
        }
    }
}



template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix& mat) const {
    if (m_nCols != mat.m_nRows) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
        return Matrix<T>();
    }

    Matrix<T> result(m_nRows, mat.m_nCols);
    for (int i = 1; i <= m_nRows; i++) {
        for (int j = 1; j <= mat.m_nCols; j++) {
            T sum = 0;
            for (int k = 1; k <= m_nCols; k++) {
                sum += this->operator()(i, k) * mat(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

template <class T>
CVector<T> Matrix<T>::operator*(const CVector<T>& vec) const {
    if (m_nCols != vec.GetSize()) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
        return CVector<T>();
    }

    CVector<T> result(m_nRows);
    for (int i = 1; i <= m_nRows; i++) {
        T sum = static_cast<T>(0);
        for (int j = 1; j <= m_nCols; j++) {
            sum += m_pData(i)(j) * vec(j);
        }
        result(i) = sum;
    }
    return result;
}

template <class T>
Matrix<T> Matrix<T>::operator*(const T& val) const {
    Matrix<T> result(*this);
    for (int i = 1; i <= m_nRows; i++) {
        for (int j = 1; j <= m_nCols; j++) {
            result(i, j) *= val;
        }
    }
    return result;
}

template <class T>
void Matrix<T>::operator*=(const CVector<T>& vec) {
    if (m_nCols != vec.GetSize()) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
        return;
    }
    CVector<T> result(m_nRows);
    for (int i = 1; i <= m_nRows; i++) {
        T sum = static_cast<T>(0);
        for (int j = 1; j <= m_nCols; j++) {
            sum += m_pData(i)(j) * vec(j);
        }
        result(i) = sum;
    }
    Matrix<T> final_mat(result);
    *this = std::move(final_mat);
}

template <class T>
void Matrix<T>::operator*=(const Matrix& mat) {
    if (m_nCols != mat.m_nRows) {
        ErrorHandler(Error::INCOMPATIBLEMATRICES);
        return;
    }

    Matrix<T> result(m_nRows, mat.m_nCols);

    for (int i = 1; i <= m_nRows; i++) {
        for (int j = 1; j <= mat.m_nCols; j++) {
            T sum = static_cast<T>(0);
            for (int k = 1; k <= m_nCols; k++) {
                sum += m_pData(i)(k) * mat(k, j);
            }
            result(i, j) = sum;
        }
    }

    *this = std::move(result);
}

template <class T>
void Matrix<T>::operator*=(const T& val) {
    for (int i = 1; i <= m_nRows; i++) {
        for (int j = 1; j <= m_nCols; j++) {
            this->operator()(i, j) *= val;
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::operator/(const T& val) const {
    if (val == 0) {
        ErrorHandler(Error::DIVISIONBYZERO);
        return Matrix<T>();
    }

    Matrix<T> result(*this);
    for (int i = 1; i <= m_nRows; i++) {
        for (int j = 1; j <= m_nCols; j++) {
            result(i, j) /= val;
        }
    }
    return result;
}


template <class T>
void Matrix<T>::operator/=(const T& val) {
    if (val == 0) {
        ErrorHandler(Error::DIVISIONBYZERO);
        return;
    }

    for (int i = 1; i <= m_nRows; i++) {
        for (int j = 1; j <= m_nCols; j++) {
            this->operator()(i, j) /= val;
        }
    }
}

template <class T>
bool Matrix<T>::operator==(const Matrix& mat) const {
    if (m_nRows != mat.m_nRows || m_nCols != mat.m_nCols) {
        return false;
    }

    for (int i = 1; i <= m_nRows; i++) {
        for (int j = 1; j <= m_nCols; j++) {
            if (this->operator()(i, j) != mat(i, j)) {
                return false;
            }
        }
    }
    return true;
}


template <class T>
bool Matrix<T>::operator!=(const Matrix& mat) const{
    return !(*this == mat);
}


template <class T>
void Matrix<T>::ErrorHandler(Error err) const {

    switch (err) {
        case Error::INCOMPATIBLEMATRICES:
            throw std::runtime_error("Matrix operation cannot take place due to incompatible matrices.");
        case Error::INVALIDINDEX:
            throw std::out_of_range("Invalid index to access matrix element.");
        case Error::INVALIDSIZE:
            throw std::invalid_argument("Invalid matrix size.");
        case Error::MEMORYALLOCATION:
            throw std::bad_alloc();
        case Error::DIVISIONBYZERO:
            throw std::runtime_error("Division by Zero.");
        case Error::NODIMMATRIX:
            throw std::runtime_error("Matrix has no dimensions.");
        case Error::BADFILEOPENING:
            throw std::runtime_error("Error: Failed to open the file.");
        case Error::BADFILEREAD:
            throw std::runtime_error("Error: Failed to read from the file.");
        case Error::BADFILEWRITE:
            throw std::runtime_error("Error: Failed to write to the file.");
        default:
            throw std::runtime_error("Unknown error occurred.");
    }
}

template <class T>
std::vector<std::vector<T>> Matrix<T>::convert_to_vector(){
    std::vector<std::vector<T>> temp(m_nRows, std::vector<double>(m_nCols));

    for (int i = 1; i <= m_nRows; i++) {
        for (int j = 1; j <= m_nCols; j++) {
            temp[i - 1][j - 1] = this->operator()(i,j);
        }
    }
    return temp;
}
template <class T> 
void Matrix<T>::plot(bool first_is_x, boost::optional<int> col) {
    Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");
    std::vector<std::vector<T>> temp(convert_to_vector());
    std::vector<std::string> column_names;
    auto plots = gp.plotGroup();

    if (col.has_value()) {
        if (first_is_x) {
            std::vector<std::vector<double>> single_column_data(2);
            single_column_data[0] = temp[0];
            single_column_data[1] = temp[col.get()];
            std::string strMessage = std::string("using 1:2 with lines title 'Col ") + std::to_string(col.get()) + "'";
            plots.add_plot1d_colmajor(single_column_data, strMessage);
        }
        else {
            std::string strMessage = std::string("with lines title 'Col ") + std::to_string(col.get()) + "'";
            plots.add_plot1d(temp[col.get()], strMessage);
        }
    }
    else if (m_nRows == 1) {
        std::string strMessage = std::string("with lines title 'Col ") + std::to_string(col.get()) + "'";
        plots.add_plot1d(temp[0], strMessage);
    }
    else {
        if (first_is_x) {
            for (size_t j = 1; j < temp.size(); j++) {
                std::vector<std::vector<double>> single_column_data(2);
                single_column_data[0] = temp[0];
                single_column_data[1] = temp[j];
                std::string strMessage = std::string("using 1:2 with lines title 'line ") + std::to_string(j) + "'";
                plots.add_plot1d_colmajor(single_column_data, strMessage);
            }
        }
        else {
            for (size_t j = 0; j < temp.size(); j++) {
                std::string strMessage = std::string("with lines title 'line ") + std::to_string(j) + "'";
                plots.add_plot1d(temp[j], strMessage);
            }
        }
    }
    gp << plots;

    std::cout << "Press enter to exit." << std::endl;
    std::cin.get();

}
template <class T>
void Matrix<T>::SwapRows(int from, int to) {
    CVector temp = std::move(this->operator()(from));
    
    this->operator()(from) = std::move(this->operator()(to));
    this->operator()(to) = std::move(temp);

}
