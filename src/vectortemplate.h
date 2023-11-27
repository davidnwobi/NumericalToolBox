// ********************************************************************************
// 
// I have adapted this class from 
// Vector Template Class
// Copyright(c), 2000-21, S. D. Rajan
//
// Implements vector as a simple container
// Indexing starts at 1

#pragma once
#include <string>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <initializer_list>
#include <memory>
#include <vector>



template <class T>
class CVector;

template <class T>
std::ostream& operator<<(std::ostream&, const CVector<T>&);


/**
 * @class CVector
 * @brief Represents a vector template class.
 * @tparam T The type of elements in the vector.
 */
template <class T>
class CVector


{
public:
    /**
        * @brief Enumerates the possible errors that can occur.
        */
    enum class Error {
        INVALIDSIZE,           /**< Invalid size */
        MEMORYALLOCATION,      /**< Memory allocation failure */
        INVALIDINDEX,          /**< Invalid index */
        INCOMPATIBLEVECTORS,    /**< Incompatible vectors */
        INVALIDOPERATION,       /**< Invalid Operation */
        BADSLICE
    };
    /**
        * @brief Default constructor. Constructs an empty vector.
        */
    CVector();

    /**
        * @brief Overloaded constructor. Constructs a vector with the specified size.
        *
        * @param size The size of the vector.
        */
    CVector(int size);

    /**
        * @brief Overloaded constructor. Constructs a vector with the specified size, initialized with a value.
        *
        * @param size The size of the vector.
        * @param iv The initial value for all elements in the vector.
        */
    CVector(int size, T iv);

    /**
        * @brief Overloaded constructor. Constructs a vector with the specified name.
        *
        * @param name The name of the vector.
        */
    CVector(const char* name);

    /**
        * @brief Overloaded constructor. Constructs a vector with the specified name and size.
        *
        * @param name The name of the vector.
        * @param size The size of the vector.
        */
    CVector(const char* name, int size);

    /**
        * @brief Overloaded constructor. Constructs a vector with the specified name, size, and initial value.
        *
        * @param name The name of the vector.
        * @param size The size of the vector.
        * @param iv The initial value for all elements in the vector.
        */
    CVector(const char* name, int size, T iv);

    /**
        * @brief Copy constructor. Constructs a vector as a copy of another vector.
        *
        * @param vec The vector to be copied.
        */
    CVector(const CVector<T>& vec);

    /**
        * @brief Move constructor. Constructs a vector from another another vector.
        *
        * @param vec The vector to be copied from.
        */
    CVector(CVector<T>&& vec) noexcept;

    CVector(std::initializer_list<T> elements);

    CVector(T* data, int start, int end);

    /**
        * @brief Destructor. Releases the dynamically allocated memory used by the vector.
        */
    ~CVector();

    /**
        * @brief Sets the size of the vector.
        *
        * @param size The size of the vector.
        */
    void SetSize(int size);

    /**
        * @brief Gets the current size of the vector.
        *
        * @return The current size of the vector.
        */
    int GetSize() const;

    /**
        * @brief Gets the name of the vector.
        *
        * @param strName A string to store the vector name.
        */
    void GetName(std::string& strName) const;

    /**
        * @brief Displays the content of the vector along with an optional message.
        *
        * @param Message The optional message to display before the vector content.
        */
    void Display(const std::string& Message) const;

    /**
        * @brief Sets the value of all elements in the vector.
        *
        * @param val The value to set for all elements.
        */
    void Set(T val);

    /**
        * @brief Sets the name of the vector.
        *
        * @param name The name of the vector.
        */
    void SetName(const std::string& name);

    /**
        * @brief Accesses the element at the specified position in the vector for both reading and writing.
        *
        * @param index The index of the element.
        * @return A reference to the element at the specified position.
        */
    T& operator()(int index);

    /**
        * @brief Accesses the element at the specified position in the vector for reading.
        *
        * @param index The index of the element.
        * @return A const reference to the element at the specified position.
        */
    T operator()(int index) const;

    CVector<T> operator()(int i, int j);

    /**
        * @brief Assigns the contents of another vector to this vector.
        *
        * @param vec The vector to be assigned.
        * @return A reference to the modified vector.
        */
    void operator=(const CVector<T>& vec);

    void operator=(CVector&& matarg) noexcept;

    void operator=(std::initializer_list<T> elements);

    /**
    * @brief Checks if this vector is to another vector
    *
    * @param vec The vector to be checked against.
    * @return A boolean
    */

    bool operator==(const CVector<T>& vec) const;

    /**
        * @brief Performs element-wise addition of two vectors.
        *
        * @param vec The vector to add.
        * @return A new vector resulting from the addition.
        */
    CVector<T> operator+(const CVector<T>& vec) const;

    /**
    * @brief Performs element-wise addition of two vectors.
    *
    * @param vec The vector to add.
    * @return A new vector resulting from the addition.
    */
    CVector<T> operator+(CVector<T>&& vec) const noexcept;

    /**
        * @brief Performs element-wise addition of a scalar value to each element in the vector.
        *
        * @param val The scalar value to add.
        * @return A new vector resulting from the addition.
        */
    CVector<T> operator+(const T& val) const;

    /**
        * @brief Adds another vector to this vector, modifying this vector.
        *
        * @param vec The vector to add.
        * @return A reference to the modified vector.
        */
    void operator+=(const CVector<T>& vec);

    /**
        * @brief Adds a scalar value to each element in the vector, modifying this vector.
        *
        * @param val The scalar value to add.
        * @return A reference to the modified vector.
        */
    void operator+=(const T& val);

    /**
        * @brief Performs element-wise subtraction of two vectors.
        *
        * @param vec The vector to subtract.
        * @return A new vector resulting from the subtraction.
        */
    CVector<T> operator-(const CVector<T>& vec) const;

    /**
        * @brief Performs element-wise subtraction of two vectors.
        *
        * @param vec The vector to subtract.
        * @return A new vector resulting from the subtraction.
        */
    CVector<T> operator-(CVector<T>&& vec) const;

    /**
        * @brief Performs element-wise subtraction of a scalar value from each element in the vector.
        *
        * @param val The scalar value to subtract.
        * @return A new vector resulting from the subtraction.
        */
    CVector<T> operator-(const T& val) const;

    /**
        * @brief Subtracts another vector from this vector, modifying this vector.
        *
        * @param vec The vector to subtract.
        * @return A reference to the modified vector.
        */
    void operator-=(const CVector<T>& vec);

    /**
        * @brief Subtracts a scalar value from each element in the vector, modifying this vector.
        *
        * @param val The scalar value to subtract.
        * @return A reference to the modified vector.
        */
    void operator-=(const T& val);

    /**
        * @brief Performs element-wise multiplication of the vector by a scalar value.
        *
        * @param val The scalar value to multiply by.
        * @return A new vector resulting from the multiplication.
        */
    CVector<T> operator*(const T& val) const;

    /**
        * @brief Performs element-wise multiplication of two vectors.
        *
        * @param vec The vector to multipy with.
        * @return A new vector resulting from the multiplication.
        */
    CVector<T> operator*(const CVector<T>& vec) const;

    /**
        * @brief Performs element-wise multiplication of two vectors.
        *
        * @param vec The vector to multipy with.
        * @return A new vector resulting from the multiplication.
        */
    CVector<T> operator*(CVector<T>&& vec) const;

    /**
        * @brief Performs element-wise multiplication of the vector by a scalar value, modifying this vector.
        *
        * @param val The scalar value to multiply by.
        * @return A reference to the modified vector.
        */
    void operator*=(const T& val);

    /**
    * @brief Multipy this vector's element by another vector's element, modifying this vector.
    *
    * @param vec The vector to multipy with.
    * @return A reference to the modified vector.
    */
    void operator*=(const CVector<T>& vec);

    /**
        * @brief Performs element-wise division of the vector by a scalar value.
        *
        * @param val The scalar value to divide by.
        * @return A new vector resulting from the division.
        */
    CVector<T> operator/(const T& val) const;

    /**
        * @brief Performs element-wise division of two vectors.
        *
        * @param vec The vector to divide by.
        * @return A new vector resulting from the division.
        */
    CVector<T> operator/(const CVector<T>& vec) const;

    /**
        * @brief Performs element-wise division of two vectors.
        *
        * @param vec The vector to divide by.
        * @return A new vector resulting from the division.
        */
    CVector<T> operator/(CVector<T>&& vec) const noexcept;

    /**
        * @brief Performs element-wise division of the vector by a scalar value, modifying this vector.
        *
        * @param val The scalar value to divide by.
        * @return A reference to the modified vector.
        */
    void operator/=(const T& val);


    /**
    * @brief divide this vector's element by another vector's element, modifying this vector.
    *
    * @param vec The vector to divide by.
    * @return A reference to the modified vector.
    */
    void operator/=(const CVector<T>& vec);

    /**
    * @brief Outputs the vector to the output stream.
    *
    * @tparam T The type of elements in the vector.
    * @param os The output stream.
    * @param vec The vector to output.
    * @return The output stream.
    */
    friend std::ostream& operator<< <T>(std::ostream& os, const CVector<T>& vec);

    /**
        * @brief Releases the dynamically allocated memory used by the vector.
        */
    void Release();

private:
    int m_nRows;             /**< Number of rows in the vector */
    std::vector<T> m_pCells;             /**< Address where the vector of type T is stored */
    std::string m_strName;   /**< Vector name */
    bool is_mem_view{false};
    /**
     * @brief Handles the specified error.
     *
     * @param err The error to handle.
     */
    void ErrorHandler(CVector<T>::Error err) const;


};

template <class T>
std::ostream& operator<<(std::ostream& os, const CVector<T>& vector) {
    os << "[ ";
    os << std::fixed << std::setprecision(4) << std::setw(8);
    for (int i = 1; i <= vector.m_nRows; i++) {
        
        os << vector(i) << " ";
    }
    os << "]";
    return os;
}


// =============== definitions ===========================================

template <class T>
CVector<T>::CVector()
// ---------------------------------------------------------------------------
// Function: default ctor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nRows = 0;
    is_mem_view = false;
}

template <class T>
CVector<T>::CVector(int size)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    size of the vector
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nRows = 0;
    SetSize(size);
    is_mem_view = false;

}

template <class T>
CVector<T>::CVector(int size, T iv)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    size of the vector, initial value
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nRows = 0;
    SetSize(size);
    Set(iv);
    is_mem_view = false;
}

template <class T>
CVector<T>::CVector(const char* strName)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    vector name
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nRows = 0; 
    m_strName = strName;
    is_mem_view = false;
}

template <class T>
CVector<T>::CVector(const char* strName, int size)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    vector name and size
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nRows = 0;
    m_strName = strName;
    SetSize(size);
    is_mem_view = false;
}

template <class T>
CVector<T>::CVector(const char* strName, int size, T iv)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    vector name, size and initial value
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nRows = 0;
    m_strName = strName;
    SetSize(size);
    if (!isnan(iv)) {
        Set(iv);
    }
    is_mem_view = false;

}



template <class T>
CVector<T>::CVector(const CVector<T>& A)
// ---------------------------------------------------------------------------
// Function: copy ctor
// Input:    vector
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nRows = A.GetSize();
    m_strName = A.m_strName;
    is_mem_view = false;
    SetSize(m_nRows);
    for (int i = 1; i <= m_nRows; i++)
        m_pCells[i] = A.m_pCells[i];
}

template <class T>
CVector<T>::CVector(CVector<T>&& vec) noexcept {
    m_pCells = std::move(vec.m_pCells);
    m_nRows = vec.m_nRows;
    m_strName = std::move(vec.m_strName);
    is_mem_view = vec.is_mem_view;
    vec.m_nRows = 0;
}

template <class T>
CVector<T>::CVector(T* data, int start, int end) {
    try {
        if (start >= end) {
            ErrorHandler(Error::BADSLICE);
        }
        m_pCells = data + (start - 1);
        m_nRows = end - start;
        is_mem_view = true;
    }
    catch (std::exception e) {
        throw;
    }


    // This constructor creates a view of the row.
    // It doesn't allocate any new memory; instead, it points to the appropriate row.
}
template <class T>
CVector<T>::CVector(std::initializer_list<T> elements) {
    

    Release();
    SetSize(elements.size());

    size_t i = 1;
    for (const T& element : elements) {
        m_pCells[i++] = element;
    }
}


template <class T>
void CVector<T>::SetSize(int nR)
// ---------------------------------------------------------------------------
// Function: dynamically allocates memory
// Input:    vector size
// Output:   none
// ---------------------------------------------------------------------------
{
    // check whether NR is legal
    if (nR <= 0)
        ErrorHandler(CVector<T>::Error::INVALIDSIZE);
    Release();
    try
    {
        m_pCells = std::vector<T>(nR + 1);
    }
    catch (std::bad_alloc)
    {
        ErrorHandler(CVector<T>::Error::MEMORYALLOCATION);
    }
    m_nRows = nR;
#ifdef __ARRAYMETER
    m_dAllocated += static_cast<double>(sizeof(T) * (nR + 1));
#endif
}

template <class T>
CVector<T>::~CVector()
// ---------------------------------------------------------------------------
// Function: dtor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // deallocate storage
    Release();
}

template <class T>
void CVector<T>::Release()
// ---------------------------------------------------------------------------
// Function: dynamically deallocates memory
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
;
}

// =============== member functions ===========================================
template <class T>
void CVector<T>::SetName(const std::string& strName)
// ---------------------------------------------------------------------------
// Function: sets the name of the vector
// Input:    vector name
// Output:   none
// ---------------------------------------------------------------------------
{
    m_strName = strName;
}

template <class T>
void CVector<T>::GetName(std::string& strName) const
// ---------------------------------------------------------------------------
// Function: gets the name of the vector
// Input:    string to hold vector name
// Output:   vector name
// ---------------------------------------------------------------------------
{
    strName = m_strName;
}

template <class T>
void CVector<T>::Set(T dV)
// ---------------------------------------------------------------------------
// Function: sets the value of all the elements in the vector to the
//           specified value
// Input:    specified value
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i = 1; i <= m_nRows; i++)
        m_pCells[i] = dV;
}

template <class T>
int CVector<T>::GetSize() const
// ---------------------------------------------------------------------------
// Function: gets the size of the vector
// Input:    none
// Output:   returns the size
// ---------------------------------------------------------------------------
{
    return m_nRows;
}

// ==================== Overloaded Operators ========================

template <class T>
T& CVector<T>::operator() (int nR)  // T& is reference
// ---------------------------------------------------------------------------
// Function: overloaded () operator to access vector contents
//           carries out bound checking
// Input:    index or location
// Output:   value at the specified index
// ---------------------------------------------------------------------------
{
    // row-column reference in bounds?
    if (nR <= 0 || nR > m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INVALIDINDEX);
        return m_pCells[1];
    }
    else
        return m_pCells[nR];
}

template <class T>
T CVector<T>::operator() (int nR) const
// ---------------------------------------------------------------------------
// Function: overloaded () operator to access vector contents
// Input:    index or location
// Output:   value at the specified index
// ---------------------------------------------------------------------------
{
    // row-column reference in bounds?
    if (nR <= 0 || nR > m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INVALIDINDEX);
        return m_pCells[1];
    }
    else
        return m_pCells[nR];
}


template <class T>
CVector<T> CVector<T>::operator() (int i, int j) {
    if ((j - 1 > m_nRows) || (i < 1)) {
        ErrorHandler(Error::BADSLICE);
    }

    CVector<T> temp(m_pCells, i, j);
    return temp;
}

template <class T>
void CVector<T>::operator= (const CVector& matarg)
// ---------------------------------------------------------------------------
// Function: overloaded = operator 
// Input:    vector to use as rvalue
// Output:   modified values
// ---------------------------------------------------------------------------
{
    // check whether vector is assigned to itself
    if (this != &matarg)
    {
        // compatible vectors?
        if (!is_mem_view) {
            Release();
            SetSize(matarg.GetSize());
        }
        m_strName = matarg.m_strName;
        for (int i = 1; i <= matarg.m_nRows; i++)
            m_pCells[i] = matarg.m_pCells[i];
    }

}


template <class T>
void CVector<T>::operator= (CVector&& vec) noexcept
// ---------------------------------------------------------------------------
// Function: overloaded = operator 
// Input:    vector to use as rvalue
// Output:   modified values
// ---------------------------------------------------------------------------
{
    // compatible vectors?
    if (this == &vec) {
        return;
    }

    if (!is_mem_view) {
        Release();
        m_pCells = std::move(vec.m_pCells);
        m_nRows = vec.m_nRows;
        m_strName = std::move(vec.m_strName);
        is_mem_view = vec.is_mem_view;
        vec.m_nRows = 0;
    }
    else {
        if (vec.GetSize() != m_nRows) {
            ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
            return;
        }
        m_strName = std::move(vec.m_strName);
        for (int i = 1; i <= vec.m_nRows; i++)
            m_pCells[i] = vec.m_pCells[i];
    }

}
template <class T>
void CVector<T>::operator=(std::initializer_list<T> elements) {
    if (is_mem_view && (elements.size() != m_nRows)) {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        return;
    }
    if (!is_mem_view) {
        Release();
        SetSize(elements.size());
    }
    size_t i = 1;
    for (const T& element : elements) {
        m_pCells[i++] = element;
    }
}


template <class T>
bool CVector<T>::operator==(const CVector<T>& vec) const {
    if (m_nRows != vec.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        return false;
    }
    for (int i = 1; i <= vec.GetSize(); i++) {
        if (m_pCells[i] != vec.m_pCells[i]) {
            return false;
        }
    }
    return true;
}


template <class T>
CVector<T> CVector<T>::operator+(const CVector<T>& matarg) const {
    if (m_nRows != matarg.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        // Return a default-constructed CVector with zero elements
        return  CVector<T>(); // Or simply "return {};" in C++11 and later
    }

    CVector<T> result(*this); // Create a copy of the current vector
    for (int i = 1; i <= matarg.m_nRows; i++)
        result.m_pCells[i] += matarg.m_pCells[i];

    return result;
}

template <class T>
CVector<T> CVector<T>::operator+(CVector<T>&& matarg) const noexcept {
    if (m_nRows != matarg.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        // Return a default-constructed CVector with zero elements
        return  CVector<T>(); // Or simply "return {};" in C++11 and later
    }

    CVector<T> result(std::move(matarg)); // Create a copy of the current vector
    for (int i = 1; i <= result.m_nRows; i++)
        result.m_pCells[i] += m_pCells[i];

    return result;
}


template <class T>
CVector<T> CVector<T>::operator+(const T& val) const {
    CVector<T> result(*this);
    for (int i = 1; i <= m_nRows; i++)
        result.m_pCells[i] += val;
    return result;

}

template <class T>
void CVector<T>::operator+=(const CVector<T>& matarg) {
    if (m_nRows != matarg.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        return;
    }
    for (int i = 1; i <= matarg.m_nRows; i++)
        m_pCells[i] += matarg.m_pCells[i];
}

template <class T>
void CVector<T>::operator+=(const T& val) {

    for (int i = 1; i <= m_nRows; i++)
        m_pCells[i] += val;

}

template <class T>
CVector<T> CVector<T>::operator-(const CVector<T>& matarg) const {
    if (m_nRows != matarg.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        return CVector<T>();
    }
    CVector<T> result(*this);
    for (int i = 1; i <= result.m_nRows; i++)
        result.m_pCells[i] -= matarg.m_pCells[i];
    return result;

}

template <class T>
CVector<T> CVector<T>::operator-(CVector<T>&& matarg) const {
    if (m_nRows != matarg.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        return CVector<T>();
    }
    CVector<T> result(matarg);
    for (int i = 1; i <= result.m_nRows; i++)
        result.m_pCells[i] -= m_pCells[i];
    return result;

}

template <class T>
CVector<T> CVector<T>::operator-(const T& val) const {
    CVector result (*this);
    for (int i = 1; i <= m_nRows; i++)
        result.m_pCells[i] -= val;
    return result;

}

template <class T>
void CVector<T>::operator-=(const CVector<T>& matarg) {
    if (m_nRows != matarg.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        return;
    }
    for (int i = 1; i <= matarg.m_nRows; i++)
        m_pCells[i] -= matarg.m_pCells[i];

}

template <class T>
void CVector<T>::operator-=(const T& val) {

    for (int i = 1; i <= m_nRows; i++)
        m_pCells[i] -= val;
}

template <class T>
CVector<T> CVector<T>::operator*(const T& val) const {
    CVector<T> result(*this);
    for (int i = 1; i <= m_nRows; i++)
        result.m_pCells[i] *= val;
    return result;

}

template <class T>
void CVector<T>::operator*=(const T& val) {

    for (int i = 1; i <= m_nRows; i++)
        m_pCells[i] *= val;
}


template <class T>
CVector<T> CVector<T>::operator/(const T& val) const {
    if (val == 0) {
        ErrorHandler(Error::INVALIDOPERATION);
        return CVector<T>();
    }
    CVector<T> result(*this);
    for (int i = 1; i <= m_nRows; i++)
        result.m_pCells[i] /= val;
    return result;

}


template <class T>
void CVector<T>::operator/=(const T& val) {
    if (val == 0) {
        ErrorHandler(Error::INVALIDOPERATION);
        return;
    }
    for (int i = 1; i <= m_nRows; i++)
        m_pCells[i] /= val;
}

template <class T>
CVector<T> CVector<T>::operator*(const CVector<T>& matarg) const {
    if (m_nRows != matarg.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        return CVector<T>();
    }
    CVector<T> result(*this);
    for (int i = 1; i <= matarg.m_nRows; i++)
        result.m_pCells[i] *= matarg.m_pCells[i];
    return result;

}

template <class T>
CVector<T> CVector<T>::operator*(CVector<T>&& matarg) const {
    if (m_nRows != matarg.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        return  CVector<T>();
    }
    CVector<T> result(std::move(matarg));
    for (int i = 1; i <= result.m_nRows; i++)
        result.m_pCells[i] *= m_pCells[i];
    return result;

}

template <class T>
void CVector<T>::operator*=(const CVector<T>& matarg) {
    if (m_nRows != matarg.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        return;
    }
    for (int i = 1; i <= matarg.m_nRows; i++)
        m_pCells[i] *= matarg.m_pCells[i];
}

template <class T>
CVector<T> CVector<T>::operator/(const CVector<T>& matarg) const {
    if (m_nRows != matarg.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        return  CVector<T>();
    }
    CVector<T> result(*this);
    for (int i = 1; i <= matarg.m_nRows; i++) {
        if (matarg.m_pCells[i] == 0) {
            ErrorHandler(Error::INVALIDOPERATION);
        }
        result.m_pCells[i] /= matarg.m_pCells[i];
    }

    return result;

}

template <class T>
CVector<T> CVector<T>::operator/(CVector<T>&& matarg) const noexcept {
    if (m_nRows != matarg.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        return CVector<T>();
    }
    CVector<T> result(std::move(matarg));
    for (int i = 1; i <= result.m_nRows; i++) {
        if (result.m_pCells[i] == 0) {
            ErrorHandler(Error::INVALIDOPERATION);
        }
        result.m_pCells[i] = m_pCells[i]/result.m_pCells[i];
    }

    return result;

}



template <class T>
void CVector<T>::operator/=(const CVector<T>& matarg) {
    if (m_nRows != matarg.m_nRows)
    {
        ErrorHandler(CVector<T>::Error::INCOMPATIBLEVECTORS);
        return;
    }
    for (int i = 1; i <= matarg.m_nRows; i++) {
        if (matarg.m_pCells[i] == 0) {
            ErrorHandler(Error::INVALIDOPERATION);
        }
        m_pCells[i] /= matarg.m_pCells[i];
    }

}

// ==================== Error Handler ========================
template <class T>
void CVector<T>::ErrorHandler(Error err) const
// ---------------------------------------------------------------------------
// Function: channels error message
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    if (err == CVector<T>::Error::INCOMPATIBLEVECTORS)
    {
        throw std::invalid_argument("Vector operation cannot take place due to incompatible vectors.");
    }
    else if (err == CVector<T>::Error::INVALIDINDEX)
    {
        throw std::invalid_argument("Invalid index to access vector element.");
    }
    else if (err == CVector<T>::Error::INVALIDSIZE)
    {
        throw std::invalid_argument("Invalid vector size to create vector.");
    }
    else if (err == CVector<T>::Error::MEMORYALLOCATION)
    {
        throw std::invalid_argument("Vector cannot be created. Memory allocation error.");
    }
    else if (err == CVector<T>::Error::BADSLICE)
    {
        throw std::invalid_argument("Error occured while slicing");
    }
    else
    {
        throw std::invalid_argument("Invalid operation");
    }
}

template <class T>
void CVector<T>::Display(const std::string& strHeader) const {
    // display specified header
    std::cout << strHeader << "\n";
    // show all values one value per line
    if (m_nRows < 1)
    {
        std::cout << "[]\n";
        return;
    }
    std::cout << "[ ";
    for (int i = 1; i <= m_nRows; i++)
    {
        // value-based At function invoked here
        std::cout << m_pCells[i] << " ";
    }
    std::cout << "]" << "\n";
}
