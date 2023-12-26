#pragma once
#include <string>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <initializer_list>
#include <vector>

template <class T>
class CVector;


template <class T>
std::ostream& operator<<(std::ostream&, const CVector<T>&);

template <typename T>
class CVector
{
private:
   std::vector<T> data;
   
public:
    enum class Error {
        INVALIDSIZE,         
        MEMORYALLOCATION,     
        INVALIDINDEX,          
        INCOMPATIBLEVECTORS,    
        INVALIDOPERATION,   
        BADSLICE
    };

    CVector() = default;
    CVector(int);
    CVector(int, const T&);
    CVector(std::initializer_list<T>) noexcept;
    CVector(const CVector&);
    CVector(CVector&&) noexcept;
    CVector(const std::vector<T>&);

    CVector& operator=(const CVector&);
    CVector& operator=(CVector&&) noexcept;
    CVector& operator=(std::initializer_list<T>) noexcept;

    T& operator()(int);
    const T& operator()(int) const;

    bool operator==(const CVector<T>&) const;
    bool operator!=(const CVector<T>&) const;

    CVector operator+(const CVector&) const;
    CVector<T> operator+(const T&) const;

    CVector operator-(const CVector&) const;
    CVector operator-(const T&) const;

    CVector operator*(const CVector&) const;
    CVector operator*(const T&) const;

    CVector operator/(const CVector&) const;
    CVector operator/(const T&) const;
    
    CVector& operator+=(const CVector&);
    CVector& operator+=(const T&);

    CVector& operator-=(const CVector&);
    CVector& operator-=(const T&);

    CVector& operator*=(const CVector&);
    CVector& operator*=(const T&);

    CVector& operator/=(const CVector&);
    CVector& operator/=(const T&);

    friend std::ostream& operator<< <T>(std::ostream&, const CVector<T>& vec);
    void ErrorHandler(Error) const;
    int GetSize() const; 
    void SetSize(int); 
    void Set(int); 
    
    void size() const;
};

template <typename T>
void CVector<T>::ErrorHandler(Error err) const
{
    switch (err)
    {
    case Error::INVALIDSIZE:
        throw std::invalid_argument("Invalid size");
        break;
    case Error::INVALIDINDEX:
        throw std::invalid_argument("Invalid index");
        break;
    case Error::INCOMPATIBLEVECTORS:
        throw std::invalid_argument("Incompatible vectors");
        break;
    case Error::INVALIDOPERATION:
        throw std::invalid_argument("Invalid operation");
        break;
    case Error::BADSLICE:
        throw std::invalid_argument("Bad slice");
        break;
    default:
        throw std::invalid_argument("Unknown error");
        break;
    }
}

template <typename T>
CVector<T>::CVector(int size)
{
    if (size < 0)
    {
        ErrorHandler(Error::INVALIDSIZE);
    }
    SetSize(size);
}

template <typename T>
CVector<T>::CVector(int size, const T& value)
{
    if (size < 0)
    {
        ErrorHandler(Error::INVALIDSIZE);
    }
    SetSize(size);
    Set(value);
}

template <typename T>
CVector<T>::CVector(std::initializer_list<T> list) noexcept
{
    SetSize(list.size());
    std::copy(list.begin(), list.end(), data.begin());
}

template <typename T>
CVector<T>::CVector(const CVector& other)
{
    SetSize(other.GetSize());
    std::copy(other.data.begin(), other.data.end(), data.begin());
}

template <typename T>
CVector<T>::CVector(CVector&& other) noexcept
{
    data = std::move(other.data);
}

template <typename T>
CVector<T>::CVector(const std::vector<T>& other)
{
    SetSize(other.size());
    std::copy(other.begin(), other.end(), data.begin());
}

template <typename T>
CVector<T>& CVector<T>::operator=(const CVector& other)
{
    if (this == &other)
        return *this;
    SetSize(other.GetSize());
    std::copy(other.data.begin(), other.data.end(), data.begin());
    return *this;
}

template <typename T>
CVector<T>& CVector<T>::operator=(CVector&& other) noexcept
{
    if (this == &other)
        return *this;
    data = std::move(other.data);
    return *this;
}

template <typename T>
CVector<T>& CVector<T>::operator=(std::initializer_list<T> list) noexcept
{
    SetSize(list.size());
    std::copy(list.begin(), list.end(), data.begin());
    return *this;
}

template <typename T>
T& CVector<T>::operator()(int index)
{
    if (index < 1 || index > data.size())
    {
        ErrorHandler(Error::INVALIDINDEX);
    }
    return data[index-1];
}

template <typename T>
const T& CVector<T>::operator()(int index) const
{
    if (index < 1 || index > data.size())
    {
        ErrorHandler(Error::INVALIDINDEX);
    }
    return data[index-1];
}

template <typename T>
bool CVector<T>::operator==(const CVector& other) const
{
    if (GetSize() != other.GetSize())
        return false;
    std::equal(data.begin(), data.end(), other.data.begin());
    return true;
}

template <typename T>
bool CVector<T>::operator!=(const CVector& other) const
{
    return !(*this == other);
}

template <typename T>
CVector<T> CVector<T>::operator+(const CVector& other) const
{
    if (GetSize() != other.GetSize())
    {
        ErrorHandler(Error::INCOMPATIBLEVECTORS);
    }
    CVector<T> result(GetSize());
    std::transform(data.begin(), data.end(), other.data.begin(), result.data.begin(), std::plus<T>());
    return result;
}

template <typename T>
CVector<T> CVector<T>::operator+(const T& value) const
{
    CVector<T> result(GetSize());
    std::transform(data.begin(), data.end(), result.data.begin(), [value](T x) { return x + value; });
    return result;
}

template <typename T>
CVector<T> CVector<T>::operator-(const CVector& other) const
{
    if (GetSize() != other.GetSize())   
    {
        ErrorHandler(Error::INCOMPATIBLEVECTORS);
    }
    CVector<T> result(GetSize());
    std::transform(data.begin(), data.end(), other.data.begin(), result.data.begin(), std::minus<T>());
    return result;
}


template <typename T>
CVector<T> CVector<T>::operator-(const T& value) const
{
    CVector<T> result(GetSize());
    std::transform(data.begin(), data.end(), result.data.begin(), [value](T x) { return x - value; });
    return result;
}

template <typename T>
CVector<T> CVector<T>::operator*(const CVector& other) const
{
    if (GetSize() != other.GetSize())
    {
        ErrorHandler(Error::INCOMPATIBLEVECTORS);
    }
    CVector<T> result(GetSize());
    std::transform(data.begin(), data.end(), other.data.begin(), result.data.begin(), std::multiplies<T>());
    return result;
}


template <typename T>
CVector<T> CVector<T>::operator*(const T& value) const
{
    CVector<T> result(GetSize());
    std::transform(data.begin(), data.end(), result.data.begin(), [value](T x) { return x * value; });
    return result;
}

template <typename T>
CVector<T> CVector<T>::operator/(const CVector& other) const
{
    if (GetSize() != other.GetSize())
    {
        ErrorHandler(Error::INCOMPATIBLEVECTORS);
    }
    CVector<T> result(GetSize());
    std::transform(data.begin(), data.end(), other.data.begin(), result.data.begin(), [this](const T& lhs, const T& rhs) { if (rhs == 0) { ErrorHandler(Error::INVALIDOPERATION); } return lhs / rhs; }); 
    return result;
}

template <typename T>
CVector<T> CVector<T>::operator/(const T& value) const
{
    CVector<T> result(GetSize());
    std::transform(data.begin(), data.end(), result.data.begin(), [this, value](T x) { if (value == 0) { ErrorHandler(Error::INVALIDOPERATION); } return x / value; }); 
    return result;
}

template <typename T>
CVector<T>& CVector<T>::operator+=(const CVector& other)
{
    if (GetSize() != other.GetSize())
    {
        ErrorHandler(Error::INCOMPATIBLEVECTORS);
    }
    std::transform(data.begin(), data.end(), other.data.begin(), data.begin(), std::plus<T>());
    return *this;
}   

template <typename T>
CVector<T>& CVector<T>::operator+=(const T& value)
{
    std::transform(data.begin(), data.end(), data.begin(), [value](T x) { return x + value; });
    return *this;
}

template <typename T>
CVector<T>& CVector<T>::operator-=(const CVector& other)
{
    if (GetSize() != other.GetSize())
    {
        ErrorHandler(Error::INCOMPATIBLEVECTORS);
    }
    std::transform(data.begin(), data.end(), other.data.begin(), data.begin(), std::minus<T>());
    return *this;
}

template <typename T>
CVector<T>& CVector<T>::operator-=(const T& value)
{
    std::transform(data.begin(), data.end(), data.begin(), [value](T x) { return x - value; });
    return *this;
}

template <typename T>
CVector<T>& CVector<T>::operator*=(const CVector& other)
{
    if (GetSize() != other.GetSize())
    {
        ErrorHandler(Error::INCOMPATIBLEVECTORS);
    }
    std::transform(data.begin(), data.end(), other.data.begin(), data.begin(), std::multiplies<T>());
    return *this;
}

template <typename T>
CVector<T>& CVector<T>::operator*=(const T& value)
{
    std::transform(data.begin(), data.end(), data.begin(), [value](T x) { return x * value; });
    return *this;
}

template <typename T>
CVector<T>& CVector<T>::operator/=(const CVector& other)
{
    if (GetSize() != other.GetSize())
    {
        ErrorHandler(Error::INCOMPATIBLEVECTORS);
    }
    std::transform(data.begin(), data.end(), other.data.begin(), data.begin(), [this](const T& lhs, const T& rhs) { if (rhs == 0) { ErrorHandler(Error::INVALIDOPERATION); } return lhs / rhs; });
    return *this;
}

template <typename T>
CVector<T>& CVector<T>::operator/=(const T& value)
{
    std::transform(data.begin(), data.end(), data.begin(), [this, value](T x) { if (value == 0) { ErrorHandler(Error::INVALIDOPERATION); } return x / value; });    
    return *this;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const CVector<T>& vec)
{
    os << "[ ";
    for (int i = 0; i < vec.data.size(); i++)
    {
        os << vec.data[i] << " ";
    }
    os << "]";
    return os;
}

template <typename T>
int CVector<T>::GetSize() const
{
    return data.size();
}

template <typename T>
void CVector<T>::SetSize(int size)
{
    if (size < 0)
    {
        ErrorHandler(Error::INVALIDSIZE);
    }
    data.resize(size);
}

template <typename T>
void CVector<T>::Set(int val)
{
    std::fill(data.begin(), data.end(), val);
}

template <typename T>
void CVector<T>::size() const
{
    return data.size();
}