#ifndef VECTOR_3_H
#define VECTOR_3_H

#include "utility/logging/logging.h"
#include "utility/Property.h"

#include "utility/math/VectorBase.h"

template<class T>
class Vector3 : public VectorBase<T, 3>
{
public:
	Vector3();
	Vector3(const T& x, const T& y, const T& z);
	Vector3(const VectorBase<T, 3>& vector);
	Vector3(const Vector3<T>& vector);
	virtual ~Vector3();

	T getValue(const unsigned int index) const;
	void setValue(const unsigned int index, const T& value);

	T operator[](const unsigned int index);

	Property<T> x;
	Property<T> y;
	Property<T> z;

	Vector3<T> normalize();
	Vector3<T> normalized() const;

	template<class U>
	void operator=(const Vector3<U>& other);

protected:
	static const unsigned int m_xIndex = 0;
	static const unsigned int m_yIndex = 1;
	static const unsigned int m_zIndex = 2;
};

template<class T>
Vector3<T>::Vector3()
	: VectorBase<T, 3>()
	, x(&VectorBase<T, 3>::m_values[m_xIndex])
	, y(&VectorBase<T, 3>::m_values[m_yIndex])
	, z(&VectorBase<T, 3>::m_values[m_zIndex])
{
	setValue(m_xIndex, 0);
	setValue(m_yIndex, 0);
	setValue(m_zIndex, 0);
}

template<class T>
Vector3<T>::Vector3(const T& x, const T& y, const T& z)
	: VectorBase<T, 3>()
	, x(&VectorBase<T, 3>::m_values[m_xIndex])
	, y(&VectorBase<T, 3>::m_values[m_yIndex])
	, z(&VectorBase<T, 3>::m_values[m_zIndex])
{
	setValue(m_xIndex, x);
	setValue(m_yIndex, y);
	setValue(m_zIndex, z);
}

template<class T>
Vector3<T>::Vector3(const VectorBase<T, 3>& vector)
	: VectorBase<T, 3>(vector)
	, x(&VectorBase<T, 3>::m_values[m_xIndex])
	, y(&VectorBase<T, 3>::m_values[m_yIndex])
	, z(&VectorBase<T, 3>::m_values[m_zIndex])
{
}

template<class T>
Vector3<T>::Vector3(const Vector3<T>& vector)
	: VectorBase<T, 3>(vector)
	, x(&VectorBase<T, 3>::m_values[m_xIndex])
	, y(&VectorBase<T, 3>::m_values[m_yIndex])
	, z(&VectorBase<T, 3>::m_values[m_zIndex])
{
}

template<class T>
Vector3<T>::~Vector3()
{
}

template<class T>
T Vector3<T>::getValue(const unsigned int index) const
{
	try
	{
		return VectorBase<T, 3>::getValue(index);
	}
	catch (std::exception& e)
	{
		LOG_ERROR(e.what());
		return 0;
	}
}

template<class T>
void Vector3<T>::setValue(const unsigned int index, const T& value)
{
	try
	{
		VectorBase<T, 3>::setValue(index, value);
	}
	catch (std::exception& e)
	{
		LOG_ERROR(e.what());
	}
}

template<class T>
T Vector3<T>::operator[](const unsigned int index)
{
	try
	{
		return VectorBase<T, 3>::getValue(index);
	}
	catch (std::exception& e)
	{
		LOG_ERROR(e.what());
		return 0;
	}
}

template<class T>
Vector3<T> Vector3<T>::normalize()
{
	return VectorBase<T, 3>::normalize();
}

template<class T>
Vector3<T> Vector3<T>::normalized() const
{
	return VectorBase<T, 3>::normalized();
}

template<class T>
template<class U>
void Vector3<T>::operator=(const Vector3<U>& other)
{
	this->assign(other);
	x = Property<T>(&VectorBase<T, 3>::m_values[m_xIndex]);
	y = Property<T>(&VectorBase<T, 3>::m_values[m_yIndex]);
	z = Property<T>(&VectorBase<T, 3>::m_values[m_zIndex]);
}

typedef Vector3<float> Vec3f;
typedef Vector3<int> Vec3i;

#endif // VECTOR_3_H