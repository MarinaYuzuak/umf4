#ifndef _VECTOR_HPP
#define _VECTOR_HPP

#include <vector>
typedef std::vector<double> V_n;

class Vector : public V_n {
public:

	Vector(int n = 0) : V_n(n) {}
	Vector(const V_n& v) : V_n(v) {}
	Vector(int n, double value) : V_n(n, value) {}
	Vector operator+(const Vector& v) const;
	Vector operator-(const Vector& v) const;
	Vector operator*(double value) const;
	double operator*(const Vector& v) const;
	Vector operator-() const;
	Vector& operator+=(const Vector& v);
	Vector& operator-=(const Vector& v);

	friend Vector operator*(double value, const Vector& v) { return v * value; }
};

Vector Vector::operator+(const Vector& v) const
{
	Vector result(v.size());
	for (int i = 0; i < v.size(); i++)
		result[i] = (*this)[i] + v[i];
	return result;
}

Vector Vector::operator-(const Vector& v) const
{
	Vector result(v.size());
	for (int i = 0; i < v.size(); i++)
		result[i] = (*this)[i] - v[i];
	return result;
}

Vector Vector::operator*(double value) const
{
	Vector result((*this).size());
	for (int i = 0; i < (*this).size(); i++)
		result[i] = value * (*this)[i];
	return result;
}

double Vector::operator*(const Vector& v) const
{
	double result = 0;
	for (int i = 0; i < (*this).size(); i++)
		result += v[i] * (*this)[i];
	return result;
}

Vector Vector::operator-() const
{
	Vector result((*this).size());
	for (int i = 0; i < (*this).size(); i++)
		result[i] = -(*this)[i];
	return result;
}

Vector& Vector::operator+=(const Vector& v)
{
	for (int i = 0; i < (*this).size(); i++)
		(*this)[i] += v[i];
	return *this;
}

Vector& Vector::operator-=(const Vector& v)
{
	for (int i = 0; i < (*this).size(); i++)
		(*this)[i] -= v[i];
	return *this;
}

#endif // !_VECTOR_HPP
