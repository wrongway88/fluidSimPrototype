#include "Boundary.h"

Boundary::Boundary()
	: m_position(0.0f, 0.0f)
	, m_normal(0.0f, -1.0f)
	, m_width(0.0f)
{
}

Boundary::~Boundary()
{
}

Vec2f Boundary::getPosition() const
{
	return m_position;
}

Vec2f Boundary::getNormal() const
{
	return m_normal;
}

float Boundary::getWidth() const
{
	return m_width;
}

void Boundary::setPosition(const Vec2f& position)
{
	m_position = position;
}

void Boundary::setNormal(const Vec2f& normal)
{
	m_normal = normal.normalized();
}

void Boundary::setWidth(const float width)
{
	m_width = width;
}

float Boundary::getDistance(const Vec2f& position) const
{
	Vec2f centerToPos;
	centerToPos = position - m_position;

	return m_normal.dotProduct(centerToPos);
}

bool Boundary::checkHit(const Vec2f& position) const
{
	Vec2f centerToPos;
	centerToPos = position - m_position;

	Vec2f tangent(-m_normal.y, m_normal.x);
	float distance = tangent.dotProduct(centerToPos);

	if(std::abs(distance) < m_width*0.5f)
	{
		return true;
	}
	return false;
}
