#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "utility/math/Vector2.h"

class Boundary
{
public:
	Boundary();
	~Boundary();

	Vec2f getPosition() const;
	Vec2f getNormal() const;
	float getWidth() const;

	void setPosition(const Vec2f& position);
	void setNormal(const Vec2f& normal);
	void setWidth(const float width);

	float getDistance(const Vec2f& position) const;
	bool checkHit(const Vec2f& position) const;

private:
	Vec2f m_position;
	Vec2f m_normal;

	float m_width;
};

#endif // BOUNDARY_H