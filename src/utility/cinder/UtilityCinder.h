#ifndef UTILITY_CINDER_H
#define UTILITY_CINDER_H

#include "cinder/Vector.h"

#include "utility/math/Vector2.h"
#include "utility/math/Vector3.h"

class UtilityCinder
{
public:
	static ci::Vec2f vectorUtiliteaToCinder(const Vec2f& vector);
	static Vec2f vectorCinderToUtilitea(const ci::Vec2f& vector);

	static ci::Vec2i vectorUtiliteaToCinder(const Vec2i& vector);
	static Vec2i vectorCinderToUtilitea(const ci::Vec2i& vector);

	static ci::Vec3f vectorUtiliteaToCinder(const Vec3f& vector);
	static Vec3f vectorCinderToUtilitea(const ci::Vec3f& vector);

	static ci::Vec3i vectorUtiliteaToCinder(const Vec3i& vector);
	static Vec3i vectorCinderToUtilitea(const ci::Vec3i& vector);
};

#endif // UTILITY_CINDER_H