#include "UtilityCinder.h"

ci::Vec2f UtilityCinder::vectorUtiliteaToCinder(const Vec2f& vector)
{
	return ci::Vec2f(vector.x, vector.y);
}

Vec2f UtilityCinder::vectorCinderToUtilitea(const ci::Vec2f& vector)
{
	return Vec2f(vector.x, vector.y);
}

ci::Vec2i UtilityCinder::vectorUtiliteaToCinder(const Vec2i& vector)
{
	return ci::Vec2i(vector.x, vector.y);
}

Vec2i UtilityCinder::vectorCinderToUtilitea(const ci::Vec2i& vector)
{
	return Vec2i(vector.x, vector.y);
}

ci::Vec3f UtilityCinder::vectorUtiliteaToCinder(const Vec3f& vector)
{
	return ci::Vec3f(vector.x, vector.y, vector.z);
}

Vec3f UtilityCinder::vectorCinderToUtilitea(const ci::Vec3f& vector)
{
	return Vec3f(vector.x, vector.y, vector.z);
}

ci::Vec3i UtilityCinder::vectorUtiliteaToCinder(const Vec3i& vector)
{
	return ci::Vec3i(vector.x, vector.y, vector.z);
}

Vec3i UtilityCinder::vectorCinderToUtilitea(const ci::Vec3i& vector)
{
	return Vec3i(vector.x, vector.y, vector.z);
}
