#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <map>

#include "utility/math/Vector2.h"
#include "utility/math/Vector3.h"

class Particle
{
public:
	Particle(const Vec2f& position, const float mass, const float density);
	~Particle();

	void setPosition(const Vec2f& position);
	void setVelocity(const Vec2f& velocity);
	void setForce(const Vec2f& force);

	void setColor(const Vec3f& color);

	void setMass(const float mass);
	void setDensity(const float density);
	void setPressure(const float pressure);
	void setExtraFloat(const std::string& key, const float value);
	void setExtraVector(const std::string& key, const Vec2f& value);

	void addPathPoint(const Vec2f& point);
	void clearPathPoints();

	Vec2f getPosition() const;
	Vec2f getVelocity() const;
	Vec2f getForce() const;

	Vec3f getColor() const;

	float getMass() const;
	float getDensity() const;
	float getPressure() const;
	float getExtraFloat(const std::string& key) const;
	Vec2f getExtraVector(const std::string& key) const;

	std::vector<Vec2f> getPathPoints() const;

	bool operator == (const Particle& other);
	bool operator != (const Particle& other);

private:
	Vec2f m_position;
	Vec2f m_velocity;
	Vec2f m_force;

	Vec3f m_color;

	float m_mass;
	float m_density;
	float m_pressure;

	std::vector<Vec2f> m_pathPoints;
	std::map<std::string, float> m_extraFloats;
	std::map<std::string, Vec2f> m_extraVectors;
};

#endif // PARTICLE_H