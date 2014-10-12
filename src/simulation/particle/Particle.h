#ifndef PARTICLE_H
#define PARTICLE_H

#include "cinder/Vector.h"

class Particle
{
public:
	Particle(const ci::Vec2f& position, const float mass, const float density);
	~Particle();

	void setPosition(const ci::Vec2f& position);
	void setVelocity(const ci::Vec2f& velocity);
	void setForce(const ci::Vec2f& force);

	void setColor(const ci::Vec3f& color);

	void setMass(const float mass);
	void setDensity(const float density);
	void setPressure(const float pressure);

	ci::Vec2f getPosition() const;
	ci::Vec2f getVelocity() const;
	ci::Vec2f getForce() const;

	ci::Vec3f getColor() const;

	float getMass() const;
	float getDensity() const;
	float getPressure() const;

	bool operator == (const Particle& other);
	bool operator != (const Particle& other);

private:
	ci::Vec2f m_position;
	ci::Vec2f m_velocity;
	ci::Vec2f m_force;

	ci::Vec3f m_color;

	float m_mass;
	float m_density;
	float m_pressure;
};

#endif // PARTICLE_H