#include "Particle.h"

Particle::Particle(const ci::Vec2f& position, const float mass, const float density):
	m_position(position),
	m_velocity(ci::Vec2f::zero()),
	m_force(ci::Vec2f::zero()),
	m_color(ci::Vec3f::zero()),
	m_mass(mass),
	m_density(density),
	m_pressure(0.0f)
{
}

Particle::~Particle()
{
}

void Particle::setPosition(const ci::Vec2f& position)
{
	m_position = position;
}

void Particle::setVelocity(const ci::Vec2f& velocity)
{
	m_velocity = velocity;
}

void Particle::setForce(const ci::Vec2f& force)
{
	m_force = force;
}

void Particle::setColor(const ci::Vec3f& color)
{
	m_color = color;
}

void Particle::setMass(const float mass)
{
	m_mass = mass;
}

void Particle::setDensity(const float density)
{
	m_density = density;
}

void Particle::setPressure(const float pressure)
{
	m_pressure = pressure;
}

ci::Vec2f Particle::getPosition() const
{
	return m_position;
}

ci::Vec2f Particle::getVelocity() const
{
	return m_velocity;
}

ci::Vec2f Particle::getForce() const
{
	return m_force;
}

ci::Vec3f Particle::getColor() const
{
	return m_color;
}

float Particle::getMass() const
{
	return m_mass;
}

float Particle::getDensity() const
{
	return m_density;
}

float Particle::getPressure() const
{
	return m_pressure;
}

bool Particle::operator == (const Particle& other)
{
	if(m_position != other.getPosition())
	{
		return false;
	}
	if(m_velocity != other.getVelocity())
	{
		return false;
	}

	if(m_force != other.getForce())
	{
		return false;
	}

	if(m_color != other.getColor())
	{
		return false;
	}

	if(m_mass != other.getMass())
	{
		return false;
	}

	if(m_density != other.getDensity())
	{
		return false;
	}

	if(m_pressure != other.getPressure())
	{
		return false;
	}

	return true;
}

bool Particle::operator != (const Particle& other)
{
	return !(*this == other);
}
