#include "Particle.h"

Particle::Particle(const Vec2f& position, const float mass, const float density)
	: m_position(position)
	, m_velocity(Vec2f(0.0f, 0.0f))
	, m_force(Vec2f(0.0f, 0.0f))
	, m_color(Vec3f(0.0f, 0.0f, 0.0f))
	, m_mass(mass)
	, m_density(density)
	, m_pressure(0.0f)
	, m_pathPoints()
	, m_extraFloats()
	, m_extraVectors()
{
}

Particle::~Particle()
{
}

void Particle::setPosition(const Vec2f& position)
{
	m_position = position;
}

void Particle::setVelocity(const Vec2f& velocity)
{
	m_velocity = velocity;
}

void Particle::setForce(const Vec2f& force)
{
	m_force = force;
}

void Particle::setColor(const Vec3f& color)
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

void Particle::setExtraFloat(const std::string& key, const float value)
{
	m_extraFloats[key] = value;
}

void Particle::setExtraVector(const std::string& key, const Vec2f& value)
{
	m_extraVectors[key] = value;
}

void Particle::addPathPoint(const Vec2f& point)
{
	if(m_pathPoints.size() > 10)
	{
		m_pathPoints.erase(m_pathPoints.begin());
	}

	m_pathPoints.push_back(point);
}

void Particle::clearPathPoints()
{
	m_pathPoints.clear();
}

Vec2f Particle::getPosition() const
{
	return m_position;
}

Vec2f Particle::getVelocity() const
{
	return m_velocity;
}

Vec2f Particle::getForce() const
{
	return m_force;
}

Vec3f Particle::getColor() const
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

float Particle::getExtraFloat(const std::string& key) const
{
	return m_extraFloats.at(key);
}

Vec2f Particle::getExtraVector(const std::string& key) const
{
	return m_extraVectors.at(key);
}

std::vector<Vec2f> Particle::getPathPoints() const
{
	return m_pathPoints;
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
