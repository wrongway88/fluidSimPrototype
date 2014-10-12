#include "particle.h"

namespace VortonSim
{
	Particle::Particle()
		: m_position(0.0f, 0.0f, 0.0f)
		, m_velocity(0.0f, 0.0f, 0.0f)
		, m_orientation(0.0f, 0.0f, 0.0f)
		, m_angularVelocity(0.0f, 0.0f, 0.0f)
		, m_mass(0.0f)
		, m_size(0.0f)
		, m_birthTime(0.0f)
	{
	}

	Particle::~Particle()
	{
	}
}