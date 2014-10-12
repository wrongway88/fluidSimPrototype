#ifndef V_PARTICLE_H
#define V_PARTICLE_H

#include "cinder/Vector.h"

namespace VortonSim
{
	class Particle
	{
	public:
		Particle();
		~Particle();

	private:
		ci::Vec3f m_position;
		ci::Vec3f m_velocity;
		ci::Vec3f m_orientation;
		ci::Vec3f m_angularVelocity;

		float m_mass;
		float m_size;
		int m_birthTime;
	};
}

#endif // PARTICLE_H