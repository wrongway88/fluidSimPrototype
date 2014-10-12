#include "Simulation.h"

#include <cmath>
#include <thread>

#include "cinder/app/AppBasic.h"
#include "cinder/gl/gl.h"

Simulation::Simulation():
	m_particles()
{

}

Simulation::~Simulation()
{

}

void Simulation::setup()
{
	m_densityCoefficient = 0.00001f;
	m_h = 90.0f;
	m_restDensity = 1000.0f;
	m_pressureStiffness = 100.0f;
	m_pressureGradientCoefficient = -0.0001f;
	m_viscosity = 0.1f;
	m_viscosityCoefficient = 0.1f;
	m_wallStiffness = 100.0f;

	m_particleMass = 0.1f;

	m_gravity = ci::Vec2f(0.0f, 9.81f);
	m_gravity *= 30.0f;

	m_params = ci::params::InterfaceGl::create(ci::app::getWindow(), "App parameters", ci::app::toPixels( ci::Vec2i( 300, 400)));
	//m_params->addParam("density coef.", &m_densityCoefficient, "min=0.0f max=1000.0f");
	m_params->addParam("h", &m_h);
	m_params->addParam("rest density", &m_restDensity, "min=0.0f max=1000.0f");
	m_params->addParam("pressure stiffness", &m_pressureStiffness, "min=0.0f max=1000.0f");
	//m_params->addParam("pressure grad. coef.", &m_pressureGradientCoefficient, "min=-1000.0f max=1000.0f");
	m_params->addParam("viscosity", &m_viscosity, "min=0.0f max=1000.0f");
	m_params->addParam("wall stiffness", &m_wallStiffness, "min=0.0f max=1000.0f");
	m_params->addParam("particle mass", &m_particleMass, "min=0.0f max=1000.0f");
	//m_params->addParam("gravity", &m_gravity, "");
	m_params->addButton("reset sim", std::bind(&Simulation::reset, this));

	// normal = [x,y], offset in -normal to [0,0] = [z]
	m_walls.push_back(ci::Vec3f(0.0f, -1.0f, 1000.0f));
	//m_walls.push_back(ci::Vec3f(0.0f, 1.0f, 0.0f));
	m_walls.push_back(ci::Vec3f(1.0f, 0.0f, 0.0f));
	m_walls.push_back(ci::Vec3f(-1.0f, 0.0f, 1920.0f));

	reset();
}

void Simulation::update(const float deltaTime)
{
	unsigned int threadCount = 8;
	unsigned int threadParticleCount = m_particles.size() / threadCount;

	// update density
	std::vector<std::thread> dThreads;
	for(unsigned int i = 0; i < threadCount; i++)
	{
		unsigned int startIdx = i * threadParticleCount;
		unsigned int endIdx = (i+1) * threadParticleCount;

		if(i == threadCount-1)
		{
			endIdx = m_particles.size();
		}

		dThreads.push_back(std::thread(std::bind(&Simulation::updateDensity, this, startIdx, endIdx)));
	}
	for(unsigned int i = 0; i < threadCount; i++)
	{
		dThreads[i].join();
	}

	// update forces
	std::vector<std::thread> fThreads;
	for(unsigned int i = 0; i < threadCount; i++)
	{
		unsigned int startIdx = i * threadParticleCount;
		unsigned int endIdx = (i+1) * threadParticleCount;

		if(i == threadCount-1)
		{
			endIdx = m_particles.size();
		}

		fThreads.push_back(std::thread(std::bind(&Simulation::updateForces, this, startIdx, endIdx)));
	}
	for(unsigned int i = 0; i < threadCount; i++)
	{
		fThreads[i].join();
	}

	// reposition
	std::vector<std::thread> rThreads;
	for(unsigned int i = 0; i < threadCount; i++)
	{
		unsigned int startIdx = i * threadParticleCount;
		unsigned int endIdx = (i+1) * threadParticleCount;

		if(i == threadCount-1)
		{
			endIdx = m_particles.size();
		}

		rThreads.push_back(std::thread(std::bind(&Simulation::reposition, this, deltaTime, startIdx, endIdx)));
	}
	for(unsigned int i = 0; i < threadCount; i++)
	{
		rThreads[i].join();
	}
}

void Simulation::draw()
{
	for(unsigned int i = 0; i < m_particles.size(); i++)
	{
		float density = m_particles[i].getDensity();
		if(density > m_maxDensity)
		{
			m_maxDensity = density * 0.8f;
		}

		float densityFactor = density / m_maxDensity;

		ci::Vec3f color = m_particles[i].getColor();
		color *= std::max(0.2f, densityFactor);
		ci::gl::color(color.x, color.y, color.z);
		ci::gl::drawSolidCircle(m_particles[i].getPosition(), m_h * 0.2f);
	}
}

void Simulation::mouseDown(const ci::app::MouseEvent& event)
{
	Particle p(event.getPos(), m_particleMass * 0.1f, 1.0f);
	p.setColor(ci::Vec3f(0.3f, 0.8f, 0.0f));
	m_particles.push_back(p);
}

void Simulation::drawParams()
{
	m_params->draw();
}

void Simulation::updateDensity(const unsigned int startIndex, const unsigned int endIndex)
{
	float h_squared = m_h * m_h;

	for(unsigned int i = startIndex; i < endIndex; i++)
	{
		ci::Vec2f position = m_particles[i].getPosition();
		std::vector<Particle> neighbours = getNeighbours(m_particles[i], m_particles);
		float density = 0.0f;

		for(unsigned int j = 0; j < neighbours.size(); j++)
		{
			ci::Vec2f neighbourPosition = neighbours[j].getPosition();
			ci::Vec2f difference = neighbourPosition - position;
			float r_squared = difference.dot(difference);

			if(r_squared < h_squared)
			{
				density += m_densityCoefficient * (h_squared - r_squared) * (h_squared - r_squared) * (h_squared - r_squared);
			}
		}

		m_particles[i].setDensity(density);
	}
}

void Simulation::updateForces(const unsigned int startIndex, const unsigned int endIndex)
{
	float h_squared = m_h * m_h;

	for(unsigned int i = startIndex; i < endIndex; i++)
	{
		if(m_particles[i].getDensity() > 0.0f)
		{
			ci::Vec2f position = m_particles[i].getPosition();
			ci::Vec2f velocity = m_particles[i].getVelocity();
			float density = m_particles[i].getDensity();
			float bar = density / m_restDensity;
			float foo = bar * bar * bar; //std::pow(bar, 3.0f);
			float pressure = m_pressureStiffness * std::max(foo - 1.0f, 0.0f);

			std::vector<Particle> neighbours = getNeighbours(m_particles[i], m_particles);

			ci::Vec2f acceleration(0.0f, 0.0f);

			for(unsigned int j = 0; j < neighbours.size(); j++)
			{
				ci::Vec2f neighbourPosition = neighbours[j].getPosition();
				ci::Vec2f difference = neighbourPosition - position;
				float r_squared = difference.dot(difference);

				if(r_squared < h_squared)
				{
					ci::Vec2f neighbourVelocity = neighbours[j].getVelocity();
					float neighbourDensity = neighbours[j].getDensity();
					float neighbourPressure = m_pressureStiffness * std::max(std::pow(neighbourDensity / m_restDensity, 3.0f) - 1.0f, 0.0f);

					float r = std::sqrt(r_squared);

					float averagePressure = 0.5f * (neighbourPressure + pressure);
					float pressureGradient = m_pressureGradientCoefficient * averagePressure / neighbourDensity * (m_h-r) * (m_h-r) / r;

					acceleration += difference * pressureGradient;

					ci::Vec2f velocityDifference = neighbourVelocity - velocity;
					float viscosity = m_viscosityCoefficient / neighbourDensity * (m_h - r);

					acceleration += velocityDifference * viscosity;
				}
			}

			m_particles[i].setForce(acceleration / density);
		}
	}
}

void Simulation::reposition(const float deltaTime, const unsigned int startIndex, const unsigned int endIndex)
{
	for(unsigned int i = startIndex; i < endIndex; i++)
	{
		ci::Vec2f position = m_particles[i].getPosition();
		ci::Vec2f velocity = m_particles[i].getVelocity();
		ci::Vec2f acceleration = m_particles[i].getForce();

		// walls
		for(unsigned int j = 0; j < m_walls.size(); j++)
		{
			ci::Vec2f normal(m_walls[j].x, m_walls[j].y);
			float offset = m_walls[j].z;
			float distance = position.dot(normal) + offset;
			acceleration += std::min(distance, 0.0f) * -m_wallStiffness * normal;
		}

		acceleration += m_gravity;
		velocity += acceleration * deltaTime;

		// friction
		velocity *= 0.99f;

		position += velocity * deltaTime;

		m_particles[i].setVelocity(velocity);
		m_particles[i].setPosition(position);
	}
}

std::vector<Particle> Simulation::getNeighbours(Particle& particle, std::vector<Particle>& particles)
{
	std::vector<Particle> result;

	for(unsigned int i = 0; i < particles.size(); i++)
	{
		if(particles[i] != particle)
		{
			if(distanceKernel(particle.getPosition().distance(particles[i].getPosition())))
			{
				result.push_back(particles[i]);
			}
		}
	}

	return result;
}

float Simulation::distanceKernel(const float distance)
{
	if(distance > m_h)
		return 0.0f;

	//float a = 1.0f / (std::pow(3.14159f, 3.0f/2.0f) * std::pow(m_h, 3.0f));
	//float b = (std::pow(distance, 2.0f)) / (std::pow(m_h, 2.0f));

	float a = 315.0f / (3.14159f * std::pow(m_h, 4.0f));
	float b = (std::pow(distance, 2.0f)) - (std::pow(m_h, 2.0f));

	return std::pow(a, b);
}

float Simulation::distanceKernelSpiky(const float distance)
{
	if(distance > m_h)
		return 0.0f;

	float a = (distance * -1.0f) * ( 45.0f / (3.14159f * std::pow(m_h, 6.0f) * distance) ) * std::pow((m_h - distance), 2.0f);

	return a;
}

void Simulation::reset()
{
	m_particles.clear();

	float w = 1920.0f * 0.5f;
	float h = 1080.0f * 0.5f;

	unsigned int hp = 20;
	unsigned int vp = 20;

	int hOffset = 0;

	for(unsigned int i = 0; i < hp; i++)
	{
		for(unsigned int j = 0; j < vp; j++)
		{
			Particle p(ci::Vec2f(i * (w/hp) + hOffset, 1080.0f - (j * (h/vp))), m_particleMass, 1.0f);
			p.setColor(ci::Vec3f(0.0f, 0.3f, 0.9f));
			m_particles.push_back(p);

			//hOffset += 50;
		}

		hOffset = 0;
	}

	m_densityCoefficient = m_particleMass * 315.0f / (64.0f * 3.14159f * std::pow(m_h, 2.0f));
	m_pressureGradientCoefficient = m_particleMass * -45.0f / (3.14159f * pow(m_h, 1.3f));
	m_viscosityCoefficient = m_particleMass * m_viscosity * 45.0f / (3.14159f * pow(m_h, 1.3f));
	m_maxDensity = 0.0f;
}
