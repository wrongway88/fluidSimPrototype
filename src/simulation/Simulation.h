#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>

#include "cinder/params/Params.h"
#include "cinder/app/App.h"

#include "particle/Particle.h"

#include "ISimulation.h"

class Simulation : public ISimulation
{
public:
	Simulation();
	~Simulation();

	virtual void setup();
	virtual void update(const float deltaTime);
	virtual void draw();

	virtual void mouseDown(const ci::app::MouseEvent& event);

	virtual void drawParams();

private:
	void updateDensity(const unsigned int startIndex, const unsigned int endIndex);
	void updateForces(const unsigned int startIndex, const unsigned int endIndex);
	void reposition(const float deltaTime, const unsigned int startIndex, const unsigned int endIndex);

	std::vector<Particle> getNeighbours(Particle& particle, std::vector<Particle>& particles);
	float distanceKernel(const float distance);
	float distanceKernelSpiky(const float distance);

	void reset();

	std::vector<Particle> m_particles;

	ci::params::InterfaceGlRef m_params;

	float m_densityCoefficient;
	float m_h;
	float m_restDensity;
	float m_pressureStiffness;
	float m_pressureGradientCoefficient;
	float m_viscosity;
	float m_viscosityCoefficient;
	float m_wallStiffness;

	float m_particleMass;

	float m_maxDensity;

	ci::Vec2f m_gravity;

	std::vector<ci::Vec3f> m_walls;
};

#endif // SIMULATION_H