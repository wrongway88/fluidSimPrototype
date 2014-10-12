#ifndef SIMULATION_PBF_H
#define SIMULATION_PBF_H

#include <vector>

#include "cinder/app/App.h"

class Particle;

class SimulationPBF
{
public:
	SimulationPBF();
	~SimulationPBF();

	virtual void setup();
	virtual void update(const float deltaTime);
	virtual void draw();

	virtual void mouseDown(const ci::app::MouseEvent& event);

	virtual void drawParams();

private:
	std::vector<Particle> getNeighbours(Particle& particle, std::vector<Particle>& particles);
	ci::Vec2f getConstraintGradient(const float restDensity, Particle& particle, std::vector<Particle>& neighbours);

	float distanceKernel(const float distance, const float h);
	ci::Vec2f distanceKernelGradient(const float distance, const float h);

	std::vector<Particle> m_particles;

	float m_h;
};

#endif SIMULATION_PBF_H