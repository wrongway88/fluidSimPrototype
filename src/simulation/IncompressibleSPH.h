#ifndef INCOMPRESSIBLE_SPH_H
#define INCOMPRESSIBLE_SPH_H

#include <vector>

#include "ISimulation.h"

#include "particle/Particle.h"
#include "particle/Boundary.h"

#include "cinder/params/Params.h"

class IncompressibleSPH : public ISimulation
{
public:
	IncompressibleSPH();
	~IncompressibleSPH();

	virtual void setup(float containerWidth = 1920.0f, float containerHeight = 1080.0f);
	virtual void update(const float deltaTime);
	virtual void draw();

	virtual void mouseDown(const ci::app::MouseEvent& event);

	virtual void drawParams();

	virtual void reset();

private:
	void initParticles();
	void initBoundaries();

	void updateParticleDensity();

	void calculatePressureForces(const float deltaTime);
	void calculateViscosity(const float deltaTime);
	void calculateBoundaryForces(const float deltaTime);

	float distanceKernel(const float distance);
	Vec2f distanceKernelGradient(const Vec2f& pos0, const Vec2f& pos1);
	float distanceKernelViscosity(const float distance);

	float m_containerWidth;
	float m_containerHeight;

	int m_verticleParticleCount;
	int m_horizontalParticleCount;

	float m_h;
	float m_visosityFactor;
	float m_scale;
	float m_particleMass;

	float m_restDensity;

	std::vector<Particle> m_particles;
	std::vector<Boundary> m_boundaries;

	ci::params::InterfaceGlRef m_params;

	Vec2f m_bodyForce;
	ci::Vec3f m_bodyForceParamsProxy; // Params can't handle ci::Vec2f ...
};

#endif // INCOMPRESSIBLE_SPH_H