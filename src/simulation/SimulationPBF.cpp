#include "SimulationPBF.h"

#include "particle/Particle.h"

SimulationPBF::SimulationPBF()
	: m_h(90.0f)
{
}

SimulationPBF::~SimulationPBF()
{
}

void SimulationPBF::setup()
{
}

void SimulationPBF::update(const float deltaTime)
{
	for(unsigned int i = 0; i < m_particles.size(); i++)
	{

	}
}

void SimulationPBF::draw()
{
}

void SimulationPBF::mouseDown(const ci::app::MouseEvent& event)
{
}

void SimulationPBF::drawParams()
{
}

std::vector<Particle> SimulationPBF::getNeighbours(Particle& particle, std::vector<Particle>& particles)
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

ci::Vec2f SimulationPBF::getConstraintGradient(const float restDensity, Particle& particle, std::vector<Particle>& neighbours)
{
	ci::Vec2f result(0.0f, 0.0f);

	for(unsigned int i = 0; i < neighbours.size(); i++)
	{
		if(particle != neighbours[i])
		{
			result += distanceKernelGradient((particle.getPosition() - neighbours[i].getPosition()).length(), m_h);
		}
	}

	result *= 1.0f / restDensity;

	return result;
}

float SimulationPBF::distanceKernel(const float distance, const float h)
{
	if(distance > h)
		return 0.0f;

	float a = 315.0f / (3.14159f * std::pow(h, 4.0f));
	float b = (std::pow(distance, 2.0f)) - (std::pow(h, 2.0f));

	return std::pow(a, b);
}

ci::Vec2f distanceKernelGradient(const float distance, const float h)
{
	float exponent = distance*distance - h*h;
	float pi = 3.1415926f;
	float h_4 = std::pow(h, 4.0f);

	// d W/d distance
	float a = 2*distance*std::powf(1.0f/h_4, exponent);
	float b = std::powf(315.0f/pi, exponent);
	float c = std::log(315.0f/(pi * h_4));

	
	// d W/d h
	float d = -2.0f * std::pow((315.0f / pi), exponent);
	float e = std::pow(1.0f/h_4, exponent);
	float f = 2.0f * distance * distance + h*h*(std::log(315.0f / pi) - 2.0f) + h*h * std::log(1.0f/h_4);


	return ci::Vec2f(a*b*c, d*e*f*(1.0f/h));
}