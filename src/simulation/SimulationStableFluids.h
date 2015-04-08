#ifndef SIMULATION_STABLE_FLUIDS_H
#define SIMULATION_STABLE_FLUIDS_H

#include "utility/math/Vector2.h"

#include "ISimulation.h"

/**
 * @note: After the paper 'Real-Time Fluid Dynamics for Games'
 */
class SimulationStableFluids : public ISimulation
{
public:
	SimulationStableFluids();
	~SimulationStableFluids();

	virtual void setup(float containerWidth = 1920.0f, float containerHeight = 1080.0f);
	virtual void update(const float deltaTime);
	virtual void draw();

	virtual void mouseDown(const ci::app::MouseEvent& event);

	virtual void drawParams();

	virtual void reset();

private:
	void setupCells();

	void densityStep(const float deltaTime);

	void velocityStep(const float deltaTime);

	template<class T>
	void advect(const float deltaTime, std::vector<std::vector<T>>& vector, std::vector<std::vector<T>>& previousVector, const std::vector<std::vector<Vec2f>>& velocities);
	template<class T>
	void diffuse(const float deltaTime, std::vector<std::vector<T>>& vector, std::vector<std::vector<T>>& previousVector);
	template<class T>
	void addSource(std::vector<std::vector<T>>& target, const std::vector<std::vector<T>>& source, const float deltaTime);

	void setBoundaries(std::vector<std::vector<float>>& vector);
	void setBoundaries(std::vector<std::vector<Vec2f>>& vector);

	void project();

	void resetArray(std::vector<std::vector<float>>& array);
	void resetArray(std::vector<std::vector<Vec2f>>& array);

	Vec2i m_gridSize;
	Vec2f m_containerSize;

	std::vector<std::vector<Vec2f>> m_velocities;
	std::vector<std::vector<Vec2f>> m_previousVelocities;

	std::vector<std::vector<float>> m_densities;
	std::vector<std::vector<float>> m_previousDensities;

	float m_diffusionCoefficient;
};

#endif // SIMULATION_STABLE_FLUIDS_H