#include "SimulationStableFluids.h"

#include "../utility/math/Vector2.h"

#include "cinder/gl/gl.h"

SimulationStableFluids::SimulationStableFluids()
	: m_gridSize(0, 0)
	, m_containerSize(0.0f, 0.0f)
	, m_velocities()
	, m_previousVelocities()
	, m_densities()
	, m_previousDensities()
	, m_diffusionCoefficient(1.0f)
{
}

SimulationStableFluids::~SimulationStableFluids()
{
}

void SimulationStableFluids::setup(float containerWidth, float containerHeight)
{
	m_containerSize.x = containerWidth;
	m_containerSize.y = containerHeight;

	m_gridSize.x = 10;
	m_gridSize.y = 10;

	setupCells();
}

void SimulationStableFluids::update(const float deltaTime)
{
	velocityStep(deltaTime);
	densityStep(deltaTime);

	resetArray(m_previousDensities);
	resetArray(m_previousVelocities);
}

void SimulationStableFluids::draw()
{
	Vec2f cellSize(0.0f, 0.0f);
	cellSize.x = m_containerSize.x / (m_gridSize.x+2);
	cellSize.y = m_containerSize.y / (m_gridSize.y+2);

	ci::gl::color(ci::ColorA8u(255, 255, 255, 255));

	for(unsigned int i = 0; i < m_gridSize.x+2; i++)
	{
		ci::gl::drawLine(ci::Vec2f(i*cellSize.x, 0.0f), ci::Vec2f(i*cellSize.x, m_containerSize.y));
	}

	for(unsigned int i = 0; i < m_gridSize.y+2; i++)
	{
		ci::gl::drawLine(ci::Vec2f(0.0f, i*cellSize.y), ci::Vec2f(m_containerSize.x, i*cellSize.y));
	}

	float baseRadius = 2.0f;
	for(unsigned int i = 0; i < m_densities.size(); i++)
	{
		for(unsigned int j = 0; j < m_densities[0].size(); j++)
		{
			float density = m_densities[i][j];

			if(density > 0.0f)
			{
				ci::gl::color(ci::ColorA8u(200, 100, 100, 255));
			}
			else
			{
				ci::gl::color(ci::ColorA8u(100, 100, 200, 255));
			}

			float radius = baseRadius + std::abs(density);
			ci::Vec2f pos(i*cellSize.x + cellSize.x*0.5f, j*cellSize.y + cellSize.y*0.5f);
			ci::gl::drawSolidCircle(pos, radius);

			ci::gl::color(ci::ColorA8u(220, 80, 80, 255));
			ci::Vec2f dir(m_velocities[i][j].x, m_velocities[i][j].y);
			//dir.normalize();
			dir *= 100.0f;

			ci::gl::drawLine(pos, pos + dir);

			dir.x = m_previousVelocities[i][j].x;
			dir.y = m_previousVelocities[i][j].y;
			//dir.normalize();
			dir *= 100.0f;

			ci::gl::color(ci::ColorA8u(140, 100, 100, 255));
			ci::gl::drawLine(pos, pos + dir);
		}
	}
}

void SimulationStableFluids::mouseDown(const ci::app::MouseEvent& event)
{
	Vec2f position(0.0f, 0.0f);
	position.x = event.getPos().x;
	position.y = event.getPos().y;

	Vec2f cell(0, 0);
	Vec2f cellSize;
	cellSize.x = m_containerSize.x / (m_gridSize.x+2);
	cellSize.y = m_containerSize.y / (m_gridSize.y+2);
	cell.x = position.x / cellSize.x;
	cell.y = position.y / cellSize.y;
	cell.x = (int)cell.x;
	cell.y = (int)cell.y;

	Vec2f cellCenter;
	cellCenter.x = cellSize.x * cell.x;
	cellCenter.y = cellSize.y * cell.y;
	cellCenter += (cellSize * 0.5f);
	Vec2f velocity = position - cellCenter;

	m_previousVelocities[cell.x][cell.y] = velocity;
}

void SimulationStableFluids::drawParams()
{
}

void SimulationStableFluids::reset()
{
	setupCells();
}

void SimulationStableFluids::setupCells()
{
	std::srand(std::time(0));
	float maxDensity = 5.0f;
	
	m_densities.clear();
	m_previousDensities.clear();
	m_velocities.clear();
	m_previousVelocities.clear();

	// +2 to add rows/colls for border cells
	for(unsigned int i = 0; i < m_gridSize.x+2; i++)
	{
		m_velocities.push_back(std::vector<Vec2f>());
		m_previousVelocities.push_back(std::vector<Vec2f>());

		m_densities.push_back(std::vector<float>());
		m_previousDensities.push_back(std::vector<float>());

		for(unsigned int j = 0; j < m_gridSize.y+2; j++)
		{
			m_velocities[i].push_back(Vec2f(0.0f, 0.0f));
			m_previousVelocities[i].push_back(Vec2f(0.0f, 0.0f));

			float randDensity = ((float)std::rand() / (float)RAND_MAX) * maxDensity;
			randDensity = randDensity - maxDensity*0.5f;
			m_densities[i].push_back(randDensity);
			m_previousDensities[i].push_back(0.0f);
		}
	}
}

void SimulationStableFluids::densityStep(const float deltaTime)
{
	addSource(m_densities, m_previousDensities, deltaTime);

	std::swap(m_densities, m_previousDensities);
	diffuse(deltaTime, m_densities, m_previousDensities);

	std::swap(m_densities, m_previousDensities);
	advect(deltaTime, m_densities, m_previousDensities, m_velocities);
}

void SimulationStableFluids::velocityStep(const float deltaTime)
{
	addSource(m_velocities, m_previousVelocities, deltaTime);

	std::swap(m_velocities, m_previousVelocities);
	diffuse(deltaTime, m_velocities, m_previousVelocities);
	project();

	std::swap(m_velocities, m_previousVelocities);
	advect(deltaTime, m_velocities, m_previousVelocities, m_velocities);
	project();
}

template<class T>
void SimulationStableFluids::advect(const float deltaTime, std::vector<std::vector<T>>& vector, std::vector<std::vector<T>>& previousVector, const std::vector<std::vector<Vec2f>>& velocities)
{
	float dt0x = deltaTime * m_gridSize.x;
	float dt0y = deltaTime * m_gridSize.y;
	Vec2f pos(0.0f, 0.0f);
	for(unsigned int i = 1; i <= m_gridSize.x; i++)
	{
		for(unsigned int j = 1; j <= m_gridSize.y; j++)
		{
			pos.x = (float)i - dt0x * velocities[i][j].x;
			pos.y = (float)j - dt0y * velocities[i][j].y;

			if(pos.x < 0.5f)
			{
				pos.x = 0.5f;
			}
			if(pos.x > m_gridSize.x + 0.5f)
			{
				pos.x = m_gridSize.x + 0.5f;
			}

			int i0 = (int)pos.x;
			int i1 = i0 + 1;

			if(pos.y < 0.5f)
			{
				pos.y = 0.5f;
			}
			if(pos.y > m_gridSize.y + 0.5f)
			{
				pos.y = m_gridSize.y + 0.5f;
			}

			int j0 = (int)pos.y;
			int j1 = j0 + 1;

			float s1 = pos.x - (float)i0;
			float s0 = 1.0f - (float)s1;

			float t1 = pos.y - (float)j0;
			float t0 = 1.0f - (float)t1;

			vector[i][j] = (previousVector[i0][j0] * t0 + previousVector[i0][j1] * t1) * s0
				+ (previousVector[i1][j0] * t0 + previousVector[i1][j1] * t1) * s1;
		}
	}

	setBoundaries(vector);
}

template<class T>
void SimulationStableFluids::diffuse(const float deltaTime, std::vector<std::vector<T>>& vector, std::vector<std::vector<T>>& previousVector)
{
	float a = deltaTime * m_gridSize.x * m_gridSize.y * m_diffusionCoefficient;

	for(unsigned int k = 0; k < 20; k++)
	{
		for(unsigned int i = 1; i <= m_gridSize.x; i++)
		{
			for(unsigned int j = 1; j <= m_gridSize.y; j++)
			{
				vector[i][j] = (previousVector[i][j] + 
					(vector[i-1][j] + vector[i+1][j] + 
					vector[i][j-1] + vector[i][j+1]) * a) / (1.0f + 4.0f * a);
			}
		}
	}

	setBoundaries(vector);
}

template<class T>
void SimulationStableFluids::addSource(std::vector<std::vector<T>>& target, const std::vector<std::vector<T>>& source, const float deltaTime)
{
	if(target.size() <= 0)
	{
		LOG_ERROR_STREAM(<< "Target array is empty");

		return;
	}
	if(target.size() != source.size())
	{
		LOG_ERROR_STREAM(<< "Source has different size than velocities");

		return;
	}
	else if(target[0].size() != source[0].size())
	{
		LOG_ERROR_STREAM(<< "Source[0] has different size than velocities");

		return;
	}

	for(unsigned int i = 0; i < target.size(); i++)
	{
		for(unsigned int j = 0; j < target[0].size(); j++)
		{
			target[i][j] += source[i][j] * deltaTime;
		}
	}
}

void SimulationStableFluids::setBoundaries(std::vector<std::vector<float>>& vector)
{
	for(unsigned int i = 1; i <= m_gridSize.x; i++)
	{
		vector[0][i] = vector[1][i];
		vector[m_gridSize.x+1][i] = vector[m_gridSize.x][i];
		vector[i][0] = vector[i][1];
		vector[i][m_gridSize.y+1] = vector[i][m_gridSize.y];
	}

	vector[0][0] = 0.5f * (vector[1][0] + vector[0][1]);
	vector[0][m_gridSize.y+1] = 0.5f * (vector[1][m_gridSize.y+1] + vector[0][m_gridSize.y]);
	vector[m_gridSize.x+1][0] = 0.5f * (vector[m_gridSize.x][0] + vector[m_gridSize.x+1][1]);
	vector[m_gridSize.x+1][m_gridSize.y+1] = 0.5f * (vector[m_gridSize.x][m_gridSize.y+1] + vector[m_gridSize.x+1][m_gridSize.y]);
}

void SimulationStableFluids::setBoundaries(std::vector<std::vector<Vec2f>>& vector)
{
	for(unsigned int i = 1; i <= m_gridSize.x; i++)
	{
		vector[i][0] = vector[i][1];
		vector[i][0].y = vector[i][0].y * -1.0f;
		vector[i][m_gridSize.y+1] = vector[i][m_gridSize.y];
		vector[i][m_gridSize.y+1].y = vector[i][m_gridSize.y+1].y * -1.0f;
	}
	for(unsigned int i = 1; i <= m_gridSize.y; i++)
	{
		vector[0][i] = vector[1][i];
		vector[0][i].x = vector[0][i].x * -1.0f;
		vector[m_gridSize.x+1][i] = vector[m_gridSize.x][i];
		vector[m_gridSize.x+1][i].x = vector[m_gridSize.x+1][i].x * -1.0f;
	}

	vector[0][0] = (vector[1][0] + vector[0][1]) * 0.5f;
	vector[0][m_gridSize.y+1] = (vector[1][m_gridSize.y+1] + vector[0][m_gridSize.y]) * 0.5f;
	vector[m_gridSize.x+1][0] = (vector[m_gridSize.x][0] + vector[m_gridSize.x+1][1]) * 0.5f;
	vector[m_gridSize.x+1][m_gridSize.y+1] = (vector[m_gridSize.x][m_gridSize.y+1] + vector[m_gridSize.x+1][m_gridSize.y]) * 0.5f;
}

void SimulationStableFluids::project()
{
	float hx = 1.0f / m_gridSize.x;
	float hy = 1.0f / m_gridSize.y;

	std::vector<std::vector<float>> p;
	std::vector<std::vector<float>> div;

	for(unsigned int i = 0; i < m_velocities.size(); i++)
	{
		p.push_back(std::vector<float>());
		div.push_back(std::vector<float>());

		for(unsigned int j = 0; j < m_velocities[0].size(); j++)
		{
			p[i].push_back(0.0f);
			div[i].push_back(0.0f);
		}
	}

	for(unsigned int i = 1; i <= m_gridSize.x; i++)
	{
		for(unsigned int j = 1; j <= m_gridSize.y; j++)
		{
			div[i][j] = -0.5f * ((m_velocities[i+1][j].x - m_velocities[i-1][j].x)*hx + (m_velocities[i][j+1].y - m_velocities[i][j-1].y)*hy);
			p[i][j] = 0.0f;
		}
	}
	setBoundaries(p);
	setBoundaries(div);

	for(unsigned int k = 0; k < 20; k++)
	{
		for(unsigned int i = 1; i <= m_gridSize.x; i++)
		{
			for(unsigned int j = 1; j <= m_gridSize.y; j++)
			{
				p[i][j] = (div[i][j] + p[i-1][j] + p[i+1][j] + p[i][j-1] + p[i][j+1]) / 4.0f;
			}
		}
		setBoundaries(p);
	}

	for(unsigned int i = 1; i <= m_gridSize.x; i++)
	{
		for(unsigned int j = 1; j <= m_gridSize.y; j++)
		{
			m_velocities[i][j].x = m_velocities[i][j].x - 0.5f * (p[i+1][j] - p[i-1][j]) / 0.5f;
			m_velocities[i][j].y = m_velocities[i][j].x - 0.5f * (p[i][j+1] - p[i][j-1]) / 0.5f;
		}
	}
	setBoundaries(m_velocities);
}

void SimulationStableFluids::resetArray(std::vector<std::vector<float>>& array)
{
	for(unsigned int i = 0; i < array.size(); i++)
	{
		for(unsigned int j = 0; j < array[i].size(); j++)
		{
			array[i][j] = 0.0f;
		}
	}
}

void SimulationStableFluids::resetArray(std::vector<std::vector<Vec2f>>& array)
{
	for(unsigned int i = 0; i < array.size(); i++)
	{
		for(unsigned int j = 0; j < array[i].size(); j++)
		{
			array[i][j] = Vec2f(0.0f, 0.0f);
		}
	}
}
