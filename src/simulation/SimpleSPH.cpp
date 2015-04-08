#include "SimpleSPH.h"

#include "cinder/gl/gl.h"

#include "utility/cinder/UtilityCinder.h"
#include "utility/logging/logging.h"

#include "utility/Defines.h"

SimpleSPH::SimpleSPH()
	: m_containerWidth(0.0f)
	, m_containerHeight(0.0f)
	, m_verticleParticleCount(5)
	, m_horizontalParticleCount(5)
	, m_particles()
	, m_h(4.0f)
	, m_visosityFactor(1.0f)
	, m_scale(1.0f / 20.0f)
	, m_restDensity(80.0f)
	, m_particleMass(676.0f) // assuming a sphere the particle would have ~0.52 m�, air has ~1300 gram per m�, therefore each particle gets a mass of 676 // i forgot why I assumed 0.52 m�, well damn...
	, m_bodyForce(0.0f, 0.05f)
	, m_bodyForceParamsProxy(0.0f, 0.05f, 0.0f)
{
}

SimpleSPH::~SimpleSPH()
{
}

void SimpleSPH::setup(float containerWidth, float containerHeight)
{
	m_containerWidth = containerWidth;
	m_containerHeight = containerHeight;

	m_params = ci::params::InterfaceGl::create(ci::app::getWindow(), "SimpleSPH", ci::app::toPixels( ci::Vec2i( 300, 400)));
	m_params->addParam("Rest Density", &m_restDensity);
	m_params->addParam("H", &m_h);
	m_params->addParam("Mass", &m_particleMass);
	m_params->addParam("Visosity", &m_visosityFactor);
	m_params->addParam("Vertical particles", &m_verticleParticleCount);
	m_params->addParam("Horizontal particles", &m_horizontalParticleCount);
	m_params->addParam("Body Force", &m_bodyForceParamsProxy);

	initParticles();
	initBoundaries();
}

void SimpleSPH::update(const float deltaTime)
{
	//LOG_WARNING_STREAM(<< "update");

	m_bodyForce.x = m_bodyForceParamsProxy.x;
	m_bodyForce.y = m_bodyForceParamsProxy.y;

	updateParticleDensity();
	calculatePressureForces(deltaTime);
	calculateViscosity(deltaTime);
	calculateBoundaryForces(deltaTime);

	for(unsigned int i = 0; i < m_particles.size(); i++)
	{
		m_particles[i].setVelocity(m_particles[i].getVelocity() + m_particles[i].getForce() + m_bodyForce * deltaTime);
		m_particles[i].setPosition(m_particles[i].getPosition() + m_particles[i].getVelocity() * deltaTime);

		if(m_particles[i].getPathPoints().size() > 0)
		{
			Vec2f lastPathPoint = m_particles[i].getPathPoints()[m_particles[i].getPathPoints().size()-1];
		
			if((m_particles[i].getPosition() - lastPathPoint).getLengthSquared() > m_h * m_h * 0.01f)
			{
				m_particles[i].addPathPoint(m_particles[i].getPosition());
			}
		}
		else
		{
			m_particles[i].addPathPoint(m_particles[i].getPosition());
		}

		// LOG_INFO_STREAM(<< "position: " << m_particles[i].getPosition());
	}
}

void SimpleSPH::draw()
{
	for(unsigned int i = 0; i < m_particles.size(); i++)
	{
		ci::Vec2f pos = UtilityCinder::vectorUtiliteaToCinder(m_particles[i].getPosition());
		
		/*ci::gl::color(0.6f, 0.3f, 0.05f, 0.5f);
		ci::gl::drawSolidCircle(pos * 20.0f, m_h / m_scale);
		*/
		float pressure = m_particles[i].getPressure();

		if(pressure > 0.0f)
		{
			ci::gl::color(0.6f, 0.1f, 0.01f, 0.5f);
		}
		else
		{
			ci::gl::color(0.01f, 0.1f, 0.6f, 0.5f);
		}

		pressure = std::abs(pressure);
		ci::gl::drawSolidCircle(pos / m_scale, std::max(1.0f, std::sqrt((pressure * 5.0f))));
	}

	ci::gl::color(0.0f, 0.0f, 0.0f, 1.0f);
	for(unsigned int i = 0; i < m_particles.size(); i++)
	{
		ci::gl::begin(GL_LINE_STRIP);
		for(unsigned int j = 0; j < m_particles[i].getPathPoints().size(); j++)
		{
			ci::gl::vertex(UtilityCinder::vectorUtiliteaToCinder(m_particles[i].getPathPoints()[j]) / m_scale);
		}
		ci::gl::end();
	}

	
	for(unsigned int i = 0; i < m_particles.size(); i++)
	{
		ci::Vec2f pos = UtilityCinder::vectorUtiliteaToCinder(m_particles[i].getPosition());

		ci::gl::color(0.0f, 0.0f, 0.0f, 1.0f);
		ci::gl::drawSolidCircle(pos / m_scale, 1.0f);

		ci::gl::color(0.0f, 0.0f, 0.0f, 0.3f);
		ci::gl::drawStrokedCircle(pos / m_scale, m_h / m_scale);
	}

	ci::gl::color(0.0f, 0.7f, 0.0f, 1.0f);
	for(unsigned int i = 0; i < m_particles.size(); i++)
	{
		ci::Vec2f pos = UtilityCinder::vectorUtiliteaToCinder(m_particles[i].getPosition());
		ci::Vec2f force = UtilityCinder::vectorUtiliteaToCinder(m_particles[i].getForce());

		ci::gl::drawLine(pos  / m_scale, (pos + (force / force.length())) / m_scale);
	}

	for(unsigned int i = 0; i < m_boundaries.size(); i++)
	{
		ci::Vec2f lineStart;
		ci::Vec2f lineEnd;

		ci::Vec2f normalStart;
		ci::Vec2f normalEnd;

		Vec2f lineDirection(-m_boundaries[i].getNormal().y, m_boundaries[i].getNormal().x);
		Vec2f tmpLineStart = m_boundaries[i].getPosition() + lineDirection * m_boundaries[i].getWidth() * 0.5f;
		Vec2f tmpLineEnd = m_boundaries[i].getPosition() + lineDirection * m_boundaries[i].getWidth() * -0.5f;

		Vec2f tmpLineCenter = m_boundaries[i].getPosition();
		Vec2f tmpNormalEnd = tmpLineCenter + (m_boundaries[i].getNormal());

		normalStart.set(tmpLineCenter.x, tmpLineCenter.y);
		normalEnd.set(tmpNormalEnd.x, tmpNormalEnd.y);

		lineStart.set(tmpLineStart.x, tmpLineStart.y);
		lineEnd.set(tmpLineEnd.x, tmpLineEnd.y);

		ci::gl::color(0.7f, 0.3f, 0.0f, 1.0f);
		ci::gl::drawLine(lineStart / m_scale, lineEnd / m_scale);

		ci::gl::color(0.8f, 0.1f, 0.0f, 1.0f);
		ci::gl::drawLine(normalStart / m_scale, normalEnd / m_scale);
	}
}

void SimpleSPH::mouseDown(const ci::app::MouseEvent& event)
{
}

void SimpleSPH::drawParams()
{
	m_params->draw();
}

void SimpleSPH::reset()
{
	initParticles();
}

void SimpleSPH::initParticles()
{
	int horizontalParticleCount = m_horizontalParticleCount;
	int verticalParticleCount = m_verticleParticleCount;

	m_particles.clear();

	float xOffset = m_h * 0.2f;
	float yOffset = m_h * 0.2f;

	float xPos = (m_containerWidth * 0.5f) - (horizontalParticleCount * 0.5f * (xOffset / m_scale));
	float yPos = (m_containerHeight * 0.5f) - (verticalParticleCount * 0.5f * (yOffset / m_scale));

	xPos *= m_scale;
	yPos *= m_scale;

	float initYPos = yPos;

	for(unsigned int x = 0; x < horizontalParticleCount; x++)
	{
		xPos += xOffset;
		yPos = initYPos;

		for(unsigned int y = 0; y < verticalParticleCount; y++)
		{
			yPos += yOffset;

			Particle newParticle(Vec2f(xPos, yPos), m_particleMass, 0.0f);

			m_particles.push_back(newParticle);
		}
	}
}

void SimpleSPH::initBoundaries()
{
	float heightOffset = -50.0f;

	Boundary boundary0;
	boundary0.setPosition(Vec2f(m_containerWidth/2.0f + 242.0f, m_containerHeight + heightOffset) * m_scale);
	boundary0.setWidth(200.0f * m_scale);

	Boundary boundary1;
	boundary1.setPosition(Vec2f(m_containerWidth/2.0f + 342.0f, m_containerHeight + heightOffset - 100.0f) * m_scale);
	boundary1.setNormal(Vec2f(-1.0f, 0.0f));
	boundary1.setWidth(200.0f * m_scale);

	Boundary boundary2;
	boundary2.setPosition(Vec2f(m_containerWidth/2.0f - 342.0f, m_containerHeight + heightOffset - 100.0f) * m_scale);
	boundary2.setNormal(Vec2f(1.0f, 0.0f));
	boundary2.setWidth(200.0f * m_scale);

	Boundary boundary3;
	boundary3.setPosition(Vec2f(m_containerWidth/2.0f - 71.0f, m_containerHeight + heightOffset - 71.0f) * m_scale);
	boundary3.setNormal(Vec2f(-1.0f, -1.0f));
	boundary3.setWidth(200.0f * m_scale);

	Boundary boundary4;
	boundary4.setPosition(Vec2f(m_containerWidth/2.0f + 71.0f, m_containerHeight + heightOffset - 71.0f) * m_scale);
	boundary4.setNormal(Vec2f(1.0f, -1.0f));
	boundary4.setWidth(200.0f * m_scale);

	Boundary boundary5;
	boundary5.setPosition(Vec2f(m_containerWidth/2.0f - 242.0f, m_containerHeight + heightOffset) * m_scale);
	boundary5.setWidth(200.0f * m_scale);

	m_boundaries.push_back(boundary0);
	m_boundaries.push_back(boundary1);
	m_boundaries.push_back(boundary2);

	m_boundaries.push_back(boundary3);
	m_boundaries.push_back(boundary4);
	m_boundaries.push_back(boundary5);
}

void SimpleSPH::updateParticleDensity()
{
	for(unsigned int i = 0; i < m_particles.size(); i++)
	{
		Vec2f pos = m_particles[i].getPosition();

		float density = 0.0f;

		for(unsigned int j = 0; j < m_particles.size(); j++)
		{
			if(i == j)
			{
				continue;
			}

			Vec2f neighbourPos = m_particles[j].getPosition();

			float distance = (pos - neighbourPos).getLength();
			float w = distanceKernel(distance);
			if(w > 0.0f)
			{
				density += m_particles[j].getMass() * w;
			}
		}

		m_particles[i].setDensity(density);
		//m_particles[i].setPressure(std::abs((density - m_restDensity)));
		//m_particles[i].setPressure((density / m_restDensity) - 1.0f);
		m_particles[i].setPressure(density - m_restDensity);

		// LOG_INFO_STREAM(<< "particle " << i << ": density - " << density << " | pressure - " << density - m_restDensity);
	}
}

void SimpleSPH::calculatePressureForces(const float deltaTime)
{
	for(unsigned int i = 0; i < m_particles.size(); i++)
	{
		Vec2f pos = m_particles[i].getPosition();

		Vec2f pressureGradient(0.0f, 0.0f);

		for(unsigned int j = 0; j < m_particles.size(); j++)
		{
			if(i == j)
			{
				continue;
			}

			Vec2f neighbourPos = m_particles[j].getPosition();

			float mass = m_particles[j].getMass();

			float factor = mass * ((m_particles[i].getPressure() / std::pow(m_particles[i].getDensity(), 2.0f)) + (m_particles[j].getPressure() / std::pow(m_particles[j].getDensity(), 2.0f)));
			Vec2f gradient = distanceKernelGradient(pos, neighbourPos) * factor;

			pressureGradient -= gradient;
		}

		// LOG_INFO_STREAM(<< pressureGradient);

		m_particles[i].setForce(pressureGradient * deltaTime);
	}
}

void SimpleSPH::calculateViscosity(const float deltaTime)
{
	for(unsigned int i = 0; i < m_particles.size(); i++)
	{
		Vec2f pos = m_particles[i].getPosition();
		Vec2f vel = m_particles[i].getVelocity();

		Vec2f viscosityForce(0.0f, 0.0f);

		for(unsigned int j = 0; j < m_particles.size(); j++)
		{
			if(i == j)
			{
				continue;
			}

			Vec2f neighbourPos = m_particles[j].getPosition();
			Vec2f neighbourVel = m_particles[j].getVelocity();

			float w = distanceKernelViscosity((pos - neighbourPos).getLength());
			Vec2f dir = (neighbourVel - vel) / m_particles[j].getDensity();

			viscosityForce += dir * m_particles[j].getMass() * w;
		}

		viscosityForce *= m_visosityFactor / m_particles[i].getDensity();

		m_particles[i].setForce((m_particles[i].getForce() + viscosityForce) * deltaTime);
	}
}

void SimpleSPH::calculateBoundaryForces(const float deltaTime)
{
	for(unsigned int i = 0; i < m_particles.size(); i++)
	{
		for(unsigned int j = 0; j < m_boundaries.size(); j++)
		{
			if(m_boundaries[j].checkHit(m_particles[i].getPosition()))
			{
				float distance = m_boundaries[j].getDistance(m_particles[i].getPosition());
				if(distance > m_h)
				{
					continue;
				}

				float borderParticleMass = 999.9f;
			
				float forceMagnitude = borderParticleMass / (borderParticleMass + m_particles[i].getMass());
				if(forceMagnitude != 0.0f)
				{
					distance = 1.0f / distance;

					distance = distance * distance * distance;
				}
				else
				{
					distance = 99999.9f;
				}

				distance = std::min(distance, m_particles[i].getMass() * m_particles[i].getVelocity().getLength() * deltaTime);

				Vec2f force = m_boundaries[j].getNormal() * distance * forceMagnitude;

				m_particles[i].setForce(m_particles[i].getForce() + (force * deltaTime));
			}
		}
	}
}

float SimpleSPH::distanceKernel(const float distance)
{
	if(distance > m_h)
	{
		return 0.0f;
	}

	float a = 315.0f / (64.0f * PI * std::pow(m_h, 9.0f));
	float b = std::pow((std::pow(m_h, 2.0f)) - (std::pow(distance, 2.0f)), 3.0f);

	return a*b;
}

Vec2f SimpleSPH::distanceKernelGradient(const Vec2f& pos0, const Vec2f& pos1)
{
	float distance = (pos0 - pos1).getLength();

	if(distance > m_h)
	{
		return Vec2f(0.0f, 0.0f);
	}

	Vec2f gradient(0.0f, 0.0f);

	float a = -45.0f / (PI * std::pow(m_h, 6.0f));
	float b = std::pow(m_h - distance, 2.0f);
	Vec2f direction = (pos0 - pos1) / distance;

	return direction * a * b;
}

float SimpleSPH::distanceKernelViscosity(const float distance)
{
	if(distance > m_h)
	{
		return 0.0f;
	}

	float a = 45.0f / (PI * std::pow(m_h, 6.0f));
	float b = std::pow(m_h - distance, 2.0f);

	return a * b;
}
