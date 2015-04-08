#include "fluidSim.h"

#include "utility/logging/logging.h"
#include "utility/logging/LogManagerImplementation.h"
#include "utility/logging/ConsoleLogger.h"
#include "utility/logging/FileLogger.h"

void FluidSim::prepareSettings(Settings* settings)
{
	FILE* f;
	AllocConsole();
	freopen_s(&f,"CON","w",stdout);

	settings->setWindowSize(1200.0f, 700.0f);
}

void FluidSim::setup()
{
	m_params = ci::params::InterfaceGl::create(ci::app::getWindow(), "Global", ci::app::toPixels( ci::Vec2i( 300, 400)));
	m_params->addParam("Fixed delta time", &m_fixedDeltaTime);
	m_params->addButton("Pause", std::bind(&FluidSim::togglePauseSim, this));
	m_params->addButton("Toggle fixed dT", std::bind(&FluidSim::toggleFixedDeltaTime, this));
	m_params->addButton("Reset", std::bind(&FluidSim::resetSimulation, this));

	std::shared_ptr<ConsoleLogger> consoleLogger = std::make_shared<ConsoleLogger>();
	consoleLogger->setLogLevel(Logger::LOG_INFOS | Logger::LOG_WARNINGS | Logger::LOG_ERRORS);
	//LogManager::getInstance()->addLogger(consoleLogger);
	//LogManager::getInstance()->addLogger(std::make_shared<FileLogger>());

	m_lastTime = 0.0f;
	m_fixedDeltaTime = 0.1f;

	m_simulation.setup(1200.0f, 700.0f);

	ci::gl::enableAlphaBlending(true);

	m_pauseSim = true;
	m_useFixedDeltaTime = true;
}

void FluidSim::mouseDown(ci::app::MouseEvent event)
{
	m_simulation.mouseDown(event);
}

void FluidSim::mouseDrag(ci::app::MouseEvent event)
{}

void FluidSim::keyDown(ci::app::KeyEvent event)
{
	if(event.getChar() == 'p')
	{
		togglePauseSim();
	}
	if(event.getCode() == ci::app::KeyEvent::KEY_ESCAPE)
	{
		quit();
	}
}

void FluidSim::update()
{
	float currentTime = ci::app::getElapsedSeconds();
	float deltaTime = currentTime - m_lastTime;
	m_lastTime = currentTime;
	m_fps = 1.0f / deltaTime;

	if(m_useFixedDeltaTime)
	{
		deltaTime = m_fixedDeltaTime;
	}

	if(m_pauseSim == false)
	{
		m_simulation.update(deltaTime);
	}
}

void FluidSim::draw()
{
	ci::gl::clear(ci::Color(0.2f, 0.2f, 0.2f));

	m_simulation.draw();
	m_params->draw();
	m_simulation.drawParams();

	std::stringstream fpsMessage;
	fpsMessage << "FPS: " << m_fps;
	ci::gl::drawString(fpsMessage.str(), ci::Vec2f(10.0f, 10.0f), ci::Color(0.0f, 0.8f, 0.2f));

	if(m_pauseSim)
	{
		ci::gl::drawString("Pause", ci::Vec2f(10.0f, 25.0f), ci::Color(1.0f, 0.0f, 0.0f));
	}

	std::string deltaTimeMessage = "";
	if(m_useFixedDeltaTime)
	{
		deltaTimeMessage = "Fixed delta time";
	}
	else
	{
		deltaTimeMessage = "Variable delta time";
	}
	ci::gl::drawString(deltaTimeMessage, ci::Vec2f(10.0f, 40.0f), ci::Color(0.0f, 0.8f, 0.2f));
}

void FluidSim::togglePauseSim()
{
	m_pauseSim = !m_pauseSim;
}

void FluidSim::resetSimulation()
{
	m_simulation.reset();
}

void FluidSim::toggleFixedDeltaTime()
{
	m_useFixedDeltaTime = !m_useFixedDeltaTime;
}

// This line tells Cinder to actually create the application
CINDER_APP_BASIC( FluidSim, ci::app::RendererGl )
