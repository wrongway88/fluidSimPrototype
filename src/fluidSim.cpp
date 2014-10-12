#include "cinder/app/AppBasic.h"
#include <list>

#include "Simulation/Simulation.h"

#include "cinder/Vector.h"

using namespace ci;
using namespace ci::app;
using namespace std;

// We'll create a new Cinder Application by deriving from the AppBasic class
class FluidSim : public AppBasic
{
public:
	void prepareSettings(Settings* settings);
	void setup();
	void mouseDown(MouseEvent event);
	void mouseDrag(MouseEvent event);
	void keyDown(KeyEvent event);
	void update();
	void draw();

private:
	// This will maintain a list of points which we will draw line segments between
	list<Vec2f> m_points;

	Simulation m_simulation;

	float m_lastTime;
};

void FluidSim::prepareSettings(Settings* settings)
{
	FILE* f;
	AllocConsole();
	freopen_s(&f,"CON","w",stdout);

	settings->setWindowSize(1920.0f, 1080.0f);
}

void FluidSim::setup()
{
	m_lastTime = 0.0f;

	m_simulation.setup();
}

void FluidSim::mouseDown(MouseEvent event)
{
	m_simulation.mouseDown(event);
}

void FluidSim::mouseDrag(MouseEvent event)
{
	//m_points.push_back(event.getPos());
}

void FluidSim::keyDown(KeyEvent event)
{
	if( event.getChar() == 'f' )
	{
		setFullScreen(!isFullScreen());
	}
	if(event.getCode() == KeyEvent::KEY_ESCAPE)
	{
		quit();
	}
}

void FluidSim::update()
{
	float currentTime = ci::app::getElapsedSeconds();
	float deltaTime = currentTime - m_lastTime;
	m_lastTime = currentTime;

	m_simulation.update(deltaTime);
}

void FluidSim::draw()
{
	gl::clear(Color(0.4f, 0.4f, 0.7f));

	gl::color(1.0f, 0.5f, 0.25f);
	gl::begin(GL_LINE_STRIP);

	for(auto pointIter = m_points.begin(); pointIter != m_points.end(); ++pointIter)
	{
		gl::vertex( *pointIter );
	}

	gl::end();

	m_simulation.draw();
	m_simulation.drawParams();
}

// This line tells Cinder to actually create the application
CINDER_APP_BASIC( FluidSim, RendererGl )