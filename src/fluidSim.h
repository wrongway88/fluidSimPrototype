#ifndef FLUID_SIM_H
#define FLUID_SIM_H

#include <list>

#include "cinder/Vector.h"
#include "cinder/params/Params.h"
#include "cinder/app/AppBasic.h"

#include "Simulation/Simulation.h"
#include "Simulation/SimpleSPH.h"
#include "Simulation/IncompressibleSPH.h"
#include "Simulation/SimulationStableFluids.h"
#include "simulation/IISPH.h"

class FluidSim : public ci::app::AppBasic
{
public:
	void prepareSettings(ci::app::AppBasic::Settings* settings);
	void setup();
	void mouseDown(ci::app::MouseEvent event);
	void mouseDrag(ci::app::MouseEvent event);
	void keyDown(ci::app::KeyEvent event);
	void update();
	void draw();

	void togglePauseSim();
	void toggleFixedDeltaTime();
	void resetSimulation();

private:
	IISPH m_simulation;

	float m_lastTime;
	float m_fixedDeltaTime;
	float m_fps;

	bool m_pauseSim;
	bool m_useFixedDeltaTime;

	ci::params::InterfaceGlRef m_params;
};

#endif // FLUID_SIM_H