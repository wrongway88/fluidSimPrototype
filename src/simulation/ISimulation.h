#ifndef I_SIMULATION_H
#define I_SIMULATION_H

#include "cinder/app/App.h"

class ISimulation
{
public:
	ISimulation();
	~ISimulation();

	virtual void setup() = 0;
	virtual void update(const float deltaTime) = 0;
	virtual void draw() = 0;

	virtual void mouseDown(const ci::app::MouseEvent& event) = 0;

	virtual void drawParams() = 0;
};

#endif I_SIMULATION_H