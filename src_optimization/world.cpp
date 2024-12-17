#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				// boolean
	saveData = m_inputData.GetBoolOpt("saveData");			// boolean
	
	// Physical parameters
	RodLength = m_inputData.GetScalarOpt("RodLength");      // meter
    helixradius = m_inputData.GetScalarOpt("helixradius");  // meter
    helixpitch = m_inputData.GetScalarOpt("helixpitch");    // meter
    rodRadius = m_inputData.GetScalarOpt("rodRadius");      // meter

    gVector = m_inputData.GetVecOpt("gVector");             // m/s^2
    maxIter = m_inputData.GetIntOpt("maxIter");             // maximum number of iterations
	numVertices = m_inputData.GetIntOpt("numVertices");     // int_num
	youngM = m_inputData.GetScalarOpt("youngM");            // Pa
	Poisson = m_inputData.GetScalarOpt("Poisson");          // dimensionless
	deltaTime = m_inputData.GetScalarOpt("deltaTime");      // seconds
	totalTime= m_inputData.GetScalarOpt("totalTime");       // seconds
	tol = m_inputData.GetScalarOpt("tol");                  // small number like 10e-7
	stol = m_inputData.GetScalarOpt("stol");				// small number, e.g. 0.1%
	density = m_inputData.GetScalarOpt("density");          // kg/m^3
	viscosity = m_inputData.GetScalarOpt("viscosity");      // viscosity in Pa-s

	baVector = m_inputData.GetVecOpt("baVector");           // magnetic field
	brVector = m_inputData.GetVecOpt("brVector");           // magnetic field
	muZero = m_inputData.GetScalarOpt("muZero");            // magnetic field

	dBar = m_inputData.GetScalarOpt("dBar");
	stiffness = m_inputData.GetScalarOpt("stiffness");

	scaleRender = m_inputData.GetScalarOpt("scaleRender");

	speed = m_inputData.GetScalarOpt("speed");

	mu = m_inputData.GetScalarOpt("mu");
	epsilonV = m_inputData.GetScalarOpt("epsilonV");

	stopTime = m_inputData.GetScalarOpt("stopTime");
	caseOpt = m_inputData.GetIntOpt("caseOpt");

	deltaBaStep = m_inputData.GetScalarOpt("deltaBaStep");

	timeWait = m_inputData.GetScalarOpt("timeWait");
	
	shearM = youngM/(2.0*(1.0+Poisson));					// shear modulus
	
	// Viscous drag coefficients using Resistive Force Theory
	eta_per = 4.0*M_PI*viscosity/( log(2*helixpitch/rodRadius) + 0.5);
    eta_par = 2.0*M_PI*viscosity/( log(2*helixpitch/rodRadius) - 0.5);

    outputFlag = 0;

    totalIter = 0;

    moment(0) = 0.0;
    moment(1) = 0.0;
    moment(2) = 0.0;

    totalStep = 0;

    //deltaBaStep = 2.0;
    //timeWait = 2.0;

    gradBa.setZero(3,3);
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}
	
	int systemRet = system("mkdir datafiles"); //make the directory

	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	ostringstream name;
	name.precision(6);
	name << fixed;
    name << "datafiles/simDER";
    name << "_deltaBaStep_" << deltaBaStep;
    name << ".txt";

    outfile.open(name.str().c_str());
    outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	outfile.close();
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	
	//if (timeStep == Nstep)
	{
		/*
		for (int i = 0; i < rod->nv; i++)
		{
			Vector3d xCurrent = rod->getVertex(i);

			outfile << xCurrent(0) << " " << xCurrent(1) << " " << xCurrent(2) << endl;
		}
		*/

		//m_externalContactForce->computeFc();
		//outfile << currentTime << " " << m_externalContactForce->reactionForce.norm() << endl;


		//outfile << m_magneticForce->baVector(0) << " " << m_magneticForce->baVector(1) << " " << m_magneticForce->baVector(2) << " " << m_externalContactForce->reactionForce.norm() << endl;
	}

	if (outputFlag == 1)
	{
		for (int i = 0; i < rod->nv; i++)
		{
			Vector3d xCurrent = rod->getVertex(i);

			//outfile << stopTime - deltaBaStep << " " << xCurrent(0) << " " << xCurrent(1) << " " << xCurrent(2) << endl;
		}
		outfile << stopTime - deltaBaStep << " " << m_magneticForce->baVector(0) << " " << m_magneticForce->baVector(1) << " " << m_magneticForce->baVector(2) << endl;
		outputFlag = 0;
	}

	
	//m_externalContactForce->computeFc();
	//outfile << currentTime << " " << m_externalContactForce->reactionForce.norm() << endl;


}

void world::setRodStepper()
{
	// Set up geometry
	rodGeometry();	

	// Create the rod 
	rod = new elasticRod(vertices, vertices, density, rodRadius, deltaTime,
		youngM, shearM, RodLength, caseOpt);

	// Find out the tolerance, e.g. how small is enough?
	characteristicForce = M_PI * pow(rodRadius ,4)/4.0 * youngM / pow(RodLength, 2);
	forceTol = tol * characteristicForce;
	
	// Set up boundary condition
	rodBoundaryCondition();
	
	// setup the rod so that all the relevant variables are populated
	rod->setup();
	// End of rod setup
	
	// set up the time stepper
	stepper = new timeStepper(*rod);
	totalForce = stepper->getForce();

	// declare the forces
	m_stretchForce = new elasticStretchingForce(*rod, *stepper);
	m_bendingForce = new elasticBendingForce(*rod, *stepper);
	m_twistingForce = new elasticTwistingForce(*rod, *stepper);
	m_inertialForce = new inertialForce(*rod, *stepper);
	m_gravityForce = new externalGravityForce(*rod, *stepper, gVector);
	m_magneticForce = new externalMagneticForce(*rod, *stepper, baVector, brVector, muZero);
	m_dampingForce = new dampingForce(*rod, *stepper, viscosity, eta_per, eta_par);
	m_externalContactForce = new externalContactForce(*rod, *stepper, dBar, stiffness, mu, epsilonV);

	Nstep = totalTime/deltaTime;

	// Allocate every thing to prepare for the first iteration
	rod->updateTimeStep();
	
	timeStep = 0;
	currentTime = 0.0;

	xInitial = rod->getVertex(rod->nv - 1);
}

// Setup geometry
void world::rodGeometry()
{

    vertices = MatrixXd(numVertices, 3);

    double deltaL = RodLength / (numVertices - 1);

    for (int i = 0; i < numVertices; i++)
    {
    	vertices(i, 0) = deltaL * i - RodLength;
        vertices(i, 1) = 0.0;
        vertices(i, 2) = 0.0;
    }

}

void world::rodBoundaryCondition()
{
	// Apply boundary condition
	rod->setVertexBoundaryCondition(rod->getVertex(0), 0);
	rod->setVertexBoundaryCondition(rod->getVertex(1), 1);
	rod->setThetaBoundaryCondition(0.0, 0);
}
	

void world::updateTimeStep()
{
	if (currentTime < stopTime - deltaTime)
	{
		Vector3d x1 = rod->getVertex(0);
		x1(0) = x1(0) + speed * deltaTime;
		rod->setVertexBoundaryCondition(x1, 0);

		Vector3d x2 = rod->getVertex(1);
		x2(0) = x2(0) + speed * deltaTime;
		rod->setVertexBoundaryCondition(x2, 1);
	}

	//m_externalContactForce->computeFc();
	//constraintForce_0 = m_externalContactForce->reactionForce;
	//disOpt = m_externalContactForce->deltaEnd;
	//cout << " deltaEnd " << disOpt << endl;

	if ( currentTime > stopTime + timeWait - 1.0 && currentTime < stopTime + timeWait - 1.0 + 2 * deltaTime)
	{
		//rod->clearMap();

		//rod->setVertexBoundaryCondition(rod->getVertex(0), 0);
		//rod->setVertexBoundaryCondition(rod->getVertex(1), 1);
		//rod->setThetaBoundaryCondition(0.0, 0);

		//rod->setVertexBoundaryCondition(rod->getVertex(rod->nv - 1), rod->nv - 1);

		//rod->reSetupMap();
		//stepper->resetMap();
		//totalForce = stepper->getForce();

		//m_externalContactForce->endNodeRadius = m_externalContactForce->normRadius;

		//flag = 1;

		totalStep = 0;
	}



	
	// solve for Ba
	if (currentTime > stopTime + timeWait)
	{
		/*
		computeReactionForce();
		constraintForce_0(0) = reactionForce(4 * (rod->nv-1) + 0);
		constraintForce_0(1) = reactionForce(4 * (rod->nv-1) + 1);
		constraintForce_0(2) = reactionForce(4 * (rod->nv-1) + 2);

		Vector3d x1 = rod->getVertexOld(rod->nv - 1);
		Vector3d x2 = rod->getVertexOld(rod->nv - 2);
		Vector3d xTTTT = (x1 - x2) / (x1 - x2).norm();
		constraintForce_0 = constraintForce_0 - xTTTT * constraintForce_0.dot(xTTTT);
		*/

		m_externalContactForce->computeFc();
		//disOpt = m_externalContactForce->deltaEnd;
		constraintForce_0 = m_externalContactForce->reactionForce;
		gradBa.setZero(3,3);

		stepper->setZero();
		m_magneticForce->baVector = m_magneticForce->baVector_old;
		m_magneticForce->baVector(0) = m_magneticForce->baVector(0) + 1e-4;
		{
			double normf = forceTol * 10.0;	
			double normf0 = 0;
	
			bool solved = false;
	
			iter = 0;

			// Start with a trial solution for our solution x
			rod->updateGuess(); // x = x0 + u * dt

			while (solved == false)
			{
				rod->prepareForIteration();
		
				stepper->setZero();

				// Compute the forces and the jacobians
				m_inertialForce->computeFi();
				m_inertialForce->computeJi();
			
				m_stretchForce->computeFs();
				m_stretchForce->computeJs();
			
				m_bendingForce->computeFb();
				m_bendingForce->computeJb();
		
				m_twistingForce->computeFt();
				m_twistingForce->computeJt();

				m_magneticForce->computeFm();
				m_magneticForce->computeJm();
		
				m_dampingForce->computeFd();
				m_dampingForce->computeJd();

				m_externalContactForce->computeFc();
		
				// Compute norm of the force equations.
				normf = 0;
				for (int i = 0; i < rod->uncons; i++)
				{
					normf += totalForce[i] * totalForce[i];
				}

				normf = sqrt(normf);
				if (iter == 0) 
				{
					normf0 = normf;
				}
		
				if (normf <= forceTol)
				{
					solved = true;
				}
				else if(iter > 0 && normf <= normf0 * stol)
				{
					solved = true;
				}
		
				if (solved == false)
				{
					stepper->integrator(); // Solve equations of motion
					rod->updateNewtonX(totalForce); // new q = old q + Delta q
					iter++;
				}

				if (iter > maxIter)
				{
					//cout << "Error. Could not converge. Exiting.\n";
					break;
				}
			}
		}
		//computeReactionForce();
		//constraintForce_1(0) = reactionForce(4 * (rod->nv-1) + 0);
		//constraintForce_1(1) = reactionForce(4 * (rod->nv-1) + 1);
		//constraintForce_1(2) = reactionForce(4 * (rod->nv-1) + 2);
		//constraintForce_1 = constraintForce_1 - xTTTT * constraintForce_1.dot(xTTTT);
		m_externalContactForce->computeFc();
		constraintForce_1 = m_externalContactForce->reactionForce;
		//disNew = m_externalContactForce->deltaEnd;
		gradBa.col(0) = (constraintForce_1 - constraintForce_0) / 1e-4;



		m_magneticForce->baVector = m_magneticForce->baVector_old;
		m_magneticForce->baVector(1) = m_magneticForce->baVector(1) + 1e-4;
		{
			double normf = forceTol * 10.0;	
			double normf0 = 0;
	
			bool solved = false;
	
			iter = 0;

			// Start with a trial solution for our solution x
			rod->updateGuess(); // x = x0 + u * dt

			while (solved == false)
			{
				rod->prepareForIteration();
		
				stepper->setZero();

				// Compute the forces and the jacobians
				m_inertialForce->computeFi();
				m_inertialForce->computeJi();
			
				m_stretchForce->computeFs();
				m_stretchForce->computeJs();
			
				m_bendingForce->computeFb();
				m_bendingForce->computeJb();
		
				m_twistingForce->computeFt();
				m_twistingForce->computeJt();

				m_magneticForce->computeFm();
				m_magneticForce->computeJm();
		
				m_dampingForce->computeFd();
				m_dampingForce->computeJd();

				m_externalContactForce->computeFc();
		
				// Compute norm of the force equations.
				normf = 0;
				for (int i = 0; i < rod->uncons; i++)
				{
					normf += totalForce[i] * totalForce[i];
				}

				normf = sqrt(normf);
				if (iter == 0) 
				{
					normf0 = normf;
				}
		
				if (normf <= forceTol)
				{
					solved = true;
				}
				else if(iter > 0 && normf <= normf0 * stol)
				{
					solved = true;
				}
		
				if (solved == false)
				{
					stepper->integrator(); // Solve equations of motion
					rod->updateNewtonX(totalForce); // new q = old q + Delta q
					iter++;
				}

				if (iter > maxIter)
				{
					//cout << "Error. Could not converge. Exiting.\n";
					break;
				}
			}
		}
		//computeReactionForce();
		//constraintForce_1(0) = reactionForce(4 * (rod->nv-1) + 0);
		//constraintForce_1(1) = reactionForce(4 * (rod->nv-1) + 1);
		//constraintForce_1(2) = reactionForce(4 * (rod->nv-1) + 2);
		//constraintForce_1 = constraintForce_1 - xTTTT * constraintForce_1.dot(xTTTT);
		m_externalContactForce->computeFc();
		constraintForce_1 = m_externalContactForce->reactionForce;
		//disNew = m_externalContactForce->deltaEnd;
		gradBa.col(1) = (constraintForce_1 - constraintForce_0) / 1e-4;




		m_magneticForce->baVector = m_magneticForce->baVector_old;
		m_magneticForce->baVector(2) = m_magneticForce->baVector(2) + 1e-4;
		{
			double normf = forceTol * 10.0;	
			double normf0 = 0;
	
			bool solved = false;
	
			iter = 0;

			// Start with a trial solution for our solution x
			rod->updateGuess(); // x = x0 + u * dt

			while (solved == false)
			{
				rod->prepareForIteration();
		
				stepper->setZero();

				// Compute the forces and the jacobians
				m_inertialForce->computeFi();
				m_inertialForce->computeJi();
			
				m_stretchForce->computeFs();
				m_stretchForce->computeJs();
			
				m_bendingForce->computeFb();
				m_bendingForce->computeJb();
		
				m_twistingForce->computeFt();
				m_twistingForce->computeJt();

				m_magneticForce->computeFm();
				m_magneticForce->computeJm();
		
				m_dampingForce->computeFd();
				m_dampingForce->computeJd();

				m_externalContactForce->computeFc();
		
				// Compute norm of the force equations.
				normf = 0;
				for (int i = 0; i < rod->uncons; i++)
				{
					normf += totalForce[i] * totalForce[i];
				}

				normf = sqrt(normf);
				if (iter == 0) 
				{
					normf0 = normf;
				}
		
				if (normf <= forceTol)
				{
					solved = true;
				}
				else if(iter > 0 && normf <= normf0 * stol)
				{
					solved = true;
				}
		
				if (solved == false)
				{
					stepper->integrator(); // Solve equations of motion
					rod->updateNewtonX(totalForce); // new q = old q + Delta q
					iter++;
				}

				if (iter > maxIter)
				{
					//cout << "Error. Could not converge. Exiting.\n";
					break;
				}
			}
		}
		//computeReactionForce();
		//constraintForce_1(0) = reactionForce(4 * (rod->nv-1) + 0);
		//constraintForce_1(1) = reactionForce(4 * (rod->nv-1) + 1);
		//constraintForce_1(2) = reactionForce(4 * (rod->nv-1) + 2);
		//constraintForce_1 = constraintForce_1 - xTTTT * constraintForce_1.dot(xTTTT);
		m_externalContactForce->computeFc();
		constraintForce_1 = m_externalContactForce->reactionForce;
		//disNew = m_externalContactForce->deltaEnd;
		gradBa.col(2) = (constraintForce_1 - constraintForce_0) / 1e-4;
	}

	if (currentTime > stopTime + timeWait && constraintForce_0.norm() > 1e-7)
	{
		Eigen::Vector3d xxxxx = gradBa.colPivHouseholderQr().solve(constraintForce_0);
		
		xxxxx = xxxxx / xxxxx.norm();
		//gradBa = gradBa / gradBa.norm();

		//if (totalStep < 1000)
		{
			moment = 0.9 * moment + epsilonV * xxxxx;
			//moment = epsilonV * gradBa;
			m_magneticForce->baVector = m_magneticForce->baVector_old - moment;
		}
		//if (totalStep > 1000)
		{
			//moment = 0.9 * moment + 0.1 * epsilonV * gradBa;
		}

		//if (rod->u.norm() < 1e-2)
		{
			//gradBa = gradBa / gradBa.norm();
			//moment = epsilonV * disOpt * gradBa - mu * epsilonV * (disOpt - disOpt_old) * gradBa;
			//m_magneticForce->baVector = m_magneticForce->baVector_old - moment;
		}


		m_magneticForce->baVector_old = m_magneticForce->baVector;
		totalStep = totalStep + 1;

		//disOpt_old = disOpt;
	}

	if (currentTime > stopTime + timeWait)
	{
		//m_externalContactForce->computeFc();
		//cout << "mag field: " << m_magneticForce->baVector.transpose() << endl;
		//cout << "force: " << constraintForce_0.norm() << endl;
		//cout << "step " << totalStep << endl;
	}

	int exitSolve = 0;

	if (totalStep > 30000)
	{
		exitSolve = 1;
	}

	if (constraintForce_0.norm() < 1e-7)
	{
		exitSolve = 1;
	}
	
	if (currentTime > stopTime + timeWait && exitSolve == 1)
	{
		//timeStep = Nstep - 1;
		currentTime = stopTime;
		stopTime = stopTime + deltaBaStep;

		//rod->clearMap();

		//rod->setVertexBoundaryCondition(rod->getVertex(0), 0);
		//rod->setVertexBoundaryCondition(rod->getVertex(1), 1);
		//rod->setThetaBoundaryCondition(0.0, 0);

		//rod->reSetupMap();
		//stepper->resetMap();
		//totalForce = stepper->getForce();

		m_externalContactForce->endNodeRadius = 0.1 * m_externalContactForce->normRadius;

		outputFlag = 1;

		moment(0) = 0.0;
		moment(1) = 0.0;
		moment(2) = 0.0;

		totalStep = 0;
	}


	double normf = forceTol * 10.0;	
	double normf0 = 0;
	
	bool solved = false;
	
	iter = 0;

	// Start with a trial solution for our solution x
	rod->updateGuess(); // x = x0 + u * dt
		
	while (solved == false)
	{
		rod->prepareForIteration();
		
		stepper->setZero();

		// Compute the forces and the jacobians
		m_inertialForce->computeFi();
		m_inertialForce->computeJi();
			
		m_stretchForce->computeFs();
		m_stretchForce->computeJs();
			
		m_bendingForce->computeFb();
		m_bendingForce->computeJb();
		
		m_twistingForce->computeFt();
		m_twistingForce->computeJt();

		m_gravityForce->computeFg();
		m_gravityForce->computeJg();

		m_magneticForce->computeFm();
		m_magneticForce->computeJm();
		
		m_dampingForce->computeFd();
		m_dampingForce->computeJd();

		m_externalContactForce->computeFc();
		
		// Compute norm of the force equations.
		normf = 0;
		for (int i = 0; i < rod->uncons; i++)
		{
			normf += totalForce[i] * totalForce[i];
		}

		normf = sqrt(normf);
		if (iter == 0) normf0 = normf;
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}
		
		if (solved == false)
		{
			stepper->integrator(); // Solve equations of motion
			rod->updateNewtonX(totalForce); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			//cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}
	
	rod->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << " iter=" << iter << endl;
	}

	currentTime += deltaTime;
		
	timeStep++;
	
	if (solved == false)
	{
		//timeStep = Nstep; // we are exiting
	}
}

int world::simulationRunning()
{
	//return 1;
	
	if (stopTime <= totalTime)
	{
		return 1;
	}
	else
	{
		return 0;
	}
	
}

int world::numPoints()
{
	return rod->nv;
}

double world::getScaledCoordinate(int i)
{
	return rod->x[i] * scaleRender;
}

Vector3d world::getX(int i)
{
	Vector3d xCurrent1 = rod->getVertex(i);
	Vector3d xCurrent2 = rod->getVertex(i);

	return ((xCurrent1 + xCurrent2) / 2) * scaleRender;
}

Vector3d world::getM1(int i)
{
	Vector3d xCurrent1 = rod->getVertex(i);
	Vector3d xCurrent2 = rod->getVertex(i);

	Vector3d m1 = rod->m1.row(i);

	return ((xCurrent1 + xCurrent2) / 2 + m1 * RodLength/50) * scaleRender;
}

double world::getCurrentTime()
{
	return currentTime;
}

double world::getTotalTime()
{
	return totalTime;
}

double world::getTubeRadius()
{
	return rod->tubeRadius * scaleRender;
}

int world::getTubeNv()
{
	return rod->tubeNv;
}

Vector3d world::getTubeNode(int i)
{
	Vector3d xCurrent = rod->tubeNode.row(i);

	return xCurrent * scaleRender;
}

void world::computeReactionForce()
{
	reactionForce = VectorXd::Zero(rod->ndof);

	m_stretchForce->computeFs();
	m_bendingForce->computeFb();
	m_twistingForce->computeFt();
	m_magneticForce->computeFm();
	m_externalContactForce->computeFc();

	reactionForce = - m_magneticForce->forceVec - m_externalContactForce->forceVec - m_bendingForce->forceVec - m_twistingForce->forceVec - m_stretchForce->forceVec;
}
