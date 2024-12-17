#include "externalMagneticForce.h"

externalMagneticForce::externalMagneticForce(elasticRod &m_rod, timeStepper &m_stepper, 
	Vector3d m_bAVector, Vector3d m_bRVector, double m_muZero)
{
	rod = &m_rod;
	stepper = &m_stepper;

	baVector_ref = m_bAVector;
	brVector_ref = m_bRVector;
	muZero = m_muZero;

	Id3<<1,0,0,
         0,1,0,
         0,0,1;

    force.setZero(7, 1);
    jacob.setZero(7, 7);

    baVector(0) = 0.0;
    baVector(1) = 0.0;
    baVector(2) = 0.0;

    baVector_old = baVector;

    forceVec = VectorXd::Zero(rod->ndof);
}

externalMagneticForce::~externalMagneticForce()
{
	;
}

void externalMagneticForce::computeFm()
{
	forceVec = VectorXd::Zero(rod->ndof);

	for (int i = 0; i < rod->ne; i++)
	{
		m1_current = rod->m1_old.row(i);
		m2_current = rod->m2_old.row(i);
		m3_current = rod->tangent_old.row(i);

		m1_start = rod->m1_initial.row(i);
		m2_start = rod->m2_initial.row(i);
		m3_start = rod->tangent_initial.row(i);

		x1 = rod->getVertexOld(i); 
		x2 = rod->getVertexOld(i+1); 

		edge = (x2 - x1).norm();

		// numerical test
		/*
		m1_start(0) = 0.8682;
		m1_start(1) = 0.0;
		m1_start(2) = 0.4961;

		m2_start(0) = 0.1195;
		m2_start(1) = 0.9706;
		m2_start(2) = -0.2090;

		m3_start(0) = -0.4815;
		m3_start(1) = 0.2408;
		m3_start(2) = 0.8427;

		m1_current(0) = 0.9355;
		m1_current(1) = 0.1275;
		m1_current(2) = 0.3295;

		m2_current(0) = -0.2939;
		m2_current(1) = 0.7985;
		m2_current(2) = 0.5254;

		m3_current(0) = -0.1961;
		m3_current(1) = -0.5883;
		m3_current(2) = 0.7845;

		edge = 0.5099;

		brVector(0) =  0.3;
		brVector(1) =  0.9;
		brVector(2) = -0.1;

		baVector(0) = -0.1;
		baVector(1) =  0.2;
		baVector(2) =  0.4;

		*/

		gradientBa.setZero(3, 3);

		// we only want the magnetic at the tip
		if (i >= 95)
		{
			brVector = brVector_ref;
		}
		else
		{
			brVector.setZero(3, 1);
		}

		Mag = (m1_start.dot(brVector) * m1_current + m2_start.dot(brVector) * m2_current + m3_start.dot(brVector) * m3_current);

		dm3de = ( Id3 - m3_current * m3_current.transpose() ) / edge;
		dm1de = - ( m3_current * m1_current.transpose() ) / edge;
		dm2de = - ( m3_current * m2_current.transpose() ) / edge;

		dm1dtheta =  m2_current;
		dm2dtheta = -m1_current;
		dm3dtheta.setZero(3,1);

		dMde = m1_start.dot(brVector) * dm1de + m2_start.dot(brVector) * dm2de + m3_start.dot(brVector) * dm3de;
		dMdtheta = m1_start.dot(brVector) * dm1dtheta + m2_start.dot(brVector) * dm2dtheta + m3_start.dot(brVector) * dm3dtheta;

		dEde = dMde.transpose() * baVector;
		dEdtheta = dMdtheta.dot(baVector);

		//cout << i << " " << dEdtheta << endl;

		force.setZero(7, 1);

		force.segment(0, 3) = - dEde + (gradientBa * Mag) / 2;
		force(3) = dEdtheta;
		force.segment(4, 3) = dEde + (gradientBa * Mag) / 2;

		force = - edge * ( force * rod->crossSectionalArea / muZero);

		for (int k = 0; k < 7; k++)
		{
			int ind = 4 * i + k;
			stepper->addForce(ind, -force[k]); // subtracting elastic force

			forceVec(ind) = forceVec(ind) + force[k];
		}

		/*

		tempM3 = dm3de * currentBa;
		tempM1 = dm1de * currentBa;
		tempM2 = dm2de * currentBa;

		tempMTheta = dMdtheta.dot(currentBa);

		d2m1dtheta2 = -m1_current;
		d2m2dtheta2 = -m2_current;
		d2m3dtheta2.setZero(3,1);

		d2Temp3de2 =  - (((dm3de * currentBr) * m3_current.transpose() + m3_current.dot(currentBr) * dm3de) * edge + m3_current * (currentBr - m3_current.dot(currentBr) * m3_current).transpose() ) / ((edge) * (edge));
		d2Temp1de2 = -  ( dm1de.transpose() * currentBr * (edge) - m3_current * m1_current.dot(currentBr) ) * m3_current.transpose() / ((edge) * (edge)) - m1_current.dot(currentBr)/(edge) * dm3de;
		d2Temp2de2 = -  ( dm2de.transpose() * currentBr * (edge) - m3_current * m2_current.dot(currentBr) ) * m3_current.transpose() / ((edge) * (edge)) - m2_current.dot(currentBr)/(edge) * dm3de;

		d2Ede2 = m1_start.dot(brVector) * d2Temp1de2 + m2_start.dot(brVector) * d2Temp2de2 + m3_start.dot(brVector) * d2Temp3de2;

		d2Mdtheta2 = m1_start.dot(brVector) * d2m1dtheta2 + m2_start.dot(brVector) * d2m2dtheta2 + m3_start.dot(brVector) * d2m3dtheta2;

		d2Edtheta2 = d2Mdtheta2.dot(currentBr);

		d2Ededtheta = (m1_start.dot(brVector) * (- (m3_current * m2_current.transpose()) * dm3de ) + m2_start.dot(brVector) * ((m3_current * m1_current.transpose()) * dm3de) ) * currentBr;

		*/

		/*
		dEde.setZero(3,1);
		dEdtheta = 0.0;
		d2Ede2.setZero(3,3);
		d2Ededtheta.setZero(3,1);
		d2Edtheta2 = 0.0;
		*/

		/*

		jacob.setZero(7, 7);

		jacob.block(0,0,3,3) = d2Ede2 - 2 * dMde * gradientBr;
		jacob.block(4,4,3,3) = d2Ede2 + 2 * dMde * gradientBr;
		jacob.block(0,3,3,3) =-d2Ede2 - 2 * dMde * gradientBr;
		jacob.block(3,0,3,3) =-d2Ede2 + 2 * dMde * gradientBr;

		jacob.col(3).segment(0,3) =-d2Ededtheta;
		jacob.col(3).segment(4,3) = d2Ededtheta;

		jacob.row(3).segment(0,3) =-d2Ededtheta;
		jacob.row(3).segment(4,3) = d2Ededtheta;

		jacob(3, 3) = d2Edtheta2;

		jacob = - edge * ( jacob * rod->crossSectionalArea / muZero );

		*/

		/*

		cout << " force " << endl;
		cout << force << endl;
		cout << " jacob " << endl;
		cout << jacob << endl;

		for (int j = 0; j < 7; j++)
		{
			for (int k = 0; k < 7; k++)
			{
				int ind1 = 4 * i + j;
				int ind2 = 4 * i + k;

				//stepper->addJacobian(ind1, ind2, - jacob(j, k));
			}
		}

		*/

	}
}

void externalMagneticForce::computeJm()
{
	;
}