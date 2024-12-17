#ifndef EXTERNALCONTACTFORCE_H
#define EXTERNALCONTACTFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class externalContactForce
{
public:
	externalContactForce(elasticRod &m_rod, timeStepper &m_stepper, 
        double m_dBar, double m_stiffness, double m_mu, double m_epsilonV);
	~externalContactForce();

	void computeFc();
	void computeJc();
    void computeFcact();

    Vector3d reactionForce;

    VectorXd forceVec;

    double endNodeRadius;
    double normRadius;

    double deltaEnd;
    double deltaEnd2;

    int endIndex;
    
private:
	elasticRod *rod;
	timeStepper *stepper;

    double radius;
    double dbar;
    double stiffness;

    double mu;
    double epsilonV;

    Vector3d x_0;
    Vector3d x_1;
    Vector3d x_2;

    Vector3d computeContactForce(double x0, double y0, double z0, double x1, double y1, double z1, 
  double x2, double y2, double z2, double radius, double dbar);
    Matrix3d computeContactJacobian(double x0, double y0, double z0, double x1, double y1, double z1, 
  double x2, double y2, double z2, double radius, double dbar);
    Vector3d ListVec(double a1, double a2, double a3);
    Matrix3d ListMat(Vector3d a1, Vector3d a2, Vector3d a3);
};

#endif
