#include "externalContactForce.h"

externalContactForce::externalContactForce(elasticRod &m_rod, timeStepper &m_stepper, 
        double m_dBar, double m_stiffness, double m_mu, double m_epsilonV)
{
	rod = &m_rod;
	stepper = &m_stepper;

	dbar = m_dBar;
	stiffness = m_stiffness; 
	radius = rod->tubeRadius;

  mu = m_mu;
  epsilonV = m_epsilonV;

  normRadius = radius;
  endNodeRadius = radius;

  reactionForce = VectorXd::Zero(6);

  //cout << radius << endl;
}

externalContactForce::~externalContactForce()
{
	;
}

void externalContactForce::computeFc()
{

  /*
  for (int i = 0; i < rod->nv; i++)
  {
    Vector3d xCurrent = rod->getVertex(i);

    if (xCurrent(0) < 0.0)
    {
      double yPosHeight = xCurrent(1);
      double forceC = - 10 * stiffness * yPosHeight;
      double jacoC = 10 * stiffness;

      int ind = 4 * i + 1;
      stepper->addForce(ind, - forceC);
      stepper->addJacobian(ind, ind, jacoC);

    }
  }
  */

  reactionForce = VectorXd::Zero(6);

	for (int i = 0; i < rod->nv; i++)
	{
		Vector3d xCurrent = rod->getVertex(i);

		double minDistance = 100000.00;
		int contactIndex = 0;

		for (int kk = 0; kk < rod->tubeNv - 1; kk++)
		{
			Vector3d tubeNode1 = rod->tubeNode.row(kk);
			Vector3d tubeNode2 = rod->tubeNode.row(kk+1);

			double d1 = (xCurrent - tubeNode1).norm();
			double d2 = (xCurrent - tubeNode2).norm();

			double xCurrentDis = (d2 + d1) / 2;

			if (xCurrentDis < minDistance)
			{
				minDistance = xCurrentDis;
				contactIndex = kk;
			}
		}

		x_0 = xCurrent;
		x_1 = rod->tubeNode.row(contactIndex);
		x_2 = rod->tubeNode.row(contactIndex+1);


		Vector3d e1 = x_0 - x_1;
		Vector3d e2 = x_0 - x_2;
		Vector3d e3 = x_1 - x_2;

    //double angleTheta1 = e1.dot(e1) + e3.dot(e3) - e2.dot(e2);
    //double angleTheta2 = e2.dot(e2) + e3.dot(e3) - e1.dot(e1);

		Vector3d e1e2 = e1.cross(e2);

    if (contactIndex < 190)
    {
      radius = normRadius * 0.1;
    }
    else
    {
      radius = normRadius;
    }

    //cout << radius << endl;

		double dis = radius - e1e2.norm() / e3.norm();

		//cout << dis << endl;

		if (dis < dbar)
		{

			Vector3d force;

      force = computeContactForce(x_0(0), x_0(1), x_0(2), x_1(0), x_1(1), x_1(2), 
          x_2(0), x_2(1), x_2(2), radius, dbar);

			force = - force * stiffness;

			for (int ii = 0; ii < 3; ii++)
			{
				int ind = 4 * i + ii;
				stepper->addForce(ind, - force(ii));
			}

      if (i == rod->nv-2)
      {
        //reactionForce(0) = force(0);
        //reactionForce(1) = force(1);
        //reactionForce(2) = force(2);
      }
      if (i == rod->nv-1)
      {
        reactionForce(0) = force(0);
        reactionForce(1) = force(1);
        reactionForce(2) = force(2);
      }




			Matrix3d jaco;

      jaco = computeContactJacobian(x_0(0), x_0(1), x_0(2), x_1(0), x_1(1), x_1(2), 
        x_2(0), x_2(1), x_2(2), radius, dbar);

			jaco = jaco * stiffness;

			for (int ii = 0; ii < 3; ii++)
			{
				for (int jj = 0; jj < 3; jj++)
				{
					int ind1 = 4 * i + ii;
					int ind2 = 4 * i + jj;

					stepper->addJacobian(ind1, ind2, jaco(ii, jj));
				}
			}
		}

	}
}


Vector3d externalContactForce::computeContactForce(double x0, double y0, double z0, double x1, double y1, double z1, 
  double x2, double y2, double z2, double radius, double dbar)
{
  Vector3d vecResult;

  vecResult = ListVec(-((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
     (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
         pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
         pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
       (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
       sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
    ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
       (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
         2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
       log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
              pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
              pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
            sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
     (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
         pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
         pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
       sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
    ((dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
       (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
         2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
       log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
              pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
              pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
            sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
     (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
         pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
         pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
       sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))),
   -((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2)))/
     (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
         pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
         pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
       (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
       sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
    ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
       (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
         2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
       log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
              pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
              pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
            sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
     (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
         pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
         pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
       sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
    ((dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
       (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
         2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
       log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
              pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
              pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
            sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
     (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
         pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
         pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
       sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))),
   -((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))))/
     (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
         pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
         pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
       (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
       sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
    ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
         2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
       (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
       log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
              pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
              pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
            sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
     (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
         pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
         pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
       sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
    ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
         2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
       (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
       log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
              pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
              pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
            sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
     (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
         pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
         pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
       sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))));

  return vecResult;
}

Matrix3d externalContactForce::computeContactJacobian(double x0, double y0, double z0, double x1, double y1, double z1, 
  double x2, double y2, double z2, double radius, double dbar)
{
  Matrix3d matResult;

  matResult = ListMat(ListVec(-((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),
               2) + pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
              pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
            sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
         pow(2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
           2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2),2))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        pow(2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2),2))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        pow(2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2),2))/
      (4.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        pow(radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)),2)*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        pow(2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2),2))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*pow(y1 - y2,2) + 2*pow(-z1 + z2,2)))/
      (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     (pow(2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2),2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        pow(2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2),2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        pow(2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2),2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*pow(y1 - y2,2) + 2*pow(-z1 + z2,2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*pow(y1 - y2,2) + 2*pow(-z1 + z2,2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))),
    -(((-x1 + x2)*(y1 - y2)*(-dbar + radius - 
            sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
          (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))))/
        (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
          (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))) - 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (4.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        pow(radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)),2)*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((-x1 + x2)*(y1 - y2)*(-dbar + radius - 
          sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-x1 + x2)*(y1 - y2)*(dbar - radius + 
          sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))),
    -(((x1 - x2)*(-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + 
                 (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
          (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*(-z1 + z2))/
        (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
          (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (4.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        pow(radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)),2)*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((x1 - x2)*(-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + 
               (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*(-z1 + z2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((x1 - x2)*(dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),
              2) + pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*(-z1 + z2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))),
   ListVec(-(((-x1 + x2)*(y1 - y2)*(-dbar + radius - 
            sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
          (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))))/
        (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
          (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))) - 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (4.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        pow(radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)),2)*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((-x1 + x2)*(y1 - y2)*(-dbar + radius - 
          sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-x1 + x2)*(y1 - y2)*(dbar - radius + 
          sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))),
    -((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
              pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
              pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
            sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
         pow(2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
           2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2),2))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        pow(2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2),2))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        pow(2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2),2))/
      (4.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        pow(radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)),2)*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        pow(2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2),2))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*pow(-x1 + x2,2) + 2*pow(z1 - z2,2)))/
      (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     (pow(2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2),2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        pow(2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2),2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        pow(2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2),2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*pow(-x1 + x2,2) + 2*pow(z1 - z2,2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*pow(-x1 + x2,2) + 2*pow(z1 - z2,2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))),
    -((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
           2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
         (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
              pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
              pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
            sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
         (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
           2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2)))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2)))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2)))/
      (4.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        pow(radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)),2)*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2)))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-y1 + y2)*(-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + 
               (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*(z1 - z2))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((-y1 + y2)*(-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + 
               (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*(z1 - z2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-y1 + y2)*(dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + 
               (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*(z1 - z2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))),
   ListVec(-(((x1 - x2)*(-dbar + radius - 
            sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
          (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*(-z1 + z2))/
        (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
            pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
            pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
          (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
          sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (4.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        pow(radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)),2)*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2)))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((x1 - x2)*(-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + 
               (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*(-z1 + z2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((x1 - x2)*(dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),
              2) + pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*(-z1 + z2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2))*(y1 - y2) + 
          2*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2))*(-z1 + z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))),
    -((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
           2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
         (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
              pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
              pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
            sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
         (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
           2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2)))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2)))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2)))/
      (4.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        pow(radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)),2)*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2)))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-y1 + y2)*(-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + 
               (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*(z1 - z2))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (2*(-x1 + x2)*(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2)) + 
          2*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2))*(z1 - z2))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((-y1 + y2)*(-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + 
               (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*(z1 - z2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((-y1 + y2)*(dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + 
               (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*(z1 - z2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))),
    -(pow(2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
           2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)),2)*
         (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
              pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
              pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
            sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     (pow(2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)),2)*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     (pow(2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)),2)*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))))/
      (4.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        pow(radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)),2)*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     (pow(2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)),2)*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*pow(x1 - x2,2) + 2*pow(-y1 + y2,2))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))))/
      (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     (pow(2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)),2)*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        (pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     (pow(2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)),2)*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     ((2*pow(x1 - x2,2) + 2*pow(-y1 + y2,2))*
        (-dbar + radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) + 
     (pow(2*(x1 - x2)*((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2)) + 
          2*(-y1 + y2)*(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2)),2)*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (4.*pow(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2),1.5)*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2))) - 
     ((2*pow(x1 - x2,2) + 2*pow(-y1 + y2,2))*
        (dbar - radius + sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
             pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
             pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
           sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))*
        log((radius - sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
               pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
               pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))/
             sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))/dbar))/
      (2.*sqrt(pow(-((x0 - x2)*(y0 - y1)) + (x0 - x1)*(y0 - y2),2) + 
          pow((x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2),2) + 
          pow(-((y0 - y2)*(z0 - z1)) + (y0 - y1)*(z0 - z2),2))*
        sqrt(pow(x1 - x2,2) + pow(y1 - y2,2) + pow(z1 - z2,2)))));

  return matResult;
}

Vector3d externalContactForce::ListVec(double a1, double a2, double a3)
{
  Vector3d vecResult;

  vecResult(0) = a1;
  vecResult(1) = a2;
  vecResult(2) = a3;
  
  return vecResult;
}

Matrix3d externalContactForce::ListMat(Vector3d a1, Vector3d a2, Vector3d a3)
{
  Matrix3d matResult;

  matResult.col(0) = a1;
  matResult.col(1) = a2;
  matResult.col(2) = a3;
  
  return matResult;
}