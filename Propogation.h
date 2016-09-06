#ifndef PROPOGATION_H
#define PROPOGATION_H

#include "Eigen/Dense"
#include <iostream>
#include <time.h>	
//#include <unsupported/Eigen/Polynomials>

using namespace Eigen;
//SanjeevaReddy.pdf



MatrixXd rotateElasticConstants(float phi)
{
	// phi in radians
	// Coordinate system followed: Rokhlin et al. y-z plane on the screen, x coming out towards you.
	// Note that positive rotation about x means you stare down the positive x axis and then rotate clockwise.
	// This function rotates elastic constants by an angle phi about X (crystal) axis only.

	//phi = degtorad(phi);

	MatrixXd a(3, 3);
	a(0, 0) = 1;
	a(0, 1) = 0;
	a(0, 2) = 0;
	a(1, 0) = 0;
	a(1, 1) = cos(phi);
	a(1, 2) = sin(phi);
	a(2, 0) = 0;
	a(2, 1) = -sin(phi);
	a(2, 2) = cos(phi);

	MatrixXd M(6, 6);

	// Initialize M

	// Top left quarter.
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			M(i, j) = a(i, j)*a(i, j);
		}
	}

	// Bottom left quarter
	for (int i = 3; i < 6; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			M(i, j) = 1;
			for (int k = 0; k < 3; ++k)
			{
				if (k != (i - 3))	M(i, j) *= a(k, j);
			}
		}
	}

	// Top right quarter
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 3; j < 6; ++j)
		{
			M(i, j) = 2;
			for (int k = 0; k < 3; ++k)
			{
				if (k != (j - 3))	M(i, j) *= a(i, k);
			}
		}
	}

	// Bottom right quarter
	M(3, 3) = a(1, 1)*a(2, 2) + a(2, 1)*a(1, 2);
	M(3, 4) = a(1, 0)*a(2, 2) + a(2, 0)*a(1, 2);
	M(3, 5) = a(1, 0)*a(2, 1) + a(2, 0)*a(1, 1);

	M(4, 3) = a(0, 1)*a(2, 2) + a(2, 1)*a(0, 2);
	M(4, 4) = a(0, 0)*a(2, 2) + a(2, 0)*a(0, 2);
	M(4, 5) = a(0, 0)*a(2, 1) + a(2, 0)*a(0, 1);

	M(5, 3) = a(0, 1)*a(1, 2) + a(1, 1)*a(0, 2);
	M(5, 4) = a(0, 0)*a(1, 2) + a(1, 0)*a(0, 2);
	M(5, 5) = a(0, 0)*a(1, 1) + a(1, 0)*a(0, 1);


	
	//Phew. Generated M. Now, to do the rotation. Easy part, actually, thanks to awesome matrix libraries.
	return M*c*M.transpose();
}

void calculateAnisotropicInitial(ray &r, float current_grain_orientation = 0.0f)
{


	//    Note: As usual, whatever is passed is in degrees.
	current_grain_orientation = degtorad(current_grain_orientation);



	//    Assume that the direction information contained in the ray object is always with respect to the global axes. Yeah, because plotting will be easier in that case.
	//    printf("Current Direction : %f\n",radtodeg(r.global_direction));
	//    printf("The y : %f and the z: %f\n",(*r).init_position.y, (*r).init_position.z);
	//Rotate all the things. Boundary, ray, etc.
	//Basically the ray.
	//	(*r).local_direction = degtorad(90) - abs((*r).global_direction);	//Direction with respect to the vertical z-axis. Counter-clockwise is positive, and clockwise is negative.
	//	printf("Incident angle: %f\n", radtodeg((*r).local_direction));
	float n2 = cos(r.global_direction);
	float n3 = sin(r.global_direction);
	//	incident_angle = r.global_direction;
	//	incident_velocity = r.phase_velocity;
	//cout<<"Ray direction: "<<radtodeg(r.global_direction)<<"   ";
	MatrixXd CLower(6, 6);

	//Rotate the elastic constants in the medium.

	CLower = rotateElasticConstants(-current_grain_orientation); //NEGATIVE HERE!!!! - Because the angle of grains is in negative direction
																 //cout<<" Phase velocity: "<<r.phase_velocity<<endl;
	double m2o = n2 / r.phase_velocity;
	double m3o = n3 / r.phase_velocity;
	//cout<<"m3o = "<<m3o<<"     ";

	//cout<<"m2o = "<<m2o<<"m30 = "<<m3o<<" "<<current_grain_orientation<<endl;

	double m2 = m2o;//m2 component for old ray and new ray.

	double m3 = m3o;


	if (prev_grain_orientation != current_grain_orientation) 
	{ //This line needs to be commented every time ray path through weld is plotted. 
		first_check = 0;
	}

	//Going to check if zroots and lageur are necessary


	if (!first_check)
	{
		//    To find the m3 component now. Let m3 = B. Lower medium, coz we're trying to find the refracted ray and not the reflected one.
		MatrixXd a(3, 3), b(3, 3), d(3, 3);
		a << CLower(5, 5)*m2o*m2o - DENSITY, CLower(5, 1)*m2o*m2o, CLower(5, 3)*m2o*m2o,
			CLower(1, 5)*m2o*m2o, CLower(1, 1)*m2o*m2o - DENSITY, CLower(1, 3)*m2o*m2o,
			CLower(3, 5)*m2o*m2o, CLower(3, 1)*m2o*m2o, CLower(3, 3)*m2o*m2o - DENSITY;

		b << (CLower(4, 5) + CLower(5, 4))*m2o, (CLower(4, 1) + CLower(5, 3))*m2o, (CLower(4, 3) + CLower(5, 2))*m2o,
			(CLower(3, 5) + CLower(1, 4))*m2o, (CLower(3, 1) + CLower(1, 3))*m2o, (CLower(3, 3) + CLower(1, 2))*m2o,
			(CLower(2, 5) + CLower(3, 4))*m2o, (CLower(2, 1) + CLower(3, 3))*m2o, (CLower(2, 3) + CLower(3, 2))*m2o;

		d << CLower(4, 4), CLower(4, 3), CLower(4, 2),
			CLower(3, 4), CLower(3, 3), CLower(3, 2),
			CLower(2, 4), CLower(2, 3), CLower(2, 2);

		vector < complex <double> > sexticCoeff;
		Matrix3d tempMatrix;

		//    Power 0
		tempMatrix = a;
		sexticCoeff.push_back(tempMatrix.determinant());

		//    Power 1
		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[1] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[1] += tempMatrix.determinant();

		//    Power 2
		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[2] += tempMatrix.determinant();

		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[2] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[2] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[2] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[2] += tempMatrix.determinant();

		//    Power 3
		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[3] += tempMatrix.determinant();

		//    Power 4
		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff[4] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[4] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff[4] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[4] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[4] += tempMatrix.determinant();

		//    Power 5
		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff[5] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[5] += tempMatrix.determinant();

		//    Power 6
		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		//cout<< "A = " << a << endl << "B = " << b << endl << "D = "<< d << endl;



		vector < complex <double> > rootsOfSextic;
		rootsOfSextic.resize(6);
		zroots(sexticCoeff, rootsOfSextic, 1);

		///   Found six values of B (m3 components; m2 = m2o) for the lower medium. //FOR THE LOWER MEDIUM!!!
		//printf("\nThe roots of the sextic equation are as follows:\n");
	//    MatrixXd checkMatrix(3,3);

		//int findLongPhaseVelIndex = 0;  // Assuming the Phase velocity is largest for the long wave.

		double newPhaseVel = 0;
		int typeOfRay = 0;
		vector <double> slownessSolutions;
		for (int i = 0; i < 6; ++i)
		{
			if (!isOpposite(real(rootsOfSextic[i]), m3o))
			{ //Roots should have same sign as m3o. //CHECKING ALL SIX ROOTS
				slownessSolutions.push_back(real(rootsOfSextic[i]));
			}
		}


		sort(slownessSolutions.begin(), slownessSolutions.end(), std::greater<double>()); //Hence slownessSolutions contains the three feasible slowness values (sorted) for lower medium 
		int numsol = slownessSolutions.size();
	
		//cout<<"\nNumber of Solutions: "<<numsol<<endl;
		if (r.type == 0)
		{
			if (numsol < 3)
			{
				cout << "\nERROR!!! LONG DOES NOT EXIST" << endl;

			}
			else
			{
				m3 = slownessSolutions.at(0); vel_id = 0;
			}
		}
		else if (r.type == 1)
		{
			switch (numsol)
			{
			case 1: m3 = slownessSolutions.at(0); vel_id = 2; break;
			case 2: m3 = slownessSolutions.at(0); vel_id = 1; break;
			case 3: m3 = slownessSolutions.at(1); vel_id = 1; break;
			case 4:
				if (slownessSolutions.at(0) > -0.0001)
				{
					cout << endl << "ONE: " << endl;
					m3 = slownessSolutions.at(2); vel_id = 1;
				}
				else
				{
					cout << endl << "TWO: " << endl;
					m3 = slownessSolutions.at(3); vel_id = 2;
				}
				break;
			}
		}
		else if (r.type == 2)
		{
			switch (numsol)
			{
			case 2: m3 = slownessSolutions.at(1); vel_id = 2; break;
			case 3: m3 = slownessSolutions.at(2); vel_id = 2; break; //IMPORTANT LINE. REQUIRES FREQUENT CHANGING. 
			case 4:
				if (slownessSolutions.at(0) > -0.0001)
				{
					cout << endl << "ONE: " << endl;
					m3 = slownessSolutions.at(3); vel_id = 2;
				}
				else
				{
					cout << endl << "TWO: " << endl;
					m3 = slownessSolutions.at(2); vel_id = 1;
				}
				break;
			}
		}

		first_check = 1;


	
		//m3 = slownessSolutions.at(r.type);
		//if(abs(slownessSolutions.at(1)-slownessSolutions.at(2))<0.00001)
		//{
		//	cout<<"Yes";
		//	m3 = slownessSolutions.at(0);
		//}
		//else
		//{
		//	//Nothing
		//}
		//first_check = 1;
		//cout<<" m3 = "<<m3<<" "<<r.type<<" "<<slownessSolutions.at(0)<<" "<<slownessSolutions.at(1)<<" "<<slownessSolutions.at(2)<<endl;
	}
	else
	{
		m3 = m3o;
		//cout<<"here";
	}

	//cout<<"m3 = "<<m3<<endl;

	//At this point, m2 and m3 in the new medium are obtained properly!!
	double slownessMagnitude = sqrt(m2*m2 + m3*m3);
	double newn2 = m2 / slownessMagnitude;
	double newn3 = m3 / slownessMagnitude;

	//GET NEW PHASE VELOCITY DIRECTION

	r.global_direction = getAbsoluteAngle(newn2, newn3);


	//cout<<"New ray direction (phase): "<<radtodeg(r.global_direction);
	double newVel = 1 / slownessMagnitude;
	//            if(1/slownessMagnitude > newLongPhaseVel){
	//                newLongPhaseVel = 1/slownessMagnitude;
	//                findLongPhaseVelIndex = i;
	//            }

	//GET PHASE VELOCITY MAGNITUDE USING KNOWN PHASE VELOCITY DIRECTION IN THE NEW MEDIUM BY SOLVING CHRISTOFFEL'S EQUATION
	MatrixXd newChristoffelEquation(3, 3);
	newChristoffelEquation(0, 0) = CLower(5, 5)*newn2*newn2 + CLower(4, 4)*newn3*newn3;
	newChristoffelEquation(0, 1) = (CLower(3, 5) + CLower(1, 4))*newn2*newn3;
	newChristoffelEquation(0, 2) = CLower(3, 5)*newn2*newn2 + CLower(2, 4)*newn3*newn3;
	newChristoffelEquation(1, 0) = (CLower(1, 4) + CLower(3, 5))*newn2*newn3;
	newChristoffelEquation(1, 1) = CLower(1, 1)*newn2*newn2 + CLower(3, 3)*newn3*newn3;
	newChristoffelEquation(1, 2) = (CLower(1, 2) + CLower(3, 3))*newn2*newn3;
	newChristoffelEquation(2, 0) = CLower(3, 5)*newn2*newn2 + CLower(2, 4)*newn3*newn3;
	newChristoffelEquation(2, 1) = (CLower(3, 3) + CLower(1, 2))*newn2*newn3;
	newChristoffelEquation(2, 2) = CLower(3, 3)*newn2*newn2 + CLower(2, 2)*newn3*newn3;
	newChristoffelEquation = newChristoffelEquation / DENSITY;
	EigenSolver<MatrixXd> es;
	es.compute(newChristoffelEquation, 1);

	//cout<<christoffelEquation<<endl;

	
	//	printf("The velocities2 can be calculated as follows:\n");
	//  Choose the current ray type.
	//cout<<es.eigenvalues()<<endl;

	int phaseVelIndex = 0;
	double phaseVel = 0;
	vector <double> velocities;
	for (int i = 0; i < 3; ++i)
	{
		if (imag(es.eigenvalues()[i]) != 0) {
			//			printf("%d - imaginary eigenvalue.\n",i);
		}
		else if (real(es.eigenvalues()[i]) < 0) {
			//			printf("%d - negative. Therefore, V is imaginary.\n", i);
		}
		else {
			//			printf("%d - velocity = %f\n", i, sqrt(real(es.eigenvalues()[i])));
			velocities.push_back(sqrt(real(es.eigenvalues()[i])));
		}
	}
	sort(velocities.begin(), velocities.end(), std::greater<double>());  //Sort in reverse order so that 0 is long velocity, and so on.
	

	
	//if(shear_wide_angle){ //ALL THIS VALID ONLY FOR SHEAR1 WAVES.
	//	phaseVelIndex = 0;
	//}
	cout<<"\nChecking velocities: "<<velocities[0]<<" "<<velocities[1]<<" "<<velocities[2]<<" ";

	if (r.type == 0)
		phaseVelIndex = 0;
	else if (r.type == 1)
	{
		phaseVelIndex = 1;
		if (r.global_direction<205.4116 || r.global_direction>334.5884)
			phaseVelIndex = 2;
	}
	else if (r.type == 2)
	{
		phaseVelIndex = 2;
		if (r.global_direction<205.4116 || r.global_direction>334.5884)
			phaseVelIndex = 1;
	}

	phaseVel = velocities[phaseVelIndex];

	//cout<<"   New Phase Velocity: "<<phaseVel<<endl;
	//phaseVel = sqrt(real(es.eigenvalues()[0])); 
	//The while loop and the above line essentially pick out ONE of the three phase velocities obtained from the equation (eigenvalues). 
	//Which phase velocity is picked depends on the ray type. Longitudinal -> Highest velocity; SV -> Second; SH -> Slowest;
	//	refracted_angle = r.global_direction;
	r.phase_velocity = phaseVel; //important line: phase velocity of the ray in the new medium is set here.    
								 //	refracted_velocity = r.phase_velocity;

								 //Calculation of Transmission coefficient happens here

								 //	transmission_coefficient = transmission_coefficient * (2 * refracted_velocity * cos(incident_angle)) / ((refracted_velocity*cos(refracted_angle)) + (incident_velocity * cos(incident_angle)));
								 //	cout << "\nChecking Transmission Coefficient:" << transmission_coefficient << endl;
								 //GET GROUP VELOCITY MAGNITUDE AND DIRECTION
	double groupVel[4] = { 0 };
	double m[4] = { 0, 0, m2, m3 };
	double P[4] = { 0 };
	for (int i = 1; i < 4; ++i)
	{
		P[i] = real(es.eigenvectors()(i - 1, phaseVelIndex));
		//cout << "P - " << i << " "<< P[i] << endl;
	}
	for (int i = 1; i < 4; ++i) {
		groupVel[i] = 0;
		for (int l = 1; l < 4; ++l) {
			for (int j = 1; j < 4; ++j) {
				for (int k = 1; k < 4; k++) {
					//cout<<i<<" "<<1/DENSITY*VoigtC(C, i, j, k, l)*m[l]*P[j]*P[k]<<endl;
					groupVel[i] += 1 / DENSITY*VoigtC(CLower, i, j, k, l)*m[l] * P[j] * P[k];
				}
			}
		}
		//cout << "Group vel - " << i << " "<< groupVel[i] << endl;
	}
	r.group_vel_global_direction = getAbsoluteAngle(groupVel[2], groupVel[3]);
	//cout<<"Group Vel Direction =  "<<180/PI * r.group_vel_global_direction<<endl;
	r.group_velocity = sqrt(groupVel[2] * groupVel[2] + groupVel[3] * groupVel[3]);
	prev_grain_orientation = current_grain_orientation;
	//cout << "Group velocity magnitude : " << sqrt(groupVel[2]*groupVel[2] + groupVel[3]*groupVel[3]) << endl;
	//cout << "Direction of Group vel is : " << radtodeg(r.group_vel_global_direction) << endl;


}

void calculateAnisotropicReflection(ray &r, float current_grain_orientation = 0.0f)
{
	//    Note: As usual, whatever is passed is in degrees.
	current_grain_orientation = degtorad(current_grain_orientation);
	//    Assume that the direction information contained in the ray object is always with respect to the global axes. Yeah, because plotting will be easier in that case.
	//    printf("Current Direction : %f\n",radtodeg(r.global_direction));
	//    printf("The y : %f and the z: %f\n",(*r).init_position.y, (*r).init_position.z);
	//Rotate all the things. Boundary, ray, etc.
	//Basically the ray.
	//	(*r).local_direction = degtorad(90) - abs((*r).global_direction);	//Direction with respect to the vertical z-axis. Counter-clockwise is positive, and clockwise is negative.
	//	printf("Incident angle: %f\n", radtodeg((*r).local_direction));
	//	cout<<"\nREFLECTION: \nANGLE: "<<radtodeg(r.global_direction)<<endl;
	float n2 = cos(r.global_direction);
	float n3 = sin(r.global_direction);
	//	incident_angle = r.global_direction;
	//	incident_velocity = r.phase_velocity;
	//cout<<"Ray direction: "<<radtodeg(r.global_direction)<<"   ";
	MatrixXd CLower(6, 6);

	//Rotate the elastic constants in the medium.

	CLower = rotateElasticConstants(-current_grain_orientation); //NEGATIVE HERE!!!!   
																 //	cout<<"\nPhase velocity: "<<r.phase_velocity<<endl;
	double m2o = n2 / r.phase_velocity;
	double m3o = n3 / r.phase_velocity;
	//cout<<"m2o = "<<m2o<<" m3o = "<<m3o<<"     ";

	//cout<<"m2o = "<<m2o<<"m30 = "<<m3o<<" "<<current_grain_orientation<<endl;

	double m2 = m2o;//m2 component for old ray and new ray.
	double m3 = m3o;

	//if(prev_grain_orientation!=current_grain_orientation){ //This line needs to be commented every time ray path through weld is plotted. 

	first_check = 0;
	//}

	if (!first_check)
	{
		//    To find the m3 component now. Let m3 = B. Lower medium, coz we're trying to find the refracted ray and not the reflected one.
		MatrixXd a(3, 3), b(3, 3), d(3, 3);
		a << CLower(5, 5)*m2o*m2o - DENSITY, CLower(5, 1)*m2o*m2o, CLower(5, 3)*m2o*m2o,
			CLower(1, 5)*m2o*m2o, CLower(1, 1)*m2o*m2o - DENSITY, CLower(1, 3)*m2o*m2o,
			CLower(3, 5)*m2o*m2o, CLower(3, 1)*m2o*m2o, CLower(3, 3)*m2o*m2o - DENSITY;

		b << (CLower(4, 5) + CLower(5, 4))*m2o, (CLower(4, 1) + CLower(5, 3))*m2o, (CLower(4, 3) + CLower(5, 2))*m2o,
			(CLower(3, 5) + CLower(1, 4))*m2o, (CLower(3, 1) + CLower(1, 3))*m2o, (CLower(3, 3) + CLower(1, 2))*m2o,
			(CLower(2, 5) + CLower(3, 4))*m2o, (CLower(2, 1) + CLower(3, 3))*m2o, (CLower(2, 3) + CLower(3, 2))*m2o;

		d << CLower(4, 4), CLower(4, 3), CLower(4, 2),
			CLower(3, 4), CLower(3, 3), CLower(3, 2),
			CLower(2, 4), CLower(2, 3), CLower(2, 2);

		vector < complex <double> > sexticCoeff;
		Matrix3d tempMatrix;

		//    Power 0
		tempMatrix = a;
		sexticCoeff.push_back(tempMatrix.determinant());

		//    Power 1
		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[1] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[1] += tempMatrix.determinant();

		//    Power 2
		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[2] += tempMatrix.determinant();

		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[2] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[2] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[2] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[2] += tempMatrix.determinant();

		//    Power 3
		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[3] += tempMatrix.determinant();

		//    Power 4
		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff[4] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[4] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff[4] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[4] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[4] += tempMatrix.determinant();

		//    Power 5
		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff[5] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[5] += tempMatrix.determinant();

		//    Power 6
		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		//cout<< "A = " << a << endl << "B = " << b << endl << "D = "<< d << endl;

		//	printf("\nThe Co-efficients of the sextic equation are as follows:\n");
		//	for (int i=0; i<sexticCoeff.size(); ++i) {
		//	cout<<sexticCoeff[i]<<endl;
		//	}

		vector < complex <double> > rootsOfSextic;
		rootsOfSextic.resize(6);
		zroots(sexticCoeff, rootsOfSextic, 1);

		//Found six values of B (m3 components; m2 = m2o) for the lower medium. //FOR THE LOWER MEDIUM!!!

	

		//    MatrixXd checkMatrix(3,3);

		//int findLongPhaseVelIndex = 0;  // Assuming the Phase velocity is largest for the long wave.
		double newPhaseVel = 0;
		int typeOfRay = 0;
		vector <double> slownessSolutions;
		for (int i = 0; i < 6; ++i) {
			if (isOpposite(real(rootsOfSextic[i]), m3o))
			{ //Roots should have same sign as m3o. //CHECKING ALL SIX ROOTS
				slownessSolutions.push_back(real(rootsOfSextic[i]));
			}
		}

		sort(slownessSolutions.begin(), slownessSolutions.end(), std::less<double>()); //Hence slownessSolutions contains the three feasible slowness values (sorted) for lower medium 
		int numsol = slownessSolutions.size();
		//cout << "\nSlowness solutions: " << slownessSolutions << endl;
		//	cout<<"\n\nNumber of Solutions: "<<numsol<<endl;
		if (r.type == 0)
		{
			if (numsol < 3)
			{
				cout << "ERROR!!! LONG DOES NOT EXIST" << endl;
			}
			else
			{
				m3 = slownessSolutions.at(0); vel_id = 0;
			}
		}
		else if (r.type == 1)
		{
			switch (numsol)
			{
			case 1: m3 = slownessSolutions.at(0); vel_id = 2; break;
			case 2: m3 = slownessSolutions.at(0); vel_id = 1; break;
			case 3: m3 = slownessSolutions.at(1); vel_id = 1; break;
			case 4:
				if (slownessSolutions.at(0) > -0.0001)
				{
					cout << endl << "ONE: " << endl;
					m3 = slownessSolutions.at(2); vel_id = 1;
				}
				else {
					cout << endl << "TWO: " << endl;
					m3 = slownessSolutions.at(3); vel_id = 2;
				}
				break;
			}
		}
		else if (r.type == 2)
		{
			switch (numsol)
			{
			case 2: m3 = slownessSolutions.at(1); vel_id = 2; break;
			case 3: m3 = slownessSolutions.at(2); vel_id = 2; break; //IMPORTANT LINE. REQUIRES FREQUENT CHANGING. 
			case 4:
				if (slownessSolutions.at(0) > -0.0001)
				{
					cout << endl << "ONE: " << endl;
					m3 = slownessSolutions.at(3); vel_id = 2;
				}
				else {
					cout << endl << "TWO: " << endl;
					m3 = slownessSolutions.at(2); vel_id = 1;
				}
				break;
			}
		}

		first_check = 1;

	}
	else
	{
		m3 = m3o;
	}


	//cout<<"m3 = "<<m3<<endl;

	//At this point, m2 and m3 in the new medium are obtained properly!!
	double slownessMagnitude = sqrt(m2*m2 + m3*m3);
	double newn2 = m2 / slownessMagnitude;
	double newn3 = m3 / slownessMagnitude;

	//GET NEW PHASE VELOCITY DIRECTION
	r.global_direction = getAbsoluteAngle(newn2, newn3);
	double newVel = 1 / slownessMagnitude;
	//            if(1/slownessMagnitude > newLongPhaseVel){
	//                newLongPhaseVel = 1/slownessMagnitude;
	//                findLongPhaseVelIndex = i;
	//            }

	//GET PHASE VELOCITY MAGNITUDE USING KNOWN PHASE VELOCITY DIRECTION IN THE NEW MEDIUM BY SOLVING CHRISTOFFEL'S EQUATION
	MatrixXd newChristoffelEquation(3, 3);
	newChristoffelEquation(0, 0) = CLower(5, 5)*newn2*newn2 + CLower(4, 4)*newn3*newn3;
	newChristoffelEquation(0, 1) = (CLower(3, 5) + CLower(1, 4))*newn2*newn3;
	newChristoffelEquation(0, 2) = CLower(3, 5)*newn2*newn2 + CLower(2, 4)*newn3*newn3;
	newChristoffelEquation(1, 0) = (CLower(1, 4) + CLower(3, 5))*newn2*newn3;
	newChristoffelEquation(1, 1) = CLower(1, 1)*newn2*newn2 + CLower(3, 3)*newn3*newn3;
	newChristoffelEquation(1, 2) = (CLower(1, 2) + CLower(3, 3))*newn2*newn3;
	newChristoffelEquation(2, 0) = CLower(3, 5)*newn2*newn2 + CLower(2, 4)*newn3*newn3;
	newChristoffelEquation(2, 1) = (CLower(3, 3) + CLower(1, 2))*newn2*newn3;
	newChristoffelEquation(2, 2) = CLower(3, 3)*newn2*newn2 + CLower(2, 2)*newn3*newn3;
	newChristoffelEquation = newChristoffelEquation / DENSITY;
	EigenSolver<MatrixXd> es;
	es.compute(newChristoffelEquation, 1);

	
	int phaseVelIndex = 0;
	double phaseVel = 0;
	vector <double> velocities;
	for (int i = 0; i < 3; ++i)
	{
		if (imag(es.eigenvalues()[i]) != 0)
		{
			//			printf("%d - imaginary eigenvalue.\n",i);
		}
		else if (real(es.eigenvalues()[i]) < 0)
		{
			//			printf("%d - negative. Therefore, V is imaginary.\n", i);
		}
		else
		{
			//			printf("%d - velocity = %f\n", i, sqrt(real(es.eigenvalues()[i])));
			velocities.push_back(sqrt(real(es.eigenvalues()[i])));
		}
	}
	sort(velocities.begin(), velocities.end(), std::greater<double>());  //Sort in reverse order so that 0 is long velocity, and so on.

																		 //	cout<<"\nVelocities: "<<velocities[0]<<"\t"<<velocities[1]<<"\t"<<velocities[2]<<endl;
																		 //    Associating index with velocities. Good thing to do?
	while (velocities[vel_id] != sqrt(real(es.eigenvalues()[phaseVelIndex])))
	{
		phaseVelIndex++;
	}
	//if(shear_wide_angle){ //ALL THIS VALID ONLY FOR SHEAR1 WAVES.
	//	phaseVelIndex = 0;
	//}
	//cout<<endl<<velocities[0]<<" "<<velocities[1]<<" "<<velocities[2]<<" ";
	phaseVel = sqrt(real(es.eigenvalues()[phaseVelIndex]));

	//phaseVel = sqrt(real(es.eigenvalues()[0])); 
	//The while loop and the above line essentially pick out ONE of the three phase velocities obtained from the equation (eigenvalues). 
	//Which phase velocity is picked depends on the ray type. Longitudinal -> Highest velocity; SV -> Second; SH -> Slowest;
	//reflected_angle = r.global_direction;

	r.phase_velocity = phaseVel; //important line: phase velocity of the ray in the new medium is set here.    

								 //reflected_velocity = r.phase_velocity;
								 //incident_angle = degtorad(270) - incident_angle;
								 //reflected_angle = reflected_angle - degtorad(90);
								 //Calculation of Transmission coefficient happens here
								 //cout << "\nChecking angles in reflection function: " << radtodeg(incident_angle) << "\t " << radtodeg(reflected_angle)<<endl;
								 //transmission_coefficient = transmission_coefficient * ((refracted_velocity*cos(refracted_angle)) - (incident_velocity*cos(incident_angle))) / ((incident_velocity * cos(incident_angle)) + (refracted_velocity*cos(refracted_angle)));
								 //cout << "\nChecking Transmission Coefficient in reflection :" << transmission_coefficient << endl;

								 //GET GROUP VELOCITY MAGNITUDE AND DIRECTION
	double groupVel[4] = { 0 };
	double m[4] = { 0, 0, m2, m3 };
	double P[4] = { 0 };
	for (int i = 1; i < 4; ++i)
	{
		P[i] = real(es.eigenvectors()(i - 1, phaseVelIndex));
		//cout << "P - " << i << " "<< P[i] << endl;
	}
	for (int i = 1; i < 4; ++i)
	{
		groupVel[i] = 0;
		for (int l = 1; l < 4; ++l)
		{
			for (int j = 1; j < 4; ++j)
			{
				for (int k = 1; k < 4; k++)
				{
					//cout<<i<<" "<<1/DENSITY*VoigtC(C, i, j, k, l)*m[l]*P[j]*P[k]<<endl;
					groupVel[i] += 1 / DENSITY*VoigtC(CLower, i, j, k, l)*m[l] * P[j] * P[k];
				}
			}
		}
		//cout << "Group vel - " << i << " "<< groupVel[i] << endl;
	}
	r.group_vel_global_direction = getAbsoluteAngle(groupVel[2], groupVel[3]);
	//cout << "Group Vel Direction =  " << 180 / PI * r.group_vel_global_direction << endl;
	r.group_velocity = sqrt(groupVel[2] * groupVel[2] + groupVel[3] * groupVel[3]);
	prev_grain_orientation = current_grain_orientation;
	//cout << "Group velocity magnitude : " << sqrt(groupVel[2]*groupVel[2] + groupVel[3]*groupVel[3]) << endl;
	//cout << "Direction of Group vel is : " << radtodeg(r.group_vel_global_direction) << endl;
}




void calculateAnisotropic(ray &r, float upper_grain_orientation = 0.0f, float lower_grain_orientation = 0.0f)
{
	//Note: As usual, whatever is passed is in degrees.
	
//	ofstream file2; // For Validation by plotting graphs
	ofstream file1,file2;
	
	file2.open("file2.txt", std::ofstream::app);
	file1.open("file1.txt", std::ofstream::app);

	upper_grain_orientation = degtorad(upper_grain_orientation);
	lower_grain_orientation = degtorad(lower_grain_orientation);

	float incident_angle1;
	incident_angle1 = r.global_direction;
	float n2 = cos(r.global_direction);
	float n3 = sin(r.global_direction);

	//	incident_angle = r.global_direction;
	//	incident_velocity = r.phase_velocity;

	MatrixXd CUpper(6, 6);
	MatrixXd CLower(6, 6);

	//Rotate the elastic constants in the medium.
	CUpper = rotateElasticConstants(-upper_grain_orientation); //NEGATIVE HERE!!!!
	CLower = rotateElasticConstants(-lower_grain_orientation); //NEGATIVE HERE!!!!

															   //UPPER PORTION!!
	MatrixXd christoffelEquation(3, 3);
	christoffelEquation(0, 0) = CUpper(5, 5)*n2*n2 + CUpper(4, 4)*n3*n3;
	christoffelEquation(0, 1) = (CUpper(3, 5) + CUpper(1, 4))*n2*n3;
	christoffelEquation(0, 2) = CUpper(3, 5)*n2*n2 + CUpper(2, 4)*n3*n3;
	christoffelEquation(1, 0) = (CUpper(1, 4) + CUpper(3, 5))*n2*n3;
	christoffelEquation(1, 1) = CUpper(1, 1)*n2*n2 + CUpper(3, 3)*n3*n3;
	christoffelEquation(1, 2) = (CUpper(1, 2) + CUpper(3, 3))*n2*n3;
	christoffelEquation(2, 0) = CUpper(3, 5)*n2*n2 + CUpper(2, 4)*n3*n3;
	christoffelEquation(2, 1) = (CUpper(3, 3) + CUpper(1, 2))*n2*n3;
	christoffelEquation(2, 2) = CUpper(3, 3)*n2*n2 + CUpper(2, 2)*n3*n3;
	christoffelEquation = christoffelEquation / DENSITY;
	EigenSolver<MatrixXd> es;
	es.compute(christoffelEquation, 1);

	int phaseVelIndex = 0;
	float phaseVel = 0.0f;
	vector <double> velocities;

	for (int i = 0; i < 3; ++i) 
	{
		if (imag(es.eigenvalues()[i]) != 0) 
		{
			//printf("%d - imaginary eigenvalue.\n",i);
		}
		else if (real(es.eigenvalues()[i]) < 0) 
		{
			//printf("%d - negative. Therefore, V is imaginary.\n", i);
		}
		else 
		{
			//printf("%d - velocity = %f\n", i, sqrt(real(newes.eigenvalues()[i])));
			velocities.push_back(sqrt(real(es.eigenvalues()[i])));
		}
	}

	sort(velocities.begin(), velocities.end(), std::greater<double>());  //Sort in reverse order so that 0 is long velocity, and so on.

	file1 << velocities[0] << "\t" << velocities[1] << "\t" << velocities[2]<<endl;

	if (r.type == 0)
	{
		vel_id = 0;
		while (velocities[vel_id] != sqrt(real(es.eigenvalues()[phaseVelIndex]))) 
		{
			phaseVelIndex++;
		}
	}
	else if (r.type == 1)
	{
		for (int i = 0; i < 3; ++i)
		{
			if ((real(es.eigenvectors()(1, i)) == 0) && (real(es.eigenvectors()(2, i)) == 0))
			{
				phaseVelIndex = i;
			}
		}

	}
	else if (r.type == 2)
	{
		for (int i = 0; i < 3; ++i)
		{
			if (!((real(es.eigenvectors()(1, i)) == 0) && (real(es.eigenvectors()(2, i)) == 0)) && (velocities[0] != sqrt(real(es.eigenvalues()[i]))))
			{
				phaseVelIndex = i;
			}
		}

	}

	////BELOW FORMULATION IS VALID ONLY FOR SHEAR1 WAVES. BASICALLY, IN THAT ANGLE RANGE, SHEAR1 WAVES HAVE LARGET VELOCITY (corresponding to phasevelindex of 0)
	//
	if(((radtodeg(r.global_direction)<205.4116)	||	(radtodeg(r.global_direction)>334.5884))	&&	(r.type==1))
	{
		phaseVelIndex = 0; 
	}
	else if (((radtodeg(r.global_direction)<205.4116) || (radtodeg(r.global_direction)>334.5884)) && (r.type == 2))
	{
		phaseVelIndex = 1;
	}

//	cout << "\nChecking: " << sqrt(real(es.eigenvalues()[phaseVelIndex])) << "\t"<< phaseVelIndex<< endl;
	phaseVel = sqrt(real(es.eigenvalues()[phaseVelIndex]));

	double m2o = n2 / phaseVel;
	double m3o = n3 / phaseVel;

	double m2 = m2o;//m2 component for old ray and new ray.
	double m3;

//	if (upper_grain_orientation != lower_grain_orientation)
//	{
		first_check = 0;
	//}

	if (!first_check)
	{
		
		//  To find the m3 component now. Let m3 = B. Lower medium, coz we're trying to find the refracted ray and not the reflected one.
		MatrixXd a(3, 3), b(3, 3), d(3, 3);
		a << CLower(5, 5)*m2o*m2o - DENSITY, CLower(5, 1)*m2o*m2o, CLower(5, 3)*m2o*m2o,
			CLower(1, 5)*m2o*m2o, CLower(1, 1)*m2o*m2o - DENSITY, CLower(1, 3)*m2o*m2o,
			CLower(3, 5)*m2o*m2o, CLower(3, 1)*m2o*m2o, CLower(3, 3)*m2o*m2o - DENSITY;

		b << (CLower(4, 5) + CLower(5, 4))*m2o, (CLower(4, 1) + CLower(5, 3))*m2o, (CLower(4, 3) + CLower(5, 2))*m2o,
			(CLower(3, 5) + CLower(1, 4))*m2o, (CLower(3, 1) + CLower(1, 3))*m2o, (CLower(3, 3) + CLower(1, 2))*m2o,
			(CLower(2, 5) + CLower(3, 4))*m2o, (CLower(2, 1) + CLower(3, 3))*m2o, (CLower(2, 3) + CLower(3, 2))*m2o;

		d << CLower(4, 4), CLower(4, 3), CLower(4, 2),
			CLower(3, 4), CLower(3, 3), CLower(3, 2),
			CLower(2, 4), CLower(2, 3), CLower(2, 2);

		vector < complex <double> > sexticCoeff;
		Matrix3d tempMatrix;

		//    Power 0
		tempMatrix = a;
		sexticCoeff.push_back(tempMatrix.determinant());

		//    Power 1
		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[1] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[1] += tempMatrix.determinant();

		//    Power 2
		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[2] += tempMatrix.determinant();

		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[2] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[2] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[2] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[2] += tempMatrix.determinant();

		//    Power 3
		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[3] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[3] += tempMatrix.determinant();

		//    Power 4
		tempMatrix.col(0) = a.col(0);   //A1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff[4] += tempMatrix.determinant();

		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[4] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = a.col(1);   //A2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff[4] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[4] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = a.col(2);   //A3
		sexticCoeff[4] += tempMatrix.determinant();

		//    Power 5
		tempMatrix.col(0) = b.col(0);   //B1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = b.col(1);   //B2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff[5] += tempMatrix.determinant();

		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = b.col(2);   //B3
		sexticCoeff[5] += tempMatrix.determinant();

		//    Power 6
		tempMatrix.col(0) = d.col(0);   //D1
		tempMatrix.col(1) = d.col(1);   //D2
		tempMatrix.col(2) = d.col(2);   //D3
		sexticCoeff.push_back(tempMatrix.determinant());

		vector < complex <double> > rootsOfSextic;
		rootsOfSextic.resize(6);
		zroots(sexticCoeff, rootsOfSextic, 1);

		double newPhaseVel = 0;
		int typeOfRay = 0;
		vector <double> slownessSolutions;
		for (int i = 0; i < 6; ++i)
		{
			if (!isOpposite(real(rootsOfSextic[i]), m3o))
			{ //Roots should have same sign as m3o. //CHECKING ALL SIX ROOTS
				slownessSolutions.push_back(real(rootsOfSextic[i]));
			}
		}

		sort(slownessSolutions.begin(), slownessSolutions.end(), std::greater<double>()); //Hence slownessSolutions contains the three feasible slowness values (sorted) for lower medium 
		int numsol = slownessSolutions.size();
		//cout<<"!! "<<numsol;
		/*
		if(r.type ==0)
		{
			m3 = slownessSolutions.at(0); 
			vel_id = 0;
		}
		else if(r.type ==1)
		{
			m3 = slownessSolutions.at(1); 
			vel_id = 1;
		}
		*/
		
		if (r.type == 0)
		{
			if (numsol < 3)
			{
				cout << "ERROR!!! LONG DOES NOT EXIST" << endl;
			}
			else
			{
				m3 = slownessSolutions.at(0); vel_id = 0;
			}
		}
		else if (r.type == 1)
		{

			switch (numsol)
			{
			case 1: m3 = slownessSolutions.at(0); vel_id = 2; break;
			case 2: m3 = slownessSolutions.at(0); vel_id = 1; break;
			case 3: m3 = slownessSolutions.at(1); vel_id = 1; break;
			case 4:
				if (slownessSolutions.at(0) > -0.0001)
				{
					//cout<<endl<<"ONE"<<endl;
					m3 = slownessSolutions.at(2); vel_id = 1;
				}
				else {
					//cout<<endl<<"TWO"<<endl;
					m3 = slownessSolutions.at(3); vel_id = 2;
				}
				break;
			}

			if (numsol == 3)
			{
				if (abs(m2o) > 0.0003)
				{
					m3 = slownessSolutions.at(0);
					vel_id = 2;
					//cout<<"CHECK11";
				}
			}

		}
		else if (r.type == 2)
		{
			switch (numsol)
			{
			case 2:
				m3 = slownessSolutions.at(1);
				vel_id = 2;

				break;
			case 3: m3 = slownessSolutions.at(2); vel_id = 2; break;
			case 4:
				if (slownessSolutions.at(0) > -0.0001)
				{
					//cout<<" ONE ";
					m3 = slownessSolutions.at(3); vel_id = 2;
				}
				else {
					//cout<<endl<<"TWO"<<endl;
					m3 = slownessSolutions.at(2); vel_id = 1;
				}
				break;
			}
			if (numsol == 3)
			{
				if (abs(m3 - m3o) > 0.00005)
				{
					for (int i = 0; i < 3; ++i)
					{
						if (abs(slownessSolutions.at(i) - m3o) < 0.00005)
						{
							m3 = slownessSolutions.at(i);
							vel_id = 1;
							//cout<<" CHECK ";
						}
					}
				}

			}
		}
		
		
	}
	else
	{
		m3 = m3o;
	}
//	cout << "\nChecking 2: " << vel_id << endl;
	//At this point, m2 and m3 in the new medium are obtained properly!!
	double slownessMagnitude = sqrt(m2*m2 + m3*m3);
	double newn2 = m2 / slownessMagnitude;
	double newn3 = m3 / slownessMagnitude;

	//GET NEW PHASE VELOCITY DIRECTION
	r.global_direction = getAbsoluteAngle(newn2, newn3);
	//	cout<<"New ray direction (phase): "<<radtodeg(r.global_direction);
	double newVel = 1 / slownessMagnitude;
	//GET PHASE VELOCITY MAGNITUDE USING KNOWN PHASE VELOCITY DIRECTION IN THE NEW MEDIUM BY SOLVING CHRISTOFFEL'S EQUATION
	MatrixXd newChristoffelEquation(3, 3);
	newChristoffelEquation(0, 0) = CLower(5, 5)*newn2*newn2 + CLower(4, 4)*newn3*newn3;
	newChristoffelEquation(0, 1) = (CLower(3, 5) + CLower(1, 4))*newn2*newn3;
	newChristoffelEquation(0, 2) = CLower(3, 5)*newn2*newn2 + CLower(2, 4)*newn3*newn3;
	newChristoffelEquation(1, 0) = (CLower(1, 4) + CLower(3, 5))*newn2*newn3;
	newChristoffelEquation(1, 1) = CLower(1, 1)*newn2*newn2 + CLower(3, 3)*newn3*newn3;
	newChristoffelEquation(1, 2) = (CLower(1, 2) + CLower(3, 3))*newn2*newn3;
	newChristoffelEquation(2, 0) = CLower(3, 5)*newn2*newn2 + CLower(2, 4)*newn3*newn3;
	newChristoffelEquation(2, 1) = (CLower(3, 3) + CLower(1, 2))*newn2*newn3;
	newChristoffelEquation(2, 2) = CLower(3, 3)*newn2*newn2 + CLower(2, 2)*newn3*newn3;
	newChristoffelEquation = newChristoffelEquation / DENSITY;
	EigenSolver<MatrixXd> newes;
	newes.compute(newChristoffelEquation, 1);

	int phaseVelIndex2 = 0;
	phaseVel = 0.0f;
	vector <double> velocities2;
	for (int i = 0; i < 3; ++i)
	{
		if (imag(newes.eigenvalues()[i]) != 0) 
		{
			//			printf("%d - imaginary eigenvalue.\n",i);
		}
		else if (real(newes.eigenvalues()[i]) < 0) 
		{
			//			printf("%d - negative. Therefore, V is imaginary.\n", i);
		}
		else 
		{
			//			printf("%d - velocity = %f\n", i, sqrt(real(newes.eigenvalues()[i])));
			velocities2.push_back(sqrt(real(newes.eigenvalues()[i])));
		}
	}

	sort(velocities2.begin(), velocities2.end(), std::greater<double>());  //Sort in reverse order so that 0 is long velocity, and so on.
//	cout << "velocities2: " << velocities2[0] << " " << velocities2[1] << " " << velocities2[2] << " " << endl;
	//file2 << velocities2[0] << "\t" << velocities2[1] << "\t" << velocities2[2] << endl;

	

		////BELOW FORMULATION IS VALID ONLY FOR SHEAR1 WAVES. BASICALLY, IN THAT ANGLE RANGE, SHEAR1 WAVES HAVE LARGET VELOCITY (corresponding to phasevelindex of 0)
		//
		if (((radtodeg(r.global_direction)<205.4116) || (radtodeg(r.global_direction)>334.5884)) && (r.type == 1))
		{
			phaseVelIndex2 = 0;
		}
		else if (((radtodeg(r.global_direction)<205.4116) || (radtodeg(r.global_direction)>334.5884)) && (r.type == 2))
		{
			phaseVelIndex2 = 1;
		}


	phaseVel = sqrt(real(newes.eigenvalues()[phaseVelIndex2]));

	//The while loop and the above line essentially pick out ONE of the three phase velocities2 obtained from the equation (eigenvalues). 
	//Which phase velocity is picked depends on the ray type. Longitudinal -> Highest velocity; SV -> Second; SH -> Slowest;
	//refracted_angle = r.global_direction;
	r.phase_velocity = phaseVel; //important line: phase velocity of the ray in the new medium is set here.    
    //refracted_velocity = r.phase_velocity;
    //cout << "\nChecking angles in refraction function: " << radtodeg(incident_angle) << "\t " << radtodeg(refracted_angle) << endl;

    //Calculation of Transmission coefficient happens here

	//transmission_coefficient = transmission_coefficient * (2 * incident_velocity * cos(incident_angle)) / ((refracted_velocity*cos(refracted_angle)) + (incident_velocity * cos(incident_angle)));
    //cout << "\nChecking Transmission Coefficient in weld: " << transmission_coefficient << endl;

    //cout<<"\nNew pv = "<<phaseVel<<endl;

    //GET GROUP VELOCITY MAGNITUDE AND DIRECTION
	
	double groupVel[4] = { 0 };
	double m[4] = { 0, 0, m2, m3 };
	double P[4] = { 0 };
	for (int i = 1; i < 4; ++i) 
	{
		P[i] = real(newes.eigenvectors()(i - 1, phaseVelIndex2));
	}
	for (int i = 1; i < 4; ++i) 
	{
		groupVel[i] = 0;
		for (int l = 1; l < 4; ++l) 
		{
			for (int j = 1; j < 4; ++j) 
			{
				for (int k = 1; k < 4; ++k) 
				{
					groupVel[i] += 1 / DENSITY*VoigtC(CLower, i, j, k, l)*m[l] * P[j] * P[k];
				}
			}
		}
	}
	r.group_vel_global_direction = getAbsoluteAngle(groupVel[2], groupVel[3]);
	r.group_velocity = sqrt(groupVel[2] * groupVel[2] + groupVel[3] * groupVel[3]);
	file2 << radtodeg(r.group_vel_global_direction) << endl;

 //   file2 << radtodeg(incident_angle1) <<"\t"<<  r.phase_velocity<<"\t" << radtodeg(r.global_direction) << "\t" << radtodeg(r.group_vel_global_direction) << endl;
	file1.close();
	file2.close();


}







ray plotRayWeld(ray r, float eta, float T)
{
	/*
	ofstream file1, file2;

	file1.open("file1.txt", std::ofstream::trunc);
	file2.open("file2.txt", std::ofstream::trunc);
	file1.close();
	file2.close();
	*/
	//XR1.open("XR1.txt");

	vector <vertex_type> rP;
	cout << endl << "\n\t\t\t\t\t\t\tINSIDE WELD\n" << endl;
	double current_crystal_orientation, next_crystal_orientation, boundary_orientation;
	//r.current_position = r.init_position;
	//r.global_direction = r.init_global_direction;
	//r.phase_velocity = r.init_phase_velocity;

	//r.group_vel_global_direction = 0; //Why was this required?
	vertex_type v;

	//glBegin(GL_LINE_STRIP);
	//glVertex3f(r.current_position.x, r.current_position.y, r.current_position.z); //Mark only first point
	rP.clear();
	rP.push_back(r.current_position);
	cout << "\nRay Parameters inside Weld Before Entry:";
	cout << "\tY: " << r.current_position.y << "\tZ: " << r.current_position.z << "\tDirection: " << radtodeg(r.global_direction) << "\tVelocity: " << r.phase_velocity;
	ray_heading_temp = r.global_direction;
	//ray_heading_temp_group = r.group_vel_global_direction; //groupvelglodir was not initialised thus far. So used glodir because ray is JUST entering anisotropic region 

	//Calling master() function to initialize the slowness values for all 3601 angles
	//	Master();


	ray_heading_temp_group = r.global_direction;
	switch (weld_entry)
	{
	case 0: boundary_orientation = 0; break;
	case 1: boundary_orientation = degtorad(90 - A / 2); break;
	case 2: boundary_orientation = degtorad(90 + A / 2); break;
	}
	/*
	cout << "\n\nWELD ENTRY: ";
	if (weld_entry == 1)
		cout << "Right Side\n";
	else if (weld_entry == 2)
		cout << "Left Side\n";
	else if (weld_entry == 0)
		cout << "Middle\n";
	*/
	dummy_angle = radtodeg(r.global_direction);
	ray rotated_ray = r;
	rotated_ray.global_direction = r.global_direction - boundary_orientation;
	//rotated_ray.group_vel_global_direction = r.group_vel_global_direction - boundary_orientation;

	//In the rotated ray, only the phase velocity direction differs. (is rotated)  

	if (rotated_ray.global_direction<0)
	{
		rotated_ray.global_direction = rotated_ray.global_direction + degtorad(360);
	}
	else if ((rotated_ray.global_direction < degtorad(180)) && (rotated_ray.global_direction>0))
	{
		rotated_ray.global_direction = rotated_ray.global_direction + degtorad(180);
	}
	//	else if (rotated_ray.global_direction>degtorad(180))
	//	{
	//Do nothing
	//	}
	//incident_angle = r.global_direction - boundary_orientation + degtorad(90);    // Measured against the normal to the interface
	//incident_velocity = r.phase_velocity;


	if ((abs(rotated_ray.global_direction - degtorad(180)) < degtorad(0.03)) || (abs(rotated_ray.global_direction) < degtorad(0.03)))
	{
		//Do nothing. Ray is almost parallel to boundary => Bypasses it. 
	}
	else
	{
		//cout<<'\n'<<"CRYSTAL ORIENTATION INSIDE WELD\n";
		calculateAnisotropicInitial(rotated_ray,  radtodeg(getCrystalOrientation(rotated_ray.current_position.y, rotated_ray.current_position.z, T, eta)));
		//	calculateAnisotropic2(rotated_ray, radtodeg(getCrystalOrientation(rotated_ray.current_position.y, rotated_ray.current_position.z, T, eta)));
		//Above calculateAnisotropic2 call is the 'first' call meant for entering the weld. New phase velocity and group velocity are obtained.
	}
	//cout<<"\nChecking actual after just entry:"<<rotated_ray.current_position.y<<" "<<rotated_ray.current_position.z<<" "<<radtodeg(rotated_ray.global_direction)<<"\n";
	//cout<<endl<<"After "<<endl<<radtodeg(rotated_ray.global_direction)<<" "<<rotated_ray.phase_velocity<<" "<<radtodeg(rotated_ray.group_vel_global_direction);
	r = rotated_ray; //Copy all attributes: position, phasevel, groupvel, directions
//	cout << "\nRay Parameters After Entry - Initial: " <<"\tY: "<< r.current_position.y << "\tZ: " << r.current_position.z << "\tDirection: " << radtodeg(r.global_direction)<<"\tVelocity: " << r.phase_velocity;;
	//Then correct only the directions
	r.global_direction = rotated_ray.global_direction + boundary_orientation;
	r.group_vel_global_direction = rotated_ray.group_vel_global_direction + boundary_orientation;


	for (int i = 0; i<2; ++i)
	{
		if (abs(r.global_direction - ray_heading_temp)>degtorad(150))
		{
			r.global_direction = r.global_direction - degtorad(180);
		}
		if (abs(r.group_vel_global_direction - ray_heading_temp_group) > degtorad(50))
		{
			r.group_vel_global_direction = r.group_vel_global_direction - degtorad(180);
		}
	}

	refracted_angle = r.global_direction - boundary_orientation + degtorad(90);
	refracted_velocity = r.phase_velocity;

	//Calculation of Transmission coefficient happens here

	transmission_coefficient = transmission_coefficient * (2 * incident_velocity * cos(incident_angle)) / ((refracted_velocity*cos(refracted_angle)) + (incident_velocity * cos(incident_angle)));
	
//	cout << "\n\nChecking Transmission Coefficient after entering weld:" << transmission_coefficient << endl;



	temp_first = 0;
	tot_time_weld = 0;
	refl_check = 0;
	cout << endl << "\n\n\t\t\t\t\t\t\tSTARTING WELD LOOP\n";

	temp_count = 0;

	//file1.open("file1.txt", std::ofstream::app);

	while (isInsideWeld(rP.back()))
	{
		++temp_count;

		switch (r.type)
		{
		case 0:
			glColor3f(1, 0, 0);
			break;
		case 1:
			glColor3f(0, 1, 0);
			break;
		case 2:
			glColor3f(0, 0, 1);
			break;
		default:
			break;
		}
		//file1 << "\n" << radtodeg(r.global_direction) << "\t";
		
		//Changing incident angle and velocity after entering weld
		incident_velocity = r.phase_velocity;
		phase_diff = r.global_direction - phase_diff;
		//	incident_angle = r.global_direction;
		dummy_angle = radtodeg(r.global_direction);

		// Commenting all the lines which write the data to external text file. Please uncomment if necessary.
		//cout << endl << "y: " << r.current_position.y << "\tz: " << r.current_position.z << "\tDir: " << radtodeg(r.global_direction) << "\tGVel: " << radtodeg(r.group_vel_global_direction) << "\tPVel: " << r.phase_velocity << " ";


		//myfile<<"\ny: "<<r.current_position.y<<"\tz: "<<r.current_position.z<<"\tDir: "<<radtodeg(r.global_direction)<<"\tPVel: "<<r.phase_velocity<<"\tGVel: "<<radtodeg(r.group_vel_global_direction); //!!
		ray_heading_temp = r.global_direction;
		ray_heading_temp_group = r.group_vel_global_direction;

		current_crystal_orientation = (getCrystalOrientation(r.current_position.y, r.current_position.z, T, eta));
		//cout<<radtodeg(current_crystal_orientation);

		//myfile<<"\tCrystalAngle: "<<radtodeg(current_crystal_orientation)<<"\t";//!!
		next_crystal_orientation = (getCrystalOrientation(r.current_position.y + step_size*cos(r.group_vel_global_direction), r.current_position.z + step_size*sin(r.group_vel_global_direction), T, eta));
		//myfile<<"\tNextCrystalAngle: "<<radtodeg(next_crystal_orientation)<<"\t";//!!

		//Taking mean of orientations as the boundary orientation 
		boundary_orientation = (current_crystal_orientation + next_crystal_orientation) / 2;



		//ROTATION: Rotate clockwise if on the left half, and anti-clockwise if on the right side. Basically, add the negative of the boundary orientation.
		//Basic idea is to create a new ray which reflects the rotated scenario, calculate the change in the path and rotate it back.

		rotated_ray = r;
		rotated_ray.global_direction = r.global_direction - boundary_orientation;
		rotated_ray.group_vel_global_direction = r.group_vel_global_direction - boundary_orientation;

		//In the rotated ray, only the phase velocity direction differs. (is rotated)  

		//myfile<<"Rotated Ray:\n"<<"CrystalAngle: "<<radtodeg(current_crystal_orientation - boundary_orientation)<<"\tNextCrystalAngle: "<<radtodeg(next_crystal_orientation - boundary_orientation)<<"\t";//!!


		if (rotated_ray.global_direction<0)
		{
			rotated_ray.global_direction = rotated_ray.global_direction + degtorad(360);
		}
		else if ((rotated_ray.global_direction < degtorad(180)) && (rotated_ray.global_direction>0))
		{
			rotated_ray.global_direction = rotated_ray.global_direction + degtorad(180);
		}
		else if (rotated_ray.global_direction>degtorad(180))
		{
			//Do nothing
		}

		//myfile<<"Dir: "<<radtodeg(rotated_ray.global_direction)<<"\t";//!!

		// Can I take absolute? Yes you can
		//incident_angle = abs(r.global_direction - boundary_orientation + degtorad(90));

		if ((abs(rotated_ray.global_direction - degtorad(180)) < degtorad(0.03)) || (abs(rotated_ray.global_direction) < degtorad(0.03)))
		{
			//Do nothing. Ray is almost parallel to boundary => Bypasses it. 
		}
		else
		{
			calculateAnisotropic(rotated_ray, radtodeg(current_crystal_orientation - boundary_orientation), radtodeg(next_crystal_orientation - boundary_orientation)); // Check why we can't use anisotropic2 here
		}

	//	if(refl_check ==1)
	//	cout << "\nCHECKING:" << rotated_ray.phase_velocity << "  " << radtodeg(rotated_ray.global_direction);

		r = rotated_ray; //Copy all attributes: position, phasevel, groupvel, directions
		//Then correct only the directions
		//myfile<<radtodeg(rotated_ray.global_direction)<<"\t";//!!
		//myfile<<radtodeg(rotated_ray.group_vel_global_direction)<<"\t";//!!




		r.global_direction = rotated_ray.global_direction + boundary_orientation;
		r.group_vel_global_direction = rotated_ray.group_vel_global_direction + boundary_orientation;

		//cout<<endl<<"GV1"<<radtodeg(r.group_vel_global_direction);
		//cout<<endl<<"RayHeadingGV "<<radtodeg(ray_heading_temp_group);

		for (int i = 0; i<2; ++i)
		{
			if (abs(r.global_direction - ray_heading_temp)>degtorad(90))
			{
				r.global_direction = r.global_direction - degtorad(180);
			}
			if (abs(r.group_vel_global_direction - ray_heading_temp_group) > degtorad(90))
			{
				r.group_vel_global_direction = r.group_vel_global_direction - degtorad(180);
			}
		}
		phase_diff = r.global_direction - phase_diff;
		while ((ray_heading_temp_group - r.group_vel_global_direction) > degtorad(90))
		{
			r.group_vel_global_direction = r.group_vel_global_direction + degtorad(180);
		}

		total_time = total_time + (0.1) / r.phase_velocity;

		//refracted_angle = r.global_direction;
		// r.phase_velocity = phaseVel; important line: phase velocity of the ray in the new medium is set here.    
		refracted_velocity = r.phase_velocity;

		refracted_angle = abs(r.global_direction - boundary_orientation + degtorad(90));
		transmission_coefficient = transmission_coefficient * (2 * incident_velocity * cos(incident_angle)) / ((refracted_velocity*cos(refracted_angle)) + (incident_velocity * cos(incident_angle)));
		//	cout << "\nChecking Transmission Coefficient in weld:" << transmission_coefficient << endl;

		//	cout << "\nChecking angles in refraction function: " << radtodeg(incident_angle) << "\t " << radtodeg(refracted_angle) << endl;
		//	myfile << radtodeg(r.global_direction) << "\t" << r.phase_velocity << "\t"; //!!
		//	myfile << radtodeg(r.group_vel_global_direction) << "\n ";//!!





		//Amplitude Calculation happens here
		//Amplitude = (v2 - v1)/(v2 + v1)

		//	trans = 1 - ((v2 - v1) / (v2 + v1));

		//	cout << "\nAmplitude: " << trans << "\n";
		//cout << "\nCHECKING!!!: " << v1 << " " << v2 << "\n";


		v.x = rP.back().x;
		v.y = rP.back().y + step_size*cos(r.group_vel_global_direction);
		v.z = rP.back().z + step_size*sin(r.group_vel_global_direction);

		r.time = abs(step_size / r.group_velocity);
		tot_time_weld = tot_time_weld + r.time;
		r.current_position = v;
		rP.push_back(v);

	//	file1  << r.phase_velocity << "\t" <<  radtodeg(r.global_direction)<< "\t" <<radtodeg(r.group_vel_global_direction)  << "\n";
		

		if (((r.current_position.z <= 0) || (r.current_position.z >= PLATE_HEIGHT)))
		{
			//Perform Reflection
			//Changing incident angle and velocity after passing through grain
			incident_velocity = r.phase_velocity;
			incident_angle = r.global_direction - boundary_orientation + degtorad(90);

			phase_diff = r.global_direction - phase_diff;
			if (refl_check == 1)
			{
				break;
			}
			if (r.current_position.z >= PLATE_HEIGHT)
			{
				r.global_direction = r.global_direction + degtorad(180);
			}

			calculateAnisotropicReflection(r, radtodeg(getCrystalOrientation(r.current_position.y, r.current_position.z, T, eta)));

			if (r.current_position.z >= PLATE_HEIGHT)
			{
				r.global_direction = r.global_direction - degtorad(180);
				r.group_vel_global_direction = r.group_vel_global_direction - degtorad(180);
			}

			//Angle corrections
			if (radtodeg(r.global_direction) < 0)
				r.global_direction = r.global_direction + degtorad(360);
			if (radtodeg(r.group_vel_global_direction) < 0)
				r.group_vel_global_direction = r.group_vel_global_direction + degtorad(360);
			reflected_angle = r.global_direction - degtorad(90);

			//r.phase_velocity = phaseVel; //important line: phase velocity of the ray in the new medium is set here. 

			reflected_velocity = r.phase_velocity;
			//incident_angle = degtorad(270) - incident_angle;
			//reflected_angle = reflected_angle - degtorad(90);
			//Calculation of Transmission coefficient happens here
			cout << "\nRay Parameter After Reflection:" <<"\tY: "<< r.current_position.y << "\tZ: " << r.current_position.z << "\tDirection: " << radtodeg(r.global_direction) << "\tVelocity: " << r.phase_velocity << endl;
			//	transmission_coefficient =  ((incident_velocity*cos(incident_angle)) - (refracted_velocity*cos(refracted_angle))) / ((incident_velocity * cos(incident_angle)) + (refracted_velocity*cos(refracted_angle)));
			//	cout << "\nChecking Transmission Coefficient in reflection :" << transmission_coefficient << endl;
			phase_diff = r.global_direction - phase_diff;
			v.x = rP.back().x;
			v.y = rP.back().y + step_size*cos(r.group_vel_global_direction);
			v.z = rP.back().z + step_size*sin(r.group_vel_global_direction);
			r.time = abs(step_size / r.group_velocity);



			tot_time_weld = tot_time_weld + r.time;
			r.current_position = v;
			rP.push_back(v);
			refl_check = 1;


		}
		//	cout << "\nChecking time after weld:" << total_time << endl;
	}

	
	cout << endl << "Number of grains encountered = " << temp_count<< endl;
	cout << "\nRay Parameters Before Exiting Weld: " <<"\tY: "<< r.current_position.y << "\tZ: " << r.current_position.z << "\tDirection: " << radtodeg(r.global_direction) << "\tVelocity: " << r.phase_velocity << endl;;
	cout << endl << "\n\t\t\t\t\t\t\tEND OF LOOP";
	rayPaths.push_back(rP);
	glEnd();

//	file1.close();

	return r;

}


ray plotRayisotropic(ray r)
{
	//Find equations of both left and right weld boundaries
	cout << "\n\n\n==========================================================================================================================\n";
    cout << endl << "\n\t\t\t\t\t\tISOTROPIC PROPOGATION - BEFORE WELD" << endl;

	//float m1 = tan(degtorad(90-A/2));
	float c = -m1*D / 2; // y = mx + c at y=0 and x = D/2
	float c_ray;
	float side = 1; //1 -> Right, -1 -> Left
	vertex_type v;
	float new_direction = 0.0f;
	vector <vertex_type> rP;
	//glBegin(GL_LINE_STRIP);
	//glVertex3f(r.current_position.x, r.current_position.y, r.current_position.z); //Mark only first point
	//	rayPaths.clear();
	rP.clear();
	rP.push_back(r.current_position);
	tot_time_weld = 0;
	tot_time_iso_1 = 0;
	tot_time_iso_2 = 0;
	total_time = 0;
	phase_diff = r.global_direction;

	cout << "\nInside Wedge Parameters:" << "\ty: " << r.current_position.y << "\tz: " << r.current_position.z << "\t Direction: " << radtodeg(r.group_vel_global_direction) << "\t Velocity: " << r.phase_velocity << endl;

	//Wedge Propogation
	v.y = (PLATE_HEIGHT - r.current_position.z) / tan(r.group_vel_global_direction) + r.current_position.y;
	v.z = PLATE_HEIGHT;
	rP.push_back(v);
	//First entry into slab - 2 material interface
	calculateAnisotropicInitial(r, isotopic_angle);
	

	//Slab Propogation

	while (r.current_position.y > -PLATE_WIDTH / 2 && r.current_position.y < PLATE_WIDTH / 2)
	{
		cout << "\nBefore Parameters:" << "\ty: " << r.current_position.y << "\tz: " << r.current_position.z << "\t Direction: " << radtodeg(r.group_vel_global_direction) << "\t Velocity: " << r.phase_velocity << endl;
		if (isInsideWeld(rP.back()))
		{
	cout << endl << "Inside Weld => Anisotropy" << endl;
			break;
		}
		if (r.group_vel_global_direction > degtorad(180))
		{
			v.y = -r.current_position.z / tan(r.group_vel_global_direction) + r.current_position.y;
			v.z = 0;
		}
		else
		{
			v.y = (PLATE_HEIGHT - r.current_position.z) / tan(r.group_vel_global_direction) + r.current_position.y;
			v.z = PLATE_HEIGHT;
		}
		new_direction = degtorad(360) - r.group_vel_global_direction;

		//cout<<endl<<v.y<<" "<<v.z;
		if (isInsideWeld(v) || ((r.current_position.y*v.y) < 0))
		{
			c_ray = -r.current_position.y*tan(r.group_vel_global_direction) + r.current_position.z;

			if (r.current_position.y < 0)
				side = -1;
			else
				side = 1;

			if ((side*m1 - tan(r.group_vel_global_direction)) < (10 ^ -5))
			{
				cout << "Ray and weld almost parallel";
				break;
			}
			r.current_position.y = (c_ray - c) / (side*m1 - tan(r.group_vel_global_direction));
			r.current_position.z = side*m1*r.current_position.y + c;
			rP.push_back(r.current_position);
			total_time = total_time + sqrt((r.current_position.y - v.y)*(r.current_position.y - v.y) + (r.current_position.z - v.z)*(r.current_position.z - v.z)) / r.phase_velocity;

			if (side == 1)
				weld_entry = 1;
			else
				weld_entry = 2;
			break;

		}

		r.global_direction = new_direction;
		r.group_vel_global_direction = new_direction;
		//	r.time = sqrt((r.current_position.y - v.y)*(r.current_position.y - v.y) + (r.current_position.z - v.z)*(r.current_position.z - v.z)) / phase_velocity;
		tot_time_iso_1 = tot_time_iso_1 + r.time;

		r.current_position = v;
		//	phase_diff = r.global_direction - phase_diff;
		//cout<<endl<<r.current_position.y<<" "<<r.current_position.z<<" "<<r.current_position.x;
		rP.push_back(v);
	}
	rayPaths.push_back(rP);
   cout << "\n\nChecking time before weld: " << total_time*1000 <<" Micro Seconds"<< endl;

  cout << endl << "\n\t\t\t\t\tCOMPLETED ISOTROPIC RAY TRACING - BEFORE WELD" << endl;
  cout << "\n==========================================================================================================================\n";
	return r;
}
ray plotRayisotropicExit(ray r)
{
	//Find equations of both left and right weld boundaries
	cout << "\n==========================================================================================================================\n";
    cout << endl << "\n\t\t\t\t\t\tISOTROPIC PROPOGATION - AFTER WELD\n" << endl;
	//float m1 = tan(degtorad(90-A/2));
	float c = -m1*D / 2;
	float c_ray;
	float side = -1; //1 -> Right, -1 -> Left
	vertex_type v;
	float new_direction = 0.0f;
	//double boundary_orientation;
	ray_heading_temp = r.global_direction;
	vector <vertex_type> rP;
	rP.clear();

	rP.push_back(r.current_position);
	tot_time_iso_2 = 0;
	int k = 0;
	while (r.current_position.y >= (-PLATE_WIDTH / 2) && r.current_position.y <= (PLATE_WIDTH / 2) && r.current_position.z <= PLATE_HEIGHT && r.current_position.z >= 0.0)
	{
		
	cout << "\nRay Parameters After Exiting Weld: " << "\tY: " << r.current_position.y << "\tZ: " << r.current_position.z << "\tDirection: " << radtodeg(r.global_direction) << "\tVelocity: " << r.phase_velocity << endl;;

		/*if(isInsideWeld(rP.back()))
		{
		cout<<endl<<"Inside Weld => Anisotropy";
		break;
		}*/
		if (r.group_vel_global_direction > degtorad(180))
		{
			v.y = -r.current_position.z / tan(r.group_vel_global_direction) + r.current_position.y;
			v.z = 0;
		}
		else
		{
			/*Here there are three conditions to be checked:

			1. If the ray hits the top part of plate
			2. If the ray hits the side wall of plate
			3. If the ray exactly hits the corner

			This is checked by comparing the slope of ray with slope made by weld-ray point and the corner point
			*/

			float m_1 = radtodeg(atan((PLATE_HEIGHT - r.current_position.z) / (-PLATE_WIDTH / 2 - r.current_position.y)));
			if (m_1 < 0)
				m_1 = m_1 + 180;

			if (radtodeg(r.group_vel_global_direction) > m_1)
			{
				v.y = -PLATE_WIDTH / 2;
				v.z = tan(r.group_vel_global_direction) * (v.y - r.current_position.y);

			}
			else if (radtodeg(r.group_vel_global_direction) < m_1)
			{
				v.y = (PLATE_HEIGHT - r.current_position.z) / tan(r.group_vel_global_direction) + r.current_position.y;
				v.z = PLATE_HEIGHT;
			}

			else if (radtodeg(r.group_vel_global_direction) == m_1)
			{
				v.y = -PLATE_WIDTH / 2;
				v.z = PLATE_HEIGHT;
			}


		}
		new_direction = degtorad(360) - r.group_vel_global_direction;
		//	total_time = total_time + sqrt((r.current_position.y - v.y)*(r.current_position.y - v.y) + (r.current_position.z - v.z)*(r.current_position.z - v.z)) / r.phase_velocity;

		if (isInsideWeld(v) || ((r.current_position.y*v.y) < 0))
		{
			c_ray = -r.current_position.y*tan(r.group_vel_global_direction) + r.current_position.z;
			if (r.current_position.y < 0)
				side = -1;
			else
				side = 1;
			if ((side*m1 - tan(r.group_vel_global_direction)) < (10 ^ -5))
			{
				cout << "Ray and weld almost parallel";
				break;
			}
			r.current_position.y = (c_ray - c) / (side*m1 - tan(r.group_vel_global_direction));
			r.current_position.z = side*m1*r.current_position.y + c;
			rP.push_back(r.current_position);
			if (side == 1)
				weld_entry = 1;
			else
				weld_entry = 2;
			break;
		}
		r.global_direction = new_direction;
		r.group_vel_global_direction = new_direction;
		//	r.time = sqrt((r.current_position.y - v.y)*(r.current_position.y - v.y) + (r.current_position.z - v.z)*(r.current_position.z - v.z)) / phase_velocity;
		//	tot_time_iso_2 = tot_time_iso_2 + r.time;
		r.current_position = v;
		phase_diff = r.global_direction - phase_diff;
		total_time = total_time + sqrt((r.current_position.y - v.y)*(r.current_position.y - v.y) + (r.current_position.z - v.z)*(r.current_position.z - v.z)) / r.phase_velocity;
		rP.push_back(v);
		if (r.current_position.z == PLATE_HEIGHT || r.current_position.y == -PLATE_WIDTH / 2)
			break;
	}
	rayPaths.push_back(rP);

	r.global_direction = degtorad(360) - r.global_direction;  // This is becuase the ray angle is measured anticlockwise from horizontal axis
	// Amp calc
	amplitude = cos(r.global_direction - degtorad(90));
	amplitude = amplitude * transmission_coefficient;
   cout << "\nTime taken to travel inside weld: " << tot_time_weld * 1000 << " Micro seconds" << endl;
	cout << "\nTotal time taken to travel completely: " << total_time * 1000 << " Micro seconds" << endl;

	cout << "\nRay Parameters Final: " << "\tY: " << r.current_position.y << "\tZ: " << r.current_position.z << "\tDirection: " << radtodeg(r.global_direction) << "\tVelocity: " << r.phase_velocity << endl;;
	cout << endl << "\n\t\t\t\t\tCOMPLETED ISOTROPIC RAY TRACING - AFTER WELD" << endl;
	cout << "\n==========================================================================================================================\n";


	return r;
}


void displayWeldRayPath()
{
	for (size_t l = 0, rayl = rayPaths.size(); l < rayl; l++)
	{
		glBegin(GL_LINE_STRIP);
		switch (0) {
		case 0:
			glColor3f(1, 0, 0);
			break;
		case 1:
			glColor3f(0, 1, 0);
			break;
		case 2:
			glColor3f(0, 0, 1);
			break;
		default:
			break;
		}
		//cout<<endl<<"**********************************************";
		for (size_t i = 0, rayl1 = rayPaths[l].size(); i < rayl1; ++i)
		{
			//cout<<endl<<rayPaths[l][i].y;
			glVertex3f(0, rayPaths[l][i].y, rayPaths[l][i].z);
		}
		//cout<<endl<<"**********************************************";
		glEnd();
	}

}

void bruteforce(ray theRay_initial)
{
	ofstream graph;
	graph.open("3dGraph.txt", std::ofstream::trunc);
	graph.close();
	graph.open("3dGraph.txt", std::ofstream::app);
	ofstream odd;
	odd.open("Odd.txt", std::ofstream::trunc);
	odd.close();
	odd.open("Odd.txt", std::ofstream::app);

	//float T_temp = 0.5f, eta_temp = 0.0f;
	ray theRay_r;
	clock_t tStart = clock();
	float temp_time;
	temp_time = total_time;
	int q = 0;
	int i = 0;
	T_temp = 0.5f;
	eta_temp = 0.0f;
	//cout<<"\nChumma: "<<theRay.current_position.y<<" "<<theRay.current_position.z<<" "<<radtodeg(theRay.global_direction);
	cout << "\nTemp Time: " << temp_time<<endl;

	//myreadtext<<"\n\n\nT: "<<T_graph[i]<<"\tEta: 0-1"<<"\n";
	eta_temp = 0.0f;
	int j = 0;
	float y_diff = 0;
	while (i < 96) //96
	{
		eta_temp = -2.0f;
		int j = 0;
		while (j < 41)
		{
			
			cout << "\n\t\tCount : " << (i * 41) + j + 1;
			cout << "\n\n" << "Eta: " << eta_temp << "  T: " << T_temp;
			//			theRay_r = ray_prop(theRay_initial, eta_temp, T_temp);
			
			theRay_r = plotRayisotropic(theRay_initial);
			//cout<<"\nCHECKING(3):\nActual after iso: "<<theRay_temp.current_position.y<<" "<<theRay_temp.current_position.z<<" "<<radtodeg(theRay_temp.global_direction);
			theRay_r = plotRayWeld(theRay_r,eta_temp, T_temp);
			//cout<<"\nCHECKING(3):\nActual after weld: "<<theRay_temp.current_position.y<<" "<<theRay_temp.current_position.z<<" "<<radtodeg(theRay_temp.global_direction);
			theRay_r = plotRayisotropicExit(theRay_r);
			//cout<<"\nCHECKING(3):Actual final: "<<theRay_r.current_position.y<<" "<<theRay_r.current_position.z<<" "<<radtodeg(theRay_r.global_direction);

			if (abs(theRay_r.current_position.y) < 1.00005*abs(theRay.current_position.y) && abs(theRay_r.current_position.y) > 0.99995*abs(theRay.current_position.y))
			{ 
				if (abs(temp_time - total_time) < 0.00001) // Double Checking with time and position
				{

					T_final[q] = T_temp;
					eta_final[q] = eta_temp;
					y_final[q] = theRay_r.current_position.y;
					
				}
				q++;

			}

			
			j++;
			y_diff = (theRay_r.current_position.y - theRay.current_position.y);
			if( v == 0)
			graph << eta_temp << "\t" << T_temp << "\t" <<y_diff*y_diff << "\n";
			else if(v==1)
			odd << eta_temp << "\t" << T_temp << "\t" << y_diff*y_diff << "\n";
			v = 0;
			eta_temp = eta_temp + 0.1;
		}
		T_temp = T_temp + 0.1;
		i++;
	}


	if (q == 0)
		cout << "\n\nNot Found!";
	else
	{
		cout << "\nFound!";
		for (i = 0; i < q; i++)
		cout << "\nEta and T: " << eta_final[i] << " " << T_final[i]<<" V Value: "<<v;

	}
	cout << "\n\nTime taken: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << " seconds";
	graph.close();
	odd.close();
}

void testcpp()
{
	/*
	ray r;
	r = theRay;
	cout << "\nBefore Parameters:" << "\ty: " << r.current_position.y << "\tz: " << r.current_position.z << "\t Direction: " << radtodeg(r.group_vel_global_direction) <<"\t Velocity: "<<r.phase_velocity<< endl;
	calculateAnisotropicInitial(r, 0.0);
	cout << "\nAfter Parameters:" << "\ty: " << r.current_position.y << "\tz: " << r.current_position.z << "\t Direction: " << radtodeg(r.group_vel_global_direction) << "\t Velocity: " << r.phase_velocity << endl;

	*/
	

	ofstream file2,file1;
	file2.open("file2.txt", std::ofstream::trunc);
	file1.open("file1.txt", std::ofstream::trunc);
	file1.close();
	file2.close();
	


	for (float i = 0.0;i < 360.0;i += 5.0)
	{
		theRay.global_direction = degtorad(i);
		calculateAnisotropic(theRay, 13, 44);
	}
	
}
#endif