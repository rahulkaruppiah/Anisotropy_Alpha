#ifndef PROPAGATION1_H
#define PROPAGATION1_H

#include "Eigen/Dense"
//#include <unsupported/Eigen/Polynomials>
#include <iostream>
#include <time.h>	
#include <ppl.h>



using namespace Eigen;
//SanjeevaReddy.pdf
MatrixXd rotateElasticConstants(float phi)
{

	//phi in radians
	// Coordinate system followed: Rokhlin et al. y-z plane on the screen, x coming out towards you.
	//Note that positive rotation about x means you stare down the positive x axis and then rotate clockwise.
	//This function rotates elastic constants by an angle phi about X (crystal) axis only.

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

void Master()
{
	MatrixXd C1(6, 6);

	//Rotate the elastic constants in the medium.
	//CUpper = rotateElasticConstants(-upper_grain_orientation); //NEGATIVE HERE!!!!
	//CLower = rotateElasticConstants(-lower_grain_orientation); //NEGATIVE HERE!!!!

	MatrixXd christoffelEquation(3, 3);
	//RowVectorXd f(50);

	double n2, n3;
	vector <double> X_temp, Y_temp;
	double rayangle;
	double phi = 0.0f;
	C1 = c;

	int phaseVelIndex = 0;
	float phaseVel = 0.0f;


	vector <double> longslow;
	int y = 0;

	for (double i = 180.0; i < 360.05; i += 0.05)
	{
		vector <double> velocities;
		phi = PI / 180 * 0;
		rayangle = (double)PI / 180.0 * i;
		n2 = cos(rayangle);
		n3 = sin(rayangle);



		C1 = rotateElasticConstants(phi);

		christoffelEquation(0, 0) = C1(5, 5)*n2*n2 + C1(4, 4)*n3*n3;
		christoffelEquation(0, 1) = (C1(3, 5) + C1(1, 4))*n2*n3;
		christoffelEquation(0, 2) = C1(3, 5)*n2*n2 + C1(2, 4)*n3*n3;
		christoffelEquation(1, 0) = (C1(1, 4) + C1(3, 5))*n2*n3;
		christoffelEquation(1, 1) = C1(1, 1)*n2*n2 + C1(3, 3)*n3*n3;
		christoffelEquation(1, 2) = (C1(1, 2) + C1(3, 3))*n2*n3;
		christoffelEquation(2, 0) = C1(3, 5)*n2*n2 + C1(2, 4)*n3*n3;
		christoffelEquation(2, 1) = (C1(3, 3) + C1(1, 2))*n2*n3;
		christoffelEquation(2, 2) = C1(3, 3)*n2*n2 + C1(2, 2)*n3*n3;
		christoffelEquation = christoffelEquation / 7850; //Divide by density

		EigenSolver<MatrixXd> es;

		es.compute(christoffelEquation, 1);

		for (int j = 0; j < 3; ++j)
			velocities.push_back(sqrt(real(es.eigenvalues()[j])));


		sort(velocities.begin(), velocities.end(), std::greater<double>());  //Sort in reverse order so that 0 is long velocity, and so on.



		phaseVel = velocities[0];

		double slownessMagnitude = 1 / phaseVel;
		longslow.push_back(slownessMagnitude);
		X.push_back(slownessMagnitude*n2);
		Y.push_back(slownessMagnitude*n3);
		Y_temp.push_back(-slownessMagnitude*n3);
	}


	X_temp = X;
	std::reverse(X_temp.begin(), X_temp.end());

	X.reserve(2 * X.size());
	X.insert(X.end(), X_temp.begin(), X_temp.end());

	Y.reserve(2 * X.size());
	Y.insert(Y.end(), Y_temp.begin(), Y_temp.end());

}


void calculateAnisotropic2(ray &r, float current_grain_orientation = 0.0f)

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
	cout << "\nDirection with respect to vertical axis: " << r.global_direction << "\t " << degtorad(90) - abs(r.global_direction) << endl;
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


	if (prev_grain_orientation != current_grain_orientation){ //This line needs to be commented every time ray path through weld is plotted. 
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
		for (size_t i = 0; i < rootsOfSextic.size(); ++i)
		{
			cout << "\nRoots of sextic: " << rootsOfSextic[i] << endl;
		}

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
		for (int i = 0; i < numsol; ++i)
		{
			cout << "\nSlowness solutions from zroots: " << slownessSolutions[i] << endl;
		}
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


		cout << "m2o = " << m2o << " SS: ";
		for (size_t temp = 0; temp < slownessSolutions.size(); temp++)
		{
			cout << slownessSolutions.at(temp) << " ";
		}
		cout << endl << " m3o = " << m3o << " type =  " << r.type << " ";


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
	cout << "NewVel: " << newVel << endl;

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

	cout << "\nThe eigenvalues of A are:" << endl << (es.eigenvalues()) << endl;

	cout << "\nm3o" << m3o << endl;

	cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
	//	printf("The velocities2 can be calculated as follows:\n");
	//  Choose the current ray type.
	//cout<<es.eigenvalues()<<endl;

	int phaseVelIndex = 0;
	double phaseVel = 0;
	vector <double> velocities;
	for (int i = 0; i < 3; ++i)
	{
		if (imag(es.eigenvalues()[i]) != 0){
			//			printf("%d - imaginary eigenvalue.\n",i);
		}
		else if (real(es.eigenvalues()[i]) < 0){
			//			printf("%d - negative. Therefore, V is imaginary.\n", i);
		}
		else{
			//			printf("%d - velocity = %f\n", i, sqrt(real(es.eigenvalues()[i])));
			velocities.push_back(sqrt(real(es.eigenvalues()[i])));
		}
	}
	sort(velocities.begin(), velocities.end(), std::greater<double>());  //Sort in reverse order so that 0 is long velocity, and so on.
	cout << "\nVelocity check: " << r.phase_velocity;
	cout << "\nVELOCITIES: " << velocities[0] << " \t " << velocities[1] << " \t " << velocities[2] << endl;


	//  Associating index with velocities. Good thing to do?

	while (velocities[vel_id] != sqrt(real(es.eigenvalues()[phaseVelIndex])))
	{
		phaseVelIndex++;
	}

	//if(shear_wide_angle){ //ALL THIS VALID ONLY FOR SHEAR1 WAVES.
	//	phaseVelIndex = 0;
	//}
	//cout<<"\nChecking velocities: "<<velocities[0]<<" "<<velocities[1]<<" "<<velocities[2]<<" ";
	phaseVel = sqrt(real(es.eigenvalues()[phaseVelIndex]));

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

	cout << "\nCurrent grain orientation: " << current_grain_orientation << endl;
	cout << "\nCurrent ray velocity: " << r.phase_velocity << endl;
	cout << "\nCurrent ray orientation: " << r.global_direction << endl;
}

void calculateAnisotropic2refl(ray &r, float current_grain_orientation = 0.0f)
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

		printf("\nThe roots of the sextic equation are as follows:\n");
		for (size_t i = 0; i<rootsOfSextic.size(); ++i)
		{
			cout << rootsOfSextic[i] << endl;
		}


		//    MatrixXd checkMatrix(3,3);

		//int findLongPhaseVelIndex = 0;  // Assuming the Phase velocity is largest for the long wave.
		double newPhaseVel = 0;
		int typeOfRay = 0;
		vector <double> slownessSolutions;
		for (int i = 0; i < 6; ++i){
			if (isOpposite(real(rootsOfSextic[i]), m3o))
			{ //Roots should have same sign as m3o. //CHECKING ALL SIX ROOTS
				slownessSolutions.push_back(real(rootsOfSextic[i]));
			}
		}

		sort(slownessSolutions.begin(), slownessSolutions.end(), std::less<double>()); //Hence slownessSolutions contains the three feasible slowness values (sorted) for lower medium 
		int numsol = slownessSolutions.size();
		cout << "\nSlowness solutions :";
		for (std::vector<double>::iterator it = slownessSolutions.begin(); it != slownessSolutions.end(); ++it)
			std::cout << ' ' << *it;
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
				else{
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
				else{
					cout << endl << "TWO: " << endl;
					m3 = slownessSolutions.at(2); vel_id = 1;
				}
				break;
			}
		}

		first_check = 1;

		cout << "m2o = " << m2o << " SS: ";
		for (size_t temp = 0; temp < slownessSolutions.size(); temp++)
		{
			cout << slownessSolutions.at(temp) << " ";
		}
		cout << endl << " m3o = " << m3o << " type =  " << r.type << " ";

		cout << endl;
		cout << "vel_id = " << vel_id << endl;


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
	cout << "New ray direction (phase): " << radtodeg(r.global_direction);
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
	cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
	cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
	//	printf("The velocities2 can be calculated as follows:\n");
	//  Choose the current ray type.
	//cout<<es.eigenvalues()<<endl;

	for (int i = 0; i < 3; ++i)
	{

		cout << "The eigenvalues of AA are:" << endl << 1 / sqrt(real(es.eigenvalues()[i])) << endl;

	}
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
	cout << "Group Vel Direction =  " << 180 / PI * r.group_vel_global_direction << endl;
	r.group_velocity = sqrt(groupVel[2] * groupVel[2] + groupVel[3] * groupVel[3]);
	prev_grain_orientation = current_grain_orientation;
	//cout << "Group velocity magnitude : " << sqrt(groupVel[2]*groupVel[2] + groupVel[3]*groupVel[3]) << endl;
	//cout << "Direction of Group vel is : " << radtodeg(r.group_vel_global_direction) << endl;
}

void calculateAnisotropic4a(ray &r, float upper_grain_orientation = 0.0f, float lower_grain_orientation = 0.0f)
{
	//    Note: As usual, whatever is passed is in degrees.
	//cout << "\nChecking: 1. Upper crystal orientations:" << upper_grain_orientation << "\t2. " << lower_grain_orientation << endl;

	if (temp_count > 147 && temp_count < 150)
	{
		cout << "\nBefore solving (inside Anisotropic4a):";
		cout << "\nCurrent grain orientation: " << upper_grain_orientation << " " << lower_grain_orientation;
		cout << "\nCurrent ray velocity: " << r.phase_velocity;
		cout << "\nCurrent ray orientation: " << radtodeg(r.global_direction) << endl;
	}
	upper_grain_orientation = degtorad(upper_grain_orientation);
	lower_grain_orientation = degtorad(lower_grain_orientation);


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

	for (int i = 0; i < 3; ++i){
		if (imag(es.eigenvalues()[i]) != 0){
			//			printf("%d - imaginary eigenvalue.\n",i);
		}
		else if (real(es.eigenvalues()[i]) < 0){
			//			printf("%d - negative. Therefore, V is imaginary.\n", i);
		}
		else{
			//			printf("%d - velocity = %f\n", i, sqrt(real(newes.eigenvalues()[i])));
			velocities.push_back(sqrt(real(es.eigenvalues()[i])));
		}
	}
	//cout<<endl<<"velocities1: "<<velocities[0]<<" "<<velocities[1]<<" "<<velocities[2]<<" "<<endl;
	sort(velocities.begin(), velocities.end(), std::greater<double>());  //Sort in reverse order so that 0 is long velocity, and so on.
	//    Associating index with velocities. Good thing to do?

	if (r.type == 0)
	{
		vel_id = 0;
		while (velocities[vel_id] != sqrt(real(es.eigenvalues()[phaseVelIndex]))) {
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
	//if(((radtodeg(r.global_direction)<205) || (radtodeg(r.global_direction)>334))&&(r.type==1))
	//{
	//	phaseVelIndex = 0; 
	//}

	phaseVel = sqrt(real(es.eigenvalues()[phaseVelIndex]));

	// ORIGINAL COMMENTED PORTION
	//    std::vector < double > eigenv;
	//    for(int i=0; i<3; i++){
	//        eigenv.push_back(sqrt(real(es.eigenvalues()[i])));
	//    }
	//    sort(eigenv.begin(),eigenv.end());
	////    for(int i=0; i<eigenv.size(); i++){
	////        printf("%f\n", eigenv[i]);
	////    }
	//    float phaseVel = eigenv[2];    //Assuming this is the max. Safe, since it is sorted.
	//cout<<"n2 ="<<n2<<" n3 = "<<n3<<endl;


	double m2o = n2 / phaseVel;
	double m3o = n3 / phaseVel;

	//m2o and m3o FOUND!!

	//cout<<" m2o = "<<m2o<<"     "<<" m3o = "<<m3o<<" pv = "<<phaseVel<<endl;

	double m2 = m2o;//m2 component for old ray and new ray.
	double m3;
	if (upper_grain_orientation != lower_grain_orientation)
	{
		first_check = 0;
	}

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

		//	    printf("\nThe Co-efficients of the sextic equation are as follows:\n");
		//		for (int i=0; i<sexticCoeff.size(); ++i) 
		//	{
		//	cout<<sexticCoeff[i]<<endl;
		//	}


		vector < complex <double> > rootsOfSextic;
		rootsOfSextic.resize(6);
		zroots(sexticCoeff, rootsOfSextic, 1);

		///   Found six values of B (m3 components; m2 = m2o) for the lower medium. //FOR THE LOWER MEDIUM!!!
		//printf("\nThe roots of the sextic equation are as follows:\n");
		//for (int i=0; i<rootsOfSextic.size(); ++i) {
		//    cout<<rootsOfSextic[i]<<endl;
		//}

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
		//cout<<"!! "<<numsol;
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
				else{
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
				else{
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
		//cout<<"   "<<m3<<"   "<<vel_id<<endl;
		//	first_check = 1;


		//m3 = slownessSolutions.at(r.type);

		////cout<<slownessSolutions.size();
		////cout<<slownessSolutions.at(2);
		////cout<<"Whaaat?";
		//if(slownessSolutions.size()>2)
		//{
		//	if((abs(slownessSolutions.at(1)-slownessSolutions.at(2))<0.00001))
		//	{
		//		//cout<<"Yes";
		//		m3 = slownessSolutions.at(0);
		//	}
		//}
		//first_check = 1;
		////cout<<"abs";
	}
	else
	{
		m3 = m3o;
		l++;
		//cout<<"here";
	}

	//cout<<"m3 = "<<m3<<endl;

	//At this point, m2 and m3 in the new medium are obtained properly!!
	double slownessMagnitude = sqrt(m2*m2 + m3*m3);
	double newn2 = m2 / slownessMagnitude;
	double newn3 = m3 / slownessMagnitude;

	//GET NEW PHASE VELOCITY DIRECTION
	r.global_direction = getAbsoluteAngle(newn2, newn3);
	//	cout<<"New ray direction (phase): "<<radtodeg(r.global_direction);
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
	EigenSolver<MatrixXd> newes;
	newes.compute(newChristoffelEquation, 1);

	//cout<<christoffelEquation<<endl;
	//cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
	//cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
	//	printf("The velocities can be calculated as follows:\n");
	//  Choose the current ray type.
	//cout<<es.eigenvalues()<<endl;
	int phaseVelIndex2 = 0;
	phaseVel = 0.0f;
	vector <double> velocities2;
	for (int i = 0; i < 3; ++i)
	{
		if (imag(newes.eigenvalues()[i]) != 0){
			//			printf("%d - imaginary eigenvalue.\n",i);
		}
		else if (real(newes.eigenvalues()[i]) < 0){
			//			printf("%d - negative. Therefore, V is imaginary.\n", i);
		}
		else{
			//			printf("%d - velocity = %f\n", i, sqrt(real(newes.eigenvalues()[i])));
			velocities2.push_back(sqrt(real(newes.eigenvalues()[i])));
		}
	}
	//cout<<"velocities2: "<<velocities2[0]<<" "<<velocities2[1]<<" "<<velocities2[2]<<" "<<endl;
	sort(velocities2.begin(), velocities2.end(), std::greater<double>());  //Sort in reverse order so that 0 is long velocity, and so on.
	//    Associating index with velocities2. Good thing to do?
	while (velocities2[vel_id] != sqrt(real(newes.eigenvalues()[phaseVelIndex2]))) {
		phaseVelIndex2++;
	}


	phaseVel = sqrt(real(newes.eigenvalues()[phaseVelIndex2]));
	//	cout<<"   New Phase Velocity: "<<phaseVel<<endl;
	//phaseVel = sqrt(real(newes.eigenvalues()[0])); 





	//The while loop and the above line essentially pick out ONE of the three phase velocities2 obtained from the equation (eigenvalues). 
	//Which phase velocity is picked depends on the ray type. Longitudinal -> Highest velocity; SV -> Second; SH -> Slowest;
	//refracted_angle = r.global_direction;
	r.phase_velocity = phaseVel; //important line: phase velocity of the ray in the new medium is set here.    
	//refracted_velocity = r.phase_velocity;

	//cout << "\nChecking angles in refraction function: " << radtodeg(incident_angle) << "\t " << radtodeg(refracted_angle) << endl;

	//Calculation of Transmission coefficient happens here

	//	transmission_coefficient = transmission_coefficient * (2 * incident_velocity * cos(incident_angle)) / ((refracted_velocity*cos(refracted_angle)) + (incident_velocity * cos(incident_angle)));
	//cout << "\nChecking Transmission Coefficient in weld: " << transmission_coefficient << endl;

	//	cout<<"\nNew pv = "<<phaseVel<<endl;

	//GET GROUP VELOCITY MAGNITUDE AND DIRECTION
	double groupVel[4] = { 0 };
	double m[4] = { 0, 0, m2, m3 };
	double P[4] = { 0 };
	for (int i = 1; i < 4; ++i) {
		P[i] = real(newes.eigenvectors()(i - 1, phaseVelIndex2));
		//cout << "P - " << i << " "<< P[i] << endl;
	}
	for (int i = 1; i < 4; ++i) {
		groupVel[i] = 0;
		for (int l = 1; l < 4; ++l) {
			for (int j = 1; j < 4; ++j) {
				for (int k = 1; k < 4; ++k) {
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


	//cout << "\nAfter solving (inside Anisotropic 4a):";
	//	cout << "\nCurrent grain orientation: " << radtodeg(upper_grain_orientation) << " " << radtodeg(lower_grain_orientation);
	//	cout << "\nCurrent ray velocity: " << r.phase_velocity;
	cout << "\nCurrent ray orientation: " << radtodeg(r.global_direction) << "\t" << radtodeg(r.group_vel_global_direction) << endl;

	//prev_grain_orientation = current_grain_orientation;
	//cout << "Group velocity magnitude : " << sqrt(groupVel[2]*groupVel[2] + groupVel[3]*groupVel[3]) << endl;
	//cout << "Direction of Group vel is : " << radtodeg(r.group_vel_global_direction) << endl;

}


void RotateSlownessTop(vector <double> &XR1, vector <double> &YR1, float theta = 0.0)
{
	Matrix2d a;

	a(0, 0) = cos(theta);
	a(0, 1) = -sin(theta);
	a(1, 0) = sin(theta);
	a(1, 1) = cos(theta);

	vector <double> XR2(7202), YR2(7202);

	for (size_t i = 0; i < X.size(); ++i)
	{
		XR2[i] = (a(0, 0)*X[i] + a(0, 1)*Y[i]);
		YR2[i] = (a(1, 0)*X[i] + a(1, 1)*Y[i]);
		//	cout << "Displaying: " << XR2.at(i) << " " << YR2.at(i) << "\t" << i + 1 << endl;
	}
	size_t i = 0;
	size_t j = 0;

	//for (int q = 0; q < 74; ++q)
	//cout << YR2.at(q) << "\t" << q + 1 << endl;

	while (i < YR2.size())
	{
		if (i > 0)
			if ((XR2[i]>0) && (XR2[i - 1]<0))
				j = i;

		if (YR2[i] > (-1e-10) && YR2[i]<0.0)          //This is to eliminate the elements that are extremely small(basically zero) through the next code snippet- Get data from Anurag
			YR2[i] = 0.0;

		if (YR2[i] < (1e-10) && YR2[i]>0.0)
			YR2[i] = -YR2[i];                       //This is to eliminate the elements that are extremely small(basically zero) through the next code snippet

		if (YR2[i]<0.0)
		{
			YR2.erase(YR2.begin() + i);
			XR2.erase(XR2.begin() + i);
			--i;
		}
		++i;
	}
	size_t k = 0;
	while (j < YR2.size())
	{
		YR1[k] = (YR2[j]);
		XR1[k] = (XR2[j]);
		YR2.erase(YR2.begin() + j);
		XR2.erase(XR2.begin() + j);
		++k;
	}
	//	for (int i = 0; i < (int)YR3.size(); i++)
	//	cout << YR3.at(i) << endl;
	size_t l = 0;
	for (; k < XR1.size(); ++k)
	{
		XR1[k] = XR2[l];
		YR1[k] = YR2[l];
		l++;
	}
}


void RotateSlownessBottom(vector <double> &XR1, vector <double> &YR1, float theta)
{
	Matrix2d a;

	a(0, 0) = cos(theta);
	a(0, 1) = -sin(theta);
	a(1, 0) = sin(theta);
	a(1, 1) = cos(theta);

	vector <double> XR2(7202), YR2(7202);

	for (size_t i = 0; i < X.size(); ++i)
	{
		XR2[i] = (a(0, 0)*X[i] + a(0, 1)*Y[i]);
		YR2[i] = (a(1, 0)*X[i] + a(1, 1)*Y[i]);
		//cout << "Displaying: " << s3.at(i) << " " << s4.at(i) << "\t" << i + 1 << endl;
	}

	size_t i = 0;
	size_t j = 0;

	while (i < YR2.size())
	{
		if (i>0)
			if (XR2[i]<0 && XR2[i - 1]>0)
				j = i;

		if (abs(YR2[i]) < (1e-10))
			YR2[i] = -YR2[i];                       //This is to eliminate the elements that are extremely small(basically zero) through the next code snippet

		if (YR2[i] > 0.0)
		{
			YR2.erase(YR2.begin() + i);
			XR2.erase(XR2.begin() + i);
			i = i - 1;
		}
		++i;
	}
	size_t k = 0;



	while (j < YR2.size())
	{
		YR1[k] = (YR2[j]);
		XR1[k] = (XR2[j]);
		YR2.erase(YR2.begin() + j);
		XR2.erase(XR2.begin() + j);
		++k;
	}
	size_t l = 0;

	for (; k < XR1.size(); ++k)
	{
		XR1[k] = XR2[l];
		YR1[k] = YR2[l];
		++l;
	}

	//cout << "Checking for crashes: " << k << " " << l << " " << XR1.size();

}










void calculateAnisotropic4c(ray &r, double theta1 = 0.0, double theta2 = 0.0)
{

	theta1 = degtorad(theta1);
	theta2 = degtorad(theta2);

	vector <double> XR1(3601), YR1(3601), XR2(3601), YR2(3601);
	//	double XR1, YR1;
	float incident_angle;
	double x_i, y_i, x_o, normal_length, slowness_vector_magnitude, phase_velocity_magnitude;
	float phase_direction, slope, normal_slope;

	//	if (refl_check == 1)
	//	{
	RotateSlownessTop(XR1, YR1, theta1);
	RotateSlownessBottom(XR2, YR2, theta2);
	//	}
	//	else if (refl_check == 0)
	//	{
	//		RotateSlownessTop(XR2, YR2, theta2);
	//		RotateSlownessBottom(XR1, YR1, theta1);
	//	}

	phase_velocity_magnitude = r.phase_velocity;
	slowness_vector_magnitude = 1 / phase_velocity_magnitude;

	incident_angle = r.global_direction;

	if (incident_angle > PI)
		incident_angle = incident_angle - PI;

	//XR1 = cos(incident_angle) * slowness_vector_magnitude;
	//YR1 = sin(incident_angle) * slowness_vector_magnitude;

	//cout << "Displaying: " << XR1 << " " << YR1;
	//    The following two while loops essentially perform the tasks described below.
	//    while loop 1 = > find point of intersection of incident ray with upper slowness surface knowing angle of incidence.
	//    while loop 1 = > find x coordinate of point on lower slowness surface which is negative of previous obtained x coordinate.

	double increment = 0.0;
	size_t i = 0;

	//		for (int i = 0; i < (int)YR1.size(); i++)
	//		cout << YR1[i] << endl;

	//while loop #1

	while (i < XR1.size())
	{
		if (abs(tan(incident_angle) - (YR1[i] / XR1[i]))< (0.0005 + increment))
			break;
		++i;
		if (i == (XR2.size() - 1))
		{
			increment = increment + 0.0001;
			i = 0;
			continue;
		}
	}

	x_i = XR1[i];
	y_i = YR1[i];
	x_o = -x_i;
	//	cout << "Found Xo: " << x_o << "  " << y_i << " "<<XR2.size()<<endl;

	increment = 0.0000004;
	i = 0;

	//while loop #2
	while (i < XR2.size())
	{
		if (abs(x_o - XR2[i]) < (increment))
		{
			break;
		}

		++i;
		if (i == (XR2.size() - 1))
		{
			increment = increment + 0.0000004;
			i = 0;
			continue;
		}
	}

	//m2 and m3 in the new medium is found here
	//	cout << "Found XR2 and YR2: " << XR2[i] << "  " << YR2[i] << endl;
	double m2 = XR2[i];
	double m3 = YR2[i];
	double c = 0.0, m = 0.0, gx = 0.0, gy = 0.0, gy1 = 0.0;
	// Calculating direction of phase velocity from point of intersection.
	phase_direction = atan(YR2[i] / XR2[i]);

	if (phase_direction >= 0)
		phase_direction = phase_direction + PI;
	else
		phase_direction = (2 * PI) + phase_direction;

	// Calculating Slope of slowness surface at point of intersection and the normal to it.

	if (i == 0)
	{
		slope = atan((YR2[i + 1] - YR2[i]) / (XR2[i + 1] - XR2[i]));
	}
	else if (i == YR2.size())
	{
		slope = atan((YR2[i] - YR2[i - 1]) / (XR2[i] - XR2[i - 1]));
	}
	else
	{
		slope = atan((YR2[i + 1] - YR2[i - 1]) / (XR2[i + 1] - XR2[i - 1]));
	}


	//	slope = radtodeg(slope);

	//normal_slope = slope + (PI / 2) + PI;


	if (slope >= 0)
		normal_slope = slope - (PI / 2) + 2 * PI;
	else
		normal_slope = slope + (PI / 2) + PI;


	m = tan(slope) - tan(normal_slope);
	c = YR2[i] - tan(normal_slope)*XR2[i];
	gx = c / m;
	gy = tan(normal_slope)*gx + c;
	//	gy1 = tan(degtorad(slope))*gx;

	double group_slowness_magnitude = sqrt(((gy - m3)*(gy - m3)) + ((gx - m2)*(gx - m2)));
	double group_velocity = 1 / group_slowness_magnitude;
	//cout << "Variables: " << gx << " " << gy << " " << c << " " << m << " " << gy1<<" "<<endl;

	//	normal_slope = degtorad(normal_slope);

	normal_length = 0.2*(10 ^ (-4));
	//uble group_velocity = sqrt((g2*g2) + (g3*g3));
	//Results Displayed.



	slowness_vector_magnitude = sqrt((XR2[i] * XR2[i]) + (YR2[i] * YR2[i])); // * (10 ^ -4);
	phase_velocity_magnitude = 1 / slowness_vector_magnitude;

	r.global_direction = phase_direction;
	r.phase_velocity = phase_velocity_magnitude;
	r.group_vel_global_direction = normal_slope;
	r.group_velocity = group_velocity;


	cout << "\nPhase Velocity Direction = " << radtodeg(phase_direction) << endl;
	//cout << "\nPhase velocity Magnitude = " << phase_velocity_magnitude << endl;
	cout << "\nGroup Velocity Direction = " << radtodeg(normal_slope) << endl;
	//cout << "\nGroup Velocity Magnitude = " << 1 / group_velocity << endl;


}


void calculateAnisotropic4crefl(ray &r, double theta1 = 0.0, double theta2 = 0.0)
{
	//cout << "Theta:" << theta1 << " " << theta2;

	theta1 = degtorad(theta1);
	theta2 = degtorad(theta2);



	//	MatrixXd C1(6, 6);
	//	RowVectorXd XR1, YR1;
	//	RowVectorXd XR2, YR2;
	vector <double> XR1(3601), YR1(3601), XR2(3601), YR2(3601);

	double incident_angle;
	double x_i, y_i, x_o, normal_length, slowness_vector_magnitude, phase_velocity_magnitude;
	double phase_direction, slope, normal_slope;
	//Rotate the elastic constants in the medium.
	//CUpper = rotateElasticConstants(-upper_grain_orientation); //NEGATIVE HERE!!!!
	//CLower = rotateElasticConstants(-lower_grain_orientation); //NEGATIVE HERE!!!!

	//	MatrixXd christoffelEquation(3, 3);

	//	Matrix<double, 2, Dynamic, RowMajor> s1(2, 3601);
	//	Matrix<double, 2, Dynamic, RowMajor> s2(2, 3601);

	RotateSlownessTop(XR1, YR1, theta1);
	//	RotateSlownessBottom(XR2, YR2, theta2);

	//	XR1 = s1.row(0);
	//	YR1 = s1.row(1);

	//	XR2 = s2.row(0);
	//	YR2 = s2.row(1);

	//for (int i = 0; i < 3601;++i)
	//cout << XR1(i)<<endl;


	incident_angle = r.global_direction;

	if (incident_angle > PI)
		incident_angle = incident_angle - PI;

	//cout << "Incident angle: " << incident_angle;
	//    The following two while loops essentially perform the tasks described below.
	//    while loop 1 = > find point of intersection of incident ray with upper slowness surface knowing angle of incidence.
	//    while loop 1 = > find x coordinate of point on lower slowness surface which is negative of previous obtained x coordinate.

	double increment_i = 0.0;
	size_t i = 0;

	//		for (int i = 0; i < (int)YR1.size(); i++)
	//		cout << YR1[i] << endl;

	//while loop #1
	while (i < XR1.size())
	{
		if (abs(tan(incident_angle) - (YR1[i] / XR1[i]))< (0.0005 + increment_i))
			break;
		++i;
		if (i == (XR2.size() - 1))
		{
			increment_i = increment_i + 0.0001;
			i = 0;
			continue;
		}
	}

	cout << "Found XR1 and YR1: " << XR1[i] << "  " << YR1[i] << endl;
	XR2 = XR1;
	YR2 = YR1;
	x_i = XR1[i];
	y_i = YR1[i];
	x_o = -x_i;
	//	cout << "Found Xo: " << x_o << "  " << y_i << " "<<XR2.size()<<endl;

	double increment = 0.0000004;
	i = 0;

	//while loop #2
	while (i < XR2.size())
	{
		if (abs(x_o - XR2[i]) < (increment))
		{
			//	cout << "Ethi"<<endl;
			break;
		}

		++i;
		//	cout << "i: " << i << endl;
		if (i == (XR2.size() - 1))
		{
			increment = increment + 0.0000004;
			i = 0;
			continue;
		}
	}

	//m2 and m3 in the new medium is found here
	cout << "Found XR2 and YR2: " << XR2[i] << "  " << YR2[i] << endl;
	double m2 = XR2[i];
	double m3 = YR2[i];
	double c = 0.0, m = 0.0, gx = 0.0, gy = 0.0, gy1 = 0.0;
	// Calculating direction of phase velocity from point of intersection.
	phase_direction = atan(YR2[i] / XR2[i]);

	//Since its reflection the angle should be between under 180



	if (phase_direction < 0)
		phase_direction = (PI)+phase_direction;
	else if (phase_direction>PI)
		phase_direction = phase_direction - PI;




	// Calculating Slope of slowness surface at point of intersection and the normal to it.

	if (i == 0)
	{
		slope = atan((YR2[i + 1] - YR2[i]) / (XR2[i + 1] - XR2[i]));
	}
	else if (i == YR2.size())
	{
		slope = atan((YR2[i] - YR2[i - 1]) / (XR2[i] - XR2[i - 1]));
	}
	else
	{
		slope = atan((YR2[i + 1] - YR2[i - 1]) / (XR2[i + 1] - XR2[i - 1]));
	}


	//slope = radtodeg(slope);





	normal_slope = slope + (PI / 2);

	m = tan(slope) - tan(normal_slope);
	c = YR2[i] - tan(normal_slope)*XR2[i];
	gx = c / m;
	gy = tan(normal_slope)*gx + c;
	//gy1 = tan(degtorad(slope))*gx;
	double group_slowness_magnitude = sqrt(((gy - m3)*(gy - m3)) + ((gx - m2)*(gx - m2)));
	double group_velocity = 1 / group_slowness_magnitude;
	//cout << "Variables: " << gx << " " << gy << " " << c << " " << m << " " << gy1<<" "<<endl;

	//	normal_slope = degtorad(normal_slope);

	normal_length = 0.2*(10 ^ (-4));
	//uble group_velocity = sqrt((g2*g2) + (g3*g3));
	//Results Displayed.



	slowness_vector_magnitude = sqrt((XR2[i] * XR2[i]) + (YR2[i] * YR2[i])); // * (10 ^ -4);
	phase_velocity_magnitude = 1 / slowness_vector_magnitude;


	r.global_direction = phase_direction;
	r.phase_velocity = phase_velocity_magnitude;
	r.group_vel_global_direction = normal_slope;
	r.group_velocity = group_velocity;

	//slowness_vector_magnitude*sin(pi / 180 * phase_direction)
	cout << "\nPhase Velocity Direction = " << phase_direction << endl;
	cout << "\nPhase velocity Magnitude = " << phase_velocity_magnitude << endl;
	cout << "\nGroup Velocity Direction = " << radtodeg(normal_slope) << endl;
	cout << "\nGroup Velocity Magnitude = " << 1 / group_velocity << endl;




}



void testcpp()
{


	//Matrix<double, 2, Dynamic, RowMajor> s1(2, 37);
	Master();
	//vector <double> XR1, YR1;
	calculateAnisotropic4c(theRay, 250.0);
	//RotateSlownessTop(XR1, YR1, degtorad(150.0));

	/*
	vector<double> C (5), D;
	C = { 1, 2, 3, 4, 5 };
	for (int i = 0; i < 17;i=i+2)
	D.push_back(i);
	cout <<endl<< C.size()<<endl;
	cout << D.size()<<endl;

	C.erase(C.begin() + 2);
	D.erase(D.begin() + 2);

	cout << C.size()<<endl;
	cout << D.size()<<endl;
	C.erase(C.begin() + 2);
	cout << C.size() << endl;
	*/

	//for (int i = 0; i < XR1.size(); ++i)
	//cout << "Displaying: " << XR1.at(i) << "  " << YR1.at(i) << "\t" << i + 1 << endl;

}



ray plotRayWeld2(ray r, float eta, float T)
{
	ofstream myfile;
	myfile.open("Data.txt");
	cout << "\nActual tracing(50): " << eta << " " << T;

	vector <vertex_type> rP;
	cout << "\n********************************************************************************\n";
	cout << endl << "\n\t\t\t\tINSIDE WELD\n" << endl;
	double current_crystal_orientation, next_crystal_orientation, boundary_orientation;
	cout << "\nActual Checking - Before entry into weld: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction);
	//r.current_position = r.init_position;
	//r.global_direction = r.init_global_direction;
	//r.phase_velocity = r.init_phase_velocity;

	//r.group_vel_global_direction = 0; //Why was this required?
	vertex_type v;

	//glBegin(GL_LINE_STRIP);
	//glVertex3f(r.current_position.x, r.current_position.y, r.current_position.z); //Mark only first point
	rP.clear();
	rP.push_back(r.current_position);
	cout << "\nRay Parameters inside Weld (Initial):\n";
	cout << endl << "y: " << r.current_position.y << "\tz: " << r.current_position.z << "\tDir: " << radtodeg(r.global_direction) << "\tCry-Ang: " << radtodeg(getCrystalOrientation(r.current_position.y, r.current_position.z, T, eta)) << "\tVel: " << r.phase_velocity;
	ray_heading_temp = r.global_direction;
	//ray_heading_temp_group = r.group_vel_global_direction; //groupvelglodir was not initialised thus far. So used glodir because ray is JUST entering anisotropic region 




	ray_heading_temp_group = r.global_direction;
	switch (weld_entry)
	{
	case 0: boundary_orientation = 0; break;
	case 1: boundary_orientation = degtorad(90 - A / 2); break;
	case 2: boundary_orientation = degtorad(90 + A / 2); break;
	}
	cout << "\n\nWELD ENTRY: ";
	if (weld_entry == 1)
		cout << "Right Side\n";
	else if (weld_entry == 2)
		cout << "Left Side\n";
	else if (weld_entry == 0)
		cout << "Middle\n";

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
	incident_angle = r.global_direction - boundary_orientation + degtorad(90);    // Measured against the normal to the interface
	incident_velocity = r.phase_velocity;

	cout << endl << "\nROTATED RAY: " << endl << "Dir: " << radtodeg(rotated_ray.global_direction) << "\tVel: " << rotated_ray.phase_velocity;

	if ((abs(rotated_ray.global_direction - degtorad(180)) < degtorad(0.03)) || (abs(rotated_ray.global_direction) < degtorad(0.03)))
	{
		//Do nothing. Ray is almost parallel to boundary => Bypasses it. 
	}
	else
	{
		//cout<<'\n'<<"CRYSTAL ORIENTATION INSIDE WELD\n";
		calculateAnisotropic2(rotated_ray, radtodeg(getCrystalOrientation(rotated_ray.current_position.y, rotated_ray.current_position.z, T, eta)));
		//Above calculateAnisotropic2 call is the 'first' call meant for entering the weld. New phase velocity and group velocity are obtained.
	}
	//cout<<"\nChecking actual after just entry:"<<rotated_ray.current_position.y<<" "<<rotated_ray.current_position.z<<" "<<radtodeg(rotated_ray.global_direction)<<"\n";
	//cout<<endl<<"After "<<endl<<radtodeg(rotated_ray.global_direction)<<" "<<rotated_ray.phase_velocity<<" "<<radtodeg(rotated_ray.group_vel_global_direction);
	r = rotated_ray; //Copy all attributes: position, phasevel, groupvel, directions
	cout << "\nActual Checking - After 1st change: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction);
	//Then correct only the directions
	r.global_direction = rotated_ray.global_direction + boundary_orientation;
	r.group_vel_global_direction = rotated_ray.group_vel_global_direction + boundary_orientation;

	cout << "\nActual Checking - Before ray heading: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction) << endl;

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
	cout << "\nChecking Transmission Coefficient after entering weld:" << transmission_coefficient << endl;
	cout << "\nActual Checking - After ray heading: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction) << endl;
	cout << "\nChecking angles just after entering weld function: " << radtodeg(incident_angle) << "\t " << radtodeg(refracted_angle) << endl;



	temp_first = 0;
	tot_time_weld = 0;
	refl_check = 0;
	cout << endl << "\t\t\t\tSTARTING WELD LOOP\n";
	cout << endl << "y: " << r.current_position.y << "\tz: " << r.current_position.z << "\tDir: " << radtodeg(r.global_direction) << "\tGVel: " << radtodeg(r.group_vel_global_direction) << "\tPVel: " << r.phase_velocity << " ";
	cout << endl << "\nInside Weld?" << isInsideWeld(rP.back()) << endl;
	temp_count = 0;

	//Calling master() function to initialize the slowness values for all 3601 angles
	Master();


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

		//Changing incident angle and velocity after entering weld
		incident_velocity = r.phase_velocity;
		phase_diff = r.global_direction - phase_diff;
		//	incident_angle = r.global_direction;

		// Commenting all the lines which write the data to external text file. Please uncomment if necessary.


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
		incident_angle = abs(r.global_direction - boundary_orientation + degtorad(90));

		if ((abs(rotated_ray.global_direction - degtorad(180)) < degtorad(0.03)) || (abs(rotated_ray.global_direction) < degtorad(0.03)))
		{
			//Do nothing. Ray is almost parallel to boundary => Bypasses it. 
		}
		else
		{
			calculateAnisotropic4c(rotated_ray, radtodeg(next_crystal_orientation - boundary_orientation)); // Check why we can't use anisotropic2 here
		}

		//cout << "\nCHECKING:" << rotated_ray.phase_velocity << "  " << r.phase_velocity;
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

		//DEFECT
		if (isAtDefect(r.current_position) && temp_first == 0 && defect_enable)
		{
			if (num_defect)
			{
				r.type = abs(r.type - 1);
				num_defect = 0;
			}
			boundary_orientation = defect_angle;
			ray_heading_temp = r.global_direction;
			ray_heading_temp_group = r.group_vel_global_direction;
			rotated_ray = r;
			rotated_ray.global_direction = r.global_direction - boundary_orientation;
			//rotated_ray.group_vel_global_direction = r.group_vel_global_direction - boundary_orientation;

			//In the rotated ray, only the phase velocity direction differs. (is rotated)  

			if (rotated_ray.global_direction<0){
				rotated_ray.global_direction = rotated_ray.global_direction + degtorad(360);
				rotflag_temp = 1;
			}
			else if ((rotated_ray.global_direction < degtorad(180)) && (rotated_ray.global_direction>0)){
				rotated_ray.global_direction = rotated_ray.global_direction + degtorad(180);
				rotflag_temp = 2;
			}
			else if (rotated_ray.global_direction>degtorad(180))
			{
				rotflag_temp = 3;
				//Do nothing
			}
			cout << endl << "Defect Rotated Ray" << endl << radtodeg(rotated_ray.global_direction) << " " << radtodeg(boundary_orientation) << " " << rotated_ray.phase_velocity;
			if ((abs(rotated_ray.global_direction - degtorad(180)) < degtorad(0.03)) || (abs(rotated_ray.global_direction) < degtorad(0.03))){
				//Do nothing. Ray is almost parallel to boundary => Bypasses it. 
			}
			else
			{
				cout << "*** " << radtodeg(getCrystalOrientation(rotated_ray.current_position.y, rotated_ray.current_position.z, T, eta)) - radtodeg(boundary_orientation) << endl;
				calculateAnisotropic2refl(rotated_ray, radtodeg(getCrystalOrientation(rotated_ray.current_position.y, rotated_ray.current_position.z, T, eta)) - radtodeg(boundary_orientation));
				//Above calculateAnisotropic2 call is the 'first' call meant for entering the weld. New phase velocity and group velocity are obtained.
			}
			cout << endl << "Defect After " << endl << radtodeg(rotated_ray.global_direction) << " " << rotated_ray.phase_velocity << " " << radtodeg(rotated_ray.group_vel_global_direction);
			r = rotated_ray; //Copy all attributes: position, phasevel, groupvel, directions
			//Then correct only the directions
			r.global_direction = rotated_ray.global_direction + boundary_orientation;
			r.group_vel_global_direction = rotated_ray.group_vel_global_direction + boundary_orientation;
			//for(int i = 0;i<2; i++)
			//{
			//	if(abs(r.global_direction - ray_heading_temp)>degtorad(150)){
			switch (rotflag_temp)
			{
			case 1:
				r.global_direction = r.global_direction - degtorad(360);
				r.group_vel_global_direction = r.group_vel_global_direction - degtorad(360);
				break;
			case 2:
				r.global_direction = r.global_direction - degtorad(180);
				r.group_vel_global_direction = r.group_vel_global_direction - degtorad(180);
				break;
			case 3:
				break;
			}
			phase_diff = r.global_direction - phase_diff;
			//if(r.global_direction > degtorad(180))
			//{
			//		r.global_direction = r.global_direction - degtorad(180);
			//}
			////	}
			////	if(abs(r.group_vel_global_direction - ray_heading_temp_group)>degtorad(50)){
			//if(r.group_vel_global_direction > degtorad(180))
			//{
			//		r.group_vel_global_direction = r.group_vel_global_direction - degtorad(180);
			//}

			/*if(abs((ray_heading_temp - r.global_direction)/2) < degtorad(90))
			{
			r.global_direction = r.global_direction - degtorad(180);
			}
			if(abs((ray_heading_temp_group - r.group_vel_global_direction)/2) < degtorad(90))
			{
			r.global_direction = r.global_direction - degtorad(180);
			}*/
			//	}
			//}	


			/*
			cout<<endl<<"REFL "<<r.current_position.z<<" "<<radtodeg(r.global_direction)<<" ";

			cout<<" "<<radtodeg(r.global_direction)<<" "<<radtodeg(getCrystalOrientation(rotated_ray.current_position.y, rotated_ray.current_position.z,T,eta))<<endl;

			calculateAnisotropic2refl(r,radtodeg(getCrystalOrientation(rotated_ray.current_position.y, rotated_ray.current_position.z),T,eta));

			cout<<"AFTER: "<<radtodeg(r.global_direction)<<" "<<radtodeg(r.group_vel_global_direction)<<endl;*/

			//Angle corrections
			if (radtodeg(r.global_direction) < 0)
				r.global_direction = r.global_direction + degtorad(360);
			if (radtodeg(r.group_vel_global_direction) < 0)
				r.group_vel_global_direction = r.group_vel_global_direction + degtorad(360);

			cout << "AFTER2: " << radtodeg(r.global_direction) << " " << radtodeg(r.group_vel_global_direction) << endl;
			v.x = rP.back().x;
			v.y = rP.back().y + step_size*cos(r.group_vel_global_direction);
			v.z = rP.back().z + step_size*sin(r.group_vel_global_direction);
			r.time = abs(step_size / r.group_velocity);
			tot_time_weld = tot_time_weld + r.time;
			r.current_position = v;
			rP.push_back(v);
			temp_first = 1;
			num_defect = 1;
		}




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
			cout << endl << "REFL " << r.current_position.z << " " << radtodeg(r.global_direction) << " ";
			if (r.current_position.z >= PLATE_HEIGHT)
			{
				r.global_direction = r.global_direction + degtorad(180);
			}
			cout << " " << radtodeg(r.global_direction) << " " << radtodeg(getCrystalOrientation(r.current_position.y, r.current_position.z, T, eta)) << endl;

			calculateAnisotropic4crefl(r, radtodeg(getCrystalOrientation(r.current_position.y, r.current_position.z, T, eta)));

			cout << "AFTER: " << radtodeg(r.global_direction) << " " << radtodeg(r.group_vel_global_direction) << endl;
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
			cout << "\nChecking reflection function: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction) << " " << r.phase_velocity << endl;
			cout << "\nChecking angles in reflection function: " << radtodeg(incident_angle) << "\t " << radtodeg(reflected_angle) << endl;
			//	transmission_coefficient =  ((incident_velocity*cos(incident_angle)) - (refracted_velocity*cos(refracted_angle))) / ((incident_velocity * cos(incident_angle)) + (refracted_velocity*cos(refracted_angle)));
			//	cout << "\nChecking Transmission Coefficient in reflection :" << transmission_coefficient << endl;
			phase_diff = r.global_direction - phase_diff;
			cout << "AFTER2: " << radtodeg(r.global_direction) << " " << radtodeg(r.group_vel_global_direction) << endl;
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

	if (defect_enable)
		cout << endl << "Defect Position: " << defect_y1 << " " << defect_z1 << " " << defect_y2 << " " << defect_z2 << endl;
	cout << endl << "Number of grains encountered = " << temp_count << " " << l << endl;
	cout << "\nActual Checking - After full change: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction);
	cout << endl << "\n\t\t\t\tEND OF LOOP";
	rayPaths.push_back(rP);
	glEnd();
	myfile.close();



	return r;

}

ray plotRayWeld3(ray r, float eta, float T)
{

	cout << "\nEta and T: " << eta << " " << T;

	vector <vertex_type> rP;
	cout << "\n********************************************************************************\n";
	cout << endl << "\n\t\t\t\tINSIDE WELD\n" << endl;

	double current_crystal_orientation, next_crystal_orientation, boundary_orientation;

	cout << "\nActual Checking - Before entry into weld: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction);

	vertex_type v;

	rP.clear();
	rP.push_back(r.current_position);

	ray_heading_temp = r.global_direction;
	ray_heading_temp_group = r.global_direction;

	switch (weld_entry)
	{
	case 0: boundary_orientation = 0; break;
	case 1: boundary_orientation = degtorad(90 - A / 2); break;
	case 2: boundary_orientation = degtorad(90 + A / 2); break;
	}

	cout << "\n\nWELD ENTRY: ";
	if (weld_entry == 1)
		cout << "Right Side\n";
	else if (weld_entry == 2)
		cout << "Left Side\n";
	else if (weld_entry == 0)
		cout << "Middle\n";


	ray rotated_ray = r;

	current_crystal_orientation = (getCrystalOrientation(r.current_position.y, r.current_position.z, T, eta));

	calculateAnisotropic4c(r, 0.0, current_crystal_orientation);

	cout << "\nActual Checking - Before entry into weld: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction);



	rotated_ray.global_direction = r.global_direction - boundary_orientation;
	rotated_ray.group_vel_global_direction = r.group_vel_global_direction - boundary_orientation;

	//In the rotated ray, only the phase velocity direction differs. (is rotated)  

	if (rotated_ray.global_direction<0)
	{
		rotated_ray.global_direction = rotated_ray.global_direction + degtorad(360);
	}
	else if ((rotated_ray.global_direction < degtorad(180)) && (rotated_ray.global_direction>0))
	{
		rotated_ray.global_direction = rotated_ray.global_direction + degtorad(180);
	}

	incident_angle = r.global_direction - boundary_orientation + degtorad(90);    // Measured against the normal to the interface
	incident_velocity = r.phase_velocity;

	cout << endl << "\nROTATED RAY: " << endl << "Dir: " << radtodeg(rotated_ray.global_direction) << "\tVel: " << rotated_ray.phase_velocity;

	if ((abs(rotated_ray.global_direction - degtorad(180)) < degtorad(0.03)) || (abs(rotated_ray.global_direction) < degtorad(0.03)))
	{
		//Do nothing. Ray is almost parallel to boundary => Bypasses it. 
	}
	else
	{
		//cout<<'\n'<<"CRYSTAL ORIENTATION INSIDE WELD\n";
		calculateAnisotropic2(rotated_ray, radtodeg(getCrystalOrientation(rotated_ray.current_position.y, rotated_ray.current_position.z, T, eta)));
		//Above calculateAnisotropic2 call is the 'first' call meant for entering the weld. New phase velocity and group velocity are obtained.
	}

	r = rotated_ray; //Copy all attributes: position, phasevel, groupvel, directions

	cout << "\nActual Checking - After 1st change: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction);
	//Then correct only the directions
	r.global_direction = rotated_ray.global_direction + boundary_orientation;
	r.group_vel_global_direction = rotated_ray.group_vel_global_direction + boundary_orientation;

	cout << "\nActual Checking - Before ray heading: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction) << endl;

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
	cout << "\nChecking Transmission Coefficient after entering weld:" << transmission_coefficient << endl;
	cout << "\nActual Checking - After ray heading: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction) << endl;
	cout << "\nChecking angles just after entering weld function: " << radtodeg(incident_angle) << "\t " << radtodeg(refracted_angle) << endl;



	temp_first = 0;
	tot_time_weld = 0;
	refl_check = 0;
	cout << endl << "\t\t\t\tSTARTING WELD LOOP\n";
	cout << endl << "y: " << r.current_position.y << "\tz: " << r.current_position.z << "\tDir: " << radtodeg(r.global_direction) << "\tGVel: " << radtodeg(r.group_vel_global_direction) << "\tPVel: " << r.phase_velocity << " ";
	cout << endl << "\nInside Weld?" << isInsideWeld(rP.back()) << endl;
	temp_count = 0;

	//Calling master() function to initialize the slowness values for all 3601 angles
	Master();


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

		//Changing incident angle and velocity after entering weld
		incident_velocity = r.phase_velocity;
		phase_diff = r.global_direction - phase_diff;
		//	incident_angle = r.global_direction;

		// Commenting all the lines which write the data to external text file. Please uncomment if necessary.


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

		//cout<<"\nDir: "<<radtodeg(rotated_ray.global_direction);//!!

		// Can I take absolute? Yes you can
		incident_angle = abs(r.global_direction - boundary_orientation + degtorad(90));

		if ((abs(rotated_ray.global_direction - degtorad(180)) < degtorad(0.03)) || (abs(rotated_ray.global_direction) < degtorad(0.03)))
		{
			//Do nothing. Ray is almost parallel to boundary => Bypasses it. 
		}
		else
		{
			calculateAnisotropic4c(rotated_ray, radtodeg(current_crystal_orientation - boundary_orientation), radtodeg(next_crystal_orientation - boundary_orientation)); // Check why we can't use anisotropic2 here
		}

		//	cout << "\nCHECKING-New:" << rotated_ray.phase_velocity << "  " << r.phase_velocity;
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
			cout << endl << "REFL " << r.current_position.z << " " << radtodeg(r.global_direction) << " ";
			if (r.current_position.z >= PLATE_HEIGHT)
			{
				r.global_direction = r.global_direction + degtorad(180);
			}
			cout << " " << radtodeg(r.global_direction) << " " << radtodeg(getCrystalOrientation(r.current_position.y, r.current_position.z, T, eta)) << endl;

			calculateAnisotropic2refl(r, radtodeg(getCrystalOrientation(r.current_position.y, r.current_position.z, T, eta)));

			cout << "AFTER: " << radtodeg(r.global_direction) << " " << radtodeg(r.group_vel_global_direction) << endl;
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
			cout << "\nChecking reflection function: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction) << " " << r.phase_velocity << endl;
			cout << "\nChecking angles in reflection function: " << radtodeg(incident_angle) << "\t " << radtodeg(reflected_angle) << endl;
			//	transmission_coefficient =  ((incident_velocity*cos(incident_angle)) - (refracted_velocity*cos(refracted_angle))) / ((incident_velocity * cos(incident_angle)) + (refracted_velocity*cos(refracted_angle)));
			//	cout << "\nChecking Transmission Coefficient in reflection :" << transmission_coefficient << endl;
			phase_diff = r.global_direction - phase_diff;
			cout << "AFTER2: " << radtodeg(r.global_direction) << " " << radtodeg(r.group_vel_global_direction) << endl;
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

	cout << endl << "Number of grains encountered = " << temp_count << " " << l << endl;
	cout << "\nActual Checking - After full change: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction);
	cout << endl << "\n\t\t\t\tEND OF LOOP";
	rayPaths.push_back(rP);
	glEnd();
	//	myfile.close();



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


ray plotRayisotropic(ray r)
{
	//Find equations of both left and right weld boundaries
	cout << "\n\n\n********************************************************************************\n";
	cout << endl << "\t\tISOTROPIC PROPOGATION STARTED - BEFORE WELD" << endl;
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
	tot_time_iso_1 = 0;

	phase_diff = r.global_direction;

	while (r.current_position.y > -PLATE_WIDTH / 2 && r.current_position.y < PLATE_WIDTH / 2)
	{
		cout << "\nInitial Ray Parameters:" << endl;
		cout << "y: " << r.current_position.y << "\tz: " << r.current_position.z << "\t Direction: " << radtodeg(r.global_direction) << endl;
		if (isInsideWeld(rP.back()))
		{
			cout << endl << "Inside Weld => Anisotropy" << endl;
			break;
		}
		if (r.global_direction > degtorad(180))
		{
			v.y = -r.current_position.z / tan(r.global_direction) + r.current_position.y;
			v.z = 0;
		}
		else
		{
			v.y = (PLATE_HEIGHT - r.current_position.z) / tan(r.global_direction) + r.current_position.y;
			v.z = PLATE_HEIGHT;
		}
		new_direction = degtorad(360) - r.global_direction;

		//cout<<endl<<v.y<<" "<<v.z;
		if (isInsideWeld(v) || ((r.current_position.y*v.y) < 0))
		{
			c_ray = -r.current_position.y*tan(r.global_direction) + r.current_position.z;

			if (r.current_position.y < 0)
				side = -1;
			else
				side = 1;

			if ((side*m1 - tan(r.global_direction)) < (10 ^ -5))
			{
				cout << "Ray and weld almost parallel";
				break;
			}
			r.current_position.y = (c_ray - c) / (side*m1 - tan(r.global_direction));
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
		//	r.time = sqrt((r.current_position.y - v.y)*(r.current_position.y - v.y) + (r.current_position.z - v.z)*(r.current_position.z - v.z)) / phase_velocity;
		tot_time_iso_1 = tot_time_iso_1 + r.time;

		r.current_position = v;
		//	phase_diff = r.global_direction - phase_diff;
		//cout<<endl<<r.current_position.y<<" "<<r.current_position.z<<" "<<r.current_position.x;
		rP.push_back(v);
	}
	rayPaths.push_back(rP);
	cout << "\n\nChecking time before weld: " << total_time << endl;

	cout << endl << "\t\tCOMPLETED ISOTROPIC RAY TRACING - BEFORE WELD" << endl;
	cout << "\n********************************************************************************\n";
	return r;
}
ray plotRayisotropic_exit(ray r)
{
	//Find equations of both left and right weld boundaries
	cout << "\n********************************************************************************";
	cout << endl << "\n\t\tISOTROPIC PROPOGATION STARTED - AFTER WELD\n" << endl;
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
		k++;
		cout << "\ny: " << r.current_position.y << "\tz: " << r.current_position.z << "\tDirection: " << radtodeg(r.global_direction) << endl;

		/*if(isInsideWeld(rP.back()))
		{
		cout<<endl<<"Inside Weld => Anisotropy";
		break;
		}*/
		if (r.global_direction > degtorad(180))
		{
			v.y = -r.current_position.z / tan(r.global_direction) + r.current_position.y;
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

			cout << "\nChecking Angle after weld:" << radtodeg(r.global_direction) << " " << m_1;
			if (radtodeg(r.global_direction) > m_1)
			{
				cout << "\nCheck1";
				v.y = -PLATE_WIDTH / 2;
				v.z = tan(r.global_direction) * (v.y - r.current_position.y);

			}
			else if (radtodeg(r.global_direction) < m_1)
			{
				cout << "\nCheck2";
				v.y = (PLATE_HEIGHT - r.current_position.z) / tan(r.global_direction) + r.current_position.y;
				v.z = PLATE_HEIGHT;
			}

			else if (radtodeg(r.global_direction) == m_1)
			{
				v.y = -PLATE_WIDTH / 2;
				v.z = PLATE_HEIGHT;
			}


		}
		new_direction = degtorad(360) - r.global_direction;
		cout << endl << "Testing2: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction) << " " << isInsideWeld(v);
		//	total_time = total_time + sqrt((r.current_position.y - v.y)*(r.current_position.y - v.y) + (r.current_position.z - v.z)*(r.current_position.z - v.z)) / r.phase_velocity;

		if (isInsideWeld(v) || ((r.current_position.y*v.y) < 0))
		{
			c_ray = -r.current_position.y*tan(r.global_direction) + r.current_position.z;
			if (r.current_position.y < 0)
				side = -1;
			else
				side = 1;
			if ((side*m1 - tan(r.global_direction)) < (10 ^ -5))
			{
				cout << "Ray and weld almost parallel";
				break;
			}
			r.current_position.y = (c_ray - c) / (side*m1 - tan(r.global_direction));
			r.current_position.z = side*m1*r.current_position.y + c;
			rP.push_back(r.current_position);
			if (side == 1)
				weld_entry = 1;
			else
				weld_entry = 2;
			break;
		}
		r.global_direction = new_direction;
		//	r.time = sqrt((r.current_position.y - v.y)*(r.current_position.y - v.y) + (r.current_position.z - v.z)*(r.current_position.z - v.z)) / phase_velocity;
		//	tot_time_iso_2 = tot_time_iso_2 + r.time;
		r.current_position = v;
		phase_diff = r.global_direction - phase_diff;
		total_time = total_time + sqrt((r.current_position.y - v.y)*(r.current_position.y - v.y) + (r.current_position.z - v.z)*(r.current_position.z - v.z)) / r.phase_velocity;
		cout << endl << "TESTING1: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction);
		rP.push_back(v);
		if (r.current_position.z == PLATE_HEIGHT || r.current_position.y == -PLATE_WIDTH / 2)
			break;
	}
	rayPaths.push_back(rP);

	cout << "\nChecking number of times loop executed: " << k << endl;
	r.global_direction = degtorad(360) - r.global_direction;  // This is becuase the ray angle is measured anticlockwise from horizontal axis
	// Amp calc
	amplitude = cos(r.global_direction - degtorad(90));
	amplitude = amplitude * transmission_coefficient;
	cout << "\nTime taken to travel inside weld: " << tot_time_weld * 1000 << " Micro seconds" << endl;
	cout << "\nTotal time taken to travel completely: " << (tot_time_weld + tot_time_iso_1 + tot_time_iso_2) * 1000 << " Micro seconds" << endl;
	cout << "\nTotal time taken to travel completely2: " << total_time * 1000 << " Micro seconds" << endl;

	cout << "\nPhase Difference: " << radtodeg(phase_diff) << "\nAmplitude: " << amplitude << endl;
	cout << "\nFinal: " << r.current_position.y << " " << r.current_position.z << " " << radtodeg(r.global_direction);
	cout << endl << "\n\tCOMPLETED ISOTROPIC RAY TRACING - AFTER WELD" << endl;

	return r;
}




























#endif