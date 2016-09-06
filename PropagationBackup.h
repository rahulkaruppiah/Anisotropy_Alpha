#ifndef PROPAGATION_H
#define PROPAGATION_H

#include "Eigen/Dense"
#include <iostream>
/*
#include <algorithm>
#include <vector>
#include <math.h>
*/ 

using namespace Eigen;
//SanjeevaReddy.pdf
MatrixXd rotateElasticConstants(float phi){
	
	//phi in radians
	// Coordinate system followed: Rokhlin et al. y-z plane on the screen, x coming out towards you.
	//Note that positive rotation about x means you stare down the positive x axis and then rotate clockwise.
	//This function rotates elastic constants by an angle phi about X (crystal) axis only.

	//phi = degtorad(phi);

	MatrixXd a(3,3);
	a(0,0) = 1;
	a(0,1) = 0;
	a(0,2) = 0;
	a(1,0) = 0;
	a(1,1) = cos(phi);
	a(1,2) = sin(phi);
	a(2,0) = 0;
	a(2,1) = -sin(phi);
	a(2,2) = cos(phi);

	MatrixXd M(6,6);

	// Initialize M

	// Top left quarter.
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			M(i,j) = a(i,j)*a(i,j);
		}
	}

	// Bottom left quarter
	for (int i = 3; i < 6; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			M(i,j) = 1;
			for (int k = 0; k < 3; ++k)
			{
				if(k!=(i-3))	M(i,j)*=a(k,j);
			}
		}
	}

	// Top right quarter
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 3; j < 6; ++j)
		{
			M(i,j) = 2;
			for (int k = 0; k < 3; ++k)
			{
				if(k!=(j-3))	M(i,j)*=a(i,k);
			}
		}
	}

	// Bottom right quarter
	M(3,3) = a(1,1)*a(2,2) + a(2,1)*a(1,2);
	M(3,4) = a(1,0)*a(2,2) + a(2,0)*a(1,2);
	M(3,5) = a(1,0)*a(2,1) + a(2,0)*a(1,1);

	M(4,3) = a(0,1)*a(2,2) + a(2,1)*a(0,2);
	M(4,4) = a(0,0)*a(2,2) + a(2,0)*a(0,2);
	M(4,5) = a(0,0)*a(2,1) + a(2,0)*a(0,1);

	M(5,3) = a(0,1)*a(1,2) + a(1,1)*a(0,2);
	M(5,4) = a(0,0)*a(1,2) + a(1,0)*a(0,2);
	M(5,5) = a(0,0)*a(1,1) + a(1,0)*a(0,1);

	//Phew. Generated M. Now, to do the rotation. Easy part, actually, thanks to awesome matrix libraries.
	return M*c*M.transpose();
	
	//Old attempt
	//From SanjeevaReddy.pdf. The coordinates are screwed. He uses a different system from what I'm using now.

	// C[1][1] = c[1][1]*pow(cos(phi),4) + 2*c[1][3]*pow(cos(phi),2)*pow(sin(phi),2) + c[3][3]*sin(phi) + c[4][4]*pow(sin(2*phi),2);
	// C[1][2] = c[1][2]*pow(cos(phi),2) + c[1][3]*pow(sin(phi),2);
	// C[1][3] = c[1][1]*pow(cos(phi),2)*pow(sin(phi),2) + c[1][3]*(pow(cos(phi),4) + pow(sin(phi),4)) + c[3][3]*pow(cos(phi),2)*pow(sin(phi),2) -c[4][4]*pow(sin(2*phi),2);
	// C[1][5] = (c[1][1]*pow(cos(phi),2) - c[3][3]*pow(sin(phi),2) - c[1][3]*cos(2*phi))*0.5*sin(2*phi) - c[4][4]*cos(2*phi)*sin(2*phi);
	// C[2][3] = c[1][2]*pow(sin(phi),2) + c[1][3]*pow(cos(phi),2);
	// C[2][5] = (c[1][2] - c[1][3])*0.5*sin(2*phi);
	// C[3][3] = c[1][1]*pow(sin(phi),4) + 2*c[1][3]*pow(cos(phi),2)*pow(sin(phi),2) + c[3][3]*pow(sin(phi),4) + c[4][4]*pow(sin(2*phi),2);
	// C[3][5] = (c[1][1]*pow(sin(phi),2) + c[3][3]*pow(cos(phi),2) + c[1][3]*cos(2*phi))*0.5*sin(2*phi) + c[4][4]*cos(2*phi)*sin(2*phi);
	// C[4][4] = c[4][4]*pow(cos(phi),2) + c[6][6]*pow(sin(phi),2);
	// C[4][6] = (c[6][6] - c[4][4])*0.5*sin(2*phi);
	// C[5][5] = (c[1][1] - 2*c[1][3] + c[3][3])*pow(0.5*sin(2*phi),2) - c[4][4]*pow(cos(2*phi),2);
	// C[6][6] = c[4][4]*pow(sin(phi),2) + c[6][6]*pow(cos(phi),2);

}

void calculateAnisotropic2(ray &r, float current_grain_orientation = 0.0f){ 
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
	//cout<<"Ray direction: "<<radtodeg(r.global_direction)<<"   ";
    MatrixXd CLower(6,6);
    
    //Rotate the elastic constants in the medium.

    CLower = rotateElasticConstants(-current_grain_orientation); //NEGATIVE HERE!!!!
	//cout<<" Phase velocity: "<<r.phase_velocity<<endl;
	double m2o = n2/r.phase_velocity;
	double m3o = n3/r.phase_velocity;
	//cout<<"m3o = "<<m3o<<"     ";
	
	//cout<<"m2o = "<<m2o<<"m30 = "<<m3o<<" "<<current_grain_orientation<<endl;

	double m2 = m2o;//m2 component for old ray and new ray.
	double m3;
	//if(prev_grain_orientation!=current_grain_orientation){
	//	first_check = 0;
	//}
	//if(!first_check){
//    To find the m3 component now. Let m3 = B. Lower medium, coz we're trying to find the refracted ray and not the reflected one.
		MatrixXd a(3,3), b(3,3), d(3,3);
		a << CLower(5,5)*m2o*m2o - DENSITY, CLower(5,1)*m2o*m2o, CLower(5,3)*m2o*m2o,
			CLower(1,5)*m2o*m2o, CLower(1,1)*m2o*m2o - DENSITY, CLower(1,3)*m2o*m2o,
		CLower(3,5)*m2o*m2o, CLower(3,1)*m2o*m2o, CLower(3,3)*m2o*m2o - DENSITY;
	    
		b << (CLower(4,5) + CLower(5,4))*m2o, (CLower(4,1)+CLower(5,3))*m2o, (CLower(4,3)+CLower(5,2))*m2o,
			(CLower(3,5) + CLower(1,4))*m2o, (CLower(3,1) + CLower(1,3))*m2o, (CLower(3,3)+CLower(1,2))*m2o,
		(CLower(2,5) + CLower(3,4))*m2o, (CLower(2,1) + CLower(3,3))*m2o, (CLower(2,3) + CLower(3,2))*m2o;
	    
		d << CLower(4,4), CLower(4,3), CLower(4,2),
			CLower(3,4), CLower(3,3), CLower(3,2),
		CLower(2,4), CLower(2,3), CLower(2,2);
	    
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
	    
	 /*   printf("\nThe Co-efficients of the sextic equation are as follows:\n");
		for (int i=0; i<sexticCoeff.size(); ++i) {
			cout<<sexticCoeff[i]<<endl;
		}*/
	    
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
		for(int i = 0; i< 6; i++){
			if(!isOpposite(real(rootsOfSextic[i]), m3o)){ //Roots should have same sign as m3o. //CHECKING ALL SIX ROOTS
				slownessSolutions.push_back(real(rootsOfSextic[i]));
			}
		}

		sort(slownessSolutions.begin(),slownessSolutions.end(),std::greater<double>()); //Hence slownessSolutions contains the three feasible slowness values (sorted) for lower medium 
		m3 = slownessSolutions.at(0);
		first_check = 1;
		//cout<<"abs";
	//}else{
	//	m3 = m3o;
	//	//cout<<"here";
	//}

	//cout<<"m3 = "<<m3<<endl;
	
	//At this point, m2 and m3 in the new medium are obtained properly!!
	double slownessMagnitude = sqrt(m2*m2 + m3*m3);
	double newn2 = m2/slownessMagnitude;
    double newn3 = m3/slownessMagnitude;
    
	//GET NEW PHASE VELOCITY DIRECTION
	r.global_direction = getAbsoluteAngle(newn2,newn3);
	//cout<<"New ray direction (phase): "<<radtodeg(r.global_direction);
	double newVel = 1/slownessMagnitude;
//            if(1/slownessMagnitude > newLongPhaseVel){
//                newLongPhaseVel = 1/slownessMagnitude;
//                findLongPhaseVelIndex = i;
            //            }

	//GET PHASE VELOCITY MAGNITUDE USING KNOWN PHASE VELOCITY DIRECTION IN THE NEW MEDIUM BY SOLVING CHRISTOFFEL'S EQUATION
            MatrixXd newChristoffelEquation(3,3);
            newChristoffelEquation(0,0) = CLower(5,5)*newn2*newn2 + CLower(4,4)*newn3*newn3;
            newChristoffelEquation(0,1) = (CLower(3,5)+CLower(1,4))*newn2*newn3;
            newChristoffelEquation(0,2) = CLower(3,5)*newn2*newn2 + CLower(2,4)*newn3*newn3;
            newChristoffelEquation(1,0) = (CLower(1,4)+CLower(3,5))*newn2*newn3;
            newChristoffelEquation(1,1) = CLower(1,1)*newn2*newn2 + CLower(3,3)*newn3*newn3;
            newChristoffelEquation(1,2) = (CLower(1,2) + CLower(3,3))*newn2*newn3;
            newChristoffelEquation(2,0) = CLower(3,5)*newn2*newn2 + CLower(2,4)*newn3*newn3;
            newChristoffelEquation(2,1) = (CLower(3,3) + CLower(1,2))*newn2*newn3;
            newChristoffelEquation(2,2) = CLower(3,3)*newn2*newn2 + CLower(2,2)*newn3*newn3;
			newChristoffelEquation = newChristoffelEquation/DENSITY;
            EigenSolver<MatrixXd> es;
            es.compute(newChristoffelEquation,1);

	//cout<<christoffelEquation<<endl;
	//cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
	//cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
//	printf("The velocities can be calculated as follows:\n");
//  Choose the current ray type.
	//cout<<es.eigenvalues()<<endl;
    int phaseVelIndex = 0;
    double phaseVel = 0;
    vector <double> velocities;
    for(int i = 0 ; i < 3; i++){
		if(imag(es.eigenvalues()[i])!=0){
//			printf("%d - imaginary eigenvalue.\n",i);
		}
		else if(real(es.eigenvalues()[i])<0){
//			printf("%d - negative. Therefore, V is imaginary.\n", i);
		}
		else{
//			printf("%d - velocity = %f\n", i, sqrt(real(es.eigenvalues()[i])));
            velocities.push_back(sqrt(real(es.eigenvalues()[i])));
		}
	}
	//cout<<"velocities: "<<velocities[0]<<" "<<velocities[1]<<" "<<velocities[2]<<" "<<endl;
	sort(velocities.begin(),velocities.end(), std::greater<double>());  //Sort in reverse order so that 0 is long velocity, and so on.
//    Associating index with velocities. Good thing to do?
    while (velocities[r.type]!=sqrt(real(es.eigenvalues()[phaseVelIndex]))) {
        phaseVelIndex++;
    }
    
	phaseVel = sqrt(real(es.eigenvalues()[phaseVelIndex])); 
	//cout<<"   New Phase Velocity: "<<phaseVel<<endl;
	//phaseVel = sqrt(real(es.eigenvalues()[0])); 
	//The while loop and the above line essentially pick out ONE of the three phase velocities obtained from the equation (eigenvalues). 
	//Which phase velocity is picked depends on the ray type. Longitudinal -> Highest velocity; SV -> Second; SH -> Slowest;
	
	r.phase_velocity = phaseVel; //important line: phase velocity of the ray in the new medium is set here.    
    
	//GET GROUP VELOCITY MAGNITUDE AND DIRECTION
	double groupVel[4] = {0};
    double m[4] = {0,0,m2,m3};
    double P[4] = {0};
    for (int i = 1; i<4; ++i) {
        P[i] = real(es.eigenvectors()(i-1,phaseVelIndex));
        //cout << "P - " << i << " "<< P[i] << endl;
    }
    for (int i = 1; i<4; ++i) {
        groupVel[i] = 0;
        for (int l = 1; l < 4; ++l) {
            for (int j = 1; j < 4; ++j) {
                for (int k = 1; k< 4; k++) {
                    //cout<<i<<" "<<1/DENSITY*VoigtC(C, i, j, k, l)*m[l]*P[j]*P[k]<<endl;
					groupVel[i] += 1/DENSITY*VoigtC(CLower, i, j, k, l)*m[l]*P[j]*P[k];
                }
            }
        }
        //cout << "Group vel - " << i << " "<< groupVel[i] << endl;
    }
    r.group_vel_global_direction = getAbsoluteAngle(groupVel[2], groupVel[3]);
	//cout<<"Group Vel Direction =  "<<180/PI * r.group_vel_global_direction<<endl;
    r.group_velocity = sqrt(groupVel[2]*groupVel[2] + groupVel[3]*groupVel[3]);
	prev_grain_orientation = current_grain_orientation;
    //cout << "Group velocity magnitude : " << sqrt(groupVel[2]*groupVel[2] + groupVel[3]*groupVel[3]) << endl;
    //cout << "Direction of Group vel is : " << radtodeg(r.group_vel_global_direction) << endl;



//	for (int i=0; i<6; ++i) {
//        if(imag(rootsOfSextic[i]))
//            continue;   // For now. Ignore Imaginary solutions to the slowness vectors.
//        double m3 = real(rootsOfSextic[i]);
//        double slownessMagnitude = sqrt(m2*m2 + m3*m3);
////        checkMatrix(0,0) = CLower(5,5)*m2*m2 + CLower(4,4)*m3*m3;
////        checkMatrix(0,1) = CLower(3,5)+CLower(1,4)*m2*m3;
////        checkMatrix(0,2) = CLower(3,5)*m2*m2 + CLower(2,4)*m3*m3;
////        checkMatrix(1,0) = CLower(1,4)+CLower(3,5)*m2*m3;
////        checkMatrix(1,1) = CLower(1,1)*m2*m2 + CLower(3,3)*m3*m3;
////        checkMatrix(1,2) = CLower(1,2) + CLower(3,3)*m2*m3;
////        checkMatrix(2,0) = CLower(3,5)*m2*m2 + CLower(2,4)*m3*m3;
////        checkMatrix(2,1) = CLower(3,3) + CLower(1,2)*m2*m3;
////        checkMatrix(2,2) = CLower(3,3)*m2*m2 + CLower(2,2)*m3*m3;
////        checkMatrix = (1/DENSITY)*checkMatrix;
////        checkMatrix = checkMatrix - Matrix3d::Identity();
////        double newn2, newn3;
////        newn2 = m2/slownessMagnitude;
////        newn3 = m3/slownessMagnitude;
////        cout << "n2 = " << newn2 << " n3 = "<< newn3<< endl;
////        cout << "Quick check on direction cosines : " << sqrt(newn2*newn2 + newn3*newn3) << endl;
////        cout << checkMatrix.determinant() << endl;
//        if(!isOpposite(m3, m3o)){   //Do Transmitted waves for now. TODO: Need to do reflected waves.
////            Same direction of 3-component of slowness. TODO: Validate the assumption that for the physically valid solutions, slowness velocity 3-component should be along the same direction along with the 3-component of the group velocity direction.
////            cout << i << endl << "m2 = " <<m2 << "\t m3 = " << m3 << endl;
////            cout << " The magnitude of the velocity = " << 1/slownessMagnitude << endl;
//            typeOfRay = 0; 
//			cout<<"Slowness "<<m2<<" "<<m3<<" "<<slownessMagnitude<<endl;
//            while(m3 != slownessSolutions[typeOfRay])
//                typeOfRay++;
//            double newn2, newn3;
//            newn2 = m2/slownessMagnitude;
//            newn3 = m3/slownessMagnitude;
//            double newVel = 1/slownessMagnitude;
////            if(1/slownessMagnitude > newLongPhaseVel){
////                newLongPhaseVel = 1/slownessMagnitude;
////                findLongPhaseVelIndex = i;
//            //            }
//            MatrixXd newChristoffelEquation(3,3);
//            newChristoffelEquation(0,0) = CLower(5,5)*newn2*newn2 + CLower(4,4)*newn3*newn3;
//            newChristoffelEquation(0,1) = (CLower(3,5)+CLower(1,4))*newn2*newn3;
//            newChristoffelEquation(0,2) = CLower(3,5)*newn2*newn2 + CLower(2,4)*newn3*newn3;
//            newChristoffelEquation(1,0) = (CLower(1,4)+CLower(3,5))*newn2*newn3;
//            newChristoffelEquation(1,1) = CLower(1,1)*newn2*newn2 + CLower(3,3)*newn3*newn3;
//            newChristoffelEquation(1,2) = (CLower(1,2) + CLower(3,3))*newn2*newn3;
//            newChristoffelEquation(2,0) = CLower(3,5)*newn2*newn2 + CLower(2,4)*newn3*newn3;
//            newChristoffelEquation(2,1) = (CLower(3,3) + CLower(1,2))*newn2*newn3;
//            newChristoffelEquation(2,2) = CLower(3,3)*newn2*newn2 + CLower(2,2)*newn3*newn3;
//            EigenSolver<MatrixXd> newes;
//            newes.compute(newChristoffelEquation,1);
//            newPhaseVel = 0;
//            vector <double> velocities;
//            for(int i = 0 ; i < 3; i++){
//                if(imag(newes.eigenvalues()[i])!=0){
//                    //			printf("%d - imaginary eigenvalue.\n",i);
////                    Don't do anything
//                }
//                else if(real(newes.eigenvalues()[i])<0){
//                    //			printf("%d - negative. Therefore, V is imaginary.\n", i);
////                    Don't do anything     
//                }
//                else{
//                    //			printf("%d - velocity = %f\n", i, sqrt(real(newes.eigenvalues()[i])/DENSITY));
//                    velocities.push_back(sqrt(real(newes.eigenvalues()[i])/DENSITY));
//                }
//            }
//            sort(velocities.begin(),velocities.end(), std::greater<double>());  //Sort in reverse - longVel is first, and so on. As usual.
//            int phaseVelIndex = 0;
//            for (int i = 0; i < velocities.size(); ++i) {
//                phaseVelIndex = 0;
//                while (velocities[i]!=sqrt(real(newes.eigenvalues()[phaseVelIndex])/DENSITY)) {
//                    phaseVelIndex++;
//                }
//                double groupVel[4] = {0};
//                double m[4] = {0,0,m2,m3};
//                double P[4] = {0};
//                for (int j = 1; j<4; ++j) {
//                    P[j] = real(newes.eigenvectors()(j-1,phaseVelIndex));
//                    //        cout << "P - " << i << P[i] << endl;
//                }
//                for (int q = 1; q<4; ++q) {
//                    groupVel[q] = 0;
//                    for (int l = 1; l < 4; ++l) {
//                        for (int j = 1; j < 4; ++j) {
//                            for (int k = 1; k< 4; k++) {
//                                groupVel[q] += 1/DENSITY*VoigtC(CLower, q, j, k, l)*m[l]*P[j]*P[k];
//                            }
//                        }
//                    }
//                    //        cout << "Group vel - " << i << " "<< groupVel[i] << endl;
//                }
//                //Check if it is of current wave type. Else, create new wave.
//                if(r.type == i){
//                    r.global_direction = getAbsoluteAngle(newn2, newn3);    //Must not do this. Must create new wave here itself for other modes.
//                    r.group_vel_global_direction = getAbsoluteAngle(groupVel[2], groupVel[3]);
//                    r.group_velocity = sqrt(groupVel[2]*groupVel[2] + groupVel[3]*groupVel[3]);
//					//cout<<global_count_temp<<" * "<<r.global_direction<<endl;
//					//global_count_temp++;
//                }
//      //          else{
//      //              //Create new wave and add it to the rayList.
//      //              if( i == typeOfRay){
//      //                  ray newRay;
//      //                  newRay.init_position = r.current_position;
//      //                  newRay.init_global_direction = getAbsoluteAngle(newn2, newn3);
//      //                  newRay.global_direction = newRay.init_global_direction;
//      //                  newRay.type = i;
//      //                  rayList.push_back(newRay);
//						////cout<<global_count_temp<<endl;
//						////global_count_temp++;
//      //              }
//      //          }
//            }
//
//        }
//    }
//	rayList.push_back(r);
////    cout << "Max Vel = " << newLongPhaseVel << " Index = " << findLongPhaseVelIndex << endl;
////    double m3 = real(rootsOfSextic[findLongPhaseVelIndex]);
//	
////    cout << "The new direction of the phase velocity is : " << radtodeg(getAbsoluteAngle(newn2, newn3)) << endl;
//    
////    Should find displacement vectors along this direction.
//    //	cout << "The eigenvalues of A' are:" << endl << newes.eigenvalues() << endl;
////	cout << "The matrix of eigenvectors, V, is:" << endl << newes.eigenvectors() << endl << endl;
////	printf("The NEW velocities can be calculated as follows:\n");
//        r.phase_velocity = newPhaseVel;
////    To find group velocity direction
//    //    cout << "Group velocity magnitude : " << sqrt(groupVel[2]*groupVel[2] + groupVel[3]*groupVel[3]) << endl;
////    cout << "Direction of Group vel is : " << radtodeg(getAbsoluteAngle(groupVel[2], groupVel[3])) << endl;
//
//
////	MatrixXd christoffelEquation(3,3);
////	christoffelEquation(0,0) = C(5,5)*n2*n2 + C(4,4)*n3*n3;
////	christoffelEquation(0,1) = (C(3,5)+C(1,4))*n2*n3;
////	christoffelEquation(0,2) = C(3,5)*n2*n2 + C(2,4)*n3*n3;
////	christoffelEquation(1,0) = (C(1,4)+C(3,5))*n2*n3;
////	christoffelEquation(1,1) = C(1,1)*n2*n2 + C(3,3)*n3*n3;
////	christoffelEquation(1,2) = (C(1,2) + C(3,3))*n2*n3;
////	christoffelEquation(2,0) = C(3,5)*n2*n2 + C(2,4)*n3*n3;
////	christoffelEquation(2,1) = (C(3,3) + C(1,2))*n2*n3;
////	christoffelEquation(2,2) = C(3,3)*n2*n2 + C(2,2)*n3*n3;
////    christoffelEquation = christoffelEquation/DENSITY;
////	EigenSolver<MatrixXd> es;
////	es.compute(christoffelEquation,1);
////	//cout<<christoffelEquation<<endl;
////	//cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
////	//cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
//////	printf("The velocities can be calculated as follows:\n");
//////  Choose the current ray type.
////	cout<<es.eigenvalues()<<endl;
////    int phaseVelIndex = 0;
////    double phaseVel = 0;
////    vector <double> velocities;
////    for(int i = 0 ; i < 3; i++){
////		if(imag(es.eigenvalues()[i])!=0){
//////			printf("%d - imaginary eigenvalue.\n",i);
////		}
////		else if(real(es.eigenvalues()[i])<0){
//////			printf("%d - negative. Therefore, V is imaginary.\n", i);
////		}
////		else{
//////			printf("%d - velocity = %f\n", i, sqrt(real(es.eigenvalues()[i])));
////            velocities.push_back(sqrt(real(es.eigenvalues()[i])));
////		}
////	}
////	sort(velocities.begin(),velocities.end(), std::greater<double>());  //Sort in reverse order so that 0 is long velocity, and so on.
//////    Associating index with velocities. Good thing to do?
////    while (velocities[r.type]!=sqrt(real(es.eigenvalues()[phaseVelIndex]))) {
////        phaseVelIndex++;
////    }
////    
////	phaseVel = sqrt(real(es.eigenvalues()[phaseVelIndex])); 
////	//phaseVel = sqrt(real(es.eigenvalues()[0])); 
////	//The while loop and the above line essentially pick out ONE of the three phase velocities obtained from the equation (eigenvalues). 
////	//Which phase velocity is picked depends on the ray type. Longitudinal -> Highest velocity; SV -> Second; SH -> Slowest;
////	//cout<<endl<<" ** "<<phaseVel<<endl;
////    r.phase_velocity = phaseVel;
////    m2o = n2/phaseVel;
////    m3o = n3/phaseVel;
////	//cout<<" "<<m2o<<" "<<m3o;
////    double groupVel[4] = {0};
////    double m[4] = {0,0,m2o,m3o};
////    double P[4] = {0};
////    for (int i = 1; i<4; ++i) {
////        P[i] = real(es.eigenvectors()(i-1,phaseVelIndex));
////        //cout << "P - " << i << " "<< P[i] << endl;
////    }
////    for (int i = 1; i<4; ++i) {
////        groupVel[i] = 0;
////        for (int l = 1; l < 4; ++l) {
////            for (int j = 1; j < 4; ++j) {
////                for (int k = 1; k< 4; k++) {
////                    //cout<<i<<" "<<1/DENSITY*VoigtC(C, i, j, k, l)*m[l]*P[j]*P[k]<<endl;
////					groupVel[i] += 1/DENSITY*VoigtC(C, i, j, k, l)*m[l]*P[j]*P[k];
////                }
////            }
////        }
////        //cout << "Group vel - " << i << " "<< groupVel[i] << endl;
////    }
////    r.group_vel_global_direction = getAbsoluteAngle(groupVel[2], groupVel[3]);
////    r.group_velocity = sqrt(groupVel[2]*groupVel[2] + groupVel[3]*groupVel[3]);
////    //cout << "Group velocity magnitude : " << sqrt(groupVel[2]*groupVel[2] + groupVel[3]*groupVel[3]) << endl;
////    //cout << "Direction of Group vel is : " << radtodeg(r.group_vel_global_direction) << endl;

}


void calculateAnisotropic4(ray &r, float upper_grain_orientation = 0.0f, float lower_grain_orientation = 0.0f){ 
//    Note: As usual, whatever is passed is in degrees.
	upper_grain_orientation = degtorad(upper_grain_orientation);
	lower_grain_orientation = degtorad(lower_grain_orientation);

	float n2 = cos(r.global_direction);  
	float n3 = sin(r.global_direction);
	//cout<<"Ray direction: "<<radtodeg(r.global_direction)<<"   ";
    
	MatrixXd CUpper(6,6);
	MatrixXd CLower(6,6);
    
    //Rotate the elastic constants in the medium.
	CUpper = rotateElasticConstants(-upper_grain_orientation); //NEGATIVE HERE!!!!
    CLower = rotateElasticConstants(-lower_grain_orientation); //NEGATIVE HERE!!!!
	
	//UPPER PORTION!!
	MatrixXd christoffelEquation(3,3);
	christoffelEquation(0,0) = CUpper(5,5)*n2*n2 + CUpper(4,4)*n3*n3;
	christoffelEquation(0,1) = (CUpper(3,5)+CUpper(1,4))*n2*n3;
	christoffelEquation(0,2) = CUpper(3,5)*n2*n2 + CUpper(2,4)*n3*n3;
	christoffelEquation(1,0) = (CUpper(1,4)+CUpper(3,5))*n2*n3;
	christoffelEquation(1,1) = CUpper(1,1)*n2*n2 + CUpper(3,3)*n3*n3;
	christoffelEquation(1,2) = (CUpper(1,2) + CUpper(3,3))*n2*n3;
	christoffelEquation(2,0) = CUpper(3,5)*n2*n2 + CUpper(2,4)*n3*n3;
	christoffelEquation(2,1) = (CUpper(3,3) + CUpper(1,2))*n2*n3;
	christoffelEquation(2,2) = CUpper(3,3)*n2*n2 + CUpper(2,2)*n3*n3;
    christoffelEquation = christoffelEquation/DENSITY;
	EigenSolver<MatrixXd> es;
	es.compute(christoffelEquation,1);
//	cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
//	cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
//	printf("The velocities can be calculated as follows:\n");
//	for(int i = 0 ; i < 3; i++){
//		if(imag(es.eigenvalues()[i])!=0){
//			printf("%d - imaginary eigenvalue.\n",i);
//		}
//		else if(real(es.eigenvalues()[i])<0){
//			printf("%d - negative. Therefore, V is imaginary.\n", i);
//		}
//		else{
//			printf("%d - velocity = %f\n", i, sqrt(real(es.eigenvalues()[i])));
//		}
//	}
//    cout << "Sorted velocities:" << endl;
    std::vector < double > eigenv;
    for(int i=0; i<3; i++){
        eigenv.push_back(sqrt(real(es.eigenvalues()[i])));
    }
    sort(eigenv.begin(),eigenv.end());
//    for(int i=0; i<eigenv.size(); i++){
//        printf("%f\n", eigenv[i]);
//    }
    float phaseVel = eigenv[2];    //Assuming this is the max. Safe, since it is sorted.
    double m2o = n2/phaseVel;
    double m3o = n3/phaseVel;

	//m2o and m3o FOUND!!

	//cout<<" m2o = "<<m2o<<"     "<<" m3o = "<<m3o<<" pv = "<<phaseVel<<endl;

	double m2 = m2o;//m2 component for old ray and new ray.
	double m3;
	if(upper_grain_orientation!=lower_grain_orientation){
		first_check = 0;
	}
	if(!first_check){
//    To find the m3 component now. Let m3 = B. Lower medium, coz we're trying to find the refracted ray and not the reflected one.
		MatrixXd a(3,3), b(3,3), d(3,3);
		a << CLower(5,5)*m2o*m2o - DENSITY, CLower(5,1)*m2o*m2o, CLower(5,3)*m2o*m2o,
			CLower(1,5)*m2o*m2o, CLower(1,1)*m2o*m2o - DENSITY, CLower(1,3)*m2o*m2o,
		CLower(3,5)*m2o*m2o, CLower(3,1)*m2o*m2o, CLower(3,3)*m2o*m2o - DENSITY;
	    
		b << (CLower(4,5) + CLower(5,4))*m2o, (CLower(4,1)+CLower(5,3))*m2o, (CLower(4,3)+CLower(5,2))*m2o,
			(CLower(3,5) + CLower(1,4))*m2o, (CLower(3,1) + CLower(1,3))*m2o, (CLower(3,3)+CLower(1,2))*m2o,
		(CLower(2,5) + CLower(3,4))*m2o, (CLower(2,1) + CLower(3,3))*m2o, (CLower(2,3) + CLower(3,2))*m2o;
	    
		d << CLower(4,4), CLower(4,3), CLower(4,2),
			CLower(3,4), CLower(3,3), CLower(3,2),
		CLower(2,4), CLower(2,3), CLower(2,2);
	    
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
	    
	 /*   printf("\nThe Co-efficients of the sextic equation are as follows:\n");
		for (int i=0; i<sexticCoeff.size(); ++i) {
			cout<<sexticCoeff[i]<<endl;
		}*/
	    
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
		for(int i = 0; i< 6; i++){
			if(!isOpposite(real(rootsOfSextic[i]), m3o)){ //Roots should have same sign as m3o. //CHECKING ALL SIX ROOTS
				slownessSolutions.push_back(real(rootsOfSextic[i]));
			}
		}

		sort(slownessSolutions.begin(),slownessSolutions.end(),std::greater<double>()); //Hence slownessSolutions contains the three feasible slowness values (sorted) for lower medium 
		m3 = slownessSolutions.at(0);
		first_check = 1;
		//cout<<"abs";
	}else{
		m3 = m3o;
		//cout<<"here";
	}

	//cout<<"m3 = "<<m3<<endl;
	
	//At this point, m2 and m3 in the new medium are obtained properly!!
	double slownessMagnitude = sqrt(m2*m2 + m3*m3);
	double newn2 = m2/slownessMagnitude;
    double newn3 = m3/slownessMagnitude;
    
	//GET NEW PHASE VELOCITY DIRECTION
	r.global_direction = getAbsoluteAngle(newn2,newn3);
	//cout<<"New ray direction (phase): "<<radtodeg(r.global_direction);
	double newVel = 1/slownessMagnitude;
//            if(1/slownessMagnitude > newLongPhaseVel){
//                newLongPhaseVel = 1/slownessMagnitude;
//                findLongPhaseVelIndex = i;
            //            }

	//GET PHASE VELOCITY MAGNITUDE USING KNOWN PHASE VELOCITY DIRECTION IN THE NEW MEDIUM BY SOLVING CHRISTOFFEL'S EQUATION
            MatrixXd newChristoffelEquation(3,3);
            newChristoffelEquation(0,0) = CLower(5,5)*newn2*newn2 + CLower(4,4)*newn3*newn3;
            newChristoffelEquation(0,1) = (CLower(3,5)+CLower(1,4))*newn2*newn3;
            newChristoffelEquation(0,2) = CLower(3,5)*newn2*newn2 + CLower(2,4)*newn3*newn3;
            newChristoffelEquation(1,0) = (CLower(1,4)+CLower(3,5))*newn2*newn3;
            newChristoffelEquation(1,1) = CLower(1,1)*newn2*newn2 + CLower(3,3)*newn3*newn3;
            newChristoffelEquation(1,2) = (CLower(1,2) + CLower(3,3))*newn2*newn3;
            newChristoffelEquation(2,0) = CLower(3,5)*newn2*newn2 + CLower(2,4)*newn3*newn3;
            newChristoffelEquation(2,1) = (CLower(3,3) + CLower(1,2))*newn2*newn3;
            newChristoffelEquation(2,2) = CLower(3,3)*newn2*newn2 + CLower(2,2)*newn3*newn3;
			newChristoffelEquation = newChristoffelEquation/DENSITY;
            EigenSolver<MatrixXd> newes;
            newes.compute(newChristoffelEquation,1);

	//cout<<christoffelEquation<<endl;
	//cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
	//cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
//	printf("The velocities can be calculated as follows:\n");
//  Choose the current ray type.
	//cout<<es.eigenvalues()<<endl;
    int phaseVelIndex = 0;
    phaseVel = 0.0f;
    vector <double> velocities;
    for(int i = 0 ; i < 3; i++){
		if(imag(newes.eigenvalues()[i])!=0){
//			printf("%d - imaginary eigenvalue.\n",i);
		}
		else if(real(newes.eigenvalues()[i])<0){
//			printf("%d - negative. Therefore, V is imaginary.\n", i);
		}
		else{
//			printf("%d - velocity = %f\n", i, sqrt(real(newes.eigenvalues()[i])));
            velocities.push_back(sqrt(real(newes.eigenvalues()[i])));
		}
	}
	//cout<<"velocities: "<<velocities[0]<<" "<<velocities[1]<<" "<<velocities[2]<<" "<<endl;
	sort(velocities.begin(),velocities.end(), std::greater<double>());  //Sort in reverse order so that 0 is long velocity, and so on.
//    Associating index with velocities. Good thing to do?
    while (velocities[r.type]!=sqrt(real(newes.eigenvalues()[phaseVelIndex]))) {
        phaseVelIndex++;
    }
    
	phaseVel = sqrt(real(newes.eigenvalues()[phaseVelIndex])); 
	//cout<<"   New Phase Velocity: "<<phaseVel<<endl;
	//phaseVel = sqrt(real(newes.eigenvalues()[0])); 
	//The while loop and the above line essentially pick out ONE of the three phase velocities obtained from the equation (eigenvalues). 
	//Which phase velocity is picked depends on the ray type. Longitudinal -> Highest velocity; SV -> Second; SH -> Slowest;
	
	r.phase_velocity = phaseVel; //important line: phase velocity of the ray in the new medium is set here.    
    
	//cout<<" m2 = "<<m2<<" m3 = "<<m3<<" pv = "<<phaseVel<<endl;

	//GET GROUP VELOCITY MAGNITUDE AND DIRECTION
	double groupVel[4] = {0};
    double m[4] = {0,0,m2,m3};
    double P[4] = {0};
    for (int i = 1; i<4; ++i) {
        P[i] = real(newes.eigenvectors()(i-1,phaseVelIndex));
        //cout << "P - " << i << " "<< P[i] << endl;
    }
    for (int i = 1; i<4; ++i) {
        groupVel[i] = 0;
        for (int l = 1; l < 4; ++l) {
            for (int j = 1; j < 4; ++j) {
                for (int k = 1; k< 4; k++) {
                    //cout<<i<<" "<<1/DENSITY*VoigtC(C, i, j, k, l)*m[l]*P[j]*P[k]<<endl;
					groupVel[i] += 1/DENSITY*VoigtC(CLower, i, j, k, l)*m[l]*P[j]*P[k];
                }
            }
        }
        //cout << "Group vel - " << i << " "<< groupVel[i] << endl;
    }
    r.group_vel_global_direction = getAbsoluteAngle(groupVel[2], groupVel[3]);
	//cout<<"Group Vel Direction =  "<<180/PI * r.group_vel_global_direction<<endl;
    r.group_velocity = sqrt(groupVel[2]*groupVel[2] + groupVel[3]*groupVel[3]);
	//prev_grain_orientation = current_grain_orientation;
    //cout << "Group velocity magnitude : " << sqrt(groupVel[2]*groupVel[2] + groupVel[3]*groupVel[3]) << endl;
    //cout << "Direction of Group vel is : " << radtodeg(r.group_vel_global_direction) << endl;

}


void plotSinglePlate(ray &r, float current_grain_orientation = 0.0f, double min_y=-PLATE_WIDTH/2, double max_y=PLATE_WIDTH/2, double min_z=0, double max_z=PLATE_HEIGHT){
    glBegin(GL_LINE_STRIP);
    
	cout<<"!!!Before plate: "<<radtodeg(r.global_direction)<<"   "<<r.phase_velocity<<"   "<<radtodeg(r.group_vel_global_direction)<<endl;
	
	while(r.current_position.y >= min_y && r.current_position.y <= max_y && r.current_position.z >= min_z && r.current_position.z <= max_z){
        switch (r.type) {
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
        glVertex3f(0, r.current_position.y, r.current_position.z);
        calculateAnisotropic2(r, current_grain_orientation);
        r.current_position.y += 1.0*cos(r.group_vel_global_direction);
        r.current_position.z += 1.0*sin(r.group_vel_global_direction);
		
        glVertex3f(0, r.current_position.y, r.current_position.z);
 //       drawPhaseFront(r.current_position.y, r.current_position.z, radtodeg(r.global_direction));
//        printf("Phase vel direction: %f; Group vel direction: %f",radtodeg(r.global_direction),radtodeg(r.group_vel_global_direction));
    }
	cout<<"After plate: "<<radtodeg(r.global_direction)<<"   "<<r.phase_velocity<<"   "<<radtodeg(r.group_vel_global_direction)<<endl;
    glEnd();
}

void plotMultiplePlates(ray &r){
    
	//This function is called for every display() routine. Therefore, we need to configure ray starting conditions each time. 
	r.current_position = r.init_position;
    r.global_direction = r.init_global_direction;
	//r.current_direction = r.init_global_direction;
	r.phase_velocity = r.init_phase_velocity;
	cout<<radtodeg(r.group_vel_global_direction)<<endl;
	first_check = 0;
	
	
    int i = 0, j=0;
	do{
		cout<<"\n\nMOVING TO PLATE "<<i+1<<endl<<endl;
		plotSinglePlate(r,multiplePlateElasticOrientation[i],-PLATE_WIDTH/2,PLATE_WIDTH/2,-i*PLATE_HEIGHT, PLATE_HEIGHT - i*PLATE_HEIGHT);
		++i;
	 }while(i<NUM_EXTRA_PLATES+1);
}




void plotRayWeld2(ray r){
	ofstream myfile;
	myfile.open("Data.txt");
	

    vector <vertex_type> rP;
    cout<<endl<<"Entering the weld";
    double current_crystal_orientation, next_crystal_orientation, boundary_orientation;
    
	//r.current_position = r.init_position;
 //   r.global_direction = r.init_global_direction;
	//r.phase_velocity = r.init_phase_velocity;
	//
    //r.group_vel_global_direction = 0; //Why was this required?
    
	vertex_type v;
    glBegin(GL_LINE_STRIP);
    //glVertex3f(r.current_position.x, r.current_position.y, r.current_position.z); //Mark only first point
    rP.clear();
    rP.push_back(r.current_position);
	
	cout<<endl<<r.current_position.y<<" "<<r.current_position.z<<" "<<radtodeg(r.global_direction)<<" "<<radtodeg(getCrystalOrientation(r.current_position.y, r.current_position.z))<<" "<<r.phase_velocity;
	ray_heading_temp = r.global_direction;
	ray_heading_temp_group = r.group_vel_global_direction;

	switch(weld_entry){
		case 0: boundary_orientation = 0; break;
		case 1: boundary_orientation = degtorad(90-A/2); break;
		case 2: boundary_orientation = degtorad(90+A/2); break;
	}
	cout<<endl<<weld_entry;

	ray rotated_ray = r;
    rotated_ray.global_direction = r.global_direction - boundary_orientation;
	cout<<endl<<"Before: "<<rotated_ray.global_direction<<" "<<radtodeg(boundary_orientation);
	//rotated_ray.group_vel_global_direction = r.group_vel_global_direction - boundary_orientation;
	
	//In the rotated ray, only the phase velocity direction differs. (is rotated)  
	
	if(rotated_ray.global_direction<0){
		rotated_ray.global_direction = rotated_ray.global_direction + degtorad(360);
	}else if((rotated_ray.global_direction<degtorad(180)) && (rotated_ray.global_direction>0)){
		rotated_ray.global_direction = rotated_ray.global_direction + degtorad(180);
	}else if(rotated_ray.global_direction>degtorad(180)){
		//Do nothing
	}
	cout<<endl<<"Before: "<<radtodeg(rotated_ray.global_direction)<<" "<<rotated_ray.phase_velocity;

	if((abs(rotated_ray.global_direction - degtorad(180))<degtorad(0.03))||(abs(rotated_ray.global_direction)<degtorad(0.03))){
		//Do nothing. Ray is almost parallel to boundary => Bypasses it. 
	}else{
		calculateAnisotropic2(rotated_ray,radtodeg(getCrystalOrientation(rotated_ray.current_position.y, rotated_ray.current_position.z)));
		//Above calculateAnisotropic2 call is the 'first' call meant for entering the weld. New phase velocity and group velocity are obtained.
	}
	cout<<endl<<"After: "<<radtodeg(rotated_ray.global_direction)<<" "<<rotated_ray.phase_velocity<<" "<<rotated_ray.group_vel_global_direction;
	r = rotated_ray; //Copy all attributes: position, phasevel, groupvel, directions
	//Then correct only the directions
	r.global_direction = rotated_ray.global_direction + boundary_orientation;
	if(abs(r.global_direction - ray_heading_temp)>degtorad(150)){
		r.global_direction = r.global_direction - degtorad(180);
	}
    r.group_vel_global_direction = rotated_ray.group_vel_global_direction + boundary_orientation;
	if(abs(r.group_vel_global_direction - ray_heading_temp_group)>degtorad(150)){
		r.group_vel_global_direction = r.group_vel_global_direction - degtorad(180);
	}


	cout<<endl<<"Realigned inside weld";
	cout<<endl<<r.current_position.y<<" "<<r.current_position.z<<" "<<radtodeg(r.global_direction)<<" "<<radtodeg(getCrystalOrientation(r.current_position.y, r.current_position.z))<<" "<<r.phase_velocity<<" "<<radtodeg(r.group_vel_global_direction);
	//cout<<endl<<"#### "<<isInsideWeld(rP.back());
	while(isInsideWeld(rP.back())){
		cout<<"Check";
		switch (r.type) {
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
		cout<<endl<<r.current_position.y<<" "<<r.current_position.z<<" "; 
		myfile<<r.current_position.y<<"\t"<<r.current_position.z<<"\t"<<radtodeg(r.global_direction)<<"\t"<<r.phase_velocity<<"\t"<<radtodeg(r.group_vel_global_direction)<<"\t"; //!!
		ray_heading_temp = r.global_direction;
		ray_heading_temp_group = r.group_vel_global_direction;

        current_crystal_orientation = (getCrystalOrientation(r.current_position.y,r.current_position.z));
		cout<<radtodeg(current_crystal_orientation);
		myfile<<radtodeg(current_crystal_orientation)<<"\t";//!!
        next_crystal_orientation = (getCrystalOrientation(r.current_position.y+step_size*cos(r.group_vel_global_direction),r.current_position.z+step_size*sin(r.group_vel_global_direction)));
		myfile<<radtodeg(next_crystal_orientation)<<"\t";//!!
        boundary_orientation = (current_crystal_orientation + next_crystal_orientation)/2;
//        cout<<endl<<"Current Position:"<<r.current_position.y<<" "<<r.current_position.z;
//        cout<<endl<<"Boundary orientation: "<<radtodeg(boundary_orientation);
//        cout<<endl<<"Ray - global direction group: "<<radtodeg(r.group_vel_global_direction);
        //ROTATION: Rotate clockwise if on the left half, and anti-clockwise if on the right side. Basically, add the negative of the boundary orientation.
        //Basic idea is to create a new ray which reflects the rotated scenario, calculate the change in the path and rotate it back.

        rotated_ray = r;
        rotated_ray.global_direction = r.global_direction - boundary_orientation;
		rotated_ray.group_vel_global_direction = r.group_vel_global_direction - boundary_orientation;
		//In the rotated ray, only the phase velocity direction differs. (is rotated)  
		myfile<<radtodeg(current_crystal_orientation - boundary_orientation)<<"\t"<<radtodeg(next_crystal_orientation - boundary_orientation)<<"\t";//!!
		myfile<<radtodeg(rotated_ray.global_direction)<<"\t";//!!

		if(rotated_ray.global_direction<0){
			rotated_ray.global_direction = rotated_ray.global_direction + degtorad(360);
		}else if((rotated_ray.global_direction<degtorad(180)) && (rotated_ray.global_direction>0)){
			rotated_ray.global_direction = rotated_ray.global_direction + degtorad(180);
		}else if(rotated_ray.global_direction>degtorad(180)){
			//Do nothing
		}
		myfile<<radtodeg(rotated_ray.global_direction)<<"\t";//!!
		if((abs(rotated_ray.global_direction - degtorad(180))<degtorad(0.03))||(abs(rotated_ray.global_direction)<degtorad(0.03))){
			//Do nothing. Ray is almost parallel to boundary => Bypasses it. 
		}else{
			//cout<<endl<<"Check1";
			calculateAnisotropic4(rotated_ray, radtodeg(current_crystal_orientation - boundary_orientation), radtodeg(next_crystal_orientation - boundary_orientation));
		}
		r = rotated_ray; //Copy all attributes: position, phasevel, groupvel, directions
		//Then correct only the directions
		myfile<<radtodeg(rotated_ray.global_direction)<<"\t";//!!
		myfile<<radtodeg(rotated_ray.group_vel_global_direction)<<"\t";//!!
		r.global_direction = rotated_ray.global_direction + boundary_orientation;
		r.group_vel_global_direction = rotated_ray.group_vel_global_direction + boundary_orientation;
		//cout<<endl<<"RayHeadingGV "<<radtodeg(ray_heading_temp_group);
		for(int i=0; i<1; i++){
			if(abs(r.global_direction - ray_heading_temp)>degtorad(90)){
				r.global_direction = r.global_direction - degtorad(180);
			}
			//cout<<endl<<"*! "<<radtodeg(r.group_vel_global_direction);
			//cout<<endl<<"** "<<radtodeg(abs(r.group_vel_global_direction - ray_heading_temp_group));
			if(abs(r.group_vel_global_direction - ray_heading_temp_group)>degtorad(90)){
				r.group_vel_global_direction = r.group_vel_global_direction - degtorad(180);
				//cout<<endl<<"$! "<<radtodeg(r.group_vel_global_direction);
			}
		}
		//cout<<endl<<"GV "<<radtodeg(r.group_vel_global_direction);
		myfile<<radtodeg(r.global_direction)<<"\t"<<r.phase_velocity<<"\t"; //!!
		myfile<<radtodeg(r.group_vel_global_direction)<<"\n ";//!!
        //r.group_velocity = rotated_ray.group_velocity;

        //cout<<endl<<"Global ray params: "<<radtodeg(r.global_direction)<<" "<<radtodeg(r.group_vel_global_direction)<<" "<<r.group_velocity;
//        cout<<endl<<"Current position in: "<<r.current_position.x<<" "<<r.current_position.y<<" "<<r.current_position.z;
        //Step forward
//        r.current_position.y += 1.0*cos(r.group_vel_global_direction);
//        r.current_position.z += 1.0*sin(r.group_vel_global_direction);
        v.x = rP.back().x;
        v.y = rP.back().y + step_size*cos(r.group_vel_global_direction);
        v.z = rP.back().z + step_size*sin(r.group_vel_global_direction);


        r.current_position = v;
		//cout<<"$$$$ "<<v.x<<" "<<v.y<<" "<<v.z;
        
		
		//cout<<endl<<"Current position: "<<v.x<<" "<<v.y<<" "<<v.z;
        
		
		
		//glVertex3f(r.current_position.x, r.current_position.y, r.current_position.z);
        rP.push_back(v);
    }
    cout<<endl<<"This is the end...";
    rayPaths.push_back(rP);
    glEnd();
	myfile.close();
}

void displayWeldRayPath(){
    for (int l = 0; l<rayPaths.size(); l++) {
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
        for (int i = 0; i<rayPaths[l].size(); ++i) {
			//cout<<endl<<rayPaths[l][i].y;
            glVertex3f(0,rayPaths[l][i].y,rayPaths[l][i].z);
        }
		//cout<<endl<<"**********************************************";
        glEnd();
    }
}


ray plotRayisotropic(ray r)
{
	//Find equations of both left and right weld boundaries
	cout<<endl<<"START"<<endl;
	//float m1 = tan(degtorad(90-A/2));
	float c = -m1*D/2;
	float c_ray; 
	float side = 1; //1 -> Right, -1 -> Left
	vertex_type v;
	float new_direction;
	vector <vertex_type> rP;
	//glBegin(GL_LINE_STRIP);
    //glVertex3f(r.current_position.x, r.current_position.y, r.current_position.z); //Mark only first point
    rP.clear();
    rP.push_back(r.current_position);

	while(r.current_position.y>-PLATE_WIDTH/2 && r.current_position.y<PLATE_WIDTH/2)
	{
		cout<<endl<<r.current_position.y<<" "<<r.current_position.z<<" "<<r.global_direction;

		if(isInsideWeld(rP.back()))
		{
			cout<<endl<<"Inside Weld => Anisotropy";
			break;
		}
		if(r.global_direction>degtorad(180))
		{
			v.y = -r.current_position.z/tan(r.global_direction) + r.current_position.y;
			v.z = 0;
		}
		else
		{
			v.y = (PLATE_HEIGHT-r.current_position.z)/tan(r.global_direction) + r.current_position.y;
			v.z = PLATE_HEIGHT;
		}
		new_direction = degtorad(360) - r.global_direction;
		if(isInsideWeld(v) || ((r.current_position.y*v.y)<0))
		{
			c_ray = -r.current_position.y*tan(r.global_direction) + r.current_position.z;
			if(r.current_position.y<0)
				side = -1;
			else 
				side = 1;
			if((side*m1 - tan(r.global_direction))<(10^-5))
			{
				cout<<"Ray and weld almost parallel";
				break;
			}
			r.current_position.y = (c_ray - c)/(side*m1 - tan(r.global_direction));
			r.current_position.z = side*m1*r.current_position.y + c;
			rP.push_back(r.current_position);
			if(side==1)
				weld_entry = 1;
			else
				weld_entry = 2;
			break;
		}	
		r.global_direction = new_direction; 
		r.current_position = v;
		rP.push_back(v);
	}
	rayPaths.push_back(rP);
	cout<<endl<<r.current_position.y<<" "<<r.current_position.z<<" "<<r.global_direction;
	cout<<endl<<"Completed Isotropic Ray Tracing"<<endl;
	return r;
}
#endif	