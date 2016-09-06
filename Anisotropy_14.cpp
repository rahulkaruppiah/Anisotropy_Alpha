#include <stdlib.h>
#include <stdio.h>
#include "GLUT/glut.h"
#include <math.h>
#include <vector>
#include <cstring>
#include "Eigen/Dense"
#include <fstream>
#include <iostream>

#include "Geometry.h"
#include "Definitions.h"
//#include "Old-Propagation.h"     Contains a lot of old functions - Keeping it just for future reference
#include "Propogation.h"
//#include "Algorithm.h"           Algorithms for reverse tracing

#include "Graph.h"
#include <time.h>


using namespace std;
using namespace Eigen;


// The width and height of your window, change it as you like
int screen_width=1366;
int screen_height=768;

int x_translate = -100;
int y_translate = 0;
int z_translate = -5;

// Absolute rotation values (0-359 degrees) and rotiation increments for each frame
double rotation_x=0, rotation_x_increment=0.0;
double rotation_y=0, rotation_y_increment=0.0;
double rotation_z=0, rotation_z_increment=0.0;

float mouse_x, mouse_y;

// Flag for rendering as lines or filled polygons
int filling=1; //0=OFF 1=ON

void IntroDisplay()
{
	cout<<"==========================================================================================================================\n";
	cout<<"\n\t\t\t\t\tRAY PROPOGATION IN ANISOTROPIC WELDS"<<"\n";
	cout<<"\n==========================================================================================================================\n\n\n";
	
	//Elastic constants
	initElasticConstants();
	cout << endl << "Elastc Constants (Transversely Isotropic) have been initialized\n" << endl;
	displayC(c);
}





void init(void)
{
	glClearColor(1.0, 1.0, 1.0, 0.0); // This clears the background color to white

	glShadeModel(GL_SMOOTH); // Type of shading for the polygons

	// Viewport transformation
	glViewport(0,0,screen_width,screen_height);

	// Projection transformation
	glMatrixMode(GL_PROJECTION); // Specifies which matrix stack is the target for matrix operations
	glLoadIdentity(); // We initialize the projection matrix as identity
	gluPerspective(90.0f,(GLfloat)screen_width/(GLfloat)screen_height,1.0f,1500.0f); // We define the "viewing volume"

	glEnable(GL_DEPTH_TEST); // We enable the depth test (also called z buffer)

	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL); // Polygon rasterization mode (polygon filled)


	//Trying to achieve transparency
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	//Trying to enable points
	glEnable(GL_POINT_SMOOTH);



	//Initializing probes
	/* 
	&&&
	At this point it is possible to initialise multiple probes. If necessary, a loop can be implemented in case a sensor array is required. 
	Also, some aspects need to be noted. 
	1. Angle is always measured from positive y axis. 
	2. Angle is always given positive. (=> Positive for both positive and negative x)
	3. Coordinate system followed: Rokhlin et al. y-z plane on the screen, x coming out towards you.
	4. If it so happens that the ray is directed away from the centre, the loop runs infinitely. (since the exit condition is not achieved)
	*/

	theRay.init_position.x = ray_x;
	theRay.init_position.y = ray_y; 
	theRay.init_position.z = ray_z;
	theRay.type = ray_type; 
	theRay.current_position = theRay.init_position;
	theRay.global_direction = ray_angle;
	theRay.group_vel_global_direction = ray_angle;
	theRay.init_global_direction = ray_angle;
	theRay.init_phase_velocity = phase_velocity;
	theRay.phase_velocity = theRay.init_phase_velocity;
	theRay_initial = theRay;
	
	rayList.push_back(theRay); //Newly added 250114
	cout << "\n\n==========================================================================================================================\n";
	cout<<"\n\t\t\t\t\t\tRENDERING COUNT - "<<ren_count<<"\n";
	cout << "\n==========================================================================================================================\n";


	//Weld Parameters:
	cout<<"\nWeld Parameters:\teta: "<<eta<<"\tT: "<<T<<"\n";
	cout << "\n==========================================================================================================================\n";


	//cout<<"\nDo you want to display the Shear Wave?";
	//cin>>SHEAR;

	//displayC(rotateElasticConstants(PI/180 * 45));
	//cout << endl<< rotateElasticConstants(45) << endl; //Was commented for convenience

	//solveChristoffelEquation(&theRay);
	// updatep();
	// init_probe(&newp,-100.0,200.0,68.0);
	// p.push_back(newp);
}



/**********************************************************
* resize(int,int)
*
* This routine must be called everytime we resize our window.
*********************************************************/

void resize (int width, int height)
{
	screen_width=width; // We obtain the new screen width values and store it
	screen_height=height; // Height value

	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // We clear both the color and the depth buffer so to draw the next frame
	glViewport(0,0,screen_width,screen_height); // Viewport transformation

	glMatrixMode(GL_PROJECTION); // Projection transformation
	glLoadIdentity(); // We initialize the projection matrix as identity
	gluPerspective(90.0f,(GLfloat)screen_width/(GLfloat)screen_height,1.0f,1500.0f);

	glutPostRedisplay (); // This command redraw the scene (it calls the same routine of glutDisplayFunc)
}



/**********************************************************
* keyboard(unsigned char,int,int)
* Used to handle the keyboard input (ASCII Characters)
*********************************************************/

void keyboard (unsigned char key, int x, int y)
{

	switch (key)
	{

	case ' ':
		rotation_x_increment=0;
		rotation_y_increment=0;
		rotation_z_increment=0;
		break;
	case 'r': case 'R':
		if (filling==0)
		{
			glPolygonMode (GL_FRONT_AND_BACK, GL_FILL); // Polygon rasterization mode (polygon filled)
			filling=1;
		}
		else
		{
			glPolygonMode (GL_FRONT_AND_BACK, GL_POINT); // Polygon rasterization mode (polygon outlined)
			filling=0;
		}
		break;
	case 'e': case 'E':
		dispElastic = !dispElastic;
		break;
	case 'f': case 'F':
		dispElastic2 = !dispElastic2;
		break;
	case 'd':
		eta+=0.1f;
		init();
		DUMMY = 0;
		break;
	case 'a':
		eta-=0.1f;
		init();
		DUMMY = 0;
		break;
	case 'w':
		T+=0.2f;
		init();
		DUMMY=0;
		break;
	case 'q':
		T-=0.2f;
		init();
		DUMMY=0;
		break;
	case 'z':
		ray_angle+=degtorad(5);
		init();
		DUMMY=0;
		break;
	case 'x':
		ray_angle-=degtorad(5);
		init();
		DUMMY=0;
		break;
		//case 'c':
		//	multiplePlateElasticOrientation[0]-=1; 
		//	init();
		//	break;
		//case 'v':
		//	multiplePlateElasticOrientation[0]+=1;
		//	init();
		//	break;
		//case 'b':
		//	multiplePlateElasticOrientation[1]-=1; 
		//	init();
		//	break;
		//case 'n':
		//	multiplePlateElasticOrientation[1]+=1;
		//	init();
		//	break;
		//case 'j':
		//	cin>>ray_y;
		//	init();
	case 'p':
		//cout<<"CurIniPos: "<<ray_y<<" "<<ray_z<<endl;
		//cout<<"CurIniDir: "<<radtodeg(ray_angle);
		//cout<<endl<<"NewIniPos: ";
		//cin>>ray_y>>ray_z;
		//cout<<"NewIniDir: ";
		cin>>ray_angle;
		ray_angle = degtorad(ray_angle);
		//cin>>phase_velocity;
		//cout<<"CurElasOri: "<<multiplePlateElasticOrientation[0]<<endl;
		//cout<<"NewElasOri: "; cin>>multiplePlateElasticOrientation[0];

		init();
		break;
	case 'n':
		ray_y = ray_y - 5;
		/*if(ray_y == 0){
		ray_y = ray_y - 5;
		}*/
		init();
		break;
	case 'b':
		ray_y = ray_y + 5;
		init();
		break;
	case 'm':
		//cin>>step_size;
		cin>>phase_velocity;
		init();
		break;

	case 27:
		exit(0);
		break;
	}
}



/**********************************************************
*
* keyboard_s(int,int,int)
*
* Used to handle the keyboard input (not ASCII Characters)
*
*********************************************************/

void keyboard_s (int key, int x, int y)
{

	switch (key)
	{
	case GLUT_KEY_UP:
		//rotation_x_increment = rotation_x_increment +0.005;
		// 	    rotation_x_increment = rotation_x_increment + 2;
		z_translate += 1;
		break;
	case GLUT_KEY_DOWN:
		//rotation_x_increment = rotation_x_increment -0.005;
		// 	    rotation_x_increment = rotation_x_increment - 2;
		z_translate -= 1;
		break;
	case GLUT_KEY_LEFT:
		//rotation_y_increment = rotation_y_increment +0.005;
		// 	    rotation_y_increment = rotation_y_increment + 2;
		ray_y -= 1;
		init();
		DUMMY=0;
		break;
	case GLUT_KEY_RIGHT:
		//rotation_y_increment = rotation_y_increment -0.005;
		// 	    rotation_y_increment = rotation_y_increment - 2;
		ray_y  += 1;
		init();
		DUMMY=0;
		break;
	}
}

/**********************************************************
*
* mouse(int,int,int, int)
*
* Used to handle the mouse input
*
*********************************************************/

void mouse(int key, int state, int x, int y)
{
	switch(key)
	{
	case GLUT_LEFT_BUTTON:
		if(state == GLUT_UP)
		{
			mouse_x = mouse_y = 0;
		}
		else if(state == GLUT_DOWN)
		{
			mouse_x = x;
			mouse_y = y;
		}
		break;

	case 3:
		x_translate = x_translate + 2;
		glTranslatef(x_translate,0,0);
		break;

	case 4:
		x_translate = x_translate - 2;
		glTranslatef(x_translate,0,0);
		break;
	}
}

void mouseMove(int x, int y)
{
	if(mouse_x>0)
	{
		y_translate = y_translate + (x - mouse_x)*0.5f;
		mouse_x = x;
	}
	if(mouse_y>0)
	{
		z_translate = z_translate - (y - mouse_y)*0.5f;
		mouse_y = y;
	}
}


/**********************************************************
*
* display()
*
* This is our main rendering subroutine, called each frame
*
*********************************************************/

void display(void)
{
	//     int l_index;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // This clear the background color to white
	glMatrixMode(GL_MODELVIEW); // Modeling transformation
	glLoadIdentity(); // Initialize the model matrix as identity

	glRotatef(-90,0.0,1.0,0.0);	//Rotate about Y to set Y-Z plane as the one on the screen.
	glRotatef(-90,1.0,0.0,0.0); // Rotate about X to set the Y-Z as described in Rokhlin.



	rotation_x = rotation_x + (rotation_x_increment - rotation_x)/50;
	rotation_y = rotation_y + (rotation_y_increment - rotation_y)/50;

	//rotation_x = rotation_x + rotation_x_increment;
	//rotation_y = rotation_y + rotation_y_increment;
	rotation_z = rotation_z + rotation_z_increment;

	if (rotation_x > 359) rotation_x = 0;
	if (rotation_y > 359) rotation_y = 0;
	if (rotation_z > 359) rotation_z = 0;

	if(rotation_x_increment > 359) rotation_x_increment = 0;
	if(rotation_y_increment > 359) rotation_y_increment = 0;
	if(rotation_z_increment > 359) rotation_z_increment = 0;

	glRotatef(rotation_x,1.0,0.0,0.0); // Rotations of the object (the model matrix is multiplied by the rotation matrices)
	glRotatef(rotation_y,0.0,1.0,0.0);
	glRotatef(rotation_z,0.0,0.0,1.0);

	//glTranslatef(x_translate,0.0,0.0);
	//glTranslatef(0.0,y_translate,0.0);
	//glTranslatef(0,0,z_translate); // We move the object 400 points forward (the model matrix is multiplied by the translation matrix)
	glTranslatef(x_translate, y_translate, z_translate);

	//Render axes.
	drawaxes();

	//Draw the pipe
	drawpipe2(PLATE_WIDTH,PLATE_HEIGHT,NUM_EXTRA_PLATES);



	//Show elastic constants
	if(dispElastic)
		showElastic(T,eta);

	if(dispElastic2)
		showElastic2();

	//plotMultiplePlates(theRay_initial);


	//WELD!!
	drawweld(); 
	drawprobe();
	if(defect_enable)
		drawdefect();
	//vertex_type temp;
	//temp.x = 0;
	//temp.y = 0.05;
	//temp.z = 25;
	//cout<<isAtDefect(temp);

	if(DUMMY<2)
	{
		

		if(DUMMY==0)
			ren_count++;


	/*
		
		theRay_temp = plotRayisotropic(theRay_initial);  // Before Weld 
		theRay_temp = plotRayWeld(theRay_temp,eta, T);   // In Weld
		theRay = plotRayisotropicExit(theRay_temp);      // After Weld
	    
		*/
		
		testcpp();
		
		//reversetracing(theRay_initial, theRay);
		//fib_alg(3,PI,2*PI);
		//graphprint(theRay_initial);
		//bruteforce(theRay_initial);
		//new_alg(theRay_initial);


		if(num_defect)
		{
			theRay_temp = plotRayisotropic(theRay_initial);
			theRay_temp = plotRayWeld(theRay_temp,eta, T);
			num_defect = 0;
		}

		if(tot_count<0)
		{
			new_angle = new_angle + degtorad(1);
			theRay_initial.global_direction = new_angle;
			theRay_initial.init_global_direction = new_angle;
			tot_count++;
		}



		if (SHEAR == 1)
		{
			DUMMY++;
			cout << "\n\nDummy7" << DUMMY;

		}
		else
			DUMMY = 2;

		
		//Graph Display
	//	plotResults();

	}

	displayWeldRayPath();
	glutSwapBuffers(); // In double buffered mode we invert the positions of the visible buffer and the writing buffer
	glEnd();
	
}



/**********************************************************
* The main routine
***********************************************************/

int main(int argc, char **argv)
{

	// Use the GLUT utility to initialize the window, to handle the input and to interact with the windows system
//	Master();
	IntroDisplay();  // Initializes the Material Constant and displays the program heading
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(screen_width,screen_height);
	glutInitWindowPosition(0,0);
	glutCreateWindow("Ultrasonic Testing");
	
	glutDisplayFunc(display);
	glutIdleFunc(display);

	glutReshapeFunc (resize);
	glutKeyboardFunc (keyboard);
	glutSpecialFunc (keyboard_s);
	glutMouseFunc(mouse);
	glutMotionFunc(mouseMove);
	init();
	
	glutMainLoop();
	
    return 0;
}
