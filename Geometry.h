#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Definitions.h"

/*********************************************************
*
*GEOMETRY DEFINITIONS
*
**********************************************************/

//To draw an arc
void drawarc(float cx, float cy, float z, float r, float start_angle, float arc_angle)
{
	glBegin(GL_LINES);
	float theta = start_angle*2*PI/360;
	float end_angle = (start_angle + arc_angle)*2*PI/360;
	for(; theta < end_angle; theta+=0.005)
	{
		glVertex3f(cx + r*cos(theta), cy + r*sin(theta), z);
		glVertex3f(cx + r*cos(theta+0.005), cy + r*sin(theta+0.005), z);
	}
	glEnd();
}

//TO DISPLAY THE GLOBAL AXES
void drawaxes()
{
	vertex_type v;
	v.x = v.y = v.z = 0;
	glBegin(GL_LINES);
	glColor3f(1.0,0.0,0.0);  //R
	glVertex3f(v.x, v.y, v.z); //X
	glVertex3f(v.x+10,v.y,v.z);
	glColor3f(0.0,1.0,0.0);  //G
	glVertex3f(v.x,v.y,v.z); //Y
	glVertex3f(v.x,v.y+10,v.z);
	glColor3f(0.0,0.0,1.0);  //B
	glVertex3f(v.x,v.y,v.z); //Z
	glVertex3f(v.x,v.y,v.z+10);
	glEnd();
}

//TO DRAW THE WELD
void drawweld()
{
	// 	glBegin(GL_LINES);
	glColor3f(0.0,0.0,0.0);
	glBegin(GL_LINES);
	glVertex3f(0,D/2,0);
	glVertex3f(0,D/2+tan(degtorad(A/2))*PLATE_HEIGHT,PLATE_HEIGHT);
	glVertex3f(0,-D/2,0);
	glVertex3f(0,-D/2-tan(degtorad(A/2))*PLATE_HEIGHT,PLATE_HEIGHT);
	glEnd();
}

void drawdefect()
{
	glColor3f(0.0,0.0,1.0);
	glBegin(GL_LINES);
	glVertex3f(0,defect_y2,defect_z2);
	glVertex3f(0,defect_y1,defect_z1);
	glEnd();
}

void drawpipe(float W, float H)
{
	//Main pipe arcs.
	glColor3f(0.0,0.0,0.0);
	//Display The pipe side lines.
	glBegin(GL_LINE_STRIP);
	glVertex3f(0,-W/2,0);
	glVertex3f(0,-W/2,H);
	// glVertex2f(-100,20);
	glVertex3f(0,W/2,H);
	glVertex3f(0,W/2,0);
	glVertex3f(0,-W/2,0);
	glEnd();
}

void drawpipe2(float W, float H, int extra_plates = 0)
{

	//Main pipe arcs.
	glColor3f(0.0,0.0,0.0);
	//Display The pipe side lines.
	glBegin(GL_LINE_STRIP);
	glVertex3f(0,-W/2,0 );
	glVertex3f(0,-W/2,H);
	// glVertex2f(-100,20);
	glVertex3f(0,W/2,H);
	glVertex3f(0,W/2,0);
	glVertex3f(0,-W/2,0);
	glEnd();

	for (int i = 1; i <= extra_plates; ++i) 
	{
		//Main pipe arcs.
		glColor3f(0.0,0.0,0.0);
		//Display The pipe side lines.
		glBegin(GL_LINE_STRIP);
		glVertex3f(0,-PLATE_WIDTH/2,-i*PLATE_HEIGHT);
		glVertex3f(0,-PLATE_WIDTH/2,PLATE_HEIGHT-i*PLATE_HEIGHT);
		glVertex3f(0,PLATE_WIDTH/2,PLATE_HEIGHT-i*PLATE_HEIGHT);
		glVertex3f(0,PLATE_WIDTH/2,-i*PLATE_HEIGHT);
		glVertex3f(0,-PLATE_WIDTH/2,-i*PLATE_HEIGHT);
		glEnd();
	}
}

void drawgraph()
{
	glColor3f(0.0,0.0,0.0);
	glBegin(GL_LINE_STRIP);
	//glColor3d(1.0f, 1.0f, 1.0f);
	glVertex3f(0, -gwidth/2, -gheight + offset);
	//glColor3d(1.0f, 0.0f, 0.0f);
	glVertex3f(0, -gwidth/2, offset);
	//glColor3d(0.0f, 1.0f, 0.0f);
	glVertex3f(0, gwidth/2, offset);
	//glColor3d(0.0f, 0.0f, 1.0f);
	glVertex3f(0, gwidth/2, -gheight + offset);
	glVertex3f(0, -gwidth/2, -gheight + offset);
	glColor4f(0.0, 0.0, 0.0, 0.25);
	glVertex3f(0, -gwidth/2, -gheight/2 + offset);
	glVertex3f(0, gwidth/2, -gheight/2 + offset);
	glEnd();
	
	glPushAttrib(GL_ENABLE_BIT); 
	glLineStipple(4, 0xAAAA);  
	glColor4f(0.0, 0.0, 0.0, 0.15);
	glEnable(GL_LINE_STIPPLE);
	glBegin(GL_LINES);
	for(int i=1; i<150;i++)
	{	
		glVertex3f(0,-gwidth/2 + 2*i, offset);
		glVertex3f(0,-gwidth/2 + 2*i, -gheight + offset);
	}

	glEnd();
	glPopAttrib();




	glColor3f(0.0, 0.0, 1.0);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glLineWidth(1.0f);
	glBegin(GL_LINE_STRIP);
	glVertex3f(0, - gwidth/2, - gheight/2 + offset);
	for(int i=0; i<135;++i)
	{	
		glVertex3f(0,pulse[i].x/4 - gwidth/2, (pulse[i].y)*scale - gheight/2 + offset);
		//glVertex3f(0,i,i);
	}
	scan_diff = han_freq * (tot_time_weld+tot_time_iso_1+tot_time_iso_2)*100;
	//cout<<"\nCheck "<<scan_diff;
	for(int i=0; i<135;++i)
	{	

		glVertex3f(0,pulse[i].x/4 - gwidth/2 + scan_diff, (pulse[i].y)*scale*amplitude - gheight/2 + offset);
		//glVertex3f(0,i,i);
	}
	glVertex3f(0,  gwidth/2, - gheight/2 + offset);



	glEnd();
	glLineWidth(1.0f);
	glDisable(GL_LINE_SMOOTH);
	glColor3f(0.0,0.0,0.0);
	outtext(0,-gwidth/2 - 25,-gheight/2 + offset,"Amplitude");
	outtext(0,-10,-gheight + offset - 5,"Number");

}

void drawprobe()
{
	//glColor4f(0.0, 0.0, 0.0, 0.15);
	//Probe 1
	glBegin(GL_LINE_STRIP);
	glVertex3f( ray_x, ray_y-(wedge_width/2), PLATE_HEIGHT);
	glVertex3f( ray_x, ray_y-(wedge_width / 2), PLATE_HEIGHT+wedge_height+(tan(wedge_angle)*wedge_width));
	glVertex3f( ray_x, ray_y+ (wedge_width / 2), PLATE_HEIGHT+wedge_height);
	glVertex3f( ray_x, ray_y+ (wedge_width / 2), PLATE_HEIGHT);
	glEnd();

	/*
	//Probe 2
	glBegin(GL_QUADS);
	glVertex3f( ray_x, -ray_y-6, ray_z);
	glVertex3f( ray_x, -ray_y-6, ray_z+6);
	glVertex3f( ray_x, -ray_y+6, ray_z+6);
	glVertex3f( ray_x, -ray_y+6, ray_z);
	glEnd();
	*/
}
void showElastic(float T, float eta)
{
	glBegin(GL_LINES);
	for(float y = 1; y<PLATE_HEIGHT; y += 2)
	{
		int maxX = D/2 + tan(degtorad(A/2))*y;
		int minX = -D/2 - tan(degtorad(A/2))*y;
		for(float x = 0; x<maxX; x+=2)
		{
			float tantheta = -1.0f*T*(D/2+y*tan(degtorad(A/2)))/pow(x,eta);
			glVertex3f(0,x,y);
			glVertex3f(0,x+1,y+tantheta);
		}
		for(float x = 0; x>minX; x-=2)
		{
			float tantheta = 1.0f*T*(D/2+y*tan(degtorad(A/2)))/pow(-x,eta);
			glVertex3f(0,x,y);
 			glVertex3f(0,x-1,y-tantheta);
		}
	}
	glVertex3f(0,0,0);
	glVertex3f(0,0,PLATE_HEIGHT);
	glEnd();
}  

void showElastic2()
{
	float tantheta = tan(degtorad(multiplePlateElasticOrientation[0]));
	glBegin(GL_LINES);
	for(float y = 1; y<PLATE_HEIGHT; y+=2)
	{
		int maxX = D/2 + tan(degtorad(A/2))*y;
		int minX = -D/2 - tan(degtorad(A/2))*y;
		for(float x = -PLATE_WIDTH/2; x<minX; x+=2)
		{
			glVertex3f(0,x,y);
			glVertex3f(0,x+1,y+tantheta);
		}
		for(float x = maxX; x<PLATE_WIDTH/2; x+=2)
		{
			glVertex3f(0,x,y);
			glVertex3f(0,x+1,y+tantheta);
		}
	}
	glEnd();
}  

void showConstThetaLines(float T, float eta)
{
	glColor3f(0,0,1);
	glBegin(GL_LINES);
	for(float y = 1; y<PLATE_HEIGHT;y+=1)
	{
		int maxX = D/2 + tan(degtorad(A))*y;
		int minX = -D/2 - tan(degtorad(A))*y;
		for(float x=0; x>minX; x-=1)
		{
			float C = (D+y*tan(degtorad(A)))/(tan(degtorad(A))*(pow(-x,eta)));
			float pointXBelow = -1.0f*pow((pow(-x,eta)*D)/C,1/eta);
			float pointXAbove = -1.0f*pow((pow(-x,eta)*D + PLATE_HEIGHT)/C,1/eta);
			glVertex2f(pointXBelow,0);
			glVertex2f(pointXAbove,PLATE_HEIGHT);
		}
	}
	glEnd();
}

#endif