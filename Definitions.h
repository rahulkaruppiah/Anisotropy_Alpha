#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <vector>
#include <iostream>
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;

#define PI 3.14159265
#define SPEED 5.893	    //mm/us; STEEL 1020 (5893m/s) Wavelength (1.47 x 10^-3)
#define FREQ 4000000	//4MHz frequency
#define WAVELENGTH SPEED/FREQ
#define PROBE_DIA 25
/**********************************************************
*
* TYPES DECLARATION
*
***********************************************************/

// Vertex type
typedef struct
{
	float x,y,z;
}vertex_type;

//Probe type
typedef struct
{
	vertex_type position;
	float absAngle;
	float localAngle;
	bool onOD;
	float path_length;
}probe;

typedef struct
{
	int type; //0 for qP, 1 for qSV, and 2 for SH -> all in Anisotropic medium.
	vertex_type init_position;
	vertex_type current_position;
	float phase_velocity;
	float init_phase_velocity;
	//float slowness;
	//float current_direction;
	float local_direction;
	//Newly added lines 240114
	float group_velocity; 
	float global_direction;
	float group_vel_global_direction;
	float init_global_direction;
	float time;
}ray;

typedef struct 
{
	float x;
	float y;
}point;

//ANGLE CONVERSIONS

float radtodeg(float angle)
{
	return angle*360/(2*PI);
}

float degtorad(float angle)
{
	return angle*2*PI/360;
}

float slopetoangle(float slope)
{
	return atan(slope); //angle returned in radians
}

void quadratic(float A, float B, float C, float* x1, float* x2)
{
	(*x1) = (-B + sqrt(B*B - 4*A*C))/2*A;
	(*x2) = (-B - sqrt(B*B - 4*A*C))/2*A;
}

void outtext(int x, int y, int z, char *string)
{
	int len, i;
	glRasterPos3f(x, y, z);
	len = (int) strlen(string);
	for (i = 0; i < len; i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, string[i]);
	}
}

//For type 308 austenitic stainless steel.
//Elastic constants

// float c[7][7] = {0};
MatrixXd c(6,6);

void initElasticConstants()
{
	//From [GDConnolly thesis - transversely isotropic material - table 3.1]
	//Also refer http://en.wikipedia.org/wiki/Hooke's_law#Anisotropic_materials

	c = MatrixXd::Zero(6,6);
	c(0,0) = 249e9;
	c(0,1) = 124e9;
	c(0,2) = 133e9;
	c(1,0) = c(0,1);
	c(1,1) = c(0,0);
	c(1,2) = c(0,2);
	c(2,0) = c(0,2);
	c(2,1) = c(0,2);
	c(2,2) = 205e9;
	c(3,3) = 125e9;
	c(4,4) = c(3,3);
	c(5,5) = 0.5f*(c(0,0) - c(0,1));

}

void displayC(MatrixXd L)
{
	/*for(int i=0;i<6;i++)
	{
	for(int j=0;j<6;j++)
	cout<<"\t"<<L(i,j); 
	cout<<endl;
	}*/
	cout<<"\n"<<L<<endl;
}
float density_austenitic = 7850;//Austenite 308 [3D ray tracing in Austenite materials, Fig. 3]
//ABOVE VALUE CHANGED FROM 8100 TO 7850 BECAUSE 7850 IS WHAT IS MENTIONED IN THE TABLE FOR TRANSVERSELY ISOTROPIC MATERIAL CORRESPONDING TO THE ABOVE VALUES OF ELASTIC CONSTANT USED.
float DENSITY = density_austenitic;



/*********************************************************
IMPORTANT DEFINITIONS
Note: Lengths are in mm.
**********************************************************/
float weld_t = 20;
float weld_h = 20;
float OD = 762;
float ID = 737;

float OR = OD/2;
float IR = ID/2;

float pipe_thickness = OR - IR;

float PLATE_HEIGHT = 80;//60
float PLATE_WIDTH = 300;
float D = 25;//4	//Width of weld at bottom
float A = 30;//43	//Full angle of weld
float D2 = 2*(D/2+PLATE_HEIGHT*tan(degtorad(A/2))); //Width of weld at top
float m1 = tan(degtorad(90-A/2));
float tolerance = 0.001f;
int DUMMY = 0;
int ren_count = 1;
int SHEAR = 0;
int temp_count = 0;
//float WELD_NORMAL = degtorad(A/2);

//float weld_angle = asin(weld_t/(2*OR));
//float weld_angle_bottom = asin(weld_t/(2*IR));

//float pipe_length = 600;
//float offset_plane = pipe_length/2;

//float cyw, cywbelow, weld_r;

//Graph Parameters

int gheight = 50;
int gwidth = 300;
int offset = -20;
int han_freq = 50; //Hz

//Wedge Parameters
float wedge_angle = degtorad(30.0);
float wedge_height = 6.0f;
float wedge_width = 20.0f;
float isotopic_angle = 0.0f;

float phase_diff = 0.0f;
int scale = 20;
float scan_diff = 0.0f;

vector <probe> p;
vector <probe> pL;
vector <probe> pR;

point pulse[2000];

bool dispElastic = 0;
bool dispElastic2 = 0;
int global_count_temp = 0;
bool first_check = 1;
ray theRay, theRay_initial, theRay_temp;
int l = 0;

vector <ray> rayList;
vector < vector <vertex_type> > rayPaths;
int refl_check = 1;
float ray_x = 0.0f;
float ray_y = 80.0f;                                                               //Ray Starting Y coordinate
float ray_z = 79.9f + wedge_height + (tan(wedge_angle)*(wedge_width/2));           //Ray Starting Z coordinate
float ray_angle = degtorad(270) - wedge_angle;         //Ray Starting Angle - Measured from positive Y axis(Horizontal axis here) 
float phase_velocity = 2730; //           //5700 -  Ray Starting initial (outer) phase velocity. - Velocity in 100% Ferrite. 2730 - Velocity in Probe.
int ray_type = 1;                        //0 -> longitudinal wave | 1 -> Shear 1 Wave | 2 -> Shear 2 Wave
float step_size = 0.1;                   //Length ofg steps with which the ray propogates - Smaller values will cause the ray to take smaller steps
float tot_time_weld = 0;
float tot_time_iso_1 = 0;
float tot_time_iso_2 = 0;
int tot_count = 0;
float new_angle = ray_angle;
int weld_entry = 0;
int vel_id = 0;
double prev_grain_orientation = 0;
float ray_heading_temp = 0.0f;
float ray_heading_temp_group = 0.0f;
float external_phase_velocity = phase_velocity;

// Necessary for the calculation of amplitude
float amplitude = 1.0f; float total_time = 0.0f;
float incident_angle = 0.0f, refracted_angle = 0.0f, reflected_angle = 0.0f ;  // Measured against the normal to the interface
float transmission_coefficient = 1.0f;                                       // Multiplied at each interface to calculate the final amplitude
float incident_velocity = 0.0f, reflected_velocity = 0.0f, refracted_velocity = 0.0f;
float dummy_angle = 0.0f;


//New definitions
//Matrix<double, 2, Dynamic, RowMajor> s(2, 7202);
//vector <double> X, Y;
//int check_master = 0;

/*************************************************************************
IMPORTANT DEFINITIONS - WELD PARAMETERS
From G.D. Conolly thesis - -5.0<eta<5.0 and 0.5<T<10.0
But for illustrutive purposes -1.0<eta<1.0 is considered
as other models will produce too little variation in elastic orientations 
**************************************************************************/

double T= 1.3f;  //1.3
double eta = 0.6f; //0.6

/*****************************************************************
IMPORTANT DEFINITIONS - WELD PARAMETERS FOR REVERSE TRACING MODEL
******************************************************************/
float eta_temp = 0.0f, T_temp = 0.5f, eta_final[100], T_final[100], y_final[100];
float eta_start = 0.0f, eta_end = 1.0f;
float T_start = 0.5f, T_end = 10.0f;
float T_graph[1200], eta_graph[1200];
float y_end[1200];

/*****************************************************************
IMPORTANT DEFINITIONS - DEFECT PARAMETERS
******************************************************************/
int defect_enable = 0;
float defect_y = 0.0f;
float defect_z = 15.0f;
float defect_angle = degtorad(80);
float defect_length = 10.0f;
float defect_y2 = defect_y + (defect_length/2)*cos(defect_angle);
float defect_z2 = defect_z + (defect_length/2)*sin(defect_angle);
float defect_y1 = defect_y - (defect_length/2)*cos(defect_angle);
float defect_z1 = defect_z - (defect_length/2)*sin(defect_angle);
float dely = defect_y2 - defect_y1;
float delz = defect_z2 - defect_z1;
float ty = 0.0f;
float tz = 0.0f;
int rotflag_temp = 0;
int num_defect = 0;            //Refers to number of times a defect is encountered. 
int defect_count = 0;
ray newray_defect[10];
int temp_first = 0;
const int NUM_EXTRA_PLATES = 0;
double multiplePlateElasticOrientation[] = {-0.004,0.004,0};
//multiplePlateElasticOrientation[0] = 0;
//multiplePlateElasticOrientation[1] = 40;
int v = 0;


double m2o = 0, m3o;


double VoigtC(MatrixXd const C, int i, int j, int k, int l)
{
	int m, n;
	if(i == 1 && j == 1)    m = 1;
	else if (i == 2 && j == 2)  m = 2;
	else if (i == 3 && j == 3)  m = 3;
	else if ((i == 2 && j == 3) || (j == 2 && i == 3))  m = 4;
	else if ((i == 1 && j == 3) || (j == 1 && i == 3))  m = 5;
	else if ((i == 2 && j == 1) || (j == 2 && i == 1))  m = 6;
	if(k == 1 && l == 1)    n = 1;
	else if (k == 2 && l == 2)  n = 2;
	else if (k == 3 && l == 3)  n = 3;
	else if ((k == 2 && l == 3) || (l == 2 && k == 3))  n = 4;
	else if ((k == 1 && l == 3) || (l == 1 && k == 3))  n = 5;
	else if ((k == 2 && l == 1) || (l == 2 && k == 1))  n = 6;
	return C(m-1,n-1);
}

double getAbsoluteAngle(double y, double z)
{
	double angle = atan2(z, y);
	if (angle < 0)
		angle += degtorad(360);
	return angle;
}



void laguer(vector < complex <double> > &a, complex<double> &x, int &its)
{
	const int MR = 8, MT = 10, MAXIT = MT*MR;
	const double EPS = numeric_limits<double>::epsilon();
	static const double frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
	complex<double> dx, x1, b, d, f, g, h, sq, gp, gm, g2;
	int m = a.size() - 1;
	for (int iter=1;iter<=MAXIT;++iter) 
	{
		b=a[m];
		double err=abs(b);
		d=f=0.0; 
		double abx=abs(x);
		its = iter;
		for (int j=m-1;j>=0;--j) 
		{
			f=x*f+d;
			d=x*d+b;
			b=x*b+a[j];
			err=abs(b)+abx*err;
			//            Efficient computation of the polynomial and its first two derivatives. f stores P00=2.
		}
		err *= EPS;
		
		//        Estimate of roundoff error in evaluating polynomial.
		if (abs(b) <= err) return;
		g=d/b;
		g2=g*g;
		h=g2-2.0*f/b;
		sq=sqrt(double(m-1)*(double(m)*h-g2));
		gp=g+sq;
		gm=g-sq;
		double abp=abs(gp);
		double abm=abs(gm);
		if (abp < abm) gp=gm;
		dx=max(abp,abm) > 0.0 ? double(m)/gp : polar(1+abx,double(iter));
		x1=x-dx;
		if (x == x1) return;
		if (iter % MT != 0) x=x1;
		else x -= frac[iter/MT]*dx;
	}
	throw("too many iterations in laguer");
}


void zroots(vector < complex <double> > &a, vector <complex<double> > &roots, const bool &polish)
{
	const double EPS=1.0e-14;
	int i,its;
	complex<double> x,b,c;
	int m=a.size()-1;
	vector < complex <double> > ad(m+1);
	for (int j=0;j<=m;j++) ad[j]=a[j];
	for (int j=m-1;j>=0;--j) 
	{
		x=0.0;
		vector< complex <double> > ad_v(j+2);
		for (int jj=0;jj<j+2;++jj) ad_v[jj]=ad[jj]; laguer(ad_v,x,its);
		if (abs(imag(x)) <= 2.0*EPS*abs(real(x)))
			x= complex<double>(real(x),0.0);
		roots[j]=x;
		b=ad[j+1];
		for (int jj=j;jj>=0;--jj) 
		{
			c=ad[jj];
			ad[jj] = b;
			b = x*b+c;
		}
	}
	if (polish)
		for (int j=0;j<m;++j)
			laguer(a,roots[j],its);
	for(int j=1;j<m;++j) 
	{
		x=roots[j];
		for (i=j-1;i>=0;--i) 
		{
			if (real(roots[i]) <= real(x)) break;
			roots[i+1]=roots[i];
		}
		roots[i+1]=x;
	}
}

bool isOpposite(double a, double b)
{
	return (((a<0)?-1:1) * ((b<0)?-1:1) <0) ? 1: 0;
}

bool isInsideWeld(vertex_type v)
{
	if ((v.z <= PLATE_HEIGHT) && (v.y >= 1.005*(-D / 2 - tan(degtorad(A / 2))*v.z)) && (v.y <= 1.005*(D / 2 + tan(degtorad(A / 2))*v.z)))
	{
		return 1;
	}
	else 
	{
		return 0;
	}
}

bool isAtDefect(vertex_type v)
{
	//cout<<endl<<dely<<" "<<delz<<" "<<defect_y1<<" "<<defect_z1;
	ty = (v.y - defect_y1)/(dely);
	tz = (v.z - defect_z1)/(delz);
	//cout<<endl<<ty<<" "<<tz<<endl;
	if(abs(ty-tz)<0.04 && (v.z > defect_z1) && (v.z < defect_z2))
	{
		cout<<"ty = "<<ty<<" tz = "<<tz;
		return 1;
	}
	return 0;	
}


double getCrystalOrientation(double y, double z, double T  , double eta  )
{
	double tantheta;
	
	if(y > 0)
	{
		tantheta = -atan(1.0f*T*(D/2+z*tan(degtorad(A/2)))/pow(y,eta)) ;  //Changed PI to 2*PI
	}
	else if(y < 0)
	{
		tantheta = atan(1.0f*T*(D/2+z*tan(degtorad(A/2)))/pow(-y,eta))  ;
	}
	else if(y == 0)
	{
		return degtorad(90);
	}

	return tantheta;
	// Returns in radians.
}

#endif