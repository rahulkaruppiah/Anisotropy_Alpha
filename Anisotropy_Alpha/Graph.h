#define _CRT_SECURE_NO_DEPRECATE 
#include <fstream>
#include<iostream>
#include <stdlib.h>
#include <stdio.h>
#include "Definitions.h"

using namespace std;









void plotResults()   //double* xData, double* yData, int dataSize
{
	FILE *gnuplotPipe, *tempDataFile;
	FILE * pFile;
	int output1[136];
	float f, g;
	ifstream myReadFile;
	pFile = fopen("PulseFinal.dat", "w+");
	int h = 0;
	myReadFile.open("Pulse.txt");
	//myWriteFile.open("PulseFinal.dat");
	float output[135];
	while (!myReadFile.eof())
	{
		myReadFile >> output[h];
		h++;
	}
	h = 0;

	while (h<135)
	{

		pulse[h].x = h + 200;
		pulse[h].y = output[h];
		fprintf(pFile, "%f %f\n", pulse[h].x, pulse[h].y);
		//cout<<pulse[h].x<<" "<<pulse[h].y<<"\n";

		//cout<<output[h]<<"\n";
		h++;
	}
	/*
	rewind(pFile);
	h = 0;

	while (h<135)
	{
	fscanf(pFile, "%f %f", &f, &g);
	printf("I have read: %f %f\n", f, g);
	h++;
	}

	*/
	myReadFile.close();
	fclose(pFile);



	char *tempDataFileName;
	char *Datafile;
	double x, y;
	int i;
	tempDataFileName = "Pulse.txt";
	Datafile = "PulseFinal.dat";
	gnuplotPipe = _popen("gnuplot", "w");

	// set axis ranges
	fprintf(gnuplotPipe, "set xrange [0:400]\n");
	fprintf(gnuplotPipe, "set yrange [-3:3]\n");

	if (gnuplotPipe)
	{
		fprintf(gnuplotPipe, "plot \"%s\" with lines, \"%s\" with lines\n", tempDataFileName, Datafile);
		//	fflush(gnuplotPipe);
		//	fprintf(gnuplotPipe, "plot \"%s\" with lines\n", Datafile);

		fflush(gnuplotPipe);


		/*
		tempDataFile = fopen(tempDataFileName, "w");
		for (i = 0; i <= dataSize; i++)
		{
		x = xData[i];
		y = yData[i];
		fprintf(tempDataFile, "%lf %lf\n", x, y);
		}
		fclose(tempDataFile);
		*/

		printf("press enter to continue...");

		//    getchar();
		remove(tempDataFileName);
		fprintf(gnuplotPipe, "exit \n");
	}
	else
	{
		printf("gnuplot not found...");
	}

}

void gplot()
{
	int i = 0;
	int nIntervals = 100;
	double intervalSize = 1.0;
	double stepSize = intervalSize / nIntervals;
	double* xData = (double*)malloc((nIntervals + 1)*sizeof(double));
	double* yData = (double*)malloc((nIntervals + 1)*sizeof(double));
	xData[0] = 0.0;
	double x0 = 0.0;
	for (i = 0; i < nIntervals; i++)
	{
		x0 = xData[i];
		xData[i + 1] = x0 + stepSize;
	}
	for (i = 0; i <= nIntervals; i++) {
		x0 = xData[i];
		yData[i] = sin(x0)*cos(10 * x0);
	}
	//plotResults(xData, yData, nIntervals);

}
