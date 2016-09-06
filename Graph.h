#ifndef GRAPH_H
#define GRAPH_H

#define _CRT_SECURE_NO_DEPRECATE 
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable:4996)

#include <fstream>
#include<iostream>
#include "Definitions.h"
#include <random>



void plotResults()   
{






	// Part 1 - Data Extraction and Inputting the resulting A Scan data into a file 

	FILE * pFile, * nFile;
	//float f, g;
	ifstream myReadFile;
	ofstream myWriteFile;
	pFile = fopen("PulseFinal.dat", "w+");
	nFile = fopen("NoisePulse.dat", "w+");
	int h = 0;
	myReadFile.open("Pulse.txt");
	float output[135];

		while (!myReadFile.eof())
	{
		myReadFile >> output[h];
		++h;
	}



	while (!myReadFile.eof())
	{
		myReadFile >> output[h];
		++h;
	}
	myReadFile.close();
	myWriteFile.open("Pulse.txt");


	//Noise
	int const Dim = 135;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution<> d;
	float a[Dim];
	int i = 0;
	float noise;
	while (i<h)
	{
		noise = d(gen);
		a[i] = output[i] + noise;
		cout << "Noise: " << noise << endl;
		fprintf(nFile, "%f %f\n", i, a[i]);
		++i;
	}
	scan_diff = han_freq * (tot_time_weld + tot_time_iso_1 + tot_time_iso_2) * 100;
	i = 0;
	while (i<h)
	{
		pulse[i].x = i + scan_diff;
	//	cout << "\nPulse.x: " << pulse[i].x;

		pulse[i].y = output[i] * amplitude;
		fprintf(pFile, "%f %f\n", pulse[i].x, pulse[i].y);
		++i;
	}
	
	fclose(pFile);
	// Part 2 - Plotting using gnuplot happens here
	FILE *gnuplotPipe;
	char *tempDataFileName, *Datafile;
	char *NoiseDataFile = "NoisePulse.dat";

	//double x, y;
	//int i;
	tempDataFileName = "Pulse.txt";
	Datafile = "PulseFinal.dat";
	gnuplotPipe = _popen("gnuplot", "w");

	//closes any active windows - does not work. Check later.
    //fprintf(gnuplotPipe, "reset\n");
	// set axis ranges
	fprintf(gnuplotPipe, "set xrange [0: \"%f\"]\n", (pulse[h - 2].x+50));
	fprintf(gnuplotPipe, "set yrange [-3:3]\n");
	fprintf(gnuplotPipe, "clear\n");

	if (gnuplotPipe)
	{
		fprintf(gnuplotPipe, "plot \"%s\" with lines, \"%s\" with lines, \"%s\" with lines\n", tempDataFileName, NoiseDataFile, Datafile);
		fflush(gnuplotPipe);
		printf("\nA-Scan has successfully plotted.\n");
		remove(tempDataFileName);
		fprintf(gnuplotPipe, "exit \n");
	}
	else
	{
		printf("\nError Loading gnuplot\n");
	}
}

#endif