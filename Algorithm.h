#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <iostream>
#include <time.h>	

ray ray_prop(ray theRay_initial, float eta_temp, float T_temp)
{
	ray theRay_r, theRay_temp_r;
	//first_check = 0;
	cout << "\nEta and T: " << eta_temp << " " << T_temp;
	theRay_temp_r = plotRayisotropicreverse(theRay_initial);
	//cout<<"\nCHECKING(1):\nReverse iso: "<<theRay_temp_r.current_position.y<<" "<<theRay_temp_r.current_position.z<<" "<<radtodeg(theRay_temp_r.global_direction);

	theRay_temp_r = plotRayWeldreverse(theRay_temp_r, eta_temp, T_temp);
	//cout<<"\nCHECKING(1):\nReverse After weld: "<<theRay_temp_r.current_position.y<<" "<<theRay_temp_r.current_position.z<<" "<<radtodeg(theRay_temp_r.global_direction);
	theRay_r = plotRayisotropicreverse_exit(theRay_temp_r);
	cout << "\nCHECKING(1):Reverse final: " << theRay_r.current_position.y << " " << theRay_r.current_position.z << " " << radtodeg(theRay_r.global_direction);
	return theRay_r;
}

void graphprint(ray theRay_initial)
{
	//float T_temp = 0.5f, eta_temp = 0.0f;
	ray theRay_r;
	clock_t tStart = clock();
	ofstream myreadtext;


	myreadtext.open("Graph_text.txt");
	myreadtext << "Y = 60 and Angle = 220: \n\n";
	//ray_y = 60.0f;
	//ray_angle = degtorad(220);
	T_temp = 0.5f;
	int i = 0;
	while (i < 96)
	{

		//myreadtext<<"\n\n\nT: "<<T_graph[i]<<"\tEta: 0-1"<<"\n";
		eta_temp = 0.0f;
		int j = 0;
		while (j < 11)
		{
			cout << "\n\t\tCount : " << (i * 11) + j;

			theRay_r = ray_prop(theRay_initial, eta_temp, T_temp);
			y_end[(i * 11) + j] = theRay_r.current_position.y;
			cout << "\nY value: " << theRay_r.current_position.y;
			eta_graph[(i * 11) + j] = eta_temp;
			T_graph[(i * 11) + j] = T_temp;
			eta_temp = eta_temp + 0.1;
			j++;
			//myreadtext<<eta_graph[j]<<","<<y_end[j]<<"\n";
		}
		T_temp = T_temp + 0.1;
		i++;
	}

	myreadtext << "T: ";
	for (int i = 0; i < 1056; i++)
		myreadtext << T_graph[i] << ", ";
	myreadtext << "\n\neta: ";
	for (int i = 0; i < 1056; i++)
		myreadtext << eta_graph[i] << ", ";
	myreadtext << "\n\nY: ";
	for (int i = 0; i < 1056; i++)
		myreadtext << y_end[i] << ", ";
	myreadtext << "\n\n";

	/*
	ray_y = 75.0f;
	ray_angle = degtorad(230);

	myreadtext<<"Y = 75 and Angle = 230: \n\n";
	T_temp = 0.5f;
	i = 0;
	while(i<96)
	{

	//myreadtext<<"\n\n\nT: "<<T_graph[i]<<"\tEta: 0-1"<<"\n";
	eta_temp = 0.0f;
	int j =0;
	while(j<11)
	{
	cout<<"\n\t\tCount : "<<(i*11)+j;
	theRay_r = ray_prop(theRay_initial, eta_temp, T_temp);
	y_end[(i*11)+j] = theRay_r.current_position.y;

	eta_graph[(i*11)+j] = eta_temp;
	T_graph[(i*11)+j] = T_temp;
	eta_temp = eta_temp + 0.1;
	j++;
	//myreadtext<<eta_graph[j]<<","<<y_end[j]<<"\n";
	}
	T_temp = T_temp + 0.1;
	i++;
	}
	myreadtext<<"Y = 65 and Angle = 235: \n\n";
	myreadtext<<"T: ";
	for(int i=0;i<1056;i++)
	myreadtext<<T_graph[i]<<", ";
	myreadtext<<"\n\neta: ";
	for(int i=0;i<1056;i++)
	myreadtext<<eta_graph[i]<<", ";
	myreadtext<<"\n\nY: ";
	for(int i=0;i<1056;i++)
	myreadtext<<y_end[i]<<", ";
	myreadtext<<"\n\n";


	while(i<100)
	{

	//myreadtext<<"Eta\tT\ty-Endpoint\n";
	myreadtext<<T_graph[i]<<"\t"<<y_end[i]<<"\n";
	i++;
	}
	i=0;
	*/
	myreadtext.close();
	cout << "\n\nTime taken: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << " seconds";
}
void bruteforce(ray theRay_initial)
{
	//float T_temp = 0.5f, eta_temp = 0.0f;
	ray theRay_r;
	clock_t tStart = clock();

	int q = 0;
	int i = 0;
	T_temp = 0.5f;
	eta_temp = 0.0f;
	//cout<<"\nChumma: "<<theRay.current_position.y<<" "<<theRay.current_position.z<<" "<<radtodeg(theRay.global_direction);



	//myreadtext<<"\n\n\nT: "<<T_graph[i]<<"\tEta: 0-1"<<"\n";
	eta_temp = 0.0f;
	int j = 0;

	while (i < 96) //96
	{
		eta_temp = -2.0f;
		int j = 0;
		while (j < 41)
		{
			cout << "\n\t\tCount : " << (i * 41) + j + 1;
//			theRay_r = ray_prop(theRay_initial, eta_temp, T_temp);
			/*
			theRay_r_temp = plotRayisotropic(theRay_initial);
			//cout<<"\nCHECKING(3):\nActual after iso: "<<theRay_temp.current_position.y<<" "<<theRay_temp.current_position.z<<" "<<radtodeg(theRay_temp.global_direction);
			theRay_r_temp = plotRayWeld2(theRay_r_temp,eta_temp, T_temp);
			//cout<<"\nCHECKING(3):\nActual after weld: "<<theRay_temp.current_position.y<<" "<<theRay_temp.current_position.z<<" "<<radtodeg(theRay_temp.global_direction);
			theRay_r = plotRayisotropic_exit(theRay_r_temp);
			cout<<"\nCHECKING(3):Actual final: "<<theRay_r.current_position.y<<" "<<theRay_r.current_position.z<<" "<<radtodeg(theRay_r.global_direction);
			*/
			if (abs(theRay_r.current_position.y) < 1.00005*abs(theRay.current_position.y) && abs(theRay_r.current_position.y) > 0.99995*abs(theRay.current_position.y))
			{
				T_final[q] = T_temp;
				eta_final[q] = eta_temp;
				y_final[q] = theRay_r.current_position.y;
				q++;
			}

			eta_temp = eta_temp + 0.1;
			j++;
			//myreadtext<<eta_graph[j]<<","<<y_end[j]<<"\n";
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
			cout << "\nEta and T: " << eta_final[i] << " " << T_final[i];

	}
	cout << "\n\nTime taken: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << " seconds";

}

#endif