#include <iostream>
#include "gnuplot.h"
using namespace std;
int main()
{
	system("D:\\Acads\\Summer-Project\\Anisotropy_Alpha\\Anisotropy_Alpha\\gnuplot\\binary\\gnuplot.exe");
	Gnuplot plot;
	
	

	system("plot sin(x)");
	system("pause");
	plot("plot cos(x)");
	system("pause");
	return 0;
}