/*
 * CubicHermiteSpline.cpp
 *
 *  Created on: Oct 15, 2018
 *      Author: leeja
 */

#include <iostream>

double* generateEquation(double p0, double p1, double v0, double v1) {
	double term0 = 2 * p0 - 2 * p1 + v0 + v1;
	double term1 = -3 * p0 + 3 * p1 - 2 * v0 - v1;
	double term2 = v0;
	double term3 = p0;

	static double equation[4];
	equation[0] = term0;
	equation[1] = term1;
	equation[2] = term2;
	equation[3] = term3;
	return equation;
}

int main() {
	double x0 = 0;
	double x1 = 100;
	double vx0 = 1;
	double vx1 = 0;

	double *xEquation = generateEquation(x0, x1, vx0, vx1);
	std::cout << xEquation[0] << " " << xEquation[1] << " " << xEquation[2] << " " << xEquation[3] << std::endl;

	double y0 = 0;
	double y1 = 100;
	double vy0 = 0;
	double vy1 = 1;

	double *yEquation = generateEquation(y0, y1, vy0, vy1);
	std::cout << yEquation[0] << " " << yEquation[1] << " " << yEquation[2] << " " << yEquation[3] << std::endl;

	return 0;
}
