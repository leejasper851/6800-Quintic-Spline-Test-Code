/*
 * QuinticHermiteSpline.cpp
 *
 *  Created on: Oct 16, 2018
 *      Author: leeja
 */

// for std::array
#include <array>
// for pow and sqrt
#include <cmath>
// for std::cout
#include <iostream>
// for std::setprecision
#include <iomanip>

// this class contains a quintic spline and robot movements related to it
class QuinticSpline {
	// this is half of the width between the left and right wheels
	static constexpr double m_HALF_DRIVE_WIDTH { 0.2 };
	// these are the positions, velocities, and accelerations of the parametric equations
	double m_px0, m_px1, m_vx0, m_vx1, m_ax0, m_ax1;
	double m_py0, m_py1, m_vy0, m_vy1, m_ay0, m_ay1;
	// these are the parametric equations and their derivatives
	std::array<double, 6> m_xEquation { };
	std::array<double, 6> m_yEquation { };
	std::array<double, 6> m_xDerEquation { };
	std::array<double, 6> m_yDerEquation { };
	std::array<double, 6> m_xDoubleDerEquation { };
	std::array<double, 6> m_yDoubleDerEquation { };

public:
	// the constructor takes in the positions, velocities, and accelerations, and generates equations based on them
	QuinticSpline(double px0, double px1, double vx0, double vx1, double ax0, double ax1, double py0, double py1, double vy0, double vy1, double ay0, double ay1) {
		m_px0 = px0;
		m_px1 = px1;
		m_vx0 = vx0;
		m_vx1 = vx1;
		m_ax0 = ax0;
		m_ax1 = ax1;
		m_py0 = py0;
		m_py1 = py1;
		m_vy0 = vy0;
		m_vy1 = vy1;
		m_ay0 = ay0;
		m_ay1 = ay1;
		generateEquations();
		generateDerEquations();
		generateDoubleDerEquations();
	}

private:
	// generates the parametric equations
	void generateEquations() {
		double term0 = -6 * m_px0 - 3 * m_vx0 - 0.5 * m_ax0 + 0.5 * m_ax1 - 3 * m_vx1 + 6 * m_px1;
		double term1 = 15 * m_px0 + 8 * m_vx0 + 1.5 * m_ax0 - m_ax1 + 7 * m_vx1 - 15 * m_px1;
		double term2 = -10 * m_px0 - 6 * m_vx0 - 1.5 * m_ax0 + 0.5 * m_ax1 - 4 * m_vx1 + 10 * m_px1;
		double term3 = 0.5 * m_ax0;
		double term4 = m_vx0;
		double term5 = m_px0;
		m_xEquation = { term0, term1, term2, term3, term4, term5 };

		term0 = -6 * m_py0 - 3 * m_vy0 - 0.5 * m_ay0 + 0.5 * m_ay1 - 3 * m_vy1 + 6 * m_py1;
		term1 = 15 * m_py0 + 8 * m_vy0 + 1.5 * m_ay0 - m_ay1 + 7 * m_vy1 - 15 * m_py1;
		term2 = -10 * m_py0 - 6 * m_vy0 - 1.5 * m_ay0 + 0.5 * m_ay1 - 4 * m_vy1 + 10 * m_py1;
		term3 = 0.5 * m_ay0;
		term4 = m_vy0;
		term5 = m_py0;
		m_yEquation = { term0, term1, term2, term3, term4, term5 };
	}

	// generates the parametric derivative equations
	void generateDerEquations() {
		double term0 = 0;
		double term1 = 5 * m_xEquation.at(0);
		double term2 = 4 * m_xEquation.at(1);
		double term3 = 3 * m_xEquation.at(2);
		double term4 = 2 * m_xEquation.at(3);
		double term5 = m_xEquation.at(4);
		m_xDerEquation = { term0, term1, term2, term3, term4, term5 };

		term0 = 0;
		term1 = 5 * m_yEquation.at(0);
		term2 = 4 * m_yEquation.at(1);
		term3 = 3 * m_yEquation.at(2);
		term4 = 2 * m_yEquation.at(3);
		term5 = m_yEquation.at(4);
		m_yDerEquation = { term0, term1, term2, term3, term4, term5 };
	}

	// generates the parametric double derivative equations
	void generateDoubleDerEquations() {
		double term0 = 0;
		double term1 = 0;
		double term2 = 4 * m_xDerEquation.at(1);
		double term3 = 3 * m_xDerEquation.at(2);
		double term4 = 2 * m_xDerEquation.at(3);
		double term5 = m_xDerEquation.at(4);
		m_xDoubleDerEquation = { term0, term1, term2, term3, term4, term5 };

		term0 = 0;
		term1 = 0;
		term2 = 4 * m_yDerEquation.at(1);
		term3 = 3 * m_yDerEquation.at(2);
		term4 = 2 * m_yDerEquation.at(3);
		term5 = m_yDerEquation.at(4);
		m_yDoubleDerEquation = { term0, term1, term2, term3, term4, term5 };
	}

	// used by getArcLen() to get each value for the Reimann sum (calculates sqrt((dx/dt)^2 + (dy/dt)^2))
	double getArcLenVal(bool isAdding, double t) {
		double xPart1 = plugInEquation(1, true, t);
		double xPart2 = 0.5 * pow(pow(m_HALF_DRIVE_WIDTH, 2) / (1 + pow(plugInEquation(1, true, t) / plugInEquation(1, false, t), 2)), -0.5);
		double xPart3 = -1 * pow(m_HALF_DRIVE_WIDTH, 2) * pow(1 + pow(plugInEquation(1, true, t) / plugInEquation(1, false, t), 2), -2);
		double xPart4 = 2 * (plugInEquation(1, true, t) / plugInEquation(1, false, t));
		double xPart5 = (plugInEquation(2, true, t) * plugInEquation(1, false, t) - plugInEquation(1, true, t) * plugInEquation(2, false, t)) / pow(plugInEquation(1, false, t), 2);
		double finalX = 0;
		if (isAdding) {
			finalX = xPart1 + xPart2 * xPart3 * xPart4 * xPart5;
		} else {
			finalX = xPart1 - xPart2 * xPart3 * xPart4 * xPart5;
		}

		double yPart1 = plugInEquation(1, false, t);
		double yPart2 = -1 * xPart5;
		double yPart3 = sqrt(pow(m_HALF_DRIVE_WIDTH, 2) / (1 + pow(plugInEquation(1, true, t) / plugInEquation(1, false, t), 2)));
		double yPart4 = -1 * (plugInEquation(1, true, t) / plugInEquation(1, false, t));
		double yPart5 = xPart2 * xPart3 * xPart4 * xPart5;
		double finalY = 0;
		if (isAdding) {
			finalY = yPart1 + (yPart2 * yPart3 + yPart4 * yPart5);
		} else {
			finalY = yPart1 - (yPart2 * yPart3 + yPart4 * yPart5);
		}

		return sqrt(pow(finalX, 2) + pow(finalY, 2));
	}

public:

	// getter for half of the width between the left and right wheels
	double getHalfDriveWidth() {
		return m_HALF_DRIVE_WIDTH;
	}

	// getter for the original parametric values
	double getParametricValue(int der, bool is0, bool isX) {
		if (der == 0) {
			if (is0) {
				if (isX) {
					return m_px0;
				}
				return m_py0;
			}
			if (isX) {
				return m_px1;
			}
			return m_py1;
		}
		if (der == 1) {
			if (is0) {
				if (isX) {
					return m_vx0;
				}
				return m_vy0;
			}
			if (isX) {
				return m_vx1;
			}
			return m_vy1;
		}
		if (is0) {
			if (isX) {
				return m_ax0;
			}
			return m_ay0;
		}
		if (isX) {
			return m_ax1;
		}
		return m_ay1;
	}

	// getter for the equations
	std::array<double, 6> getEquation(int der, bool isX) {
		if (der == 0) {
			if (isX) {
				return m_xEquation;
			}
			return m_yEquation;
		}
		if (der == 1) {
			if (isX) {
				return m_xDerEquation;
			}
			return m_yDerEquation;
		}
		if (isX) {
			return m_xDoubleDerEquation;
		}
		return m_yDoubleDerEquation;
	}

	// plug in values of t into the parametric equations
	double plugInEquation(int der, bool isX, double t) {
		std::array<double, 6> equation = getEquation(der, isX);
		double term0 = equation.at(0) * pow(t, 5);
		double term1 = equation.at(1) * pow(t, 4);
		double term2 = equation.at(2) * pow(t, 3);
		double term3 = equation.at(3) * pow(t, 2);
		double term4 = equation.at(4) * t;
		double term5 = equation.at(5);
		return term0 + term1 + term2 + term3 + term4 + term5;
	}

	// returns the slope of the spline given t
	double getSlope(double t) {
		if (t == 0) {
			return m_vy0 / m_vx0;
		}
		if (t == 1) {
			return m_vy1 / m_vx1;
		}
		return plugInEquation(1, false, t) / plugInEquation(1, true, t);
	}

	// returns the double derivative at time t (d^2(y)/dx^2)
	double getDoubleDer(double t) {
		if (t == 0) {
			return (m_vx0 * m_ay0 - m_ax0 * m_vy0) / pow(m_vx0, 3);
		}
		if (t == 1) {
			return (m_vx1 * m_ay1 - m_ax1 * m_vy1) / pow(m_vx1, 3);
		}
		return (plugInEquation(2, false, t) * plugInEquation(1, true, t) - plugInEquation(2, true, t) * plugInEquation(1, false, t))/pow(plugInEquation(1, true, t), 3);
	}

	// returns the normal slope of the spline at time t
	double getNormalSlope(double t) {
		return -1 / getSlope(t);
	}

	// gets the difference in x of the wheel positions from the original spline at time t
	double getWheelX(double t) {
		double slope = getNormalSlope(t);
		return sqrt(pow(m_HALF_DRIVE_WIDTH, 2) / (1 + pow(slope, 2)));
	}

	// gets the difference in y of the wheel positions from the original spline at time t
	double getWheelY(double t) {
		return getNormalSlope(t) * getWheelX(t);
	}

	// gets the arc length of a portion of the wheel paths (use isAdding to specify which wheel path)
	double getArcLen(bool isAdding, double t1, double t2) {
		constexpr int NUM_STEPS = 50;
		double step = (t2 - t1) / NUM_STEPS;
		if (t1 == 0) {
			t1 = step / 1000000;
			step = (t2 - t1) / NUM_STEPS;
		}
		std::array<double, NUM_STEPS+1> arcLenVals { };
		for (int i = 0; i <= NUM_STEPS; i++) {
			arcLenVals.at(i) = getArcLenVal(isAdding, t1 + i * step);
		}
		double arcLength = 0;
		for (int i = 1; i <= NUM_STEPS; i++) {
			arcLength += ((arcLenVals.at(i-1) + arcLenVals.at(i)) / 2) * step;
		}
		return arcLength;
	}
};

int main() {
	double px0 = 0;
	double px1 = 1;
	double vx0 = 1;
	double vx1 = 1;
	double ax0 = 0;
	double ax1 = 0;

	double py0 = 0;
	double py1 = 1;
	double vy0 = 0;
	double vy1 = 0;
	double ay0 = 0;
	double ay1 = 0;

	QuinticSpline quinticSpline(px0, px1, vx0, vx1, ax0, ax1, py0, py1, vy0, vy1, ay0, ay1);

	std::array<double, 6> xEquation = quinticSpline.getEquation(0, true);
	std::cout << "f(t) = " << xEquation.at(0) << "t^5 + " << xEquation.at(1) << "t^4 + " << xEquation.at(2) << "t^3 + " << xEquation.at(3) << "t^2 + " << xEquation.at(4) << "t + " << xEquation.at(5) << '\n';

	std::array<double, 6> yEquation = quinticSpline.getEquation(0, false);
	std::cout << "g(t) = " << yEquation.at(0) << "t^5 + " << yEquation.at(1) << "t^4 + " << yEquation.at(2) << "t^3 + " << yEquation.at(3) << "t^2 + " << yEquation.at(4) << "t + " << yEquation.at(5) << '\n';

	std::cout << "position 0: (" << quinticSpline.getParametricValue(0, true, true) << ", " << quinticSpline.getParametricValue(0, true, false) << ")   ";
	std::cout << "position 1: (" << quinticSpline.getParametricValue(0, false, true) << ", " << quinticSpline.getParametricValue(0, false, false) << ")\n";

	std::cout << "1st derivative 0: " << quinticSpline.getSlope(0) << "   ";
	std::cout << "1st derivative 1: " << quinticSpline.getSlope(1) << '\n';

	std::cout << "2nd derivative 0: " << quinticSpline.getDoubleDer(0) << "   ";
	std::cout << "2nd derivative 1: " << quinticSpline.getDoubleDer(1) << '\n';

	std::array<double, 6> xDerEquation = quinticSpline.getEquation(1, true);
	std::cout << "x'(t) = " << xDerEquation.at(0) << ' ' << xDerEquation.at(1) << ' ' << xDerEquation.at(2) << ' ' << xDerEquation.at(3) << ' ' << xDerEquation.at(4) << ' ' << xDerEquation.at(5) << '\n';

	std::array<double, 6> yDerEquation = quinticSpline.getEquation(1, false);
	std::cout << "y'(t) = " << yDerEquation.at(0) << ' ' << yDerEquation.at(1) << ' ' << yDerEquation.at(2) << ' ' << yDerEquation.at(3) << ' ' << yDerEquation.at(4) << ' ' << yDerEquation.at(5) << '\n';

	constexpr int NUM_STEPS = 50;

	std::array<double, NUM_STEPS> velPropL { };
	std::array<double, NUM_STEPS> velPropR { };

	double prevPowerL;
	double prevPowerR;

	double step = 1.0 / NUM_STEPS;
	for (int i = 0; i <= NUM_STEPS; i++) {
		double t = i * step;

		double xPosC = quinticSpline.plugInEquation(0, true, t);
		double yPosC = quinticSpline.plugInEquation(0, false, t);
		double xPosL = xPosC + quinticSpline.getWheelX(t);
		double yPosL = yPosC + quinticSpline.getWheelY(t);
		double xPosR = xPosC - quinticSpline.getWheelX(t);
		double yPosR = yPosC - quinticSpline.getWheelY(t);

		double lenL = quinticSpline.getArcLen(true, t-step, t);
		double lenR = quinticSpline.getArcLen(false, t-step, t);

		if (quinticSpline.getSlope(t) > 0) {
			double temp = xPosL;
			xPosL = xPosR;
			xPosR = temp;

			temp = yPosL;
			yPosL = yPosR;
			yPosR = temp;

			temp = lenL;
			lenL = lenR;
			lenR = temp;
		}

		if (quinticSpline.getSlope(t) == 0) {
			if (quinticSpline.plugInEquation(1, true, t) > 0) {
				yPosL = xPosC + quinticSpline.getHalfDriveWidth();
				yPosR = xPosC - quinticSpline.getHalfDriveWidth();
			} else {
				yPosL = xPosC - quinticSpline.getHalfDriveWidth();
				yPosR = xPosC + quinticSpline.getHalfDriveWidth();
			}
		}

		double powerL = lenL > lenR ? 1 : lenL/lenR;
		double powerR = lenR > lenL ? 1 : lenR/lenL;
		powerL /= 3.58;
		powerR /= 3.58;

		if (t > 0) {
			velPropL[static_cast<int>(t*NUM_STEPS-0.9)] = powerL;
			velPropR[static_cast<int>(t*NUM_STEPS-0.9)] = powerR;
		}

		std::cout << std::fixed << std::setprecision(3);
		std::cout << "t = " << t << ": " << quinticSpline.getSlope(t) << "   ";
		std::cout << "C: (" << xPosC << ", " << yPosC << ") ";
		std::cout << "L: (" << xPosL << ", " << yPosL << ") ";
		std::cout << "R: (" << xPosR << ", " << yPosR << ")";
		if (t != 0) {
			std::cout << "   lenL: " << lenL;
			std::cout << " lenR: " << lenR;

			double powerDiffL = powerL - prevPowerL;
			double powerDiffR = powerR - prevPowerR;
			std::cout << "   powerDiffs: " << powerDiffL << ' ' << powerDiffR;
		}

		std::cout << '\n';

		prevPowerL = powerL;
		prevPowerR = powerR;
	}

	std::cout << "     t: ";
	for (int i = 0; i < NUM_STEPS; i++) {
		std::cout << (i+1) / 50.0;
		if (i != NUM_STEPS-1) {
			std::cout << ", ";
		}
	}
	std::cout << "\npower1: ";
	for (int i = 0; i < NUM_STEPS; i++) {
		std::cout << velPropL.at(i);
		if (i != NUM_STEPS-1) {
			std::cout << ", ";
		}
	}
	std::cout << "\npower2: ";
	for (int i = 0; i < NUM_STEPS; i++) {
		std::cout << velPropR.at(i);
		if (i != NUM_STEPS-1) {
			std::cout << ", ";
		}
	}

	return 0;
}
