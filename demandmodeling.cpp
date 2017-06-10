/**
   Copyright 2017 Shawn Gilroy

   This file is part of Demand Curve Analyzer, Qt port.

   Demand Curve Analyzer is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, version 3.

   Demand Curve Analyzer is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Demand Curve Analyzer.  If not, see http://www.gnu.org/licenses/.

   The Demand Curve Analyzer is a tool to assist researchers in behavior economics.

   Email: shawn(dot)gilroy(at)temple.edu

   ====================================================================================

   ALGLIB 3.11.0 (source code generated 2017-05-11)
   Copyright (c) Sergey Bochkanov (ALGLIB project).

   >>> SOURCE LICENSE >>>
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation (www.fsf.org); either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is available at
   http://www.fsf.org/licensing/licenses
   >>> END OF LICENSE >>>

  */

#include "demandmodeling.h"
#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

void demandmodeling::SetX(const char *mString)
{
    x = mString;
}

void demandmodeling::SetY(const char *mString)
{
    y = mString;
}

void demandmodeling::SetStarts(const char *mString)
{
    c = mString;
}

void demandmodeling::SetBounds(const char *mUpperString, const char *mLowerString)
{
    bndu = mUpperString;
    bndl = mLowerString;
}

real_1d_array demandmodeling::GetParams()
{
    return c;
}

lsfitstate demandmodeling::GetState()
{
    return state;
}

ae_int_t demandmodeling::GetInfo()
{
    return info;
}

lsfitreport demandmodeling::GetReport()
{
    return rep;
}

void exponential_demand(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
	std::vector<double> *param = (std::vector<double> *) ptr;
    double k = param->at(0);
    func = log10(c[0]) + k * (exp(-c[1] * c[0] * x[0]) - 1);
}

void exponential_demand_with_k(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
    func = log10(c[1]) + c[0] * (exp(-c[2] * c[1] * x[0]) - 1);
}

void exponentiated_demand(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
	std::vector<double> *param = (std::vector<double> *) ptr;
    double k = param->at(0);
    func = c[0] * pow(10, (k * (exp(-c[1] * c[0] * x[0]) - 1)));
}

void exponentiated_demand_with_k(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
    func = c[1] * pow(10, (c[0] * (exp(-c[2] * c[1] * x[0]) - 1)));
}

void linear_demand(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
    func = log(c[2]) + (c[1] * log(x[0])) - c[0] * x[0];
}

void demandmodeling::FitLinear(const char *mStarts)
{
    SetStarts(mStarts);

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);
    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, linear_demand);

    lsfitresults(state, info, c, rep);
}

void demandmodeling::FitExponential(const char *mStarts, std::vector<double> mParams)
{
    SetStarts(mStarts);
    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);
    real_1d_array s = "[1.0, 1.0e-3]";
    lsfitsetscale(state, s);
    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, exponential_demand, NULL, &mParams);

    lsfitresults(state, info, c, rep);
}

void demandmodeling::FitExponentialWithK(const char *mStarts)
{
    SetStarts(mStarts);

    std::vector<double> mParams;
    mParams.push_back(3);

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);

    real_1d_array s = "[1.0, 1.0e-3, 1.0]";
    lsfitsetscale(state, s);

    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, exponential_demand_with_k, NULL, &mParams);

    lsfitresults(state, info, c, rep);
}

void demandmodeling::FitExponentiated(const char *mStarts, std::vector<double> mParams)
{
    SetStarts(mStarts);

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);

    real_1d_array s = "[1.0, 1.0e-3]";
    lsfitsetscale(state, s);

    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, exponentiated_demand, NULL, &mParams);

    lsfitresults(state, info, c, rep);
}

void demandmodeling::FitExponentiatedWithK(const char *mStarts)
{
    SetStarts(mStarts);

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);

    real_1d_array s = "[1.0, 1.0e-3, 1.0]";
    lsfitsetscale(state, s);

    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, exponentiated_demand_with_k);

    lsfitresults(state, info, c, rep);
}

void demandmodeling::InitializeDefaults(bool isLog)
{
	maxits = 100;
    epsx = 0.000001;
    diffstep = 0.00001;

    raisedFlag = false;

    containsZeroes = false;
    numPosValues = 0;
    reversalCount = 0;
    bounceCount = 0;
    passingMeasures = 0;
    deltaq = 0.025;
    bounce = 0.10;
    reversals = 0;
    ncons0 = 2;
}

double demandmodeling::getPbar()
{
	std::vector<double> mPrices;

	for (int i = 0; i < (int) x.rows(); i++)
	{
		mPrices.push_back(x[i][0]);
	}

	sort( mPrices.begin(), mPrices.end() );
	mPrices.erase( unique( mPrices.begin(), mPrices.end() ), mPrices.end() );

    double sum = 0;

    for (int i = 0; i < mPrices.size(); i++)
    {
        sum = sum + mPrices[i];
    }

    return sum / (double) mPrices.size();
}

double demandmodeling::getMinimumConsumption()
{
	double minNonZeroConsumption = maxrealnumber;

    for (int i = 0; i < y.length(); i++)
    {
        if (y[i] > 0 && y[i] < minNonZeroConsumption)
        {
        	minNonZeroConsumption = y[i];
        }
    }

    if (isLogarithmic)
    {
    	minNonZeroConsumption = pow(10, minNonZeroConsumption);
    }

    return minNonZeroConsumption;
}

double demandmodeling::getMaximumConsumption()
{
	double maxNonZeroConsumption = minrealnumber;

    for (int i = 0; i < y.length(); i++)
    {
        if (y[i] > 0 && y[i] > maxNonZeroConsumption)
        {
        	maxNonZeroConsumption = y[i];
        }
    }

    if (isLogarithmic)
    {
    	maxNonZeroConsumption = pow(10, maxNonZeroConsumption);
    }

    return maxNonZeroConsumption;
}

std::string demandmodeling::getBP0String()
{
    double maxNonZeroPrice = minrealnumber;

    std::string priceString = "\"NA\"";
    std::ostringstream strs;

    for (int i = 0; i < y.length(); i++)
    {
        if (y[i] <= 0 && x[i][0] > maxNonZeroPrice)
        {
            maxNonZeroPrice = x[i][0];

            strs.str("");
            strs << maxNonZeroPrice;
            priceString = strs.str();
        }
    }

    return priceString;
}

std::string demandmodeling::getBP1String()
{
    double maxNonZeroPrice = minrealnumber;

    std::string priceString = "\"NA\"";
    std::ostringstream strs;

    for (int i = 0; i < y.length(); i++)
    {
        if (y[i] > 0 && x[i][0] > maxNonZeroPrice)
        {
            maxNonZeroPrice = x[i][0];


            strs.str("");
            strs << maxNonZeroPrice;
            priceString = strs.str();
        }
    }

	return priceString;
}

std::string demandmodeling::getOmaxEString()
{
    double maxExpendNumber = minrealnumber;

    std::ostringstream strs;

    for (int i = 0; i < y.length(); i++)
    {
        if ((x[i][0] * y[i]) >= maxExpendNumber)
        {
            maxExpendNumber = (x[i][0] * y[i]);
        }
    }

    strs << maxExpendNumber;

    return strs.str();
}

std::string demandmodeling::getPmaxEString()
{
    double maxExpendNumber = minrealnumber;
    double maxPrice = 0.0;

    std::ostringstream strs;

    for (int i = 0; i < y.length(); i++)
    {
        if ((x[i][0] * y[i]) >= maxExpendNumber)
        {
        	maxExpendNumber = (x[i][0] * y[i]);
            maxPrice = x[i][0];
        }
    }

    strs << maxPrice;

    return strs.str();
}

std::string demandmodeling::getIntensityString()
{
    double minNonZeroPrice = maxrealnumber;

    std::ostringstream strs;

    for (int i = 0; i < y.length(); i++)
    {
        if (x[i][0] < minNonZeroPrice)
        {
            minNonZeroPrice = x[i][0];
            strs.str("");
            strs << y[i];
        }
    }

    return strs.str();
}

struct QPairFirstComparer
{
    template<typename T1, typename T2>
    bool operator()(const std::pair<double, double> &one, const std::pair<double, double> &two) const
    {
        return one.first < two.first;
    }
};

demandmodeling::demandmodeling()
{

}
