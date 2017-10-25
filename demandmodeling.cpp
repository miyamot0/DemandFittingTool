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

// <!-- Generic helper -->
void demandmodeling::SetModel(const char *mString)
{
	modelMode = mString;
}

// <!-- Generic helper -->
void demandmodeling::SetX(const char *mString)
{
    x = mString;
}

// <!-- Generic helper -->
void demandmodeling::SetY(const char *mString)
{
    y = mString;
    yStored = mString;

    likelyQ0 = getMaximumConsumption();

    if (modelMode.find("Exponential")  != std::string::npos)
    {
    	for (int i=0; i<y.length(); i++)
    	{
    		y[i] = log10(y[i]);
    	}
    }
    else if (modelMode.find("Linear")  != std::string::npos)
    {
    	for (int i=0; i<y.length(); i++)
    	{
    		y[i] = log(y[i]);
    	}
    }
}

// <!-- Generic helper -->
void demandmodeling::SetStarts(const char *mString)
{
    c = mString;
}

// <!-- Generic helper -->
void demandmodeling::SetBounds(const char *mUpperString, const char *mLowerString)
{
    bndu = mUpperString;
    bndl = mLowerString;
}

// <!-- Generic helper -->
real_1d_array demandmodeling::GetParams()
{
    return c;
}

// <!-- Generic helper -->
lsfitstate demandmodeling::GetState()
{
    return state;
}

// <!-- Generic helper -->
ae_int_t demandmodeling::GetInfo()
{
    return info;
}

// <!-- Generic helper -->
lsfitreport demandmodeling::GetReport()
{
    return rep;
}

// <!-- Exponential model call -->
void exponential_demand(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
	std::vector<double> *param = (std::vector<double> *) ptr;
    double k = param->at(0);
    func = log10(c[0]) + k * (exp(-c[1] * c[0] * x[0]) - 1);
}

// <!-- Exponential model call -->
void exponential_demand_with_k(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
    func = log10(c[0]) + c[2] * (exp(-c[1] * c[0] * x[0]) - 1);
}

// <!-- Exponentiated model call -->
void exponentiated_demand(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
	std::vector<double> *param = (std::vector<double> *) ptr;
    double k = param->at(0);
    func = c[0] * pow(10, (k * (exp(-c[1] * c[0] * x[0]) - 1)));
}

// <!-- Exponentiated model call -->
void exponentiated_demand_with_k(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
    func = c[0] * pow(10, (c[2] * (exp(-c[1] * c[0] * x[0]) - 1)));
}

// <!-- Linear model call -->
void linear_demand(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
{
    func = log(c[2]) + (c[1] * log(x[0])) - c[0] * x[0];
}

// <!-- Generic helper -->
int demandmodeling::SignificantDigits()
{
    double x = abs(likelyQ0);
    return (x < 10 ? 1 :
        (x < 100 ? 2 :
        (x < 1000 ? 3 :
        (x < 10000 ? 4 :
        (x < 100000 ? 5 :
        (x < 1000000 ? 6 :
        (x < 10000000 ? 7 :
        (x < 100000000 ? 8 :
        (x < 1000000000 ? 9 :
        10)))))))));
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

void demandmodeling::BuildLinearString(std::ostringstream &out)
{
	out.str("");
	out << "{";

    if ((int) GetInfo() == 2 || (int) GetInfo() == 5)
    {
        double a = GetState().c[0];
        double b = GetState().c[1];
        double L = GetState().c[2];

        double pmaxd = (1 + b)/a;
        double omaxd = (L * pow(pmaxd, b)) / exp(a * pmaxd) * pmaxd;

        double pbar = getPbar();

        out << "\"ID\":" << 1 << ",";
        out << "\"Equation\":" << "\"Linear\"" << ",";
        out << "\"BP0\":" << getBP0String() << ",";
        out << "\"BP1\":" << getBP1String() << ",";
        out << "\"Omaxe\":" << getOmaxEString() << ",";
        out << "\"Pmaxe\":" << getPmaxEString() << ",";
        out << "\"L\":" << L << ",";
        out << "\"Lse\":" << GetReport().errpar[2] << ",";
        out << "\"b\":" << b << ",";
        out << "\"bse\":" << GetReport().errpar[1] << ",";
        out << "\"a\":" << a << ",";
        out << "\"ase\":" << GetReport().errpar[0] << ",";
        out << "\"R2\":" << GetReport().r2 << ",";
        out << "\"Elasticity\":" << "\"\"" << ",";		//pbar?
        out << "\"MeanElasticity\":" << (b - (a * pbar)) << ",";
        out << "\"Intensity\":" << getIntensityString() << ",";		// get intensity
        out << "\"Omaxd\":" << omaxd << ",";
        out << "\"Pmaxd\":" << pmaxd << ",";
        out << "\"RMSError\":" << GetReport().rmserror << ",";
        out << "\"avgerror\":" << GetReport().avgerror << ",";
        out << "\"Notes\":" << GetInfo();
    }
    else
    {
        out << "\"ID\":" << 1 << ",";
        out << "\"Equation\":" << "\"Linear\"" << ",";
        out << "\"BP0\":" << getBP0String() << ",";
        out << "\"BP1\":" << getBP1String() << ",";
        out << "\"Omaxe\":" << getOmaxEString() << ",";
        out << "\"Pmaxe\":" << getPmaxEString() << ",";
        out << "\"L\":" << "\"\"" << ",";
        out << "\"Lse\":" << "\"\"" << ",";
        out << "\"b\":" << "\"\"" << ",";
        out << "\"bse\":" << "\"\"" << ",";
        out << "\"a\":" << "\"\"" << ",";
        out << "\"ase\":" << "\"\"" << ",";
        out << "\"R2\":" << "\"\"" << ",";
        out << "\"Elasticity\":" << "\"\"" << ",";
        out << "\"MeanElasticity\":" << "\"\"" << ",";
        out << "\"Intensity\":" << getIntensityString() << ",";		// get intensity
        out << "\"Omaxd\":" << "\"\"" << ",";
        out << "\"Pmaxd\":" << "\"\"" << ",";
        out << "\"RMSError\":" << "\"\"" << ",";
        out << "\"avgerror\":" << "\"\"" << ",";
        out << "\"Notes\":" << GetInfo();
    }

	out << "}";
}

// <!-- Brute force start points for exponential model -->
void demandmodeling::BuildStartString(std::ostringstream &out, double setK, bool isExponentiated)
{
	int gridSize = 40;

	double sQ0[gridSize+1];
	double sA[gridSize+1];
	double sK[gridSize+1];

	double minQ = getMinimumConsumption();
	if (!(minQ > 0.0))
	{
		minQ = 0.01;
	}
	minQ = log(minQ);

	double maxC = log(getMaximumConsumption() * 1.5);

	double minA = 1;
	double maxA = 1.09;

	double minK = log(0.5);
	double maxK = log(log(getMaximumConsumption()) + 1);

	for (int i = 0; i <= gridSize; i++)
	{
		sQ0[i] = minQ + i * ((maxC - minQ) / gridSize);
		sA[i] = minA + i * ((maxA - minA) / gridSize);
		sK[i] = minK + i * ((maxK - minK) / gridSize);

		sQ0[i] = exp(sQ0[i]);
		sA[i] = log(sA[i]);
		sK[i] = exp(sK[i]);
	}

	double ssr = 0;
	double optimSSR = maxrealnumber;

	double bestK = 1;
	double bestQ = 1;
	double bestA = 0.001;

	if (setK < 0)
	{
		for (int k = 0; k <= gridSize; k++)
		{
			for (int q = 0; q <= gridSize; q++)
			{
				for (int a = 0; a <= gridSize; a++)
				{
					if (isExponentiated)
					{
						GetSumSquaresExponentiated(ssr, sQ0[q], sA[a], sK[k]);
					}
					else
					{
						GetSumSquaresExponential(ssr, sQ0[q], sA[a], sK[k]);
					}

					if (ssr < optimSSR)
					{
						optimSSR = ssr;

						bestK = sK[k];
						bestQ = sQ0[q];
						bestA = sA[a];
					}
				}
			}
		}
	}
	else
	{
		for (int q = 0; q <= gridSize; q++)
		{
			for (int a = 0; a <= gridSize; a++)
			{
				if (isExponentiated)
				{
					GetSumSquaresExponentiated(ssr, sQ0[q], sA[a], setK);
				}
				else
				{
					GetSumSquaresExponential(ssr, sQ0[q], sA[a], setK);
				}

				if (ssr < optimSSR)
				{
					optimSSR = ssr;

					bestQ = sQ0[q];
					bestA = sA[a];
				}
			}
		}
	}

	out.clear();
	out.str("");
	out << "[";
	out << bestQ;
	out << ", ";
	out << bestA;

	if (setK < 0)
	{
		out << ", ";
		out << bestK;
	}

	out << "]";
}

// <!-- Calculate sum squared residuals for exponential model -->
void demandmodeling::GetSumSquaresExponential(double &ssr, double &q0, double &alpha, double &k)
{
	ssr = 0;

    for (int i = 0; i < y.length(); i++)
    {
    	projection = log10(q0) + k * (exp(-alpha * q0 * x[i][0]) - 1.0);
    	ssr = ssr + pow(y[i] - projection, 2);
    }
}

// <!-- Calculate sum squared residuals for exponentiated model -->
void demandmodeling::GetSumSquaresExponentiated(double &ssr, double &q0, double &alpha, double &k)
{
	ssr = 0;

    for (int i = 0; i < y.length(); i++)
    {
    	projection = q0 * pow(10.0, (k * (exp(-(alpha) * q0 * x[i][0]) - 1.0)));

    	ssr = ssr + pow(y[i] - projection, 2);
    }
}

// <!-- Fit exponential demand model -->
void demandmodeling::FitExponential(const char *mStarts, std::vector<double> mParams)
{
    SetStarts(mStarts);
    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);

    std::ostringstream mScaleString;
    mScaleString << "[1.0, 1.0e-";
    mScaleString << (scaleAssessment + (SignificantDigits() - 3));
    mScaleString << "]";

    real_1d_array s = mScaleString.str().c_str();
    lsfitsetscale(state, s);

    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, exponential_demand, NULL, &mParams);

    lsfitresults(state, info, c, rep);
}

// <!-- Fit exponential demand model with K -->
void demandmodeling::FitExponentialWithK(const char *mStarts)
{
    SetStarts(mStarts);

    std::vector<double> mParams;

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);

    std::ostringstream mScaleString;
    mScaleString << "[1.0, 1.0e-";
    mScaleString << (scaleAssessment + (SignificantDigits() - 3));
    mScaleString << ", 1.0]";

    real_1d_array s = mScaleString.str().c_str();
    lsfitsetscale(state, s);

    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, exponential_demand_with_k, NULL, &mParams);

    lsfitresults(state, info, c, rep);
}

// <!-- Construct json-readable string from results -->
void demandmodeling::BuildExponentialString(std::ostringstream &out, std::string mK, std::vector<double> params)
{
	out.str("");
	out << "{";

    if ((int) GetInfo() == 2 || (int) GetInfo() == 5)
    {
        double alpha = GetState().c[1];
        double alphase = GetReport().errpar[1];

        double k = (mK.find("fit") != std::string::npos) ? GetState().c[2] : params[0];

        std::ostringstream mKos;

        if (mK.find("fit") != std::string::npos)
        {
        	mKos << GetReport().errpar[0];
        }
        else
        {
        	mKos << "\"---\"";
        }

        std::string kse = mKos.str();

        double q0 = GetState().c[0];
        double q0se = GetReport().errpar[0];

        double pmaxd = 1/(q0 * alpha * pow(k, 1.5)) * (0.083 * k + 0.65);
        double omaxd = (pow(10, (log10(q0) + (k * (exp(-alpha * q0 * pmaxd) - 1))))) * pmaxd;

        double EV = 1/(alpha * pow(k, 1.5) * 100);

        out << "\"ID\":" << 1 << ",";
        out << "\"Equation\":" << "\"Exponential\"" << ",";
        out << "\"BP1\":" << getBP1String() << ",";
        out << "\"Omaxe\":" << getOmaxEString() << ",";
        out << "\"Pmaxe\":" << getPmaxEString() << ",";
        out << "\"Omaxd\":" << omaxd << ",";
        out << "\"Pmaxd\":" << pmaxd << ",";
        out << "\"Alpha\":" << alpha << ",";
        out << "\"Alphase\":" << alphase << ",";
        out << "\"Q0d\":" << q0 << ",";
        out << "\"Q0dse\":" << q0se << ",";
        out << "\"K\":" << k << ",";
        out << "\"Kse\":" << kse << ",";
        out << "\"EV\":" << EV << ",";
        out << "\"kMethod\":\"" << mK << "\",";
        out << "\"R2\":" << GetReport().r2 << ",";
        out << "\"RMSError\":" << GetReport().rmserror << ",";
        out << "\"avgerror\":" << GetReport().avgerror << ",";
        out << "\"Notes\":" << GetInfo();
    }
    else
    {
        out << "\"ID\":" << 1 << ",";
        out << "\"Equation\":" << "\"Exponential\"" << ",";
        out << "\"BP1\":" << getBP1String() << ",";
        out << "\"Omaxe\":" << getOmaxEString() << ",";
        out << "\"Pmaxe\":" << getPmaxEString() << ",";
        out << "\"Omaxd\":" << "\"\"" << ",";
        out << "\"Pmaxd\":" << "\"\"" << ",";
        out << "\"Alpha\":" << "\"\"" << ",";
        out << "\"Alphase\":" << "\"\"" << ",";
        out << "\"Q0d\":" << "\"\"" << ",";
        out << "\"Q0dse\":" << "\"\"" << ",";
        out << "\"K\":" << "\"\"" << ",";
        out << "\"Kse\":" << "\"\"" << ",";
        out << "\"EV\":" << "\"\"" << ",";
        out << "\"kMethod\":\"" << mK << "\",";
        out << "\"R2\":" << "\"\"" << ",";
        out << "\"RMSError\":" << "\"\"" << ",";
        out << "\"avgerror\":" << "\"\"" << ",";
        out << "\"Notes\":" << GetInfo();
    }

    out << "}";
}

// <!-- Fit exponentiated demand model -->
void demandmodeling::FitExponentiated(const char *mStarts, std::vector<double> mParams)
{
    SetStarts(mStarts);

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);

    std::ostringstream mScaleString;
    mScaleString << "[1.0, 1.0e-";
    mScaleString << (scaleAssessment + SignificantDigits());
    mScaleString << "]";

    real_1d_array s = mScaleString.str().c_str();
    lsfitsetscale(state, s);

    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, exponentiated_demand, NULL, &mParams);

    lsfitresults(state, info, c, rep);
}

// <!-- Fit exponentiated demand model with K -->
void demandmodeling::FitExponentiatedWithK(const char *mStarts)
{
    SetStarts(mStarts);

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);

    std::ostringstream mScaleString;
    mScaleString << "[1.0, 1.0e-";
    mScaleString << (scaleAssessment + SignificantDigits());
    mScaleString << ", 1.0]";

    real_1d_array s = mScaleString.str().c_str();
    lsfitsetscale(state, s);

    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, exponentiated_demand_with_k);

    lsfitresults(state, info, c, rep);
}

// <!-- Construct json-readable string from results -->
void demandmodeling::BuildExponentiatedString(std::ostringstream &out, std::string mK, std::vector<double> params)
{
	out.str("");
	out << "{";

    if ((int) GetInfo() == 2 || (int) GetInfo() == 5)
    {
        double alpha = GetState().c[1];
        double alphase = GetReport().errpar[1];

        double k = (mK.find("fit") != std::string::npos) ? GetState().c[2] : params[0];

        std::ostringstream mKos;

        if (mK.find("fit") != std::string::npos)
        {
        	mKos << GetReport().errpar[2];
        }
        else
        {
        	mKos << "\"---\"";
        }

        std::string kse = mKos.str();

        double q0 = GetState().c[0];
        double q0se = GetReport().errpar[0];

        double pmaxd = 1/(q0 * alpha * pow(k, 1.5)) * (0.083 * k + 0.65);
        double omaxd = (q0 * (pow(10,(k * (exp(-alpha * q0 * pmaxd) - 1))))) * pmaxd;

        double EV = 1/(alpha * pow(k, 1.5) * 100);

        out << "\"ID\":" << 1 << ",";
        out << "\"Equation\":" << "\"Exponentiated\"" << ",";
        out << "\"BP1\":" << getBP1String() << ",";
        out << "\"Omaxe\":" << getOmaxEString() << ",";
        out << "\"Pmaxe\":" << getPmaxEString() << ",";
        out << "\"Omaxd\":" << omaxd << ",";
        out << "\"Pmaxd\":" << pmaxd << ",";
        out << "\"Alpha\":" << alpha << ",";
        out << "\"Alphase\":" << alphase << ",";
        out << "\"Q0d\":" << q0 << ",";
        out << "\"Q0dse\":" << q0se << ",";
        out << "\"K\":" << k << ",";
        out << "\"Kse\":" << kse << ",";
        out << "\"EV\":" << EV << ",";
        out << "\"kMethod\":\"" << mK << "\",";
        out << "\"R2\":" << GetReport().r2 << ",";
        out << "\"RMSError\":" << GetReport().rmserror << ",";
        out << "\"avgerror\":" << GetReport().avgerror << ",";
        out << "\"Notes\":" << GetInfo();
    }
    else
    {
        out << "\"ID\":" << 1 << ",";
        out << "\"Equation\":" << "\"Exponentiated\"" << ",";
        out << "\"BP1\":" << getBP1String() << ",";
        out << "\"Omaxe\":" << getOmaxEString() << ",";
        out << "\"Pmaxe\":" << getPmaxEString() << ",";
        out << "\"Omaxd\":" << "\"\"" << ",";
        out << "\"Pmaxd\":" << "\"\"" << ",";
        out << "\"Alpha\":" << "\"\"" << ",";
        out << "\"Alphase\":" << "\"\"" << ",";
        out << "\"Q0d\":" << "\"\"" << ",";
        out << "\"Q0dse\":" << "\"\"" << ",";
        out << "\"K\":" << "\"\"" << ",";
        out << "\"Kse\":" << "\"\"" << ",";
        out << "\"EV\":" << "\"\"" << ",";
        out << "\"kMethod\":\"" << mK << "\",";
        out << "\"R2\":" << "\"\"" << ",";
        out << "\"RMSError\":" << "\"\"" << ",";
        out << "\"avgerror\":" << "\"\"" << ",";
        out << "\"Notes\":" << GetInfo();
    }

	out << "}";
}

// <!-- Initialize defaults -->
void demandmodeling::InitializeDefaults()
{
	maxits = 10000;
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

    likelyQ0 = -1;

    scaleAssessment = 3;
}

// <!-- Generic helper -->
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

    for (int i = 0; i < (int)mPrices.size(); i++)
    {
        sum = sum + mPrices[i];
    }

    return sum / (double) mPrices.size();
}

// <!-- Generic helper -->
double demandmodeling::getMinimumCost()
{
	double minNonZeroConsumption = maxrealnumber;

    for (int i = 0; i < x.rows(); i++)
    {
        if (x[i][0] > 0 && x[i][0] < minNonZeroConsumption)
        {
        	minNonZeroConsumption = x[i][0];
        }
    }

    return minNonZeroConsumption;
}

// <!-- Generic helper -->
double demandmodeling::getMinimumConsumption()
{
	double minNonZeroConsumption = maxrealnumber;

    for (int i = 0; i < yStored.length(); i++)
    {
        if (yStored[i] > 0 && yStored[i] < minNonZeroConsumption)
        {
        	minNonZeroConsumption = yStored[i];
        }
    }

    return minNonZeroConsumption;
}

// <!-- Generic helper -->
double demandmodeling::getMaximumConsumption()
{
	double maxNonZeroConsumption = minrealnumber;

    for (int i = 0; i < yStored.length(); i++)
    {
        if (yStored[i] > 0 && yStored[i] > maxNonZeroConsumption)
        {
        	maxNonZeroConsumption = yStored[i];
        }
    }

    return maxNonZeroConsumption;
}

// <!-- Generic helper -->
std::string demandmodeling::buildUpperBoundsFit()
{
    std::ostringstream mUpperBoundos;
    mUpperBoundos << "[+inf,+inf,";
    mUpperBoundos << (((log(getMaximumConsumption()) - log(getMinimumConsumption())) + 0.5) * 2);
    mUpperBoundos << "]";

	return mUpperBoundos.str();
}

// <!-- Generic helper -->
std::string demandmodeling::getBP0String()
{
    double maxNonZeroPrice = minrealnumber;

    std::string priceString = "\"NA\"";
    std::ostringstream strs;

    for (int i = 0; i < yStored.length(); i++)
    {
        if (yStored[i] <= 0 && x[i][0] > maxNonZeroPrice)
        {
            maxNonZeroPrice = x[i][0];

            strs.str("");
            strs << maxNonZeroPrice;
            priceString = strs.str();
        }
    }

    return priceString;
}

// <!-- Generic helper -->
std::string demandmodeling::getBP1String()
{
    double maxNonZeroPrice = minrealnumber;

    std::string priceString = "\"NA\"";
    std::ostringstream strs;

    for (int i = 0; i < yStored.length(); i++)
    {
        if (yStored[i] > 0 && x[i][0] > maxNonZeroPrice)
        {
            maxNonZeroPrice = x[i][0];


            strs.str("");
            strs << maxNonZeroPrice;
            priceString = strs.str();
        }
    }

	return priceString;
}

// <!-- Generic helper -->
std::string demandmodeling::getOmaxEString()
{
    double maxExpendNumber = minrealnumber;

    std::ostringstream strs;

    for (int i = 0; i < yStored.length(); i++)
    {
        if ((x[i][0] * yStored[i]) >= maxExpendNumber)
        {
            maxExpendNumber = (x[i][0] * yStored[i]);
        }
    }

    strs << maxExpendNumber;

    return strs.str();
}

// <!-- Generic helper -->
std::string demandmodeling::getPmaxEString()
{
    double maxExpendNumber = minrealnumber;
    double maxPrice = 0.0;

    std::ostringstream strs;

    for (int i = 0; i < yStored.length(); i++)
    {
        if ((x[i][0] * yStored[i]) >= maxExpendNumber)
        {
        	maxExpendNumber = (x[i][0] * yStored[i]);
            maxPrice = x[i][0];
        }
    }

    strs << maxPrice;

    return strs.str();
}

// <!-- Generic helper -->
std::string demandmodeling::getIntensityString()
{
    double minNonZeroPrice = maxrealnumber;

    std::ostringstream strs;

    for (int i = 0; i < yStored.length(); i++)
    {
        if (x[i][0] < minNonZeroPrice)
        {
            minNonZeroPrice = x[i][0];
            strs.str("");
            strs << yStored[i];
        }
    }

    return strs.str();
}

demandmodeling::demandmodeling()
{

}
