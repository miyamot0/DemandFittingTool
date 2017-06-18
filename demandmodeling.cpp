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

void demandmodeling::SetModel(const char *mString)
{
	modelMode = mString;
}

void demandmodeling::SetX(const char *mString)
{
    x = mString;
}

void demandmodeling::SetY(const char *mString)
{
    y = mString;
    yStored = mString;

    likelyQ0 = getMaximumConsumption();

    //mModel.find("Exponential") != std::string::npos

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

void demandmodeling::FitExponentialWithK(const char *mStarts)
{
    SetStarts(mStarts);

    std::vector<double> mParams;

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);

    // K Q0 alpha
    std::ostringstream mScaleString;
    mScaleString << "[1.0, 1.0, 1.0e-";
    mScaleString << (scaleAssessment + (SignificantDigits() - 3));
    mScaleString << "]";

    real_1d_array s = mScaleString.str().c_str();
    lsfitsetscale(state, s);

    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, exponential_demand_with_k, NULL, &mParams);

    lsfitresults(state, info, c, rep);
}

void demandmodeling::BuildExponentialString(std::ostringstream &out, std::string mK, std::vector<double> params)
{
	out.str("");
	out << "{";

    if ((int) GetInfo() == 2 || (int) GetInfo() == 5)
    {
        double alpha = (mK.find("fit") != std::string::npos) ? GetState().c[2] : GetState().c[1];
        double alphase = (mK.find("fit") != std::string::npos) ? GetReport().errpar[2] : GetReport().errpar[1];

        double k = (mK.find("fit") != std::string::npos) ? GetState().c[0] : params[0];

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

        double q0 = (mK.find("fit") != std::string::npos) ? GetState().c[1] : GetState().c[0];
        double q0se = (mK.find("fit") != std::string::npos) ? GetReport().errpar[1] : GetReport().errpar[0];

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

void demandmodeling::FitExponentiatedWithK(const char *mStarts)
{
    SetStarts(mStarts);

    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);

    std::ostringstream mScaleString;
    mScaleString << "[1.0, 1.0, 1.0e-";
    mScaleString << (scaleAssessment + SignificantDigits());
    mScaleString << "]";

    real_1d_array s = mScaleString.str().c_str();
    lsfitsetscale(state, s);

    lsfitsetcond(state, epsx, maxits);

    alglib::lsfitfit(state, exponentiated_demand_with_k);

    lsfitresults(state, info, c, rep);
}

void demandmodeling::BuildExponentiatedString(std::ostringstream &out, std::string mK, std::vector<double> params)
{
	out.str("");
	out << "{";

    if ((int) GetInfo() == 2 || (int) GetInfo() == 5)
    {
        double alpha = (mK.find("fit") != std::string::npos) ? GetState().c[2] : GetState().c[1];
        double alphase = (mK.find("fit") != std::string::npos) ? GetReport().errpar[2] : GetReport().errpar[1];

        double k = (mK.find("fit") != std::string::npos) ? GetState().c[0] : params[0];

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

        double q0 = (mK.find("fit") != std::string::npos) ? GetState().c[1] : GetState().c[0];
        double q0se = (mK.find("fit") != std::string::npos) ? GetReport().errpar[1] : GetReport().errpar[0];

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

    for (int i = 0; i < yStored.length(); i++)
    {
        if (yStored[i] > 0 && yStored[i] < minNonZeroConsumption)
        {
        	minNonZeroConsumption = yStored[i];
        }
    }

    return minNonZeroConsumption;
}

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

std::string demandmodeling::buildUpperBoundsKSet()
{
	std::ostringstream mUpperBounds;
	mUpperBounds << "[";
	mUpperBounds << getMaximumConsumption() * 2;
	mUpperBounds << ",+inf]";

	return mUpperBounds.str();
}

std::string demandmodeling::buildStartValuesKSet(double mK)
{
    std::ostringstream mStartValues;
    mStartValues << "[";
    mStartValues << mK;
    mStartValues << ", 0.0001]";

	return mStartValues.str();
}

std::string demandmodeling::buildUpperBoundsFit()
{
    std::ostringstream mUpperBoundos;
    mUpperBoundos << "[";
    mUpperBoundos << (((log10(getMaximumConsumption()) - log10(getMinimumConsumption())) + 0.5) * 2);
    mUpperBoundos << ",";
    mUpperBoundos << (getMaximumConsumption() * 2);
    mUpperBoundos << ", +inf]";

	return mUpperBoundos.str();
}

std::string demandmodeling::buildStartValuesFit()
{
    std::ostringstream mStartValues;
    mStartValues << "[";
    mStartValues << ((log10(getMaximumConsumption()) - log10(getMinimumConsumption())) + 0.5);
    mStartValues << ",";
    mStartValues << getMaximumConsumption();
    mStartValues << ", 0.0001]";

	return mStartValues.str();
}

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
