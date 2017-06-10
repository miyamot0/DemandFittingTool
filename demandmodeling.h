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

#ifndef DEMANDMODELING_H
#define DEMANDMODELING_H

#include <iostream>
#include "interpolation.h"

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"

using namespace alglib;
using namespace std;

class demandmodeling
{
public:
    demandmodeling();

    void SetX(const char *mString);
    void SetY(const char *mString);
    void SetStarts(const char *mString);
    void SetBounds(const char *mUpperString, const char *mLowerString);

    real_1d_array GetParams();
    lsfitstate GetState();
    ae_int_t GetInfo();
    lsfitreport GetReport();

    void FitLinear(const char *mStarts);
    void FitExponential(const char *mStarts, std::vector<double> mParams);
    void FitExponentialWithK(const char *mStarts);
    void FitExponentiated(const char *mStarts, std::vector<double> mParams);
    void FitExponentiatedWithK(const char *mStarts);

    double getMinimumConsumption();
    double getMaximumConsumption();

    double getPbar();

    std::string getBP0String();
    std::string getBP1String();
    std::string getOmaxEString();
    std::string getPmaxEString();
    std::string getIntensityString();

    //std::vector<std::string> GetSteinTest(std::vector<std::string> &x, std::vector<std::string> &y);

    void InitializeDefaults(bool isLog);

    bool raisedFlag;

private:
    real_2d_array x;
    real_1d_array y;
    real_1d_array c;

    real_1d_array bndl;
    real_1d_array bndu;

    ae_int_t maxits;
    ae_int_t info;
    lsfitstate state;
    lsfitreport rep;

    double epsx;
    double diffstep;

    std::vector<std::pair<double, double> > mModPoints;

    bool checkX, checkY;
    double tempX, tempY;
    bool containsZeroes;
    int numPosValues;
    int reversalCount;
    int bounceCount;
    int passingMeasures;
    double deltaq;
    double bounce;
    double reversals;
    double ncons0;

    bool isLogarithmic;

    std::string deltaQPass;
    std::string bouncePass;
    std::string reversalpass;

    double deltaQ;
    double bounceThreshold;
    double bounceScore;

};

#endif // DEMANDMODELING_H
