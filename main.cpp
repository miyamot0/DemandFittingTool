/*
 * main.cpp
 *
 *  Created on: 8 Jun 2017
 *      Author: testing
 */

#include <iostream>
#include <sstream>
#include "demandmodeling.h"

using namespace std;

int main(int argc, char *argv[])
{
	demandmodeling mFitter;

    std::string mPrices = std::string(argv[1]);
    std::string mConsumption = std::string(argv[2]);
    std::string mModel = std::string(argv[3]);
    std::string mK = std::string(argv[4]);

	/*

    std::string mPrices;
    std::string mConsumption;
    std::string mModel;
    std::string mK;

    */

    //mModel = "Exponential";
    //mModel = "Exponentiated";
    //mPrices = "[[0.0],[0.5],[1.0],[1.5],[2.0],[2.5],[3.0],[4.0],[5.0],[10.0],[15.0]]";
    //mPrices = "[[0.0],[0.5],[1.0],[1.5],[2.0],[2.5],[3.0],[4.0],[5.0]]";
    //mConsumption = "[1000,1000,1000,800,800,700,600,500,400]";
    //mConsumption = "[10,10,10,8,8,7,6,5,4]";
    //mConsumption = "[10,10,10,8,8,7,6,5,4,2,1]";
    //mK = "fit";
    //mK = "range";
    //mK = "1";
    // ~test

	bool isLog = ((mModel.find("Linear") != std::string::npos) || (mModel.find("Exponential") != std::string::npos));

	mFitter.InitializeDefaults(isLog);

    std::ostringstream out;
    out << "{";

    /*
     * Linear Model Fit conditions
     * */
    if (mModel.find("Linear") != std::string::npos)
    {
    	mFitter.SetX(mPrices.c_str());
    	mFitter.SetY(mConsumption.c_str());
        mFitter.SetBounds("[+inf,+inf,+inf]", "[-inf,-inf,-inf]");

        /*
         * Fit Linear demand model
         */
        mFitter.FitLinear("[1,1,1]");

        /*
         * Reporting results, in JSON style return
         */
        if ((int) mFitter.GetInfo() == 2 || (int) mFitter.GetInfo() == 5)
        {
            double a = mFitter.GetState().c[0];
            double b = mFitter.GetState().c[1];
            double L = mFitter.GetState().c[2];

            double pmaxd = (1 + b)/a;
            double omaxd = (L * pow(pmaxd, b)) / exp(a * pmaxd) * pmaxd;

            double pbar = mFitter.getPbar();

            out << "\"ID\":" << 1 << ",";
            out << "\"Equation\":" << "\"Linear\"" << ",";
            out << "\"BP0\":" << mFitter.getBP0String() << ",";
            out << "\"BP1\":" << mFitter.getBP1String() << ",";
            out << "\"Omaxe\":" << mFitter.getOmaxEString() << ",";
            out << "\"Pmaxe\":" << mFitter.getPmaxEString() << ",";
            out << "\"L\":" << L << ",";
            out << "\"Lse\":" << mFitter.GetReport().errpar[2] << ",";
            out << "\"b\":" << b << ",";
            out << "\"bse\":" << mFitter.GetReport().errpar[1] << ",";
            out << "\"a\":" << a << ",";
            out << "\"ase\":" << mFitter.GetReport().errpar[0] << ",";
            out << "\"R2\":" << mFitter.GetReport().r2 << ",";
            out << "\"Elasticity\":" << "\"\"" << ",";		//pbar?
            out << "\"MeanElasticity\":" << (b - (a * pbar)) << ",";
            out << "\"Intensity\":" << mFitter.getIntensityString() << ",";		// get intensity
            out << "\"Omaxd\":" << omaxd << ",";
            out << "\"Pmaxd\":" << pmaxd << ",";
            out << "\"RMSError\":" << mFitter.GetReport().rmserror << ",";
            out << "\"avgerror\":" << mFitter.GetReport().avgerror << ",";
            out << "\"Notes\":" << mFitter.GetInfo();
        }
        else
        {
            out << "\"ID\":" << 1 << ",";
            out << "\"Equation\":" << "\"Linear\"" << ",";
            out << "\"BP0\":" << mFitter.getBP0String() << ",";
            out << "\"BP1\":" << mFitter.getBP1String() << ",";
            out << "\"Omaxe\":" << mFitter.getOmaxEString() << ",";
            out << "\"Pmaxe\":" << mFitter.getPmaxEString() << ",";
            out << "\"L\":" << "\"\"" << ",";
            out << "\"Lse\":" << "\"\"" << ",";
            out << "\"b\":" << "\"\"" << ",";
            out << "\"bse\":" << "\"\"" << ",";
            out << "\"a\":" << "\"\"" << ",";
            out << "\"ase\":" << "\"\"" << ",";
            out << "\"R2\":" << "\"\"" << ",";
            out << "\"Elasticity\":" << "\"\"" << ",";
            out << "\"MeanElasticity\":" << "\"\"" << ",";
            out << "\"Intensity\":" << mFitter.getIntensityString() << ",";		// get intensity
            out << "\"Omaxd\":" << "\"\"" << ",";
            out << "\"Pmaxd\":" << "\"\"" << ",";
            out << "\"RMSError\":" << "\"\"" << ",";
            out << "\"avgerror\":" << "\"\"" << ",";
            out << "\"Notes\":" << mFitter.GetInfo();
        }
    }
    /*
     * Exponential Model Fit conditions
     */
    else if (mModel.find("Exponential") != std::string::npos)
    {
    	mFitter.SetX(mPrices.c_str());
    	mFitter.SetY(mConsumption.c_str());

    	std::vector<double> params;

        /*
         * Fit Exponential, with K as parameter
         */
    	if (mK.find("fit") != std::string::npos)
    	{
            // K, Q0, Alpha
            std::ostringstream mUpperBoundos;
            mUpperBoundos << "[";
            mUpperBoundos << ((log10(mFitter.getMaximumConsumption()) - log10(mFitter.getMinimumConsumption())) + 0.5);
            mUpperBoundos << ", +inf, +inf]";

            mFitter.SetBounds(mUpperBoundos.str().c_str(), "[0.5, 0.001, -inf]");

            std::ostringstream mStartValues;
            mStartValues << "[";
            mStartValues << ((log10(mFitter.getMaximumConsumption()) - log10(mFitter.getMinimumConsumption())) + 0.5) / 2;
            mStartValues << ",";
            mStartValues << mFitter.getMaximumConsumption();
            mStartValues << ", 0.01]";

            /*
             * Fit Exponential, with K as parameter
             */
            mFitter.FitExponentialWithK(mStartValues.str().c_str());
    	}

		/*
		 * Fit Exponential, with K as derived scaling measure
		 */
    	else if (mK.find("range") != std::string::npos)
    	{
    		mFitter.SetBounds("[+inf, +inf]", "[0.0001, -inf]");

            std::ostringstream mStartValues;
            mStartValues << "[";
            mStartValues << mFitter.getMaximumConsumption();
            mStartValues << ", 0.01]";

    		params.push_back((log10(mFitter.getMaximumConsumption()) - log10(mFitter.getMinimumConsumption())) + 0.5);

    		/*
    		 * Fit Exponential, with K as derived scaling measure
    		 */
    		mFitter.FitExponential(mStartValues.str().c_str(), params);
    	}

		/*
		 * Fit Exponential, with K as assigned value
		 */
    	else
    	{
    		mFitter.SetBounds("[+inf, +inf]", "[0.0001, -inf]");

    		double mKparse = atof(mK.c_str());

            std::ostringstream mStartValues;
            mStartValues << "[";
            mStartValues << mFitter.getMaximumConsumption();
            mStartValues << ", 0.01]";

    		params.push_back(mKparse);

    		/*
    		 * Fit Exponential, with K as assigned value
    		 */
    		mFitter.FitExponential(mStartValues.str().c_str(), params);
    	}

        /*
         * Reporting results, in JSON style return
         */
        if ((int) mFitter.GetInfo() == 2 || (int) mFitter.GetInfo() == 5)
        {
            double alpha = (mK.find("fit") != std::string::npos) ? mFitter.GetState().c[2] : mFitter.GetState().c[1];
            double alphase = (mK.find("fit") != std::string::npos) ? mFitter.GetReport().errpar[2] : mFitter.GetReport().errpar[1];

            double k = (mK.find("fit") != std::string::npos) ? mFitter.GetState().c[0] : params[0];

            std::ostringstream mKos;

            if (mK.find("fit") != std::string::npos)
            {
            	mKos << mFitter.GetReport().errpar[0];
            }
            else
            {
            	mKos << "\"---\"";
            }

            std::string kse = mKos.str();

            double q0 = (mK.find("fit") != std::string::npos) ? mFitter.GetState().c[1] : mFitter.GetState().c[0];
            double q0se = (mK.find("fit") != std::string::npos) ? mFitter.GetReport().errpar[1] : mFitter.GetReport().errpar[0];

            double pmaxd = 1/(q0 * alpha * pow(k, 1.5)) * (0.083 * k + 0.65);
            double omaxd = (pow(10, (log10(q0) + (k * (exp(-alpha * q0 * pmaxd) - 1))))) * pmaxd;

            double EV = 1/(alpha * pow(k, 1.5) * 100);

            out << "\"ID\":" << 1 << ",";
            out << "\"Equation\":" << "\"Exponential\"" << ",";
            out << "\"BP1\":" << mFitter.getBP1String() << ",";
            out << "\"Omaxe\":" << mFitter.getOmaxEString() << ",";
            out << "\"Pmaxe\":" << mFitter.getPmaxEString() << ",";
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
            out << "\"R2\":" << mFitter.GetReport().r2 << ",";
            out << "\"RMSError\":" << mFitter.GetReport().rmserror << ",";
            out << "\"avgerror\":" << mFitter.GetReport().avgerror << ",";
            out << "\"Notes\":" << mFitter.GetInfo();
        }
        else
        {
            out << "\"ID\":" << 1 << ",";
            out << "\"Equation\":" << "\"Exponential\"" << ",";
            out << "\"BP1\":" << mFitter.getBP1String() << ",";
            out << "\"Omaxe\":" << mFitter.getOmaxEString() << ",";
            out << "\"Pmaxe\":" << mFitter.getPmaxEString() << ",";
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
            out << "\"Notes\":" << mFitter.GetInfo();
        }
    }
    /**
     * Exponentiated Model Fit conditions
     */
    else if (mModel.find("Exponentiated") != std::string::npos)
    {
    	mFitter.SetX(mPrices.c_str());
    	mFitter.SetY(mConsumption.c_str());

    	std::vector<double> params;

		/*
		 * Fit Exponentiated, with K as fitted parameter
		 */
    	if (mK.find("fit") != std::string::npos)
    	{
            // K, Q0, Alpha
            std::ostringstream mUpperBoundos;
            mUpperBoundos << "[";
            mUpperBoundos << ((mFitter.getMaximumConsumption() - mFitter.getMinimumConsumption()) + 0.5);
            mUpperBoundos << ", +inf, +inf]";

            mFitter.SetBounds(mUpperBoundos.str().c_str(), "[0.5, 0.001, -inf]");

            std::ostringstream mStartValues;
            mStartValues << "[";
            mStartValues << ((mFitter.getMaximumConsumption() - mFitter.getMinimumConsumption()) + 0.5) / 2;
            mStartValues << ",";
            mStartValues << mFitter.getMaximumConsumption();
            mStartValues << ", 0.01]";

    		/*
    		 * Fit Exponentiated, with K as fitted parameter
    		 */
            mFitter.FitExponentiatedWithK(mStartValues.str().c_str());
    	}

		/*
		 * Fit Exponentiated, with K as derived scaling measure
		 */
    	else if (mK.find("range") != std::string::npos)
    	{
    		mFitter.SetBounds("[+inf, +inf]", "[0.0001, -inf]");

            std::ostringstream mStartValues;
            mStartValues << "[";
            mStartValues << mFitter.getMaximumConsumption();
            mStartValues << ", 0.01]";

    		params.push_back(((mFitter.getMaximumConsumption() - mFitter.getMinimumConsumption()) + 0.5));

    		/*
    		 * Fit Exponentiated, with K as derived scaling measure
    		 */
    		mFitter.FitExponentiated(mStartValues.str().c_str(), params);
    	}

		/*
		 * Fit Exponentiated, with K as assigned value
		 */
    	else
    	{
    		mFitter.SetBounds("[+inf, +inf]", "[0.0001, -inf]");

    		double mKparse = atof(mK.c_str());

            std::ostringstream mStartValues;
            mStartValues << "[";
            mStartValues << mFitter.getMaximumConsumption();
            mStartValues << ", 0.01]";

    		params.push_back(mKparse);

    		/*
    		 * Fit Exponentiated, with K as assigned value
    		 */
    		mFitter.FitExponentiated(mStartValues.str().c_str(), params);
    	}

        /*
         * Reporting results, in JSON style return
         */
        if ((int) mFitter.GetInfo() == 2 || (int) mFitter.GetInfo() == 5)
        {
            double alpha = (mK.find("fit") != std::string::npos) ? mFitter.GetState().c[2] : mFitter.GetState().c[1];
            double alphase = (mK.find("fit") != std::string::npos) ? mFitter.GetReport().errpar[2] : mFitter.GetReport().errpar[1];

            double k = (mK.find("fit") != std::string::npos) ? mFitter.GetState().c[0] : params[0];

            std::ostringstream mKos;

            if (mK.find("fit") != std::string::npos)
            {
            	mKos << mFitter.GetReport().errpar[0];
            }
            else
            {
            	mKos << "\"---\"";
            }

            std::string kse = mKos.str();

            double q0 = (mK.find("fit") != std::string::npos) ? mFitter.GetState().c[1] : mFitter.GetState().c[0];
            double q0se = (mK.find("fit") != std::string::npos) ? mFitter.GetReport().errpar[1] : mFitter.GetReport().errpar[0];

            double pmaxd = 1/(q0 * alpha * pow(k, 1.5)) * (0.083 * k + 0.65);
            double omaxd = (q0 * (pow(10,(k * (exp(-alpha * q0 * pmaxd) - 1))))) * pmaxd;

            double EV = 1/(alpha * pow(k, 1.5) * 100);

            out << "\"ID\":" << 1 << ",";
            out << "\"Equation\":" << "\"Exponentiated\"" << ",";
            out << "\"BP1\":" << mFitter.getBP1String() << ",";
            out << "\"Omaxe\":" << mFitter.getOmaxEString() << ",";
            out << "\"Pmaxe\":" << mFitter.getPmaxEString() << ",";
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
            out << "\"R2\":" << mFitter.GetReport().r2 << ",";
            out << "\"RMSError\":" << mFitter.GetReport().rmserror << ",";
            out << "\"avgerror\":" << mFitter.GetReport().avgerror << ",";
            out << "\"Notes\":" << mFitter.GetInfo();
        }
        else
        {
            out << "\"ID\":" << 1 << ",";
            out << "\"Equation\":" << "\"Exponentiated\"" << ",";
            out << "\"BP1\":" << mFitter.getBP1String() << ",";
            out << "\"Omaxe\":" << mFitter.getOmaxEString() << ",";
            out << "\"Pmaxe\":" << mFitter.getPmaxEString() << ",";
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
            out << "\"Notes\":" << mFitter.GetInfo();
        }
    }

    out << "}";

	cout << out.str() << endl;

	return 0;
}
