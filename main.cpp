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
    std::ostringstream out;
	std::vector<double> params;

    std::string mPrices = std::string(argv[1]);
    std::string mConsumption = std::string(argv[2]);
    std::string mModel = std::string(argv[3]);
    std::string mK = std::string(argv[4]);

    //std::string mPrices;
    //std::string mConsumption;
    //std::string mModel;
    //std::string mK;

    //mModel = "Exponential";
    //mModel = "Exponentiated";
    //mPrices = "[[0.5],[1.0],[1.5],[2.0],[2.5],[3.0],[4.0],[5.0],[10.0],[15.0]]";
    //mPrices = "[[0.0],[0.5],[1.0],[1.5],[2.0],[2.5],[3.0],[4.0],[5.0],[10.0],[15.0]]";
    //mConsumption = "[3,3,2.9030899869919438,2.9030899869919438,2.845098040014257,2.7781512503836434,2.6989700043360187,2.6020599913279625,2.3010299956639813,2]";
    //mConsumption = "[100,100,100,80,80,70,60,50,40,20,10]";
    //mConsumption = "[10,10,10,8,8,7,6,5,4,2,1]";

    //mConsumption = "[100,100,80,80,70,60,50,40,20,10]";
    //mConsumption = "[10,10,8,8,7,6,5,4,2,1]";
    //mConsumption = "[1,1,0.8,0.8,0.7,0.6,0.5,0.4,0.2,0.1]";
    //mK = "fit";
    //mK = "range";
    //mK = "1.5";

	mFitter.InitializeDefaults();
	mFitter.SetModel(mModel.c_str());
	mFitter.SetX(mPrices.c_str());
	mFitter.SetY(mConsumption.c_str());

    if (mModel.find("Linear") != std::string::npos)
    {
        /*
         * Linear Model Fit conditions
         * */

        mFitter.SetBounds("[+inf,+inf,+inf]", "[-inf,-inf,-inf]");
        mFitter.FitLinear("[1,1,1]");

        mFitter.BuildLinearString(out);
    }
    else if (mModel.find("Exponential") != std::string::npos)
    {
        /*
         * Exponential Model Fit conditions
         */

    	if (mK.find("fit") != std::string::npos)
    	{
            /*
             * Fit Exponential, with K as parameter
             */

            mFitter.SetBounds(mFitter.buildUpperBoundsFit().c_str(), "[0.5, 0.001, -inf]");
            mFitter.FitExponentialWithK(mFitter.buildStartValuesFit().c_str());
    	}
    	else if (mK.find("range") != std::string::npos)
    	{
    		/*
    		 * Fit Exponential, with K as derived scaling measure
    		 */

    		double mKvalue = log10(mFitter.getMaximumConsumption()) - log10(mFitter.getMinimumConsumption()) + 0.5;
    		params.push_back(mKvalue);

    		mFitter.SetBounds(mFitter.buildUpperBoundsKSet().c_str(), "[0.0001,-inf]");
    		mFitter.FitExponential(mFitter.buildStartValuesKSet(mKvalue).c_str(), params);
    	}
    	else
    	{
    		/*
    		 * Fit Exponential, with K as assigned value
    		 */

    		double mKvalue = atof(mK.c_str());
    		params.push_back(mKvalue);

    		mFitter.SetBounds(mFitter.buildUpperBoundsKSet().c_str(), "[0.0001,-inf]");
            mFitter.FitExponential(mFitter.buildStartValuesKSet(mKvalue).c_str(), params);
    	}

    	mFitter.BuildExponentialString(out, mK, params);
    }
    else if (mModel.find("Exponentiated") != std::string::npos)
    {
        /**
         * Exponentiated Model Fit conditions
         */

    	if (mK.find("fit") != std::string::npos)
    	{
    		/*
    		 * Fit Exponentiated, with K as fitted parameter
    		 */

    		mFitter.SetBounds(mFitter.buildUpperBoundsFit().c_str(), "[0.5, 0.001, -inf]");
            mFitter.FitExponentiatedWithK(mFitter.buildStartValuesFit().c_str());
    	}
    	else if (mK.find("range") != std::string::npos)
    	{
    		/*
    		 * Fit Exponentiated, with K as derived scaling measure
    		 */

    		double mKvalue = log10(mFitter.getMaximumConsumption()) - log10(mFitter.getMinimumConsumption()) + 0.5;
    		params.push_back(mKvalue);

    		mFitter.SetBounds(mFitter.buildUpperBoundsKSet().c_str(), "[0.0001,-inf]");
    		mFitter.FitExponentiated(mFitter.buildStartValuesKSet(mKvalue).c_str(), params);
    	}
    	else
    	{
    		/*
    		 * Fit Exponentiated, with K as assigned value
    		 */

    		double mKvalue = atof(mK.c_str());
    		params.push_back(mKvalue);

    		mFitter.SetBounds(mFitter.buildUpperBoundsKSet().c_str(), "[0.0001,-inf]");
    		mFitter.FitExponentiated(mFitter.buildStartValuesKSet(mKvalue).c_str(), params);
    	}

    	mFitter.BuildExponentiatedString(out, mK, params);
    }

	cout << out.str() << endl;

	return 0;
}
