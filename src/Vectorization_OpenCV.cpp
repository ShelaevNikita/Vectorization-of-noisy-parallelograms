#include <iostream>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

#include <fstream>
#include <cmath>

using namespace cv;
using namespace std;

class Vectorization {

private:
    double rhoStep;
    double thetaStep;

public:

    Vectorization(double rhoStepInit, double thetaStepInit)
    {
        rhoStep = rhoStepInit;
        thetaStep = thetaStepInit;
    }

    void parallelogram(vector<pair<float, float>>* inputPoints, vector<Vec3d>* outputLines) {

        Mat lines;
        vector<Point2f> pointVector;

        double rhoApproximation = 0.0;
        double distance = 0.0;

        for (int i = 0; i < inputPoints->size(); i++)
        {
            pointVector.push_back(Point2f(inputPoints->at(i).first, inputPoints->at(i).second));
            distance = hypot(inputPoints->at(i).first, inputPoints->at(i).second);
            rhoApproximation = rhoApproximation < distance ? distance : rhoApproximation;
        }

        double rhoMin = -rhoApproximation * 1.2, rhoMax = rhoApproximation * 1.2;
        double thetaMin = 0.0, thetaMax = CV_PI;
        HoughLinesPointSet(pointVector, lines, 400, 1,
            rhoMin, rhoMax, rhoStep,
            thetaMin, thetaMax, thetaStep);

        lines.copyTo(*outputLines);
    }


    void parsePoints(string path, vector<pair<float, float>>* points) {
        ifstream inputFile(path);
        pair<float, float> temp = make_pair(0.0, 0.0);

        if (!inputFile.is_open()) 
            throw ios_base::failure("Incorrect filepath");
        else 
        {
            while (inputFile >> temp.first >> temp.second) 
                points->push_back(temp);
        }
        inputFile.close();
    }


    void fourLeaders(vector<Vec3d>* houghOutput, vector<pair<double, double>>* leaders) {

        double rhoInterval = rhoStep * 30, thetaInterval = thetaStep * 600;
        double deltaRho, deltaTheta;

        leaders->push_back(make_pair(houghOutput->at(0).val[1], houghOutput->at(0).val[2]));
        int i = 1;
        bool inBunch;
        while (leaders->size() < 4 && i < houghOutput->size()) 
        {
            inBunch = false;
            for (int j = 0; j < leaders->size(); j++) 
            {
                deltaRho = abs(houghOutput->at(i).val[1] - leaders->at(j).first);
                deltaTheta = abs(houghOutput->at(i).val[2] - leaders->at(j).second);
                if (deltaRho <= rhoInterval && deltaTheta <= thetaInterval)
                {
                    inBunch = true;
                    break;
                }
            }
            if (!inBunch) 
                leaders->push_back(make_pair(houghOutput->at(i).val[1], houghOutput->at(i).val[2]));
            i++;
        }

        double k, b;
        for (i = 0; i < 4; i++)
        {
            k = -cos(leaders->at(i).second) / sin(leaders->at(i).second); // k = -cos(theta) / sin(theta)
            b = leaders->at(i).first / sin(leaders->at(i).second); // b = r / sin(theta)
            leaders->at(i) = make_pair(k, b);
        }

    }

    void pointDistribution(vector<pair<double, double>>* leaders, vector<pair<float, float>>* points,
        vector<vector<pair<double, double>>>* pointsSet) {

        double rhoStep = 0.1;

        // Note: you should tie points to leaders from bottom to top only if you confirm that 4 leaders are actually different

        double interval = rhoStep * 5;
        double distance, x, y, k, b;
        for (int i = 0; i < points->size(); i++) 
        {
            for (int j = 3; j > -1; j--) 
            {
                x = points->at(i).first;
                y = points->at(i).second;
                k = leaders->at(j).first;
                b = leaders->at(j).second;
                distance = sqrt(pow((x + k * y - k * b) / (k * k + 1) - x, 2) +
                    pow(k * (x + k * y - k * b) / (k * k + 1) + b - y, 2));
                if (distance <= interval)
                {
                    pointsSet->at(j).push_back(points->at(i));
                    break;
                }
            }
        }
    }
};

