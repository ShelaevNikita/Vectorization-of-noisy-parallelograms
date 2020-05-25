#ifndef VECORIZATION_PARALLELOGRAMS_LINESDETECTION_H
#define VECORIZATION_PARALLELOGRAMS_LINESDETECTION_H

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <utility>
#include "vector"
#include "Structures.h"

using namespace cv;
using namespace std;

class LinesDetection {
public:
    enum distributionType {
        ShortestDistance,
        LeastVotes,
        TwoClosestLines
    };

private:
    double rhoStep,
            thetaStep,
            interval;

    distributionType type;
public:

    LinesDetection(double rhoStep,
                   double thetaStep,
                   double interval,
                   distributionType newType = ShortestDistance);

    void parallelogram(vector<pair<float, float>> *inputPoints, vector<Vec3d> *outputLines);

    void fourLeaders(vector<Vec3d> *houghOutput, vector<lineABC> *leaders);

    void pointDistribution(vector<lineABC> *leaders, vector<pair<float, float>> *points,
                           vector<vector<pair<double, double>>> *pointsSet);
};

#endif //VECORIZATION_PARALLELOGRAMS_LINESDETECTION_H
