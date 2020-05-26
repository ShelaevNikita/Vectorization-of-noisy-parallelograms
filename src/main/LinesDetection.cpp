#include <iostream>
#include <cmath>

#include "LinesDetection.h"

// First distribution type: point goes ONLY to the closest by distance line.
// Second distribution type: point goes ONLY to least voted by Hough line, if "close" to it.
// Third distribution type: point goes to 2 "closest" lines.

LinesDetection::LinesDetection(double rhoStep, double thetaStep, double interval,
                               distributionType newType) {
    this->rhoStep = rhoStep;
    this->thetaStep = thetaStep;
    this->interval = interval * rhoStep;
    this->type = newType;
}

void LinesDetection::parallelogram(vector<pair<float, float>> *inputPoints, vector<Vec3d> *outputLines) {

    Mat lines;
    vector<Point2f> pointVector;

    double rhoApproximation = 0.0;
    double distance;

    for (auto &inputPoint : *inputPoints) {
        pointVector.emplace_back(inputPoint.first, inputPoint.second);

        distance = hypot(inputPoint.first, inputPoint.second);
        rhoApproximation = rhoApproximation < distance ? distance : rhoApproximation;
    }

    double rhoMin = -rhoApproximation * 1.2, rhoMax = rhoApproximation * 1.2;
    double thetaMin = 0.0, thetaMax = CV_PI - thetaStep;
    HoughLinesPointSet(pointVector, lines, 400, 1,
                       rhoMin, rhoMax, rhoStep,
                       thetaMin, thetaMax, thetaStep);

    lines.copyTo(*outputLines);
}

void LinesDetection::fourLeaders(vector<Vec3d> *houghOutput, vector<lineABC> *leaders) {

    double rhoInterval = rhoStep * 30, thetaInterval = thetaStep * 600;
    double deltaRho, deltaTheta;

    vector<pair<double, double>> leadersTemp;
    leadersTemp.emplace_back(houghOutput->at(0).val[1], houghOutput->at(0).val[2]);
    int i = 1;
    bool inBunch;
    while (leadersTemp.size() < 4 && i < houghOutput->size()) {
        inBunch = false;

        for (auto leader : leadersTemp) {
            deltaRho = abs(houghOutput->at(i).val[1] - leader.first);
            deltaTheta = abs(houghOutput->at(i).val[2] - leader.second);
            if (deltaRho <= rhoInterval && deltaTheta <= thetaInterval) {
                inBunch = true;
                break;
            }

            if (signbit(leader.first) != signbit(houghOutput->at(i).val[1])) { // signbit returns true if x if negative
                deltaRho = abs(houghOutput->at(i).val[1] - (-leader.first));
                deltaTheta = abs(abs(houghOutput->at(i).val[2] - leader.second) - CV_PI);
                if (deltaRho <= rhoInterval && deltaTheta <= thetaInterval) {
                    inBunch = true;
                    break;
                }
            }
        }
        if (!inBunch) {
            leadersTemp.emplace_back(houghOutput->at(i).val[1], houghOutput->at(i).val[2]);

        }
        i++;
    }

    double a, b, c;
    for (i = 0; i < 4; i++) {
        a = cos(leadersTemp.at(i).second); // A = cos(theta)
        b = sin(leadersTemp.at(i).second); // B = sin(theta)
        c = -leadersTemp.at(i).first; // C = -r
        leaders->push_back(lineABC{a, b, c}); // r = x*cos(theta) + y*sin(theta) -> Ax + By + C = 0
    }

}

void LinesDetection::pointDistribution(vector<lineABC> *leaders, vector<pair<float, float>> *points,
                                       vector<vector<pair<double, double>>> *pointsSet) {

    double distance, x, y, a, b, c;

    if (type == ShortestDistance) {
        for (auto &point : *points) {
            double minDistance = INT32_MAX;
            int bestLine = 0;
            for (int j = 3; j > -1; j--) {
                x = point.first;
                y = point.second;
                a = leaders->at(j).a;
                b = leaders->at(j).b;
                c = leaders->at(j).c;

                distance = abs(a * x + b * y + c) / sqrt(a * a + b * b);

                if (distance < minDistance) {
                    bestLine = j;
                    minDistance = distance;
                }
            }
            pointsSet->at(bestLine).push_back(point);
        }
    } else {
        for (auto &point : *points) {
            for (int j = 3, used = 0; j > -1 && used < type; j--) { // used depends on type in number form
                x = point.first;
                y = point.second;
                a = leaders->at(j).a;
                b = leaders->at(j).b;
                c = leaders->at(j).c;

                distance = abs(a * x + b * y + c) / sqrt(a * a + b * b);

                if (distance <= interval) {
                    pointsSet->at(j).push_back(point);
                    used++;
                }
            }
        }
    }
}



