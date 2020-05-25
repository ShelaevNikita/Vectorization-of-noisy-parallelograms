#include <iostream>

#include <opencv2/core.hpp>
#include <fstream>
#include <utility>

#include "Vectorization_OpenCV.cpp"
#include "MSD.cpp"
#include "debug.h"

#include <random>

//#include "Structures.cpp"     :It's already included in MSD.cpp or Vect...CV.cpp

using namespace cv;
using namespace std;

struct error_vectorization {
    double rhoStep,
            thetaStep,
            error_fmin,
            error_MSD,
            count;
    
    int max_counter,
            interval;
};

class Vectorization_parallelograms {

public:
    static vector<pair<double, double >>
    vectorization_parallelograms(const string &fname, error_vectorization error = {0.1,
                                                                                   CV_PI / 3600.0,
                                                                                   1.0e-10,
                                                                                   1.0e-25,
                                                                                   2.0,
                                                                                   5000,
                                                                                   5}) {
        ifstream inputFile(fname);
        vector<pair<float, float>> points;

        pair<float, float> temp = make_pair(0.0, 0.0);

        if (!inputFile.is_open())
            throw ios_base::failure("Incorrect filepath");
        else {
            while (inputFile >> temp.first >> temp.second)
                points.push_back(temp);
        }
        inputFile.close();

        return vectorization(points, error);
    }

    static vector<pair<double, double >>
    vectorization_parallelograms(vector<pair<float, float>> points,
                                 error_vectorization error = {0.1,
                                                              CV_PI / 3600.0,
                                                              1.0e-10,
                                                              1.0e-25,
                                                              2.0,
                                                              5000,
                                                              5}) {

        return vectorization(move(points), error);
    }

private:

    static double rhoStepEstimation(vector<pair<float, float>>* points) {
        random_device rd;   // non-deterministic generator
        mt19937 gen(rd()); // Mersenne-Twister generator for better random numbers than rand()
        uniform_int_distribution<> dist(0, points->size() - 1); // set distribution range to [0, n - 1]
        int point1Index = dist(gen);
        int point2Index = dist(gen);
        pair<float, float> point1;
        pair<float, float> point2;

        SAY("Random indicies: %d %d\n", point1Index, point2Index);

        float distance = 0;
        float distance1 = 1E15;
        float distance2 = 1E15;
        int closestNeighbor1 = 0;
        int closestNeighbor2 = 0;
        for (int i = 0; i < points->size(); i++) {
            if (i != point1Index) {
                distance = hypot(points->at(point1Index).first - points->at(i).first,
                    points->at(point1Index).second - points->at(i).second); // sqrt((x2-x1)^2 + (y2-y1)^2)
                if (distance < distance1) {
                    distance2 = distance1;
                    closestNeighbor2 = closestNeighbor1;

                    distance1 = distance;
                    closestNeighbor1 = i;
                }
                else if (distance < distance2) {
                    distance2 = distance;
                    closestNeighbor2 = i;
                }
            }
        }


        distance = 0;
        distance1 = 1E15;
        distance2 = 1E15;
        int closestNeighbor3 = 0;
        int closestNeighbor4 = 0;
        for (int i = 0; i < points->size(); i++) {
            if (i != point2Index) {
                distance = hypot(points->at(point2Index).first - points->at(i).first,
                    points->at(point2Index).second - points->at(i).second); // sqrt((x2-x1)^2 + (y2-y1)^2)
                if (distance < distance1) {
                    distance2 = distance1;
                    closestNeighbor4 = closestNeighbor3;

                    distance1 = distance;
                    closestNeighbor3 = i;
                }
            }
        }

        SAY("(%f, %f) (%f, %f) (%f, %f)\n", points->at(point1Index).first, points->at(point1Index).second,
            points->at(closestNeighbor1).first, points->at(closestNeighbor1).second,
            points->at(closestNeighbor2).first, points->at(closestNeighbor2).second);

        SAY("(%f, %f) (%f, %f) (%f, %f)\n", points->at(point2Index).first, points->at(point2Index).second,
            points->at(closestNeighbor3).first, points->at(closestNeighbor3).second,
            points->at(closestNeighbor4).first, points->at(closestNeighbor4).second);

        // (y1-y2)x + (x2-x1)y + (x1y2-x2y1) = 0
        // point1 - neighbor1; point1 - neighbor2; neighbor1 - neighbor2

        double line1A = (2 * points->at(point1Index).second - 2 * points->at(closestNeighbor2).second) / 3;
        double line1B = (2 * points->at(closestNeighbor2).first - 2 * points->at(point1Index).first) / 3;
        double line1C = (points->at(point1Index).first * points->at(closestNeighbor1).second -
            points->at(closestNeighbor1).first * points->at(point1Index).second +
            points->at(point1Index).first * points->at(closestNeighbor2).second -
            points->at(closestNeighbor2).first * points->at(point1Index).second +
            points->at(closestNeighbor1).first * points->at(closestNeighbor2).second - 
            points->at(closestNeighbor2).first * points->at(closestNeighbor1).second) / 3;

        double candidate1 = abs(line1A * points->at(point1Index).first + line1B * points->at(point1Index).second + line1C) /
            sqrt(line1A * line1A + line1B * line1B);


        double line2A = (2 * points->at(point2Index).second - 2 * points->at(closestNeighbor4).second) / 3;
        double line2B = (2 * points->at(closestNeighbor4).first - 2 * points->at(point2Index).first) / 3;
        double line2C = (points->at(point2Index).first * points->at(closestNeighbor3).second -
            points->at(closestNeighbor3).first * points->at(point2Index).second +
            points->at(point2Index).first * points->at(closestNeighbor4).second -
            points->at(closestNeighbor4).first * points->at(point2Index).second +
            points->at(closestNeighbor3).first * points->at(closestNeighbor4).second -
            points->at(closestNeighbor4).first * points->at(closestNeighbor3).second) / 3;

        double candidate2 = abs(line2A * points->at(point2Index).first + line2B * points->at(point2Index).second + line2C) /
            sqrt(line2A * line2A + line2B * line2B);


        return max(candidate1, candidate2);
    }

    static vector<pair<double, double >>
    vectorization(vector<pair<float, float >> points,
                  error_vectorization error) {

        error.rhoStep = 1.5 * rhoStepEstimation(&points);
        if (error.rhoStep <= 0.01) {
            error.rhoStep = 0.01;
        }
        SAY("rhoStep = %f\n", error.rhoStep);

        LineExtraction mainObject(error.rhoStep, error.thetaStep, error.interval, LineExtraction::ShortestDistance);

        // MAIN CALL of parallelogram
        vector<Vec3d> line3dFirst;
        mainObject.parallelogram(&points, &line3dFirst);

        SAY("Number of lines: %lu\n", line3dFirst.size());

        /*
        for (auto lineHough : line3dFirst) {
            SAY("votes:%d   %f %f\n", (int) lineHough.val[0], lineHough.val[1], lineHough.val[2]);
        }
        */

        // Finding 4 leaders of 1.txt and distributing its points between the 4 leaders
        vector<lineABC> leaders;
        mainObject.fourLeaders(&line3dFirst, &leaders);

        vector<vector<pair<double, double>>> pointsSet(4);
        mainObject.pointDistribution(&leaders, &points, &pointsSet);

        for (int i = 0; i < 4; i++) {
            SAY("Leader %d: %lu points, %fx + %fy + %f = 0 \n", i + 1, pointsSet[i].size(), leaders[i].a,
                leaders[i].b, leaders[i].c);
            for (auto &j : pointsSet[i])
                SAY("\t(%f, %f)\n", j.first, j.second);
        }

       SAY("\n");

         /*
            MSD - mean square deviation:
            error_fmin - error of the "Golden section"; default = 1.0e-10
            error_MSD - error in finding the best result for MSD; default = 1.0e-20
            max_counter - maximum number of iterations per cycle; default = 3000
            count - priority of a parallelogram in choosing between it and a quadrilateral; count > 0 , default = 1.0
         */

         MSD MSD_foo(error.error_fmin, error.error_MSD, error.max_counter, error.count);
         double error_result = 0.0;
         vector<pair<double, double>> result_k = MSD_foo.MSD_main(pointsSet, leaders, error_result);

         SAY("\t\t     result: \n");
         for (int i = 0; i <= 3; i++)
            SAY(" \t x = %f \t y = %f\n", result_k[i].first, result_k[i].second);
         SAY("\t\tERROR: %f", error_result);
         SAY("\n___________________________________________________________________________\n");
        return result_k;
    }
};
