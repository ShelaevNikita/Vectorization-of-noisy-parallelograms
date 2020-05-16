#include <iostream>

#include <opencv2/core.hpp>
#include <fstream>
#include <utility>

#include "Vectorization_OpenCV.cpp"
#include "MSD.cpp"
#include "debug.h"

//#include "Structures.cpp"     :It's already included in MSD.cpp or Vect...CV.cpp

using namespace cv;
using namespace std;

struct error_vectorization {
    double rhoStep,
            thetaStep,
            error_fmin,
            error_MSD,
            d;
    int max_counter,
            interval;
};

class Vectorization_parallelograms {

public:
    static vector<pair<double, double >>
    vectorization_parallelograms(const string &fname, error_vectorization error = {0.1,
                                                                                   CV_PI / 3600.0,
                                                                                   1.0e-25,
                                                                                   1.0e-50,
                                                                                   0.15,
                                                                                   3000,
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
                                                              1.0e-25,
                                                              1.0e-50,
                                                              0.15,
                                                              3000,
                                                              5}) {

        return vectorization(move(points), error);
    }

private:
    static vector<pair<double, double >>
    vectorization(vector<pair<float, float >> points,
                  error_vectorization error) {

        LineExtraction mainObject(error.rhoStep, error.thetaStep, error.interval, LineExtraction::ShortestDistance);

        // MAIN CALL of parallelogram
        vector<Vec3d> line3dFirst;
        mainObject.parallelogram(&points, &line3dFirst);

        SAY("Number of lines: %lu\n", line3dFirst.size());

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
            error_fmin - error of the "Golden section"; default = 1.0e-15
            error_MSD - error in finding the best result for MSD; default = 1.0e-25
            max_counter - maximum number of iterations per cycle; default = 2500
            d - spread of values by b (b_min, b_max) for each line; d < 0.5; default = 0.15
        */

        MSD MSD_foo(error.error_fmin, error.error_MSD, error.max_counter, error.d);

        vector<pair<double, double>> result_k;
        /*
        MSD_foo.MSD_main(pointsSet, leaders, &result_k);

        SAY("\t\t     result: \n");
        for (int i = 0; i <= 3; i++)
            SAY(" \t x = %f \t y = %f\n", result_k[i].first, result_k[i].second);
        SAY("___________________________________________________________________________\n");
        */
        return result_k;
    }
};
