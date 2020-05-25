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
                                                                                   1.0e-20,
                                                                                   1.0,
                                                                                   2500,
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
                                                              1.0e-20,
                                                              1.0,
                                                              2500,
                                                              5}) {

        return vectorization(move(points), error);
    }

private:

    struct line_result {
        double k, bmin, bmax, b, MSD;
        bool flag = true;
    };

    static void find_min_MSD(line_result (&res), vector<pair<double, double>> array_of_points) {
        double b = (res.bmin + res.bmax) / 2;
        double result = 0.0;

        double error_MSD = 1.0e-25;

        double result_last = result + 2 * error_MSD;
        int counter = 0;
        bool flag = true;
        if (res.k >= 1.0e5) 
            flag = false;
        double err = abs(result_last - result);

        int max_counter = 5000;

        while (err >= error_MSD && counter <= max_counter) {
            result = f_min(flag, res.bmin, res.bmax, res.k, b, array_of_points);
            if (b < res.bmin + 1) res.bmin -= 30;
            if (b > res.bmax - 1) res.bmax += 30;
            err = abs(result_last - result);
            result_last = result;
            counter++;
        }
        if (flag == false) 
            res.flag = false;
        res.b = b;
        res.MSD = sqrt(result);
    }

   static double foo_fmin(bool &flag, double &k, double sum, double b, vector<pair<double, double>> points) {
        int size = points.size();
        int i;
        double k1 = 0.0;
        for (i = 0; i < size; i++) 
            k1 += points[i].first * (points[i].second - b);
        k = k1 / sum;
        if (abs(k) > 1.0e5 || sum == 0) flag = false;                          // x = b
        k1 = 0.0;
        for (i = 0; i < size; i++) {
            if (flag == true) k1 += pow((points[i].second - k * points[i].first - b), 2);
            else k1 += pow((points[i].first - b), 2);
        }
        return k1;
    }

    static double f_min(bool &flag, double bmin, double bmax, double &k, double &b, vector<pair<double, double>> points) {
        double sum = 0.0;
        for (int i = 0; i < points.size(); i++)
            sum += pow(points[i].first, 2);
        double Lcenter, Rcenter, Lfcenter, Rfcenter, beg, end, BegBuf, EndBuf;
        double delta = (sqrt(5) - 1) / 2;
        beg = bmin;
        end = bmax;
        BegBuf = foo_fmin(flag, k, sum, beg, points);
        EndBuf = foo_fmin(flag, k, sum, end, points);

        int max_counter = 5000;
        double error_fmin = 1.0e-10;

        for (int i = 0; i <= max_counter; i++) {
            double g = delta * (end - beg);
            Lcenter = end - g;
            Lfcenter = foo_fmin(flag, k, sum, Lcenter, points);
            Rcenter = beg + g;
            Rfcenter = foo_fmin(flag, k, sum, Rcenter, points);
            if (Lfcenter < Rfcenter) {
                end = Rcenter;
                EndBuf = Rfcenter;
            } else if (Lfcenter > Rfcenter) {
                beg = Lcenter;
                BegBuf = Lfcenter;
            } else {
                end = Rcenter;
                EndBuf = Rfcenter;
                beg = Lcenter;
                BegBuf = Lfcenter;
            }
            double r = fabs(end - beg);
            if (r < error_fmin) {
                b = (end + beg) / 2;
                break;
            }
        }
        return (foo_fmin(flag, k, sum, b, points));
    }

    static pair<int, int> twoNeighbors(vector<pair<float, float>>* points, int midPoint) {
        float distance = 0;
        float distance1 = 1e15;
        float distance2 = 1e15;
        int closestNeighbor1 = 0;
        int closestNeighbor2 = 0;
        for (int i = 0; i < points->size(); i++) {
            if (i != midPoint) {
                distance = hypot(points->at(midPoint).first - points->at(i).first,
                    points->at(midPoint).second - points->at(i).second); // sqrt((x2-x1)^2 + (y2-y1)^2)
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

        return make_pair(closestNeighbor1, closestNeighbor2);
    }

    static double rhoStepEstimation(vector<pair<float, float>>* points) {
        random_device rd;   // non-deterministic generator
        mt19937 gen(rd()); // Mersenne-Twister generator for better random numbers than rand()
        uniform_int_distribution<> dist(0, points->size() - 1); // set distribution range to [0, n - 1]
        int midPoint[3];
        pair<int, int> neighbors[3];
        double candidates[3];

        for (int i = 0; i < 3; i++) {
            midPoint[i] = dist(gen); // choose random point from .txt as midPoint
            neighbors[i] = twoNeighbors(points, midPoint[i]); // find 2 closest to midPoint neighbors

            line_result res;
            double lineA = points->at(midPoint[i]).second - points->at(neighbors[i].first).second;
            double lineB = points->at(neighbors[i].first).first - points->at(midPoint[i]).first;
            double lineC = points->at(midPoint[i]).first * points->at(neighbors[i].first).second -
                points->at(neighbors[i].first).first * points->at(midPoint[i]).second;

            if (abs(lineA / lineB) >= 1000 && abs(lineB) <= 0.001)
                res.flag = false;    // x = b

            double bsr;
            if (res.flag == true) {
                bsr = -lineC / lineB;
                res.k = -lineA / lineB;
            }
            else {
                bsr = -lineC / lineA;
                res.k = 1.0e5;                                     // x = b
            }
            res.bmin = bsr - 15;
            res.bmax = bsr + 15;

            vector<pair<double, double>> triplet =
            { points->at(midPoint[i]), points->at(neighbors[i].first), points->at(neighbors[i].second) };
            find_min_MSD(res, triplet);

            double x[3] = { points->at(midPoint[i]).first,
                    points->at(neighbors[i].first).first,
                    points->at(neighbors[i].second).first };
            double y[3] = { points->at(midPoint[i]).second,
                    points->at(neighbors[i].first).second,
                    points->at(neighbors[i].second).second };
            double k = res.k;
            double b = res.b;
            candidates[i] = 0;
            for (int j = 0; j < 3; j++) {
                candidates[i] = max(candidates[i], 
                    sqrt(pow((x[j] + k * y[j] - k * b) / (k * k + 1) - x[j], 2) + 
                         pow(k * (x[j] + k * y[j] - k * b) / (k * k + 1) + b - y[j], 2))
                );
            }
        }

        double rhoStep = 0.0;
        for (int i = 0; i < 3; i++)
            rhoStep += candidates[i];
        return rhoStep / 3;
    }

    static vector<pair<double, double >>
    vectorization(vector<pair<float, float >> points,
                  error_vectorization error) {

        error.rhoStep = 1.2 * rhoStepEstimation(&points);
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
            max_counter - maximum number of iterations per cycle; default = 2500
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
