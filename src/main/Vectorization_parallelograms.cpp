#include <fstream>
#include <algorithm>
#include "LinesDetection.h"
#include "MSD.h"
#include "debug.h"
#include "Vectorization_parallelograms.h"

vector<pair<double, double >>
Vectorization_parallelograms::vectorization_parallelograms(const string &fname, error_vectorization error) {
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

vector<pair<double, double >>
Vectorization_parallelograms::vectorization_parallelograms(vector<pair<float, float>> points,
                                                           error_vectorization error) {

    return vectorization(move(points), error);
}

pair<int, int> Vectorization_parallelograms::twoNeighbors(vector<pair<float, float>> *points, int midPoint) {
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
            } else if (distance < distance2) {
                distance2 = distance;
                closestNeighbor2 = i;
            }
        }
    }

    return make_pair(closestNeighbor1, closestNeighbor2);
}

double Vectorization_parallelograms::rhoStepEstimation(vector<pair<float, float>> *points, error_vectorization error) {
    pair<int, int> neighbors;
    double sum_rho = 0;
    vector<double> rho;
    for (int i = 0; i < points->size(); i++) {
        int midPoint = i;
        neighbors = twoNeighbors(points, midPoint); // find 2 closest to midPoint neighbors

        MSD::line_result res;
        double lineA = points->at(midPoint).second - points->at(neighbors.first).second;
        double lineB = points->at(neighbors.first).first - points->at(midPoint).first;
        double lineC = points->at(midPoint).first * points->at(neighbors.first).second -
                       points->at(neighbors.first).first * points->at(midPoint).second;

        if (abs(lineA / lineB) >= 1000 && abs(lineB) <= 0.001)
            res.flag = false;    // x = b

        double bsr;
        if (res.flag) {
            bsr = -lineC / lineB;
            res.k = -lineA / lineB;
        } else {
            bsr = -lineC / lineA;
            res.k = 1.0e5;                                     // x = b
        }
        res.bmin = bsr - 20;
        res.bmax = bsr + 20;
        vector<pair<double, double>> triplet =
                {points->at(midPoint), points->at(neighbors.first), points->at(neighbors.second)};

        MSD msd(error.error_fmin, error.error_MSD, error.max_counter, error.count);
        msd.find_min_MSD_mono(res, triplet);

        double x[3] = {points->at(midPoint).first,
                       points->at(neighbors.first).first,
                       points->at(neighbors.second).first};
        double y[3] = {points->at(midPoint).second,
                       points->at(neighbors.first).second,
                       points->at(neighbors.second).second};
        double k = res.k;
        double b = res.b;
        rho.emplace_back(0.0);
        for (int j = 0; j < 3; j++) {
            double tmp = (x[j] + k * y[j] - k * b) / (k * k + 1);
            rho[i] = max(rho[i], sqrt(pow(tmp - x[j], 2) + pow(k * tmp + b - y[j], 2)));
        }
    }

    sort(rho.begin(), rho.end());
    int rm = 4; // rm the largest and smallest rho
    for (int i = rm; i < rho.size() - rm; i++)
        sum_rho += rho[i];

    double r = sum_rho / (double) (points->size() - 2 * rm);
    if (r <= 0.01) return 0.01;
    if (r <= 0.1) return 0.1;
    else return 1.2 * r;
}

vector<pair<double, double >>
Vectorization_parallelograms::vectorization(vector<pair<float, float >> points,
                                            error_vectorization error) {
    if (error.rhoStep == -1)
        error.rhoStep = rhoStepEstimation(&points, error);

    SAY("rhoStep = %f\n", error.rhoStep);

    LinesDetection mainObject(error.rhoStep, error.thetaStep, error.interval);

    // MAIN CALL of parallelogram
    vector<Vec3d> line3dFirst;
    mainObject.parallelogram(&points, &line3dFirst);

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

    SAY("\t\t\t\tresult: \n");
    for (int i = 0; i <= 3; i++)
        SAY(" \t x = %f \t y = %f\n", result_k[i].first, result_k[i].second);
    SAY("\t\tERROR: %f\n", error_result);
    SAY("___________________________________________________________________________\n");
    return result_k;
}
