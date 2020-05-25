//
// Created by maratdin7 on 25.05.2020.
//

#ifndef VECORIZATION_PARALLELOGRAMS_VECTORIZATION_PARALLELOGRAMS_H
#define VECORIZATION_PARALLELOGRAMS_VECTORIZATION_PARALLELOGRAMS_H

#include <opencv2/core.hpp>
#include <vector>

using namespace std;
class Vectorization_parallelograms{
public:
    struct error_vectorization {
        double rhoStep,
                thetaStep,
                error_fmin,
                error_MSD,
                count;

        int max_counter,
                interval;
    };

    static vector<pair<double, double >>
    vectorization_parallelograms(const string &fname, error_vectorization error = {0.1,
                                                                                   CV_PI / 3600.0,
                                                                                   1.0e-10,
                                                                                   1.0e-20,
                                                                                   1.0,
                                                                                   2500,
                                                                                   5});

    static vector<pair<double, double >>
    vectorization_parallelograms(vector<pair<float, float>> points,
                                 error_vectorization error = {0.1,
                                                              CV_PI / 3600.0,
                                                              1.0e-10,
                                                              1.0e-20,
                                                              1.0,
                                                              2500,
                                                              5});

private:

    struct line_result {
        double k{}, bmin{}, bmax{}, b{};
        bool flag = true;
    };

    static void find_min_MSD(line_result (&res), const vector<pair<double, double>> &array_of_points);

    static double foo_fmin(bool &flag, double &k, double sum, double b, vector<pair<double, double>> points);

    static double
    f_min(bool &flag, double bmin, double bmax, double &k, double &b, const vector<pair<double, double>> &points);

    static pair<int, int> twoNeighbors(vector<pair<float, float>> *points, int midPoint);

    static double rhoStepEstimation(vector<pair<float, float>> *points);

    static vector<pair<double, double >>
    vectorization(vector<pair<float, float >> points,
                  error_vectorization error);
};


#endif //VECORIZATION_PARALLELOGRAMS_VECTORIZATION_PARALLELOGRAMS_H
