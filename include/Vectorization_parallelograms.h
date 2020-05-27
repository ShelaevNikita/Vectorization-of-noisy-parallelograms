#ifndef VECORIZATION_PARALLELOGRAMS_VECTORIZATION_PARALLELOGRAMS_H
#define VECORIZATION_PARALLELOGRAMS_VECTORIZATION_PARALLELOGRAMS_H

#include "opencv2/core.hpp"
#include <vector>

using namespace std;
class Vectorization_parallelograms{
public:
    struct error_vectorization {
        double  rhoStep,
                thetaStep,
                error_fmin,
                error_MSD,
                count;

        int max_counter,
                interval;
    };

    static vector<pair<double, double >>
    vectorization_parallelograms(const string &fname, error_vectorization error = {-1,
                                                                                   CV_PI / 3600.0,
                                                                                   1.0e-10,
                                                                                   1.0e-20,
                                                                                   1.0,
                                                                                   2500,
                                                                                   5});

    static vector<pair<double, double >>
    vectorization_parallelograms(vector<pair<float, float>> points,
                                 error_vectorization error = {-1,
                                                              CV_PI / 3600.0,
                                                              1.0e-10,
                                                              1.0e-20,
                                                              1.0,
                                                              2500,
                                                              5});

private:

    static pair<int, int> twoNeighbors(vector<pair<float, float>> *points, int midPoint);

    static double rhoStepEstimation(vector<pair<float, float>> *points, error_vectorization error);

    static vector<pair<double, double >>
    vectorization(vector<pair<float, float >> points,
                  error_vectorization error);
};


#endif //VECORIZATION_PARALLELOGRAMS_VECTORIZATION_PARALLELOGRAMS_H
