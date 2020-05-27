#ifndef VECORIZATION_PARALLELOGRAMS_MSD_H
#define VECORIZATION_PARALLELOGRAMS_MSD_H

#include <utility>
#include <vector>
#include "structures.h"

using namespace std;

class MSD {

private:
    double error_fmin;
    double error_MSD;
    int max_counter;
    double count;

public:
    MSD(double error_fmin, double error_MSD, int max_counter, double count);


    struct pairs_of_parallel_lines {
        int pair11 = 0;
        int pair12 = 1;
        int pair21 = 2;
        int pair22 = 3;
    };
    struct line_result {
        double k{}, bmin{}, bmax{}, b{}, MSD{};
        bool flag = true;
    };
private:

    static void find_pairs(line_result res[4], pairs_of_parallel_lines &result);

    void find_min_MSD_para(line_result (&res)[4], int pair1, int pair2,
                           vector<vector<pair<double, double>>> array_of_points);

    static double foo_fmin_para(bool &flag, double &k, double sum, double b1, double b2,
                                vector<pair<double, double>> points1, vector<pair<double, double>> points2);

    double
    f_min_para(bool &flag, double bmin1, double bmax1, double bmin2, double bmax2, double &k, double &b1, double &b2,
               const vector<pair<double, double>> &points1, const vector<pair<double, double>> &points2, double sum,
               double delta);

public:

    void find_min_MSD_mono(line_result &res, const vector<pair<double, double>> &points);

private:

    static double foo_fmin_mono(bool &flag, double &k, double sum, double b, vector<pair<double, double>> points);

    double
    f_min_mono(bool &flag, double bmin, double bmax, double &k, double &b, const vector<pair<double, double>> &points,
               double sum);

    static vector<pair<double, double>> intersection(line_result res[4], int pair1);

    vector<pair<double, double>> comparison(line_result res_para[4], line_result res_mono[4],
                                            pairs_of_parallel_lines squad_pairs,
                                            vector<pair<double, double>> result_para,
                                            vector<pair<double, double>> result_mono, double &error_result);

public:
    vector<pair<double, double>> MSD_main(const vector<vector<pair<double, double>>> &points,
                                          vector<lineABC> Lines, double &error_result);

};

#endif //VECORIZATION_PARALLELOGRAMS_MSD_H
