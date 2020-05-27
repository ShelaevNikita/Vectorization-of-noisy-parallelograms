#include "debug.h"

                ;
#include <iostream>
#include <cmath>
#include "MSD.h"

MSD::MSD(double error_fmin, double error_MSD, int max_counter, double count) {
    this->error_fmin = error_fmin;
    this->error_MSD = error_MSD;
    this->max_counter = max_counter;
    this->count = count;
}


void MSD::find_pairs(line_result res[4], pairs_of_parallel_lines &result) {
    double k = res[0].k;
    double min = abs(k);
    for (int i = 1; i <= 3; i++) {
        double a = abs(k - res[i].k);
        if (a < min) {
            result.pair12 = i;
            min = a;
        }
    }
    if (result.pair12 == 2) {
        result.pair21 = 1;
        result.pair22 = 3;
    } else if (result.pair12 == 3) {
        result.pair21 = 1;
        result.pair22 = 2;
    }
}

void MSD::find_min_MSD_para(line_result (&res)[4], int pair1, int pair2,
                            vector<vector<pair<double, double>>> array_of_points) {
    double result = 0.0;
    double result_last = 2 * error_MSD;
    int counter = 0;
    bool flag = true;
    double err = abs(result_last - result);
    double b1 = (res[pair1].bmin + res[pair1].bmax) / 2;
    double b2 = (res[pair2].bmin + res[pair2].bmax) / 2;
    double k = (res[pair1].k + res[pair2].k) / 2;
    double sum = 0.0;
    for (auto &i : array_of_points[pair1])
        sum += pow(i.first, 2);
    for (auto &j : array_of_points[pair2])
        sum += pow(j.first, 2);
    double delta = (sqrt(5) - 1) / 2;
    if (k >= 1.0e5) flag = false;
    while (err >= error_MSD && counter <= max_counter) {
        result = f_min_para(flag, res[pair1].bmin, res[pair1].bmax, res[pair2].bmin,
                            res[pair2].bmax, k, b1, b2, array_of_points[pair1], array_of_points[pair2], sum, delta);
        if (b1 < res[pair1].bmin + 1) res[pair1].bmin -= 30;
        if (b1 > res[pair1].bmax - 1) res[pair1].bmax += 30;
        if (b2 < res[pair2].bmin + 1) res[pair2].bmin -= 30;
        if (b2 > res[pair2].bmax - 1) res[pair2].bmax += 30;
        err = abs(result_last - result);
        result_last = result;
        counter++;
    }
    res[pair1].b = b1;
    res[pair2].b = b2;
    res[pair1].MSD = sqrt(result);
    res[pair2].MSD = sqrt(result);
    if (flag) {
        res[pair1].k = k;
        res[pair2].k = k;
    } else {                    // x = b
        res[pair1].flag = false;
        res[pair2].flag = false;
    }
}

double MSD::foo_fmin_para(bool &flag, double &k, double sum, double b1, double b2,
                          vector<pair<double, double>> points1, vector<pair<double, double>> points2) {
    int size1 = points1.size();
    int size2 = points2.size();
    double k1 = 0.0;
    double k2 = 0.0;
    for (int i = 0; i < size1; i++)
        k1 += points1[i].first * (points1[i].second - b1);
    for (int j = 0; j < size2; j++)
        k2 += points2[j].first * (points2[j].second - b2);
    k = (k1 + k2) / sum;
    if (abs(k) > 1.0e5 || sum == 0) flag = false;                    // x = b
    k1 = 0.0;
    k2 = 0.0;
    for (int i = 0; i < size1; i++) {
        if (flag) k1 += pow((points1[i].second - k * points1[i].first - b1), 2);
        else k1 += pow((points1[i].first - b1), 2);
    }
    for (int j = 0; j < size2; j++) {
        if (flag) k2 += pow((points2[j].second - k * points2[j].first - b2), 2);
        else k2 += pow((points2[j].first - b2), 2);
    }
    return (k1 + k2);
}

double
MSD::f_min_para(bool &flag, double bmin1, double bmax1, double bmin2, double bmax2, double &k, double &b1, double &b2,
                const vector<pair<double, double>> &points1, const vector<pair<double, double>> &points2, double sum,
                double delta) {
    double Lcenter, Rcenter, Lfcenter, Rfcenter, beg, end;
    for (int j = 0; j <= 1; j++) {
        if (j == 0) {
            beg = bmin1;
            end = bmax1;
            foo_fmin_para(flag, k, sum, bmin1, b2, points1, points2);
            foo_fmin_para(flag, k, sum, bmax1, b2, points1, points2);
        } else {
            beg = bmin2;
            end = bmax2;
            foo_fmin_para(flag, k, sum, b1, bmin2, points1, points2);
            foo_fmin_para(flag, k, sum, b1, bmax2, points1, points2);
        }
        for (int i = 0; i <= max_counter; i++) {
            double g = delta * (end - beg);
            Lcenter = end - g;
            if (j == 0) Lfcenter = foo_fmin_para(flag, k, sum, Lcenter, b2, points1, points2);
            else Lfcenter = foo_fmin_para(flag, k, sum, b1, Lcenter, points1, points2);
            Rcenter = beg + g;
            if (j == 0) Rfcenter = foo_fmin_para(flag, k, sum, Rcenter, b2, points1, points2);
            else Rfcenter = foo_fmin_para(flag, k, sum, b1, Rcenter, points1, points2);
            if (Lfcenter < Rfcenter) {
                end = Rcenter;
            } else if (Lfcenter > Rfcenter) {
                beg = Lcenter;
            } else {
                end = Rcenter;
                beg = Lcenter;
            }
            double r = fabs(end - beg);
            if (r < error_fmin) {
                double f = (end + beg) / 2;
                if (j == 0) b1 = f;
                else b2 = f;
                break;
            }
        }
    }
    return (foo_fmin_para(flag, k, sum, b1, b2, points1, points2));
}

void MSD::find_min_MSD_mono(line_result &res,
                            const vector<pair<double, double>> &points) {
    double b = (res.bmin + res.bmax) / 2;
    double result = 0.0;
    double result_last = 2 * error_MSD;
    int counter = 0;
    bool flag = true;
    if (res.k >= 1.0e5) flag = false;
    double sum = 0.0;
    for (auto &i : points)
        sum += pow(i.first, 2);

    double err = abs(result_last - result);
    while (err >= error_MSD && counter <= max_counter) {
        result = f_min_mono(flag, res.bmin, res.bmax, res.k, b, points, sum);
        if (b < res.bmin + 1) res.bmin -= 30;
        if (b > res.bmax - 1) res.bmax += 30;
        err = abs(result_last - result);
        result_last = result;
        counter++;
    }
    if (!flag) res.flag = false;
    res.b = b;
    res.MSD = sqrt(result);
}

double MSD::foo_fmin_mono(bool &flag, double &k, double sum, double b, vector<pair<double, double>> points) {
    int size = points.size();
    int i;
    double k1 = 0.0;
    for (i = 0; i < size; i++)
        k1 += points[i].first * (points[i].second - b);
    k = k1 / sum;
    if (abs(k) > 1.0e5 || sum == 0) flag = false;                          // x = b
    k1 = 0.0;
    for (i = 0; i < size; i++) {
        if (flag) k1 += pow((points[i].second - k * points[i].first - b), 2);
        else k1 += pow((points[i].first - b), 2);
    }
    return k1;
}

double
MSD::f_min_mono(bool &flag, double bmin, double bmax, double &k, double &b, const vector<pair<double, double>> &points,
                double sum) {
    double Lcenter, Rcenter, Lfcenter, Rfcenter, beg, end;
    double delta = (sqrt(5) - 1) / 2;
    beg = bmin;
    end = bmax;
    foo_fmin_mono(flag, k, sum, beg, points);
    foo_fmin_mono(flag, k, sum, end, points);
    for (int i = 0; i <= max_counter; i++) {
        double g = delta * (end - beg);
        Lcenter = end - g;
        Lfcenter = foo_fmin_mono(flag, k, sum, Lcenter, points);
        Rcenter = beg + g;
        Rfcenter = foo_fmin_mono(flag, k, sum, Rcenter, points);
        if (Lfcenter < Rfcenter) {
            end = Rcenter;
        } else if (Lfcenter > Rfcenter) {
            beg = Lcenter;
        } else {
            end = Rcenter;
            beg = Lcenter;
        }
        double r = fabs(end - beg);
        if (r < error_fmin) {
            b = (end + beg) / 2;
            break;
        }
    }
    return (foo_fmin_mono(flag, k, sum, b, points));
}

vector<pair<double, double>> MSD::intersection(line_result res[4], int pair1) {
    vector<pair<double, double>> result;
    int f = 0;
    double x, y;
    for (int j = 0; j <= 1; j++) {
        double b = res[f].b;
        double k = res[f].k;
        for (int i = 1; i <= 3; i++) {
            if (i != pair1) {
                double b_i = res[i].b;
                double k_i = res[i].k;
                if (!res[i].flag) {
                    x = b_i;
                    y = k * x + b;
                } else if (!res[f].flag) {
                    x = b;
                    y = k_i * x + b_i;
                } else {
                    x = (b_i - b) / (k - k_i);
                    y = k * x + b;
                }
                result.emplace_back(x, y);
            }
        }
        f = pair1;
    }
    return result;
}

vector<pair<double, double>> MSD::comparison(line_result res_para[4], line_result res_mono[4],
                                             pairs_of_parallel_lines squad_pairs,
                                             vector<pair<double, double>> result_para,
                                             vector<pair<double, double>> result_mono, double &error_result) {
    vector<pair<double, double>> result;
    int pair_first = squad_pairs.pair11;
    int pair_second = squad_pairs.pair12;
    for (int i = 0; i <= 1; i++) {
        double compare = res_para[pair_first].MSD / (res_mono[pair_first].MSD + res_mono[pair_second].MSD);
        SAY(" Compare = MSD_para / (MSD_mono1 + MSD_mono2) = %f\n", compare);
        if (compare > count) {
            result.emplace_back(result_mono[pair_first].first, result_mono[pair_first].second);
            SAY("\t quadrangle\n");
            result.emplace_back(result_mono[pair_second].first, result_mono[pair_second].second);
            error_result += (res_mono[pair_first].MSD + res_mono[pair_second].MSD);
        } else {
            result.emplace_back(result_para[pair_first].first, result_para[pair_first].second);
            SAY("\t parallelogram\n");
            result.emplace_back(result_para[pair_second].first, result_para[pair_second].second);
            error_result += res_para[pair_first].MSD;
        }
        pair_first = squad_pairs.pair21;
        pair_second = squad_pairs.pair22;
    }
    return result;
}


vector<pair<double, double>> MSD::MSD_main(const vector<vector<pair<double, double>>> &points,
                                           vector<lineABC> Lines, double &error_result) {
    line_result res_para[4], res_mono[4];
    int i;
    for (i = 0; i <= 3; i++)
        if (abs(Lines[i].a / Lines[i].b) >= 1000 && abs(Lines[i].b) <= 0.001) res_para[i].flag = false;    // x = b
    for (i = 0; i <= 3; i++) {
        double bsr;
        if (res_para[i].flag) {
            bsr = -Lines[i].c / Lines[i].b;
            res_para[i].k = -Lines[i].a / Lines[i].b;
        } else {
            bsr = -Lines[i].c / Lines[i].a;
            res_para[i].k = 1.0e5;                                     // x = b
        }
        res_para[i].bmin = bsr - 20;
        res_para[i].bmax = bsr + 20;
        res_mono[i] = res_para[i];
    }
    pairs_of_parallel_lines squad_pairs;
    find_pairs(res_para, squad_pairs);
    find_min_MSD_para(res_para, squad_pairs.pair11, squad_pairs.pair12, points);
    find_min_MSD_para(res_para, squad_pairs.pair21, squad_pairs.pair22, points);
    for (i = 0; i <= 3; i++)
        find_min_MSD_mono(res_mono[i], points[i]);
    vector<pair<double, double>> result_para = intersection(res_para, squad_pairs.pair12);
    vector<pair<double, double>> result_mono = intersection(res_mono, squad_pairs.pair12);
    for (i = 0; i <= 3; i++)
        SAY("%fx + %f\t MSD_para = %f\n", res_para[i].k, res_para[i].b, res_para[i].MSD);
    vector<pair<double, double>> result_full =
            comparison(res_para, res_mono, squad_pairs, result_para, result_mono, error_result);
    return result_full;
}
