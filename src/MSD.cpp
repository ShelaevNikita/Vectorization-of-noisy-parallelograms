#include <iostream>
#include <opencv2/core.hpp>
#include <cmath>

using namespace cv;
using namespace std;

class MSD {

private:
    double error_fmin;
    double error_MSD;
    int max_counter;
    double d;
    double count;

public:
    MSD(double error_fmin_Init, double error_MSD_Init, int max_counter_Init, double d_Init, double count_Init) {
        error_fmin = error_fmin_Init;
        error_MSD = error_MSD_Init;
        max_counter = max_counter_Init;
        d = d_Init;
        count = count_Init;
    }

private:
    static void find_pairs(double res[4][3], pair<pair<int, int>, pair<int, int>>(&result)) {
        double k = res[0][0];
        double min = k;
        for (int i = 1; i <= 3; i++) {
            double a = abs(k - res[i][0]);
            if (a < min) {
                result.first.second = i;
                min = a;
            }
        }
        if (result.first.second == 2) {
            result.second.first = 1;
            result.second.second = 3;
        } else if (result.first.second == 3) {
            result.second.first = 1;
            result.second.second = 2;
        }
    }

    void find_min_MSD_para(double (&res)[4][3], pair<int, int> pair1, vector<vector<pair<double, double>>> array_of_points) {
        double result = 0.0;
        double result_last = result + 2 * error_MSD;
        int counter = 0;
        double err = abs(result_last - result);
        double b1 = (res[pair1.first][1] + res[pair1.first][2]) / 2;
        double b2 = (res[pair1.second][1] + res[pair1.second][2]) / 2;
        double k = (res[pair1.first][0] + res[pair1.second][0]) / 2;
        while (err >= error_MSD && counter <= max_counter) {
             result = f_min_para(res[pair1.first][1], res[pair1.first][2], res[pair1.second][1], 
                 res[pair1.second][2], k, b1, b2, array_of_points[pair1.first], array_of_points[pair1.second]);
             err = abs(result_last - result);
             result_last = result;
             counter++;
        }
        res[pair1.first][0] = k;
        res[pair1.second][0] = k;
        res[pair1.first][1] = b1;
        res[pair1.second][1] = b2;
        res[pair1.first][2] = sqrt(result);
        res[pair1.second][2] = sqrt(result);
    }

    static double foo_fmin_para(double &k, double sum, double b1, double b2,
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
        k1 = 0.0;
        k2 = 0.0;
        for (int i = 0; i < size1; i++)
            k1 += pow((points1[i].second - k * points1[i].first - b1), 2);
        for (int j = 0; j < size2; j++)
            k2 += pow((points2[j].second - k * points2[j].first - b2), 2);
        return (k1 + k2);
    }

    double f_min_para(double bmin1, double bmax1, double bmin2, double bmax2, double &k, double &b1, double &b2,
            vector<pair<double, double>> points1, vector<pair<double, double>> points2) {
        double sum = 0.0;
        double Lcenter, Rcenter, Lfcenter, Rfcenter, beg, end, BegBuf, EndBuf;
        for (int i = 0; i < points1.size(); i++)
            sum += pow(points1[i].first, 2);
        for (int j = 0; j < points2.size(); j++)
            sum += pow(points2[j].first, 2);
        double delta = (sqrt(5) - 1) / 2;
        for (int j = 0; j <= 1; j++) {
            if (j == 0) {
                beg = bmin1;
                end = bmax1;
                BegBuf = foo_fmin_para(k, sum, bmin1, b2, points1, points2);
                EndBuf = foo_fmin_para(k, sum, bmax1, b2, points1, points2);
            } else {
                beg = bmin2;
                end = bmax2;
                BegBuf = foo_fmin_para(k, sum, b1, bmin2, points1, points2);
                EndBuf = foo_fmin_para(k, sum, b1, bmax2, points1, points2);
            }
            for (int i = 0; i <= max_counter; i++) {
                double g = delta * (end - beg);
                Lcenter = end - g;
                if (j == 0) Lfcenter = foo_fmin_para(k, sum, Lcenter, b2, points1, points2);
                    else Lfcenter = foo_fmin_para(k, sum, b1, Lcenter, points1, points2);
                Rcenter = beg + g;
                if (j == 0) Rfcenter = foo_fmin_para(k, sum, Rcenter, b2, points1, points2);
                    else Rfcenter = foo_fmin_para(k, sum, b1, Rcenter, points1, points2);
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
                    double f = (end + beg) / 2;
                    if (j == 0) b1 = f;
                    else b2 = f;
                    break;
                }
            }
        }
        return (foo_fmin_para(k, sum, b1, b2, points1, points2));
    }

    void find_min_MSD_mono(double(&res)[4][3], int number,
                vector<vector<pair<double, double>>> array_of_points) {
        double b = (res[number][1] + res[number][2]) / 2;
        double result = 0.0;
        double result_last = result + 2 * error_MSD;
        int counter = 0;
        double err = abs(result_last - result);
        while (err >= error_MSD && counter <= max_counter) {
            result = f_min_mono(res[number][1], res[number][2], res[number][0], b, array_of_points[number]);
            err = abs(result_last - result);
            result_last = result;
            counter++;
        }
        res[number][1] = b;
        res[number][2] = sqrt(result);
    }

   static double foo_fmin_mono(double &k, double sum, double b, vector<pair<double, double>> points) {
        int size = points.size();
        int i;
        double k1 = 0.0;
        for (i = 0; i < size; i++) 
            k1 += points[i].first * (points[i].second - b);
        k = k1 / sum;
        k1 = 0.0;
        for (i = 0; i < size; i++)
            k1 += pow((points[i].second - k * points[i].first - b), 2);
        return k1;
    }

    double f_min_mono(double bmin, double bmax, double &k, double &b, vector<pair<double, double>> points) {
        double sum = 0.0;
        for (int i = 0; i < points.size(); i++)
            sum += pow(points[i].first, 2);
        double Lcenter, Rcenter, Lfcenter, Rfcenter, beg, end, BegBuf, EndBuf;
        double delta = (sqrt(5) - 1) / 2;
        beg = bmin;
        end = bmax;
        BegBuf = foo_fmin_mono(k, sum, beg, points);
        EndBuf = foo_fmin_mono(k, sum, end, points);
        for (int i = 0; i <= max_counter; i++) {
            double g = delta * (end - beg);
            Lcenter = end - g;
            Lfcenter = foo_fmin_mono(k, sum, Lcenter, points);
            Rcenter = beg + g;
            Rfcenter = foo_fmin_mono(k, sum, Rcenter, points);
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
        return (foo_fmin_mono(k, sum, b, points));
    }

    static vector<pair<double, double>> intersection(double res[4][3], int pair1) {
        vector<pair<double, double>> result;
        int f = 0;
        double x, y;
        for (int j = 0; j <= 1; j++) {
            for (int i = 1; i <= 3; i++) {
                if (i != pair1) {
                    x = (res[i][1] - res[f][1]) / (res[f][0] - res[i][0]);
                    y = res[f][0] * x + res[f][1];
                    result.push_back(make_pair(x, y));
                }
            }
            f = pair1;
        }
        return result;
    }

    static vector<pair<double, double>> comparison(double res_para[4][3], double res_mono[4][3], 
                pair<pair<int, int>, pair<int, int>> pair1, vector<pair<double, double>> result_para, 
                vector<pair<double, double>> result_mono, double &error_result) {
        vector<pair<double, double>> result;
        int pair_first = pair1.first.first;
        int pair_second = pair1.first.second;
        for (int i = 0; i <= 1; i++) {
           double compare = res_para[pair_first][2] / (res_mono[pair_first][2] + res_mono[pair_second][2]);
           /*!*/ cout << " Compare = MSD_para / (MSD_mono1 + MSD_mono2) = " << compare << endl; 
           if (compare > 1) {
                result.push_back(make_pair(result_mono[pair_first].first, result_mono[pair_first].second));
                /*!*/ printf("\t quadrangle\n");
                result.push_back(make_pair(result_mono[pair_second].first, result_mono[pair_second].second));
                error_result += (res_mono[pair_first][2] + res_mono[pair_second][2]);
           } else {
                result.push_back(make_pair(result_para[pair_first].first, result_para[pair_first].second));
                /*!*/ printf("\t parallelogram\n");
                result.push_back(make_pair(result_para[pair_second].first, result_para[pair_second].second));
                error_result += res_para[pair_first][2];
           }
           pair_first = pair1.second.first;
           pair_second = pair1.second.second;
        }
        return result;
    }

 public:
     vector<pair<double, double>> MSD_main(vector<vector<pair<double, double>>> Points, 
             vector<pair<double, double>> Lines, double &error_result) {
        double res_para[4][3], res_mono[4][3];
        int i;
        for (i = 0; i <= 3; i++) {
            double bsr = Lines[i].second;
            res_para[i][0] = Lines[i].first;
            double log = log10(bsr);
            int f = floor(log);
            double g = count * pow(10, -f);
            if (log - f < d) g *= 5;
            if (log - f > (1 - d)) g *= 0.5;
            res_para[i][1] = bsr - g * bsr;
            res_para[i][2] = bsr + g * bsr;
            for (int k = 0; k <= 2; k++)
                res_mono[i][k] = res_para[i][k];
        }
        /*!*/ for (i = 0; i <= 3; i++)
        /*!*/    cout << "\t k = " << res_para[i][0] << "\t bmin = " << res_para[i][1] << "\t bmax = " << res_para[i][2] << endl;
        pair<pair<int, int>, pair<int, int>> squad_of_pairs = make_pair(make_pair(0, 1), make_pair(2, 3));
        find_pairs(res_para, squad_of_pairs);
        find_min_MSD_para(res_para, squad_of_pairs.first, Points);
        find_min_MSD_para(res_para, squad_of_pairs.second, Points);
        for (i = 0; i <= 3; i++) 
            find_min_MSD_mono(res_mono, i, Points);
        vector<pair<double, double>> result_para = intersection(res_para, squad_of_pairs.first.second);
        vector<pair<double, double>> result_mono = intersection(res_mono, squad_of_pairs.first.second);
        /*!*/ for (i = 0; i <= 3; i++)
        /*!*/    cout << " k_para = " << res_para[i][0] << "\t b_para = " << res_para[i][1] <<
        /*!*/        "\t MSD_para = " << res_para[i][2] << endl;
        /*!*/ for (i = 0; i <= 3; i++)
        /*!*/   cout << "\t k_mono = " << res_mono[i][0] << "\t b_mono = " << res_mono[i][1] <<
        /*!*/       "\t MSD_mono = " << res_mono[i][2] << endl;
        vector<pair<double, double>> result_full = 
            comparison(res_para, res_mono, squad_of_pairs, result_para, result_mono, error_result);
        return result_full;
     }
};
