
#include <iostream>
#include <opencv2/core.hpp>
#include <cmath>

using namespace cv;
using namespace std;

class MSD
{

private:

    double error_fmin;
    double error_MSD;
    int max_counter;
    double d;

public:

    MSD(double error_fmin_Init, double error_MSD_Init, int max_counter_Init, double d_Init)
    {
        error_fmin = error_fmin_Init;
        error_MSD = error_MSD_Init;
        max_counter = max_counter_Init;
        d = d_Init;
    }

private:

    void find_pairs(double res[4][3], int(&result)[2][2])
    {
        double k = res[0][0];
        double min = k;
        int para2_1, para2_2;
        int para1_2 = 0;
        for (int i = 1; i <= 3; i++)
        {
            double a = abs(k - res[i][0]);
            if (a < min)
            {
                para1_2 = i;
                min = a;
            }
        }
        if (para1_2 == 1)
        {
            para2_1 = 2;
            para2_2 = 3;
        }
        else if (para1_2 == 2)
        {
            para2_1 = 1;
            para2_2 = 3;
        }
        else
        {
            para2_1 = 1;
            para2_2 = 2;
        }
        result[0][0] = 0;
        result[0][1] = para1_2;
        result[1][0] = para2_1;
        result[1][1] = para2_2;
    }

    void find_min_MSD_para(double(&res)[4][3], int para[2], vector<vector<pair<double, double>>> array_of_points)
    {
        int para_1 = para[0];
        int para_2 = para[1];
        double result = 0.0;
        double result_last = result + 2 * error_MSD;
        int counter = 0;
        double err = abs(result_last - result);
        double b1 = (res[para_1][1] + res[para_1][2]) / 2;
        double b2 = (res[para_2][1] + res[para_2][2]) / 2;
        double k = (res[para_1][0] + res[para_2][0]) / 2;
        while (err >= error_MSD && counter <= max_counter)
        {
             result = f_min_para(res[para_1][1], res[para_1][2], res[para_2][1], res[para_2][2], k, b1, b2,
                   array_of_points[para_1], array_of_points[para_2]);
             err = abs(result_last - result);
             result_last = result;
             counter++;
        }
        res[para_1][0] = k;
        res[para_2][0] = k;
        res[para_1][1] = b1;
        res[para_2][1] = b2;
        res[para_1][2] = sqrt(result);
        res[para_2][2] = sqrt(result);
    }

    double foo_fmin_para(double& k, double sum, double b1, double b2,
        vector<pair<double, double>> points1, vector<pair<double, double>> points2)
    {
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

    double f_min_para(double bmin1, double bmax1, double bmin2, double bmax2, double& k, double &b1, double &b2,
        vector<pair<double, double>> points1, vector<pair<double, double>> points2)
    {
        double sum = 0.0;
        double Lcenter, Rcenter, Lfcenter, Rfcenter, beg, end, BegBuf, EndBuf;
        for (int i = 0; i < points1.size(); i++)
            sum += pow(points1[i].first, 2);
        for (int j = 0; j < points2.size(); j++)
            sum += pow(points2[j].first, 2);
        double delta = (sqrt(5) - 1) / 2;
        for (int j = 0; j <= 1; j++)
        {
            if (j == 0)
            {
                beg = bmin1;
                end = bmax1;
                BegBuf = foo_fmin_para(k, sum, bmin1, b2, points1, points2);
                EndBuf = foo_fmin_para(k, sum, bmax1, b2, points1, points2);
            }
            else
            {
                beg = bmin2;
                end = bmax2;
                BegBuf = foo_fmin_para(k, sum, b1, bmin2, points1, points2);
                EndBuf = foo_fmin_para(k, sum, b1, bmax2, points1, points2);
            }
            for (int i = 0; i <= max_counter; i++)
            {
                double g = delta * (end - beg);
                Lcenter = end - g;
                if (j == 0) Lfcenter = foo_fmin_para(k, sum, Lcenter, b2, points1, points2);
                else Lfcenter = foo_fmin_para(k, sum, b1, Lcenter, points1, points2);
                Rcenter = beg + g;
                if (j == 0) Rfcenter = foo_fmin_para(k, sum, Rcenter, b2, points1, points2);
                else Rfcenter = foo_fmin_para(k, sum, b1, Rcenter, points1, points2);
                if (Lfcenter < Rfcenter)
                {
                    end = Rcenter;
                    EndBuf = Rfcenter;
                }
                else if (Lfcenter > Rfcenter)
                {
                    beg = Lcenter;
                    BegBuf = Lfcenter;
                }
                else
                {
                    end = Rcenter;
                    EndBuf = Rfcenter;
                    beg = Lcenter;
                    BegBuf = Lfcenter;
                }
                double r = fabs(end - beg);
                if (r < error_fmin)
                {
                    double f = (end + beg) / 2;
                    if (j == 0) b1 = f;
                    else b2 = f;
                    break;
                }
            }
        }
        return (foo_fmin_para(k, sum, b1, b2, points1, points2));
    }

    void find_min_MSD_mono(double(&res)[4][3], int number, vector<vector<pair<double, double>>> array_of_points)
    {
        double b = (res[number][1] + res[number][2]) / 2;
        double result = 0.0;
        double result_last = result + 2 * error_MSD;
        int counter = 0;
        double err = abs(result_last - result);
        while (err >= error_MSD && counter <= max_counter)
        {
            result = f_min_mono(res[number][1], res[number][2], res[number][0], b, array_of_points[number]);
            err = abs(result_last - result);
            result_last = result;
            counter++;
        }
        res[number][1] = b;
        res[number][2] = sqrt(result);
    }

    double foo_fmin_mono(double& k, double sum, double b, vector<pair<double, double>> points)
    {
        int size = points.size();
        int i = 0;
        double k1 = 0.0;
        for (i = 0; i < size; i++) 
            k1 += points[i].first * (points[i].second - b);
        k = k1 / sum;
        k1 = 0.0;
        for (i = 0; i < size; i++)
            k1 += pow((points[i].second - k * points[i].first - b), 2);
        return k1;
    }

    double f_min_mono(double bmin, double bmax, double &k, double &b, vector<pair<double, double>> points)
    {
        double sum = 0.0;
        for (int i = 0; i < points.size(); i++)
            sum += pow(points[i].first, 2);
        double Lcenter, Rcenter, Lfcenter, Rfcenter, beg, end, BegBuf, EndBuf;
        double delta = (sqrt(5) - 1) / 2;
        beg = bmin;
        end = bmax;
        BegBuf = foo_fmin_mono(k, sum, beg, points);
        EndBuf = foo_fmin_mono(k, sum, end, points);
        for (int i = 0; i <= max_counter; i++)
        {
            double g = delta * (end - beg);
            Lcenter = end - g;
            Lfcenter = foo_fmin_mono(k, sum, Lcenter, points);
            Rcenter = beg + g;
            Rfcenter = foo_fmin_mono(k, sum, Rcenter, points);
            if (Lfcenter < Rfcenter)
            {
                end = Rcenter;
                EndBuf = Rfcenter;
            }
            else if (Lfcenter > Rfcenter)
            {
                beg = Lcenter;
                BegBuf = Lfcenter;
            }
            else
            {
                end = Rcenter;
                EndBuf = Rfcenter;
                beg = Lcenter;
                BegBuf = Lfcenter;
            }
            double r = fabs(end - beg);
            if (r < error_fmin)
            {
                b = (end + beg) / 2;
                break;
            }
        }
        return (foo_fmin_mono(k, sum, b, points));
    }

    void intersection(double res[4][3], int para, double (&result)[4][2])
    {
        int counter = 0;
        int f = 0;
        for (int j = 0; j <= 1; j++)
        {
            for (int i = 1; i <= 3; i++)
            {
                if (i != para)
                {
                    double x = (res[i][1] - res[f][1]) / (res[f][0] - res[i][0]);
                    result[counter][0] = x;
                    result[counter][1] = res[f][0] * x + res[f][1];
                    counter++;
                }
            }
            f = para;
        }
        return;
    }

    double comparison(double res_para[4][3], double res_mono[4][3], int para[2][2], 
        double result_para[4][2], double result_mono[4][2], vector<pair<double, double>>* result)
    {
        double error_result = 0.0;
        for (int i = 0; i <= 1; i++)
        {
           int para_1 = para[i][0];
           int para_2 = para[i][1];
           double compare = res_para[para_1][2] / (res_mono[para_1][2] + res_mono[para_2][2]);
           cout << " Compare = MSD_para / (MSD_mono1 + MSD_mono2) = " << compare << endl;
           if (compare > 1)
           {
                result->push_back(make_pair(result_mono[para_1][0], result_mono[para_1][1]));
                printf("\t quadrangle\n");
                result->push_back(make_pair(result_mono[para_2][0], result_mono[para_2][1]));
                error_result += (res_mono[para_1][2] + res_mono[para_2][2]);
           }
           else
           {
                result->push_back(make_pair(result_para[para_1][0], result_para[para_1][1]));
                printf("\t parallelogram\n");
                result->push_back(make_pair(result_para[para_2][0], result_para[para_2][1]));
                error_result += res_para[para_1][2];
           }
        }
        return error_result;
    }

 public:

    void MSD_main(vector<vector<pair<double, double>>> Points, vector<pair<double, double>> Lines, 
        vector<pair<double, double>>* result, double &error_result)
    {
        double res_para[4][3], res_mono[4][3], result_para[4][2], result_mono[4][2];
        int i;
        for (i = 0; i <= 3; i++)
        {
            double bsr = Lines[i].second;
            res_para[i][0] = Lines[i].first;
            double log = log10(bsr);
            int f = floor(log);
            double g = 2 * pow(10, -f);
            if (log - f < d) g *= 5;
            if (log - f > (1 - d)) g *= 0.5;
            res_para[i][1] = bsr - g * bsr;
            res_para[i][2] = bsr + g * bsr;
            for (int k = 0; k <= 2; k++)
                res_mono[i][k] = res_para[i][k];
        }
        for (i = 0; i <= 3; i++)
            cout << "\t k = " << res_para[i][0] << "\t bmin = " << res_para[i][1] << "\t bmax = " << res_para[i][2] << endl;
        int squad_of_pairs[2][2];
        find_pairs(res_para, squad_of_pairs);
        for (i = 0; i <= 1; i++)
            find_min_MSD_para(res_para, squad_of_pairs[i], Points);
        for (i = 0; i <= 3; i++) 
            find_min_MSD_mono(res_mono, i, Points);
        intersection(res_para, squad_of_pairs[0][1], result_para);
        intersection(res_mono, squad_of_pairs[0][1], result_mono);
        for (i = 0; i <= 3; i++)
            cout << " k_para = " << res_para[i][0] << "\t b_para = " << res_para[i][1] << 
                "\t MSD_para = " << res_para[i][2] << endl;
        for (i = 0; i <= 3; i++)
            cout << "\t k_mono = " << res_mono[i][0] << "\t b_mono = " << res_mono[i][1] << 
                "\t MSD_mono = " << res_mono[i][2] << endl;
        vector<pair<double, double>> result_full;
        error_result = comparison(res_para, res_mono, squad_of_pairs, result_para, result_mono, &result_full);
        for (i = 0; i < 4; i++)
            result->push_back(result_full.at(i));
    }
};