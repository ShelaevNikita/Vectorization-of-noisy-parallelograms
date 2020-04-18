
#include <iostream>
#include <opencv2/core.hpp>
#include <cmath>

using namespace cv;
using namespace std;

/*
vector<pair<double, double>> P = { { 1.47164, 75.693 }, { 1.68387, 76.2508 }, { 2.94818, 76.4739 },
{ 3.79413, 76.9201 }, { 4.64168, 77.3663 }, { 5.91321, 77.8125 }, { 6.76436, 78.2586 },
{ 8.0422, 78.8164 }, { 8.04411, 78.9279 }, { 8.89507, 79.151 }, { 9.75144, 79.5972 }, { 11.0208, 79.4857 },
{ 11.4095, 78.0355 }, { 11.8096, 77.1432 }, { 12.2081, 76.2508 }, { 12.4008, 75.5815 }, { 11.5546, 75.3584 },
{ 10.2843, 74.9122 }, { 9.43571, 74.466 }, { 8.1698, 74.0198 }, { 7.32486, 73.5736 }, { 6.48152, 73.1274 },
{ 5.22202, 72.6812 }, { 4.38229, 72.235 }, { 2.71349, 72.3466 }, { 2.301, 73.3505 }, { 1.88804, 74.6891 } };
*/

class SKO
{
private:
    static void minmax(pair<double, double> line, double (&res)[3])
    {
        double bsr = line.second;
        /*
        double bsum = 0.0;
        double ksum = 0.0;
        double b11 = 0.0;
        double b12 = 0.0;
        int sk = 0;
        int sb = 0;
        for (int i = 0; i < size; i++)
        {
            ksum += line[i].first;
            bsum += line[i].second;
        }
        res[0] = ksum / size;
        bsr = bsum / size;
        for (int i = 0; i < size; i++)
        {
            double b0 = line[i].second;
            if (b0 < bsr)
            {
                b11 += b0;
                sk++;
            }
            else if (b0 > bsr)
            {
                b12 += b0;
                sb++;
            }
        }
        res[1] = b11 / sk;
        res[2] = b12 / sb;
        */
        res[0] = line.first;
        //if (size == 1)
        //{
            double log = log10(bsr);
            int f = floor(log);
            double g = pow(10, -f);
            if (log - f < 0.25) g *= 5;
            if (log - f > 0.75) g *= 0.5;
            res[1] = bsr - g * bsr;
            res[2] = bsr + g * bsr;
        //}
        return;
    }

private:
    static void find_pairs(double res[4][3], int(&result)[2][2])
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

private:
    static void find_points(vector<pair<double, double>> points, vector<pair<double, double>> lines,
        vector<pair<double, double>>* result, double error)
    {
        int size_points = points.size();
        vector<int> number;
        for (int j = 0; j < size_points; j++) number.push_back(0);
        for (int i = 0; i < lines.size(); i++)
        {
            for (int j = 0; j < size_points; j++)
            {
                double y = lines[i].first * points[j].first + lines[i].second;
                double y_error = abs(y - points[j].second);
                if (y_error < error) number[j] += 1;
            }
        }
        int max = 0;
        int count = 0;
        for (int j = 0; j < size_points; j++)
        {
            int num = number[j];
            if (num > 0) count++;
            if (num > max) max = num;
        }
        do
        {
            for (int j = 0; j < size_points; j++)
            {
                if (number[j] == max)
                {
                    number[j] = -1;
                    result->push_back(points[j]);
                    count--;
                }
            }
            max--;
        } while (count > 0 && max > 0);
        return;
    }

private:
    static void find_min_SKO_para(double(&res)[4][3], int para[2],
        vector<vector<pair<double, double>>> array_of_points, double error, double ERROR)
    {
        double b1 = 0.0;
        double b2 = 0.0;
        int para_1 = para[0];
        int para_2 = para[1];
        double result = 0.0;
        double result_last = result + 2 * ERROR;
        int counter = 0;
        int max_counter = 1500;
        double err = abs(result_last - result);
        double k = (res[para_1][0] + res[para_2][0]) / 2;
        while (err >= ERROR && counter <= max_counter)
        {
             result = f_min_para(res[para_1][1], res[para_1][2], res[para_2][1], res[para_2][2], k, b1, b2,
                   array_of_points[para_1], array_of_points[para_2], error);
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

private:
    static double foo_fmin_para(double& k, double sum, double b1, double b2,
        vector<pair<double, double>> points1, vector<pair<double, double>> points2)
    {
        int size1 = points1.size();
        int size2 = points2.size();
        double k1 = 0.0;
        double k2 = 0.0;
        for (int i = 0; i < size1; i++)
        {
            double f = points1[i].first * (points1[i].second - b1);
            k1 += f;
        }
        for (int j = 0; j < size2; j++)
        {
            double f = points2[j].first * (points2[j].second - b2);
            k2 += f;
        }
        k = (k1 + k2) / sum;
        k1 = 0.0;
        k2 = 0.0;
        for (int i = 0; i < size1; i++)
        {
            double f = points1[i].second - k * points1[i].first - b1;
            k1 += f * f;
        }
        for (int j = 0; j < size2; j++)
        {
            double f = points2[j].second - k * points2[j].first - b2;
            k2 += f * f;
        }
        return k1 + k2;
    }

private:
    static double f_min_para(double bmin1, double bmax1, double bmin2, double bmax2, double& k, double &b1, double &b2,
        vector<pair<double, double>> points1, vector<pair<double, double>> points2, double Tol)
    {
        double sum = 0.0;
        for (int i = 0; i < points1.size(); i++)
        {
            double f = points1[i].first;
            sum += f * f;
        }
        for (int j = 0; j < points2.size(); j++)
        {
            double f = points2[j].first;
            sum += f * f;
        }
        b1 = (bmax1 + bmin1) / 2;
        b2 = (bmax2 + bmin2) / 2;
        double Lcenter, Rcenter, Lfcenter, Rfcenter, beg, end, BegBuf, EndBuf;
        double delta = (sqrt(5) - 1) / 2;
        int maxcounts = 1500;
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

            for (int i = 0; i <= maxcounts; i++)
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
                if (r < 2 * Tol)
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

private:
    static void find_min_SKO_mono(double(&res)[4][3], int number,
        vector<vector<pair<double, double>>> array_of_points, double error, double ERROR)
    {
        double b = 0.0;
        double result = 0.0;
        double result_last = result + 2 * ERROR;
        int counter = 0;
        int max_counter = 1500;
        double err = abs(result_last - result);
        while (err >= ERROR && counter <= max_counter)
        {
            result = f_min_mono(res[number][1], res[number][2], res[number][0], b, array_of_points[number], error);
            err = abs(result_last - result);
            result_last = result;
            counter++;
        }
        res[number][1] = b;
        res[number][2] = sqrt(result);
    }

private:
    static double foo_fmin_mono(double& k, double sum, double b, vector<pair<double, double>> points)
    {
        int size = points.size();
        int i = 0;
        double k1 = 0.0;
        for (i = 0; i < size; i++)
        {
            double f = points[i].first * (points[i].second - b);
            k1 += f;
        }
        k = k1 / sum;
        k1 = 0.0;
        for (i = 0; i < size; i++)
        {
            double f = points[i].second - k * points[i].first - b;
            k1 += f * f;
        }
        return k1;
    }

private:
    static double f_min_mono(double bmin, double bmax, double &k, double &b, vector<pair<double, double>> points, double Tol)
    {
        double sum = 0.0;
        for (int i = 0; i < points.size(); i++)
        {
            double f = points[i].first;
            sum += f * f;
        }
        b = (bmax + bmin) / 2;
        double Lcenter, Rcenter, Lfcenter, Rfcenter, beg, end, BegBuf, EndBuf;
        double delta = (sqrt(5) - 1) / 2;
        int maxcounts = 1500;
        beg = bmin;
        end = bmax;
        BegBuf = foo_fmin_mono(k, sum, beg, points);
        EndBuf = foo_fmin_mono(k, sum, end, points);
        for (int i = 0; i <= maxcounts; i++)
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
            if (r < 2 * Tol)
            {
                b = (end + beg) / 2;
                break;
            }
        }
        return (foo_fmin_mono(k, sum, b, points));
    }

private:
    static void intersection(double res[4][3], int para, double (&result)[4][2])
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

private:
    static void comparison(double res_para[4][3], double res_mono[4][3], int para[2][2], 
        double result_para[4][2], double result_mono[4][2], double error, double (&result)[4][2])
    {
        int counter = 0;
        for (int i = 0; i <= 1; i++)
        {
            int para_1 = para[i][0];
            int para_2 = para[i][1];
            double compare = res_para[para_1][2] / (res_mono[para_1][2] + res_mono[para_2][2]);
            cout << " compare = SKO_para / (SKO_mono1 + SKO_mono2) = " << compare << endl;
            if (compare > error)
            {
                result[counter][0] = result_mono[para_1][0];
                result[counter][1] = result_mono[para_1][1];
                counter++;
                result[counter][0] = result_mono[para_2][0];
                result[counter][1] = result_mono[para_2][1];
                printf(" quadrangle\n");
            }
            else
            {
                result[counter][0] = result_para[para_1][0];
                result[counter][1] = result_para[para_1][1];
                counter++;
                result[counter][0] = result_para[para_2][0];
                result[counter][1] = result_para[para_2][1];
                printf(" parallelogram\n");
            }
            counter++;
        }
    }

public:
    static void SKO_main(vector<vector<pair<double, double>>> Points, vector<pair<double, double>> Lines, double(&result)[4][2])
    {
        // Погрешности:
        //double error_points = 0.1;
        double error_fmin = 1.0e-10;
        double ERROR_SKO = 1.0e-50;
        double error_quadrangle = 1.5;
        //
        double res_para[4][3], res_mono[4][3], result_para[4][2], result_mono[4][2];
        for (int i = 0; i <= 3; i++)
        {
            minmax(Lines[i], res_para[i]);
            res_mono[i][0] = res_para[i][0];
            res_mono[i][1] = res_para[i][1];
            res_mono[i][2] = res_para[i][2];
        }
        for (int i = 0; i <= 3; i++)
            cout << "\t k = " << res_para[i][0] << "\t bmin = " << res_para[i][1] << "\t bmax = " << res_para[i][2] << endl;
        int squad_of_pairs[2][2];
        find_pairs(res_para, squad_of_pairs);
        /*
        vector<vector<pair<double, double>>> array_of_points;
        for (int i = 0; i <= 3; i++)
        {
            vector<pair<double, double>> array;
            find_points(Points, Lines[i], &array, error_points);
            cout << "\t points.size = " << array.size() << endl;
            array_of_points.push_back(array);
        }
        */
        find_min_SKO_para(res_para, squad_of_pairs[0], Points, error_fmin, ERROR_SKO);
        find_min_SKO_para(res_para, squad_of_pairs[1], Points, error_fmin, ERROR_SKO);
        for (int i = 0; i <= 3; i++) find_min_SKO_mono(res_mono, i, Points, error_fmin, ERROR_SKO);
        intersection(res_para, squad_of_pairs[0][1], result_para);
        intersection(res_mono, squad_of_pairs[0][1], result_mono);
        for (int i = 0; i <= 3; i++)
            cout << " k_para = " << res_para[i][0] << "\t b_para = " << res_para[i][1] << 
                "\t SKO_para = " << res_para[i][2] << endl;
        for (int i = 0; i <= 3; i++)
            cout << "\t k_mono = " << res_mono[i][0] << "\t b_mono = " << res_mono[i][1] << 
                "\t SKO_mono = " << res_mono[i][2] << endl;
        comparison(res_para, res_mono, squad_of_pairs, result_para, result_mono, error_quadrangle, result);
    }
};

void line(vector<pair<double, double>> points, pair<double, double> *lines)
{
    int size = points.size();
    double k, b;
    int f = 0;
    int s = 1;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (j >= size - 2 && i == 1 && i != j)
            {
                k = (points[i].second - points[j].second) / (points[i].first - points[j].first);
                b = points[i].second - k * points[i].first;
                lines->first = k;
                lines->second = b;
            }
        }
    }
}

void lines(vector<pair<double, double>>* Lines, vector<vector<pair<double, double>>>* points)
{
    vector<pair<double, double>> P1 = { { 1.47164, 75.693 }, { 1.68387, 76.2508 }, { 2.94818, 76.4739 },
    { 3.79413, 76.9201 }, { 4.64168, 77.3663 }, { 5.91321, 77.8125 }, { 6.76436, 78.2586 },
    { 8.0422, 78.8164 }, { 8.04411, 78.9279 }, { 8.89507, 79.151 }, { 9.75144, 79.5972 } };
    vector<pair<double, double>> P3 = { { 12.4008, 75.5815 }, { 11.5546, 75.3584 }, { 10.2843, 74.9122 },
    { 9.43571, 74.466 }, { 8.1698, 74.0198 }, { 7.32486, 73.5736 }, { 6.48152, 73.1274 },
    { 5.22202, 72.6812 }, { 4.38229, 72.235 } };
    vector<pair<double, double>> P2 = { { 9.75144, 79.5972 }, { 11.0208, 79.4857 }, { 11.4095, 78.0355 },
    { 11.8096, 77.1432 }, { 12.2081, 76.2508 }, { 12.4008, 75.5815 } };
    vector<pair<double, double>> P4 = { { 2.71349, 72.3466 }, { 2.301, 73.3505 }, { 1.88804, 74.6891 }, { 1.47164, 75.693 } };
    pair<double, double> Lines_1, Lines_2, Lines_3, Lines_4;
    line(P1, &Lines_1);
    line(P2, &Lines_2);
    line(P3, &Lines_3);
    line(P4, &Lines_4);
    Lines->push_back(Lines_1);
    Lines->push_back(Lines_2);
    Lines->push_back(Lines_3);
    Lines->push_back(Lines_4);
    points->push_back(P1);
    points->push_back(P2);
    points->push_back(P3);
    points->push_back(P4);
}

int main()
{
    vector<vector<pair<double, double>>> Points;
    vector<pair<double, double>> Lines;
    lines(&Lines, &Points);                    
    double result[4][2];
    SKO::SKO_main(Points, Lines, result);
    cout<<fixed;
    cout.precision(10);
    for (int i = 0; i <= 3; i++) cout << "\t x = " << result[i][0] << "\t y = " << result[i][1] << endl;
}