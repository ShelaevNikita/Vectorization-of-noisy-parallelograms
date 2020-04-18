
#include <iostream>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <fstream>
#include <cmath>

using namespace cv;
using namespace std;

class Vectorization {

public:
    static void parallelogram(vector<pair<float, float>>* inputPoints, vector<Vec3d>* outputLines) {

        Mat lines;
        vector<Point2f> pointVector;

        double rhoApproximation = 0.0;
        double distance = 0.0;

        for (int i = 0; i < inputPoints->size(); i++)
        {
            pointVector.push_back(Point2f(inputPoints->at(i).first, inputPoints->at(i).second));

            distance = hypot(inputPoints->at(i).first, inputPoints->at(i).second);
            rhoApproximation = rhoApproximation < distance ? distance : rhoApproximation;
        }

        double rhoMin = 0.0, rhoMax = rhoApproximation * 1.2, rhoStep = 0.1; // let's add 20% to rhoApproximation
        double thetaMin = 0.0, thetaMax = CV_PI * 2.0, thetaStep = CV_PI / 7200.0; //step ~~ 0.1 degree (now 0.025 degrees)
        HoughLinesPointSet(pointVector, lines, 400, 1,
            rhoMin, rhoMax, rhoStep,
            thetaMin, thetaMax, thetaStep);

        lines.copyTo(*outputLines);
    }


    static void parsePoints(string path, vector<pair<float, float>>* points) {
        ifstream inputFile(path);
        pair<float, float> temp = make_pair(0.0, 0.0);

        if (!inputFile.is_open()) {
            throw ios_base::failure("Incorrect filepath");
        }
        else {
            while (inputFile >> temp.first >> temp.second) {
                points->push_back(temp);
            }
        }

        inputFile.close();
    }

};

class Bunch {

public:
    static void bunching(vector<Vec3d>* houghOutput, vector<vector<pair<double, double>>>* bunches) {
        double k, b;

        for (int i = 0; i < houghOutput->size(); i++) {
            k = -cos(houghOutput->at(i).val[2]) / sin(houghOutput->at(i).val[2]); // val[1] == rho, val[2] == theta
            b = houghOutput->at(i).val[1] / sin(houghOutput->at(i).val[2]);

            if (houghOutput->at(i).val[1] >= 60.0 && houghOutput->at(i).val[1] <= 71.0 &&
                houghOutput->at(i).val[2] >= 1.89 && houghOutput->at(i).val[2] <= 2.1) {

                bunches->at(0).push_back(make_pair(k, b));
            }
            else if (houghOutput->at(i).val[1] >= 27.0 && houghOutput->at(i).val[1] <= 36.5 &&
                houghOutput->at(i).val[2] >= 5.0 && houghOutput->at(i).val[2] <= 5.2) {

                bunches->at(1).push_back(make_pair(k, b));
            }
            else if (houghOutput->at(i).val[1] >= 26.0 && houghOutput->at(i).val[1] <= 41.0 &&
                houghOutput->at(i).val[2] >= 0.3 && houghOutput->at(i).val[2] <= 0.415) {

                bunches->at(2).push_back(make_pair(k, b));
            }
            else if (houghOutput->at(i).val[1] >= 55.0 && houghOutput->at(i).val[1] <= 72.0 &&
                houghOutput->at(i).val[2] >= 3.45 && houghOutput->at(i).val[2] <= 3.56) {

                bunches->at(3).push_back(make_pair(k, b));
            }
        }
    }
};

class SKO
{
private:
    static void minmax(pair<double, double> line, double(&res)[3])
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
    static double f_min_para(double bmin1, double bmax1, double bmin2, double bmax2, double& k, double& b1, double& b2,
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
    static double f_min_mono(double bmin, double bmax, double& k, double& b, vector<pair<double, double>> points, double Tol)
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
    static void intersection(double res[4][3], int para, double(&result)[4][2])
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
        double result_para[4][2], double result_mono[4][2], double error, double(&result)[4][2])
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

int main()
{
    vector<pair<float, float>> points1, points2, points3;
    Vectorization::parsePoints("1.txt", &points1);
    Vectorization::parsePoints("2.txt", &points2);
    Vectorization::parsePoints("3.txt", &points3);

    cout << points1.size() << " " << points2.size() << " " << points3.size() << endl; //Debug only!

    // MAIN CALL of parallelogram
    vector<Vec3d> line3dFirst, line3dSecond, line3dThird; //vector of generated lines in (votes, rho, theta) form
    Vectorization::parallelogram(&points1, &line3dFirst);
    Vectorization::parallelogram(&points2, &line3dSecond);
    Vectorization::parallelogram(&points3, &line3dThird);

    // basic console output to rank all the lines by their votes
    for (int i = 0; i < line3dFirst.size(); i++) {
        printf("votes:%d, rho:%.7f, theta:%.7f\n", (int)line3dFirst.at(i).val[0], line3dFirst.at(i).val[1], line3dFirst.at(i).val[2]);
    }
    printf("Number of lines: %d \n\n", line3dFirst.size());

    for (int i = 0; i < line3dSecond.size(); i++) {
        printf("votes:%d, rho:%.7f, theta:%.7f\n", (int)line3dSecond.at(i).val[0], line3dSecond.at(i).val[1], line3dSecond.at(i).val[2]);
    }
    printf("Number of lines: %d \n\n", line3dSecond.size());

    for (int i = 0; i < line3dThird.size(); i++) {
        printf("votes:%d, rho:%.7f, theta:%.7f\n", (int)line3dThird.at(i).val[0], line3dThird.at(i).val[1], line3dThird.at(i).val[2]);
    }
    printf("Number of lines: %d \n\n", line3dThird.size());

    // Plotting 1.txt set of lines
    char hough_window1[] = "1.txt through Hough";

    Mat hough_image1 = Mat::zeros(700, 700, CV_8UC3);  // basic 700x700 image

    for (size_t i = 0; i < line3dFirst.size(); i++)
    {
        double rho = line3dFirst.at(i).val[1], theta = line3dFirst.at(i).val[2];
        Point pt1, pt2;
        double a = cos(theta), b = sin(theta);
        double x0 = a * rho, y0 = b * rho;
        pt1.x = cvRound(x0 + 300 * (-b) + 300); // сдвигаем точки на 300 пикселей для центрирования наших пучков
        pt1.y = cvRound(y0 + 300 * (a) + 300);
        pt2.x = cvRound(x0 - 300 * (-b) + 300);
        pt2.y = cvRound(y0 - 300 * (a) + 300);
        line(hough_image1, pt1, pt2, Scalar(0, 0, 40 * ((int)line3dFirst.at(i).val[0] - 1)), 1, LINE_AA);
        
    }
    //imshow(hough_window1, hough_image1);
    //moveWindow(hough_window1, 0, 10);
    //waitKey(0);

    // Plotting 2.txt set of lines
    char hough_window2[] = "2.txt through Hough";

    Mat hough_image2 = Mat::zeros(700, 700, CV_8UC3);  // basic 700x700 image

    for (size_t i = 0; i < line3dSecond.size(); i++)
    {
        double rho = line3dSecond.at(i).val[1], theta = line3dSecond.at(i).val[2];
        Point pt1, pt2;
        double a = cos(theta), b = sin(theta);
        double x0 = a * rho, y0 = b * rho;
        pt1.x = cvRound(x0 + 300 * (-b) + 300); // сдвигаем точки на 300 пикселей для центрирования наших пучков
        pt1.y = cvRound(y0 + 300 * (a) + 300);
        pt2.x = cvRound(x0 - 300 * (-b) + 300);
        pt2.y = cvRound(y0 - 300 * (a) + 300);
        line(hough_image2, pt1, pt2, Scalar(0, 0, 40 * ((int)line3dSecond.at(i).val[0] - 1)), 1, LINE_AA);

    }
    //imshow(hough_window2, hough_image2);
    //moveWindow(hough_window2, 200, 60);
    //waitKey(0);

    // Plotting 3.txt set of lines
    char hough_window3[] = "3.txt through Hough";

    Mat hough_image3 = Mat::zeros(700, 700, CV_8UC3);  // basic 700x700 image

    for (size_t i = 0; i < line3dThird.size(); i++)
    {
        double rho = line3dThird.at(i).val[1], theta = line3dThird.at(i).val[2];
        Point pt1, pt2;
        double a = cos(theta), b = sin(theta);
        double x0 = a * rho, y0 = b * rho;
        pt1.x = cvRound(x0 + 300 * (-b) + 300); // сдвигаем точки на 300 пикселей для центрирования наших пучков
        pt1.y = cvRound(y0 + 300 * (a)+300);
        pt2.x = cvRound(x0 - 300 * (-b) + 300);
        pt2.y = cvRound(y0 - 300 * (a)+300);
        line(hough_image3, pt1, pt2, Scalar(0, 0, 40 * ((int)line3dThird.at(i).val[0] - 1)), 1, LINE_AA);

    }

    //imshow(hough_window3, hough_image3);
    //moveWindow(hough_window3, 400, 110);
    //waitKey(0);   // calling waitkey(0) only after all the console output!!! next waitkey is below

    // Bunching all together
    vector<vector<pair<double, double>>> bunches(4);

    Bunch::bunching(&line3dFirst, &bunches);
    
    for (int i = 0; i < bunches[0].size(); i++) {
        cout << bunches[0][i].first << " " << bunches[0][i].second << endl;
    }
    printf("First bunch contains %d lines \n\n", bunches[0].size());

    for (int i = 0; i < bunches[1].size(); i++) {
        cout << bunches[1][i].first << " " << bunches[1][i].second << endl;
    }
    printf("Second bunch contains %d lines \n\n", bunches[1].size());

    for (int i = 0; i < bunches[2].size(); i++) {
        cout << bunches[2][i].first << " " << bunches[2][i].second << endl;
    }
    printf("Third bunch contains %d lines \n\n", bunches[2].size());

    for (int i = 0; i < bunches[3].size(); i++) {
        cout << bunches[3][i].first << " " << bunches[3][i].second << endl;
    }
    printf("Fourth bunch contains %d lines \n\n", bunches[3].size());

    double result[4][2];
    // SKO_main принимает вектор векторов точек, вектор прямых и двумерный массив результата
    //SKO::SKO_main(Points, Lines, result);
    for (int i = 0; i <= 3; i++) cout << "\t x = " << result[i][0] << "\t y = " << result[i][1] << endl;
    waitKey(0); // only after all the console output
}
