
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
    static void minmax(vector<pair<double, double>> line, int size, double& k, double& bmax, double& bmin)
    {
        double bsr;
        double bsum = 0.0;
        double ksum = 0.0;
        double b11 = 0.0;
        double b12 = 0.0;
        int sk = 0;
        int sb = 0;
        for (int i = 0; i < size; i++)
        {
            ksum += line[i].first;
            sk++;
            bsum += line[i].second;
            sb++;
        }
        k = ksum / sk;
        bsr = bsum / sb;
        sk = 0;
        sb = 0;
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
        bmin = b11 / sk;
        bmax = b12 / sb;
        if (size == 1)
        {
            bmin = bsr - 0.1 * bsr;
            bmax = bsr + 0.1 * bsr;
        }
        return;
    }

private:
    static double foo_fmin_para(double& k, double sum, double b1, double b2, 
        double** points1, int size1, double** points2, int size2)
    {
        int i = 0;
        int j = 0;
        double k1 = 0.0;
        double k2 = 0.0;
        for (i = 0; i < size1; i++)
        {
            double f = points1[i][0] * (points1[i][1] - b1);
            k1 += f;
        }
        for (j = 0; j < size2; j++)
        {
            double f = points2[j][0] * (points2[j][1] - b2);
            k2 += f;
        }
        k = (k1 + k2) / sum;
        k1 = 0.0;
        k2 = 0.0;
        for (i = 0; i < size1; i++)
        {
            double f = (points1[i][1] - k * points1[i][0] - b1);
            k1 += f * f;
        }
        for (j = 0; j < size2; j++)
        {
            double f = points2[j][1] - k * points2[j][0] - b2;
            k2 += f * f;
        }
        return k1 + k2;
    }

private:
    static double foo_fmin_mono(double& k, double sum, double b1, double** points1, int size1)
    {
        int i = 0;
        double k1 = 0.0;
        for (i = 0; i < size1; i++)
        {
            double f = points1[i][0] * (points1[i][1] - b1);
            k1 += f;
        }
        k = k1 / sum;
        k1 = 0.0;
        for (i = 0; i < size1; i++)
        {
            double f = (points1[i][1] - k * points1[i][0] - b1);
            k1 += f * f;
        }
        return k1;
    }

private:
    static double f_min_mono(double& bmin1, double bmax1, double& k, double** points1, int size1, double Tol)
    {
        double sum = 0.0;
        for (int i = 0; i < size1; i++)
        {
            double f = points1[i][0];
            sum += f * f;
        }
        double b1 = abs(bmax1 - bmin1) / 2;
        double Lcenter, Rcenter, Lfcenter, Rfcenter, beg, end, BegBuf, EndBuf;
        double delta = (sqrt(5) - 1) / 2;
        int counter;
        int maxcounts = 1000;
        beg = bmin1;
        end = bmax1;
        BegBuf = foo_fmin_mono(k, sum, beg, points1, size1);
        EndBuf = foo_fmin_mono(k, sum, end, points1, size1);
        counter = 0;
        for (int i = 0; i <= maxcounts; i++)
        {
            double g = delta * (end - beg);
            Lcenter = end - g;
            Lfcenter = foo_fmin_mono(k, sum, Lcenter, points1, size1);
            Rcenter = beg + g;
            Rfcenter = foo_fmin_mono(k, sum, Rcenter, points1, size1);
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
                b1 = (end + beg) / 2;
                break;
            }
            counter++;
        }
        bmin1 = b1;
        double result = foo_fmin_mono(k, sum, b1, points1, size1);
        return result;
    }

private:
    static double f_min_para(double& bmin1, double bmax1, double& bmin2, double bmax2, double& k,
        double** points1, int size1, double** points2, int size2, double Tol)
    {
        double sum = 0.0;
        for (int i = 0; i < size1; i++)
        {
            double f = points1[i][0];
            sum += f * f;
        }
        for (int j = 0; j < size2; j++)
        {
            double f = points2[j][0];
            sum += f * f;
        }
        double b2 = abs(bmax2 - bmin2) / 2;
        double b1 = abs(bmax1 - bmin1) / 2;
        double Lcenter, Rcenter, Lfcenter, Rfcenter, beg, end, BegBuf, EndBuf;
        double delta = (sqrt(5) - 1) / 2;
        int counter;
        int maxcounts = 1000;
        for (int j = 0; j <= 1; j++)
        {
            if (j == 0)
            {
                beg = bmin1;
                end = bmax1;
                BegBuf = foo_fmin_para(k, sum, bmin1, b2, points1, size1, points2, size2);
                EndBuf = foo_fmin_para(k, sum, bmax1, b2, points1, size1, points2, size2);
            }
            else
            {
                beg = bmin2;
                end = bmax2;
                BegBuf = foo_fmin_para(k, sum, b1, bmin2, points1, size1, points2, size2);
                EndBuf = foo_fmin_para(k, sum, b1, bmax2, points1, size1, points2, size2);
            }
            counter = 0;
            for (int i = 0; i <= maxcounts; i++)
            {
                double g = delta * (end - beg);
                Lcenter = end - g;
                if (j == 0) Lfcenter = foo_fmin_para(k, sum, Lcenter, b2, points1, size1, points2, size2);
                else Lfcenter = foo_fmin_para(k, sum, b1, Lcenter, points1, size1, points2, size2);
                Rcenter = beg + g;
                if (j == 0) Rfcenter = foo_fmin_para(k, sum, Rcenter, b2, points1, size1, points2, size2);
                else Rfcenter = foo_fmin_para(k, sum, b1, Rcenter, points1, size1, points2, size2);
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
                counter++;
            }
        }
        bmin1 = b1;
        bmin2 = b2;
        double result = foo_fmin_para(k, sum, b1, b2, points1, size1, points2, size2);
        return result;
    }

private:
    static int point(vector<pair<double, double>> lines, int size_lines, vector<pair<float, float>> points, 
        int size_points, double error, int* result)
    {
        int* number = new int[size_points];
        for (int j = 0; j < size_points; j++)
        {
            result[j] = -1;
            number[j] = 0;
        }
        for (int i = 0; i < size_lines; i++)
        {
            for (int j = 0; j < size_points; j++)
            {
                double y = lines[i].first * points[j].first + lines[i].second;
                double y_error = abs(y - points[j].second);
                if (y_error < error)
                {
                    int f = number[j];
                    number[j] = f + 1;
                }
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
        int size = count;
        do
        {
            for (int j = 0; j < size_points; j++)
            {
                if (number[j] == max)
                {
                    number[j] = -1;
                    result[j] = j;
                    count--;
                }
            }
            max--;
        } while (count > 0 && max > 0);
        return size;
    }

private:
    static void array_of_points(vector<pair<float, float>> points, int size_points, int size, int* array, double** result)
    {
        for (int i = 0; i < size; i++) result[i] = new double[2];
        int j = 0;
        for (int i = 0; i < size_points; i++)
        {
            if (array[i] == i && j < size)
            {
                result[j][0] = points[i].first;
                result[j][1] = points[i].second;
                j++;
            }
        }
        return;
    }

private:
    static double intersection(double k1, double k2, double b1, double b2, double& x)
    {
        x = (b2 - b1) / (k1 - k2);
        double y = k1 * x + b1;
        return y;
    }

private:
    static void foo_main(vector<pair<float, float>> Points, int size_points,
        vector<pair<double, double>> Lines1, vector<pair<double, double>> Lines2, vector<pair<double, double>> Lines3, 
        vector<pair<double, double>> Lines4, int size_lines,
        double error_points, double error, double ERROR, double (&result)[4][2])
    {
        double res[4][3];
        minmax(Lines1, size_lines, res[0][0], res[0][1], res[0][2]);
        minmax(Lines2, size_lines, res[1][0], res[1][1], res[1][2]);
        minmax(Lines3, size_lines, res[2][0], res[2][1], res[2][2]);
        minmax(Lines4, size_lines, res[3][0], res[3][1], res[3][2]);
        double k = res[0][0];
        double min = k;
        int d = 0;
        for (int i = 1; i <= 3; i++)
        {
            double a = abs(k - res[i][0]);
            if (a < min)
            {
                d = i;
                min = a;
            }
        }
        if (d != 1)
        {
            for (int i = 0; i < 3; i++)
            {
                min = res[1][i];
                res[1][i] = res[d][i];
                res[d][i] = min;
            }
        }
        int* array = new int[size_points];
        int size1 = point(Lines1, size_lines, Points, size_points, error_points, array);
        double** points1 = new double* [size1];
        array_of_points(Points, size_points, size1, array, points1);
        int size2 = point(Lines2, size_lines, Points, size_points, error_points, array);
        double** points2 = new double* [size2];
        array_of_points(Points, size_points, size2, array, points2);
        int size3 = point(Lines3, size_lines, Points, size_points, error_points, array);
        double** points3 = new double* [size3];
        array_of_points(Points, size_points, size3, array, points3);
        int size4 = point(Lines4, size_lines, Points, size_points, error_points, array);
        double** points4 = new double* [size4];
        array_of_points(Points, size_points, size4, array, points4);
        double k12 = (res[0][0] + res[1][0]) / 2;
        double result1 = 0.0;
        double result_last = result1 + 2 * ERROR;
        int counter = 0;
        int max_counter = 1000;
        double err = abs(result_last - result1);
        while (err >= ERROR && counter <= max_counter)
        {
            if (d == 1) result1 = f_min_para(res[0][1], res[0][2], res[1][1], res[1][2], k12, points1, size1, points2, size2, error);
            else if (d == 2) result1 = f_min_para(res[0][1], res[0][2], res[1][1], res[1][2], k12, points1, size1, points3, size3, error);
            else result1 = f_min_para(res[0][1], res[0][2], res[1][1], res[1][2], k12, points1, size1, points4, size4, error);
            err = abs(result_last - result1);
            result_last = result1;
            counter++;
        }
        cout << result1 << " " << res[0][1] << " " << res[1][1] << " " << counter << " " << k12 << endl;
        double result2 = 0.0;
        result_last = result2 + 2 * ERROR;
        double k34 = (res[2][0] + res[3][0]) / 2;
        counter = 0;
        err = abs(result_last - result2);
        while (err >= ERROR && counter <= max_counter)
        {
            if (d == 1) result2 = f_min_para(res[2][1], res[2][2], res[3][1], res[3][2], k34, points3, size3, points4, size4, error);
            else if (d == 2) result2 = f_min_para(res[2][1], res[2][2], res[3][1], res[3][2], k34, points2, size2, points4, size4, error);
            else result2 = f_min_para(res[2][1], res[2][2], res[3][1], res[3][2], k34, points2, size2, points3, size3, error);
            err = abs(result_last - result2);
            result_last = result2;
            counter++;
        }
        cout << result2 << " " << res[2][1] << " " << res[3][1] << " " << counter << " " << k34 << endl;
        double x1;
        double y1 = intersection(k12, k34, res[0][1], res[2][1], x1);
        double x2;
        double y2 = intersection(k12, k34, res[0][1], res[3][1], x2);
        double x3;
        double y3 = intersection(k12, k34, res[1][1], res[2][1], x3);
        double x4;
        double y4 = intersection(k12, k34, res[1][1], res[3][1], x4);
       result[0][0] = x1;
       result[0][1] = y1;
       result[1][0] = x2;
       result[1][1] = y2;
       result[2][0] = x3;
       result[2][1] = y3;
       result[3][0] = x4;
       result[3][1] = y4;
    }

public:
    static void SKO_main(vector<pair<float, float>> Points, int size_points, vector<pair<double, double>> Line1, 
        vector<pair<double, double>> Line2, vector<pair<double, double>> Line3,
        vector<pair<double, double>> Line4, int size_lines, double (&result)[4][2])
    {
        double error_points = 0.1;
        double error_fmin = 1.0e-10;
        double ERROR = 1.0e-50;
        foo_main(Points, size_points, Line1, Line2, Line3, Line4, size_lines, error_points, error_fmin, ERROR, result);
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
    SKO::SKO_main(points1, points1.size(), bunches[0], bunches[1], bunches[2], bunches[3], 1, result);
    for (int i = 0; i <= 3; i++) cout << "\t x = " << result[i][0] << "\t y = " << result[i][1] << endl;
    waitKey(0); // only after all the console output

}
