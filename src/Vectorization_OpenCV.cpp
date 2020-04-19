
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
    static void parallelogram(vector<pair<float, float>>* inputPoints, vector<Vec3d>* outputLines,
        double rhoStep = 0.1, double thetaStep = CV_PI / 3600.0) {

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

        double rhoMin = - rhoApproximation * 1.2, rhoMax = rhoApproximation * 1.2; // let's add 20% to rhoApproximation
        double thetaMin = 0.0, thetaMax = CV_PI; // default thetaStep == (pi/180)/40 ~~ 0.025 degrees
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

    static void fourLeaders(vector<Vec3d>* houghOutput, vector<pair<double, double>>* leaders) {

        double rhoStep = 0.1, thetaStep = CV_PI / 3600.0;

        double rhoInterval = rhoStep * 30, thetaInterval = thetaStep * 600;
        double deltaRho, deltaTheta;

        leaders->push_back(make_pair(houghOutput->at(0).val[1], houghOutput->at(0).val[2]));
        int i = 1;
        bool inBunch;
        while (leaders->size() < 4 && i < houghOutput->size()) {
            inBunch = false;
            for (int j = 0; j < leaders->size(); j++) {
                deltaRho = abs(houghOutput->at(i).val[1] - leaders->at(j).first);
                deltaTheta = abs(houghOutput->at(i).val[2] - leaders->at(j).second);
                if (deltaRho <= rhoInterval && deltaTheta <= thetaInterval) {
                    inBunch = true;
                    break;
                }
            }
            if (!inBunch) {
                leaders->push_back(make_pair(houghOutput->at(i).val[1], houghOutput->at(i).val[2]));

            }
            i++;
        }

        double k, b;
        for (i = 0; i < 4; i++) {
            k = -cos(leaders->at(i).second) / sin(leaders->at(i).second); // k = -cos(theta) / sin(theta)
            b = leaders->at(i).first / sin(leaders->at(i).second); // b = r / sin(theta)
            leaders->at(i) = make_pair(k, b);
        }

    }

    static void pointDistribution(vector<pair<double, double>>* leaders, vector<pair<float, float>>* points, 
        vector <vector<pair<double, double>>>* pointsSet) {

        double rhoStep = 0.1;

        // Note: you should tie points to leaders from bottom to top only if you confirm that 4 leaders are actually different

        double interval = rhoStep * 3;
        double distance, x, y, k, b;
        for (int i = 0; i < points->size(); i++) {
            for (int j = 3; j > -1; j--) {
                x = points->at(i).first;
                y = points->at(i).second;
                k = leaders->at(j).first;
                b = leaders->at(j).second;
                distance = sqrt(pow((x + k * y - k * b) / (k * k + 1) - x, 2) + pow(k * (x + k * y - k * b) / (k * k + 1) + b - y, 2));
                if (distance <= interval) {
                    pointsSet->at(j).push_back(points->at(i));
                    break;
                }
            }
        }
    }
};


int main()
{
    
    vector<pair<float, float>> points1, points2, points3;
    Vectorization::parsePoints(".\\Materials\\1.txt", &points1);
    Vectorization::parsePoints(".\\Materials\\2.txt", &points2);
    Vectorization::parsePoints(".\\Materials\\3.txt", &points3);


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


    /*
    // Plotting 1.txt set of lines
    char hough_window1[] = "1.txt through Hough";

    Mat hough_image1 = Mat::zeros(700, 700, CV_8UC3);  // basic 700x700 image

    for (size_t i = 0; i < line3dFirst.size(); i++)
    {
        double rho = line3dFirst.at(i).val[1], theta = line3dFirst.at(i).val[2];
        Point pt1, pt2;
        double a = cos(theta), b = sin(theta);
        double x0 = a * rho, y0 = b * rho;
        pt1.x = cvRound(x0 + 200 * (-b) + 300); // сдвигаем точки на 300 пикселей для центрирования наших пучков
        pt1.y = cvRound(y0 + 200 * (a) + 300);
        pt2.x = cvRound(x0 - 200 * (-b) + 300);
        pt2.y = cvRound(y0 - 200 * (a) + 300);
        line(hough_image1, pt1, pt2, Scalar(0, 0, 40 * ((int)line3dFirst.at(i).val[0] - 1)), 1, LINE_AA);
        
    }

    imshow(hough_window1, hough_image1);
    moveWindow(hough_window1, 0, 10);
    //waitKey(0);

    */

    /*

    // Plotting 2.txt set of lines
    char hough_window2[] = "2.txt through Hough";

    Mat hough_image2 = Mat::zeros(700, 700, CV_8UC3);  // basic 700x700 image

    for (size_t i = 0; i < line3dSecond.size(); i++)
    {
        double rho = line3dSecond.at(i).val[1], theta = line3dSecond.at(i).val[2];
        Point pt1, pt2;
        double a = cos(theta), b = sin(theta);
        double x0 = a * rho, y0 = b * rho;
        pt1.x = cvRound(x0 + 200 * (-b) + 300); // сдвигаем точки на 300 пикселей для центрирования наших пучков
        pt1.y = cvRound(y0 + 200 * (a) + 300);
        pt2.x = cvRound(x0 - 200 * (-b) + 300);
        pt2.y = cvRound(y0 - 200 * (a) + 300);
        line(hough_image2, pt1, pt2, Scalar(0, 0, 40 * ((int)line3dSecond.at(i).val[0] - 1)), 1, LINE_AA);

    }

    imshow(hough_window2, hough_image2);
    moveWindow(hough_window2, 200, 60);
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
        pt1.x = cvRound(x0 + 200 * (-b) + 300); // сдвигаем точки на 300 пикселей для центрирования наших пучков
        pt1.y = cvRound(y0 + 200 * (a) + 300);
        pt2.x = cvRound(x0 - 200 * (-b) + 300);
        pt2.y = cvRound(y0 - 200 * (a) + 300);
        line(hough_image3, pt1, pt2, Scalar(0, 0, 40 * ((int)line3dThird.at(i).val[0] - 1)), 1, LINE_AA);

    }

    imshow(hough_window3, hough_image3);
    moveWindow(hough_window3, 400, 110);
    //waitKey(0);   // calling waitkey(0) only after all the console output!!! next waitkey is below

    */


    /*

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

    */


    // Finding 4 leaders of 1.txt and distribute the points between the 4 leaders
    vector<pair<double, double>> leaders1;
    Bunch::fourLeaders(&line3dFirst, &leaders1);

    vector <vector<pair<double, double>>> pointsSet1(4);
    Bunch::pointDistribution(&leaders1, &points1, &pointsSet1);

    cout << "1.txt leaders: \n";
    for (int i = 0; i < 4; i++) {
        printf("Leader %d: %d points, y = %fx + %f \n", i + 1, pointsSet1[i].size(), leaders1[i].first, leaders1[i].second);
        for (int j = 0; j < pointsSet1[i].size(); j++) {
            cout << pointsSet1[i][j].first << " " << pointsSet1[i][j].second << endl;
        }
        cout << endl;

    }

    // Plotting bunch leaders for 1.txt
    char hough_window4[] = "1.txt leaders";

    Mat hough_image4 = Mat::zeros(700, 700, CV_8UC3);  // basic 700x700 image

    for (size_t i = 0; i < 4; i++)
    {
        Point pt1, pt2;
        pt1.x = 0 + 300;
        pt1.y = cvRound(leaders1[i].first * (pt1.x - 300) + leaders1[i].second + 300);
        pt2.x = 100 + 300;
        pt2.y = cvRound(leaders1[i].first * (pt2.x - 300) + leaders1[i].second + 300);
        line(hough_image4, pt1, pt2, Scalar(0, 0, 255), 1, LINE_AA);

    }

    for (int i = 0; i < points1.size(); i++) {
        circle(hough_image4, Point(points1[i].first + 300, points1[i].second + 300), 1, Scalar(0, 255, 0), -1);
    }

    imshow(hough_window4, hough_image4);
    moveWindow(hough_window4, 600, 10);



    // Finding 4 leaders of 2.txt and distribute the points between the 4 leaders
    vector<pair<double, double>> leaders2;
    Bunch::fourLeaders(&line3dSecond, &leaders2);

    vector <vector<pair<double, double>>> pointsSet2(4);
    Bunch::pointDistribution(&leaders2, &points2, &pointsSet2);

    cout << "2.txt leaders: \n";
    for (int i = 0; i < 4; i++) {
        printf("Leader %d: %d points, y = %fx + %f \n", i + 1, pointsSet2[i].size(), leaders2[i].first, leaders2[i].second);
        for (int j = 0; j < pointsSet2[i].size(); j++) {
            cout << pointsSet2[i][j].first << " " << pointsSet2[i][j].second << endl;
        }
        cout << endl;

    }

    // Plotting bunch leaders for 2.txt
    char hough_window5[] = "2.txt leaders";

    Mat hough_image5 = Mat::zeros(700, 700, CV_8UC3);  // basic 700x700 image

    for (size_t i = 0; i < 4; i++)
    {
        Point pt1, pt2;
        pt1.x = 0 + 300;
        pt1.y = cvRound(leaders2[i].first * (pt1.x - 300) + leaders2[i].second + 300);
        pt2.x = 100 + 300;
        pt2.y = cvRound(leaders2[i].first * (pt2.x - 300) + leaders2[i].second + 300);
        line(hough_image5, pt1, pt2, Scalar(0, 0, 255), 1, LINE_AA);

    }

    for (int i = 0; i < points2.size(); i++) {
        circle(hough_image5, Point(points2[i].first + 300, points2[i].second + 300), 1, Scalar(0, 255, 0), -1);
    }

    imshow(hough_window5, hough_image5);
    moveWindow(hough_window5, 700, 10);


    // Finding 4 leaders of 3.txt and distribute the points between the 4 leaders
    vector<pair<double, double>> leaders3;
    Bunch::fourLeaders(&line3dThird, &leaders3);

    vector <vector<pair<double, double>>> pointsSet3(4);
    Bunch::pointDistribution(&leaders3, &points3, &pointsSet3);

    cout << "3.txt leaders: \n";
    for (int i = 0; i < 4; i++) {
        printf("Leader %d: %d points, y = %fx + %f \n", i + 1, pointsSet3[i].size(), leaders3[i].first, leaders3[i].second);
        for (int j = 0; j < pointsSet3[i].size(); j++) {
            cout << pointsSet3[i][j].first << " " << pointsSet3[i][j].second << endl;
        }
        cout << endl;

    }

    // Plotting bunch leaders for 3.txt
    char hough_window6[] = "3.txt leaders";

    Mat hough_image6 = Mat::zeros(700, 700, CV_8UC3);  // basic 700x700 image

    for (size_t i = 0; i < 4; i++)
    {
        Point pt1, pt2;
        pt1.x = 0 + 300;
        pt1.y = cvRound(leaders3[i].first * (pt1.x - 300) + leaders3[i].second + 300);
        pt2.x = 100 + 300;
        pt2.y = cvRound(leaders3[i].first * (pt2.x - 300) + leaders3[i].second + 300);
        line(hough_image6, pt1, pt2, Scalar(0, 0, 255), 1, LINE_AA);

    }

    for (int i = 0; i < points3.size(); i++) {
        circle(hough_image6, Point(points3[i].first + 300, points3[i].second + 300), 1, Scalar(0, 255, 0), -1);
    }

    imshow(hough_window6, hough_image6);
    moveWindow(hough_window6, 800, 10);

    waitKey(0); // only after all the console output

    return 0;
}
