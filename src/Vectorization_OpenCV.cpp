
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
        double rhoStep = 0.1, double thetaStep = CV_PI / 7200.0) {

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

        double rhoMin = 0.0, rhoMax = rhoApproximation * 1.2; // let's add 20% to rhoApproximation
        double thetaMin = 0.0, thetaMax = CV_PI * 2.0; // default thetaStep == (pi/180)/40 ~~ 0.025 degrees
        HoughLinesPointSet(pointVector, lines, 600, 1,   // 600 линий необходимо для 2.txt, иначе для него 4 пучок оказывается слишком мал.
            rhoMin, rhoMax, rhoStep,                     // Однако в этом случае в 1.txt появляется куча слишком отклоняющихся линий. 
            thetaMin, thetaMax, thetaStep);              // Но эти линии фильтруются нашими границами в if'ах и выходит нормально.

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

/*
    For 2.txt value ranges for distinguishing 4 different bunches are as follows:
              rho>=  rho<=  theta>=  theta<=
    1 bunch:  23     39     2.33     2.5
    2 bunch:  62     73     5.59     5.66
    3 bunch:  54     68     0.67     0.88
    4 bunch:  29     41     3.81     4.0

    для того чтобы получить 4 пучка для 2.txt, надо подставить эти значения в if'ы в верхнюю функцию
*/

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

    imshow(hough_window1, hough_image1);
    moveWindow(hough_window1, 0, 10);
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
        pt1.x = cvRound(x0 + 300 * (-b) + 300); // сдвигаем точки на 300 пикселей для центрирования наших пучков
        pt1.y = cvRound(y0 + 300 * (a)+300);
        pt2.x = cvRound(x0 - 300 * (-b) + 300);
        pt2.y = cvRound(y0 - 300 * (a)+300);
        line(hough_image3, pt1, pt2, Scalar(0, 0, 40 * ((int)line3dThird.at(i).val[0] - 1)), 1, LINE_AA);

    }

    imshow(hough_window3, hough_image3);
    moveWindow(hough_window3, 400, 110);
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


    waitKey(0); // only after all the console output

    return 0;
}
