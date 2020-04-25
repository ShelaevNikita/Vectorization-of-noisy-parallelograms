#include <iostream>

#include <opencv2/core.hpp>

#include "Vectorization_OpenCV.cpp"
#include "MSD.cpp"

using namespace cv;
using namespace std;

int main(int argc, char* argv[])
{
    Vectorization mainObject(0.1, CV_PI / 3600.0);

    vector<pair<float, float>> points;
    cout << argv[1] << endl;
    mainObject.parsePoints(argv[1], &points);

    // MAIN CALL of parallelogram
    vector<Vec3d> line3dFirst;
    mainObject.parallelogram(&points, &line3dFirst);

    // basic console output to rank all the lines by their votes
    for (int i = 0; i < line3dFirst.size(); i++) {
        printf("votes:%d, rho:%.7f, theta:%.7f\n", (int)line3dFirst.at(i).val[0], line3dFirst.at(i).val[1], line3dFirst.at(i).val[2]);
    }
    printf("Number of lines: %d \n\n", line3dFirst.size());


    // Finding 4 leaders of 1.txt and distributing its points between the 4 leaders
    vector<pair<double, double>> leaders;
    mainObject.fourLeaders(&line3dFirst, &leaders);

    vector<vector<pair<double, double>>> pointsSet(4);
    mainObject.pointDistribution(&leaders, &points, &pointsSet);

    cout << argv[1] << " leaders: \n";
    for (int i = 0; i < 4; i++) {
        printf("Leader %d: %d points, y = %fx + %f \n", i + 1, pointsSet[i].size(), leaders[i].first, leaders[i].second);
        for (int j = 0; j < pointsSet[i].size(); j++) {
            cout << pointsSet[i][j].first << " " << pointsSet[i][j].second << endl;
        }
    }

    printf("\n");

    // MSD - mean square deviation:
    // error_fmin - error of the "Golden section"; default = 1.0e-15
    // error_MSD - error in finding the best result for MSD; default = 1.0e-25
    // max_counter - maximum number of iterations per cycle; default = 1500
    // error_quadrangle - error of deviation of a parallelogram from a quadrilateral; required - 1.0; preferably - < 1.5
    // d - spread of values by b (b_min, b_max) for each line; d < 1; default = 0.25

    MSD MSD_foo(1.0e-15, 1.0e-25, 5000, 2.5, 0.3);

    double result[4][2];
    MSD_foo.MSD_main(pointsSet, leaders, result);

    for (int i = 0; i <= 3; i++) 
        cout << "\t x = " << result[i][0] << " \t y = " << result[i][1] << endl;

    return 0;
}