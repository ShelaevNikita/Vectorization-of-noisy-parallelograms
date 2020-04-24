#include <iostream>

#include <opencv2/core.hpp>

#include "Vectorization_OpenCV.cpp"

using namespace cv;
using namespace std;

int main(int argc, char* argv[])
{
    Vectorization mainObject(0.1, CV_PI / 3600.0);

    vector<pair<float, float>> points1;
    cout << argv[1] << endl;
    mainObject.parsePoints(argv[1], &points1);

    // MAIN CALL of parallelogram
    vector<Vec3d> line3dFirst;
    mainObject.parallelogram(&points1, &line3dFirst);

    // basic console output to rank all the lines by their votes
    for (int i = 0; i < line3dFirst.size(); i++) {
        printf("votes:%d, rho:%.7f, theta:%.7f\n", (int)line3dFirst.at(i).val[0], line3dFirst.at(i).val[1], line3dFirst.at(i).val[2]);
    }
    printf("Number of lines: %d \n\n", line3dFirst.size());


    // Finding 4 leaders of 1.txt and distributing its points between the 4 leaders
    vector<pair<double, double>> leaders1;
    mainObject.fourLeaders(&line3dFirst, &leaders1);

    vector <vector<pair<double, double>>> pointsSet1(4);
    mainObject.pointDistribution(&leaders1, &points1, &pointsSet1);

    cout << argv[1] << " leaders: \n";
    for (int i = 0; i < 4; i++) {
        printf("Leader %d: %d points, y = %fx + %f \n", i + 1, pointsSet1[i].size(), leaders1[i].first, leaders1[i].second);
        for (int j = 0; j < pointsSet1[i].size(); j++) {
            cout << pointsSet1[i][j].first << " " << pointsSet1[i][j].second << endl;
        }
        cout << endl;

    }

    return 0;
}