#include <iostream>

#include <opencv2/core.hpp>
#include <fstream>
#include <utility>

#include "Vectorization_OpenCV.cpp"
#include "MSD.cpp"

using namespace cv;
using namespace std;

struct error_vectorization {
    double rhoStep,thetaStep,
            error_fmin, error_MSD, count;
           
    int max_counter,
            interval;
};

class Vectorization_parallelograms {

public:
    static vector<vector<pair<double, double>>>
    vectorization_parallelograms(const string input, const string output, error_vectorization error = {0.1, CV_PI / 3600.0,
           1.0e-15, 1.0e-50, 2.5, 5000, 5}) {

        return vectorization(input, output, error);
    }

private:

    static vector<string> parseString_Input(string path, int& size) {
        vector<string> argv;
        ifstream inputFile(path);
        string temp = "";
        if (!inputFile.is_open())
            throw ios_base::failure("Incorrect filepath");
        else
            while (inputFile >> temp)
                argv.push_back(temp);
        inputFile.close();
        size = argv.size();
        return argv;
    }

    static void parseString_Output(string path, vector<vector<pair<double, double>>> result, int number) {
        ofstream outputFile(path);
        if (!outputFile.is_open())
            throw ios_base::failure("Incorrect filepath");
        else {
            outputFile << fixed;
            outputFile.precision(7);
            for (int j = 0; j < number; j++) {
                for (int i = 0; i <= 3; i++)
                    outputFile << result[j][i].first << " " << result[j][i].second << "\n";
                outputFile << "\n";
            }
            outputFile.close();
        }
    }

    static vector<pair<float, float>> parsePoints(string path) {
        vector<pair<float, float>> points;
        ifstream inputFile(path);
        pair<float, float> temp = make_pair(0.0, 0.0);
        if (!inputFile.is_open())
            throw ios_base::failure("Incorrect filepath");
        else
            while (inputFile >> temp.first >> temp.second)
                points.push_back(temp);
        inputFile.close();
        return points;
    }

    static vector<vector<pair<double, double>>>
    vectorization(string input, string output, error_vectorization error) {

        int number = 0;
        vector<string> argv = parseString_Input(input, number);
        vector<vector<pair<double, double>>> full_result;
        Vectorization mainObject(error.rhoStep, error.thetaStep, error.interval);

        for (int i = 0; i < number; i++) {
            vector<pair<float, float>> points = parsePoints(argv[i]);
                // MAIN CALL of parallelogram
            vector<Vec3d> line3dFirst;
            mainObject.parallelogram(&points, &line3dFirst);

            cout << "Number of lines: %lu" << line3dFirst.size() << endl;

            // Finding 4 leaders of 1.txt and distributing its points between the 4 leaders
            vector<pair<double, double>> leaders;
            mainObject.fourLeaders(&line3dFirst, &leaders);

            vector<vector<pair<double, double>>> pointsSet(4);
            mainObject.pointDistribution(&leaders, &points, &pointsSet);

            for (int i = 0; i < 4; i++) {
                printf("Leader %d: %lu points, y = %fx + %f \n", i + 1, pointsSet[i].size(), leaders[i].first,
                    leaders[i].second);
                for (auto& j : pointsSet[i])
                    printf("\t(%f, %f)\n", j.first, j.second);
            }

            printf("\n");

            /*
                MSD - mean square deviation:
                error_fmin - error of the "Golden section"; default = 1.0e-15
                error_MSD - error in finding the best result for MSD; default = 1.0e-25
                max_counter - maximum number of iterations per cycle; default = 2500
                count - spread of values by b (bmin, bmax) for each line; count >= 1; default count = 2
            */

            MSD MSD_foo(error.error_fmin, error.error_MSD, error.max_counter, error.count);
            double error_result = 0.0;
            vector<pair<double, double>> result_k = MSD_foo.MSD_main(pointsSet, leaders, error_result);

            printf("\t\t     result: \n");
            for (int i = 0; i <= 3; i++)
                printf(" \t x = %f \t y = %f\n", result_k[i].first, result_k[i].second);
            printf("\n\t\t error = %f\n", error_result);
            printf("___________________________________________________________________________\n");
            full_result.push_back(result_k);
        }
        parseString_Output(output, full_result, number);
        return full_result;
    }
};
