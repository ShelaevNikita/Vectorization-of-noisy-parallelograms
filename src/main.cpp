
#include <iostream>

#include <opencv2/core.hpp>

#include "Vectorization_OpenCV.cpp"
#include "MSD.cpp"

using namespace cv;
using namespace std;

class Main
{

private:

    static int parseString(string path, vector<string>* argv) 
    {
        ifstream inputFile(path);
        string temp = "";
        if (!inputFile.is_open())
            throw ios_base::failure("Incorrect filepath");
        else
        {
            while (inputFile >> temp)
                argv->push_back(temp);
        }
        inputFile.close();
        return (argv->size());
    }

public:

    static void main_foo(string input, vector<vector<pair<double, double>>>* result)
    {
        vector<string> argv;
        int number = parseString(input, &argv);
        for (int i = 0; i < number; i++)
        {
            Vectorization mainObject(0.1, CV_PI / 3600.0);

            vector<pair<float, float>> points;
            cout << "\n\t\t" << argv.at(i) << ":" << endl;
            mainObject.parsePoints(argv.at(i), &points);
            // MAIN CALL of parallelogram
            vector<Vec3d> line3dFirst;
            mainObject.parallelogram(&points, &line3dFirst);

            // basic console output to rank all the lines by their votes
            /*
            for (int i = 0; i < line3dFirst.size(); i++) 
                printf("votes:%d, rho:%.7f, theta:%.7f\n", (int)line3dFirst.at(i).val[0], line3dFirst.at(i).val[1],
                    line3dFirst.at(i).val[2]);
            */
            printf("Number of lines: %d \n", line3dFirst.size());

            // Finding 4 leaders of 1.txt and distributing its points between the 4 leaders
            vector<pair<double, double>> leaders;
            mainObject.fourLeaders(&line3dFirst, &leaders);

            vector<vector<pair<double, double>>> pointsSet(4);
            mainObject.pointDistribution(&leaders, &points, &pointsSet);

            cout << "\t" << argv.at(i) << " leaders: \n";
            for (int i = 0; i < 4; i++)
            {
                printf("Leader %d: %d points, y = %f * x + %f \n", i + 1, pointsSet[i].size(), leaders[i].first, leaders[i].second);
                for (int j = 0; j < pointsSet[i].size(); j++) 
                    cout << "\t" << pointsSet[i][j].first << " " << pointsSet[i][j].second << endl;
            }

            printf("\n");

            /*
                MSD - mean square deviation:
                error_fmin - error of the "Golden section"; default = 1.0e-15
                error_MSD - error in finding the best result for MSD; default = 1.0e-25
                max_counter - maximum number of iterations per cycle; default = 2500
                d - spread of values by b (b_min, b_max) for each line; d < 0.5; default = 0.15
            */

            MSD MSD_foo(1.0e-25, 1.0e-50, 3000, 0.15);

            vector<pair<double, double>> result_k;
            double error_result;
            MSD_foo.MSD_main(pointsSet, leaders, &result_k, error_result);

            cout << "\t\t error_result = " << error_result << "\n\t\t result:" << endl;
            for (int i = 0; i <= 3; i++)
                cout << " \t x = " << result_k[i].first << " \t  y = " << result_k[i].second << endl;
            result->push_back(result_k);
        }
    }
};

int main()
{
    vector<vector<pair<double, double>>> result;
    string input = ".\\Materials\\full.txt";
    cout << fixed;
    cout.precision(6);
    Main::main_foo(input, &result);
}