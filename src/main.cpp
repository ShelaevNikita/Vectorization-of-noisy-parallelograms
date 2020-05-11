#include <iostream>
#include "Vectorization_parallelograms.cpp"
using namespace std;

int main() {
    string input = ".\\Materials\\full.txt";
    string output = ".\\Materials\\result.txt";
    cout << fixed;
    cout.precision(7);
        Vectorization_parallelograms::vectorization_parallelograms(input, output);
}
