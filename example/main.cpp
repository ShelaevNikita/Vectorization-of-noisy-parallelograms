#include <iostream>
#include <Vectorization_parallelograms.h>

using namespace std;

int main(int argc, char** argv) {
    cout << fixed<<endl;
    cout.precision(6);
    for (int i = 1; i < argc; i++) {
         Vectorization_parallelograms::vectorization_parallelograms(argv[i]);
    }
}