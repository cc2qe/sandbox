#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

class histogram {

public:    
    histogram(float minValue, float maxValue, int numBins);
    ~histogram();
    // count a data value
    int count(float value);
    int writeFile(std::string filename);

private:
    int numBins;
    float minValue;
    float maxValue;
    float binWidth;
    long overBin;
    long underBin;
    
    std::vector<long> bins;

};


