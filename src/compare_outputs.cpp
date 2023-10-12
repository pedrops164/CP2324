#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>

// compile with g++ src/compare_outputs.cpp -O2 -o compare -std=c++11

bool isEqual(const std::string& strA, const std::string& strB, int n) {
    // Remove the dot and minus sign
    std::string procStrA = strA;
    std::string procStrB = strB;
    procStrA.erase(std::remove(procStrA.begin(), procStrA.end(), '.'), procStrA.end());
    procStrA.erase(std::remove(procStrA.begin(), procStrA.end(), '-'), procStrA.end());
    procStrB.erase(std::remove(procStrB.begin(), procStrB.end(), '.'), procStrB.end());
    procStrB.erase(std::remove(procStrB.begin(), procStrB.end(), '-'), procStrB.end());

    // Consider the smallest string's length if it's less than n
    int compareLength = std::min(n, static_cast<int>(std::min(procStrA.length(), procStrB.length())));

    // Compare the two strings up to the compareLength
    return procStrA.substr(0, compareLength) == procStrB.substr(0, compareLength);
}

bool filesAreEqual(std::string filePath1, std::string filePath2, int n) {
    std::ifstream file1(filePath1);
    std::ifstream file2(filePath2);

    if (!file1.is_open() || !file2.is_open()) {
        std::cerr << "Error opening one of the files." << std::endl;
        return false;
    }

    std::string line1, line2;
    // Skip the header
    std::getline(file1, line1);
    std::getline(file2, line2);

    int lineNum = 1; // As we already skipped the header

    bool allEqual = true;

    while (std::getline(file1, line1) && std::getline(file2, line2)) {
        std::istringstream iss1(line1);
        std::istringstream iss2(line2);

        std::string value1, value2;
        bool lineEqual = true;
        while (iss1 >> value1 && iss2 >> value2) {
            if (!isEqual(value1, value2, n)) {
                lineEqual = false;
                break;
            }
        }

        if (!lineEqual) {
            std::cout << "Difference detected at line " << lineNum << ":\n";
            std::cout << "File 1: " << line1 << "\n";
            std::cout << "File 2: " << line2 << "\n";
            allEqual = false;
        }

        lineNum++;
    }

    // Check if one file has more lines than the other
    if ((file1.peek() != EOF) || (file2.peek() != EOF)) {
        std::cout << "Files have a different number of lines." << std::endl;
        return false;
    }

    return allEqual;
}

int main(int argc, char* argv[]) {
    std::string filePath1, filePath2;
    int n;

    //std::cout << "Enter path of the first file: ";
    //std::cin >> filePath1;
//
    //std::cout << "Enter path of the second file: ";
    //std::cin >> filePath2;
//
    //std::cout << "Enter the accuracy value n: ";
    //std::cin >> n;

    filePath1 = argv[1];
    filePath2 = argv[2];
    n = atoi(argv[3]);

    if (filesAreEqual(filePath1, filePath2, n)) {
        std::cout << "The files are the same with accuracy " << n << "." << std::endl;
    } else {
        std::cout << "The files are different." << std::endl;
    }

    return 0;
}