//
// Created by Suman Kalyan Bera on 2020-12-26.
//

#include <iostream>
#include <fstream>      // std::ifstream
#include <string>


int main()
{
    std::ifstream files;
    files.open("../db.txt");
    if (files.is_open())
        std::cout<< "open file" <<std::endl;
    else
        std::cout<< "closed file" <<std::endl;
    std::string line;
    while (getline(files, line)) {
        std::cout << line << std::endl;
    }
    files.close();
}
