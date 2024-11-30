#include "Spray.h"
#include <chrono>
#include <iostream>

int main(int argc, char** argv) {

    if (argc < 2) {
        printf("Please, enter the name of your data file.\n");
        exit(0);
    }

    const std::string data_file_name = argv[1];

    DataFile* df = new DataFile(data_file_name);
    Spray* spray = new Spray(df);

    auto start = std::chrono::high_resolution_clock::now();
    spray->Initialize();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> init_duration = end - start;
    std::cout << "Initialization time: " << init_duration.count() << " seconds" << std::endl;

    spray->Save("spray");

    start = std::chrono::high_resolution_clock::now();
    int progress_step = 5e6;
    for (int j = 0; j < 5e7; ++j) {
        spray->Update();
        if (j % 10 == 0) {
            spray->Save("spray");
        }
        if (j % progress_step == 0) {
            std::cout << "Progress: " << (j / 5e7 * 100) << "% complete (" << j << "/" << 5e7 << " iterations)" << std::endl;
        }
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> update_duration = end - start;
    std::cout << "Update time for 5e7 iterations: " << update_duration.count() << " seconds" << std::endl;

    delete df;
    delete spray;

    system("gnuplot ../res/visu.gnu");

    return 0;
}

