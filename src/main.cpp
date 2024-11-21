#include"Spray.h"

int main(int argc, char** argv){

    if (argc < 2)
    {
        printf("Please, enter the name of your data file.\n");
        exit(0);
    }

    const std::string data_file_name = argv[1];

    DataFile* df = new DataFile(data_file_name);
    Function* fct = new Function(df);
    Spray* spray = new Spray(df,fct);

    spray->Initialize();
    spray->Display();
    spray->Save("spray");
    for (int j = 0; j < 2e7; ++j){
        spray->Update();
        spray->Display();
        spray->Save("spray");
    }

    delete df, delete fct, delete spray;

    system("gnuplot ../res/visu.gnu");

    return 0;
}