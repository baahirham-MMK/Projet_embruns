#include"Drop.h"

int main(int argc, char** argv){

    if (argc < 2)
    {
        printf("Please, enter the name of your data file.\n");
        exit(0);
    }

    const std::string data_file_name = argv[1];

    DataFile* df = new DataFile(data_file_name);
    Function* fct = new Function(df);
    Drop* drop = new Drop(df,fct);

    drop->Initialize();
    drop->Display();
    drop->Save("drop_0");

    for (int i = 0; i<5e7; ++i){
        drop->Update();
        if (i%10 == 0){
            drop->Display();
            drop->Save("drop_0");
        }
    }

    delete df, delete fct, delete drop;

    system("gnuplot ../res/visu.gnu");

    return 0;
}