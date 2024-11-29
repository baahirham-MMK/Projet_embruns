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

    int i = 0;
    while(spray->Get_t_m() < df->Get_T_f())
    {
        spray->Update();
        //spray->Display();

        spray->Save("spray");
        i += 1;
        //printf("iteration=%d\n",i);
    }

    delete df, delete fct, delete spray;

    system("gnuplot ../res/visu.gnu");

    return 0;
}