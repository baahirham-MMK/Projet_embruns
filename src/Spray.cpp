#ifndef _SPRAY_CPP
#include "Spray.h"


Spray::Spray(DataFile* df, Function* fct) :
_df(df), _fct(fct)
{

} 

Spray::~Spray(){
    for (int i = 0; i < _df->Get_N(); ++i){
        delete this->_spray[i];
    }
}

void Spray::Initialize()
{

    this->_t = 0.0;
    this->_x_p_m = 0.0;
    this->_v_p_m = 0.0;
    this->_r_p_m = 0.0;
    this->_m_p_m = 0.0;
    this->_T_p_m = 0.0;

    this->_spray.resize(_df->Get_N());

    for (int i = 0; i < _df->Get_N(); ++i) {
        this->_spray[i] = new Drop(this->_df, this->_fct);
        this->_spray[i]->Initialize();
        this->_t += this->_spray[i]->Get_t();
        this->_x_p_m += this->_spray[i]->Get_x_p();
        this->_v_p_m += this->_spray[i]->Get_v_p();
        this->_r_p_m += this->_spray[i]->Get_r_p();
        this->_m_p_m += this->_spray[i]->Get_m_p();
        this->_T_p_m += this->_spray[i]->Get_T_p();
    }
    this->_t *= (1./double(_df->Get_N()));
    this->_x_p_m *= (1./double(_df->Get_N()));
    this->_v_p_m *= (1./double(_df->Get_N()));
    this->_r_p_m *= (1./double(_df->Get_N()));
    this->_m_p_m *= (1./double(_df->Get_N()));
    this->_T_p_m *= (1./double(_df->Get_N()));
}

void Spray::Update()
{

    this->_t = 0.0;
    this->_x_p_m = 0.0;
    this->_v_p_m = 0.0;
    this->_r_p_m = 0.0;
    this->_m_p_m = 0.0;
    this->_T_p_m = 0.0;
    
    for (int i = 0; i < _df->Get_N(); ++i) {
        this->_spray[i]->Update();
        this->_t += this->_spray[i]->Get_t();
        this->_x_p_m += this->_spray[i]->Get_x_p();
        this->_v_p_m += this->_spray[i]->Get_v_p();
        this->_r_p_m += this->_spray[i]->Get_r_p();
        this->_m_p_m += this->_spray[i]->Get_m_p();
        this->_T_p_m += this->_spray[i]->Get_T_p();
    }
    this->_t *= (1./_df->Get_N());
    this->_x_p_m *= (1./_df->Get_N());
    this->_v_p_m *= (1./_df->Get_N());
    this->_r_p_m *= (1./_df->Get_N());
    this->_m_p_m *= (1./_df->Get_N());
    this->_T_p_m *= (1./_df->Get_N());
}

void Spray::Display()
{
    std::cout << "t = " << this->_t << " [s] " << std::endl;
    std::cout << "x_p = " << this->_x_p_m << " [m] " << std::endl;
    std::cout << "v_p = " << this->_v_p_m << " [m/s] " << std::endl;
    std::cout << "r_p = " << this->_r_p_m << " [m] " << std::endl;
    std::cout << "m_p = " << this->_m_p_m << " [kg] " << std::endl;
    std::cout << "T_p = " << this->_T_p_m << " [K] " << std::endl;
}

void Spray::Save(std::string n_drop)
{
    std::string n_file = "../res/" + n_drop + ".dat";
    
    std::ofstream monflux;
    monflux.open(n_file, std::ios::app);  
    if (monflux.is_open()) {
        monflux << this->_t << " " << this->_x_p_m << " " << this->_v_p_m << " " << this->_r_p_m*1e6 << " " << this->_m_p_m*1e9 << " " << this->_T_p_m - 273.15 << std::endl;
        monflux.close();
    } else {
        std::cerr << "Erreur : impossible d'ouvrir le fichier " << n_file << std::endl;
    }
}


#define _SPRAY_CPP
#endif