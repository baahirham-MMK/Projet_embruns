#ifndef _SPRAY_CPP
#include "Spray.h"


Spray::Spray(DataFile* df, Function* fct) :
_df(df), _fct(fct)
{

} 

Spray::~Spray(){
    int N = this->_spray.size();
    for (int i = 0; i < N; ++i){
        delete this->_spray[i];
    }
}

void Spray::Initialize()
{

    this->_t_m = 0.0;
    this->_x_p_m = 0.0;
    this->_v_p_m = 0.0;
    this->_r_p_m = 0.0;
    this->_m_p_m = 0.0;
    this->_T_p_m = 0.0;

    double Mtot = _fct->rho_0()*pow(_df->Get_L(),3);
    int N;

    while(_m_p_m < Mtot) 
    {
        this->_spray.push_back(new Drop(this->_df, this->_fct));
        N = this->_spray.size();
        this->_spray[N-1]->Initialize();
        this->_t_m += this->_spray[N-1]->Get_t();
        this->_x_p_m += this->_spray[N-1]->Get_x_p();
        this->_v_p_m += this->_spray[N-1]->Get_v_p();
        this->_r_p_m += this->_spray[N-1]->Get_r_p();
        this->_m_p_m += this->_spray[N-1]->Get_m_p();
        this->_T_p_m += this->_spray[N-1]->Get_T_p();
    }
    printf("Nombre de goutte : %d\n",N);
    this->_t_m *= (1./double(N));
    this->_x_p_m *= (1./double(N));
    this->_v_p_m *= (1./double(N));
    this->_r_p_m *= (1./double(N));
    this->_m_p_m *= (1./double(N));
    this->_T_p_m *= (1./double(N));
}

void Spray::Update()
{

    this->_t_m = 0.0;
    this->_x_p_m = 0.0;
    this->_v_p_m = 0.0;
    this->_r_p_m = 0.0;
    this->_m_p_m = 0.0;
    this->_T_p_m = 0.0;
    int N = this->_spray.size();
    
    for (int i = 0; i < N; ++i) 
    {
        this->_spray[i]->Update();
        this->_t_m += this->_spray[i]->Get_t();
        this->_x_p_m += this->_spray[i]->Get_x_p();
        this->_v_p_m += this->_spray[i]->Get_v_p();
        this->_r_p_m += this->_spray[i]->Get_r_p();
        this->_m_p_m += this->_spray[i]->Get_m_p();
        this->_T_p_m += this->_spray[i]->Get_T_p();
    }
    this->_t_m *= (1./double(N));
    this->_x_p_m *= (1./double(N));
    this->_v_p_m *= (1./double(N));
    this->_r_p_m *= (1./double(N));
    this->_m_p_m *= (1./double(N));
    this->_T_p_m *= (1./double(N));
}

void Spray::Display()
{
    std::cout << "t = " << this->_t_m << " [s] " << std::endl;
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
        monflux << this->_t_m << " " << this->_x_p_m << " " << this->_v_p_m << " " << this->_r_p_m*1e6 << " " << this->_m_p_m*1e9 << " " << this->_T_p_m - 273.15 << std::endl;
        monflux.close();
    } else {
        std::cerr << "Erreur : impossible d'ouvrir le fichier " << n_file << std::endl;
    }
}


#define _SPRAY_CPP
#endif