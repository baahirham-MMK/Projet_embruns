#ifndef _SPRAY_CPP
#include "Spray.h"


Spray::Spray(DataFile* df, Function* fct) :
_df(df), _fct(fct), _Ms(0.0), _QQ(0.0)
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
    this->_M0 = Mtot;
    int N;

    this->_Abstemp10 = _df->Get_T_air_celcius()+273.15;
    this->_qs10 = 0.62197 * 611.2 / _df->Get_p_0() * exp((17.67 * (this->_Abstemp10 - 273.15)) / (this->_Abstemp10 - 29.66));
    this->_rs10 = _qs10 / (1 - _qs10);
    this->_r10 = (_df->Get_Q_RH()/100.0)*this->_rs10;
    this->_q10 = this->_r10 / (1 + this->_r10);
    this->_q = this->_q10;
    this->_rhod = 1.2929*273.13/this->_Abstemp10;

    this->_HUMIDITY = 0.0;



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
    this->_Ms = 0.0;
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

    // Masse d'eau evaporée
    this->_Ms = this->_M0 - N*this->_m_p_m;
    
    
    // humidité supplementaire 
    this->_HUMIDITY = (1.0 - this->_q)/this->_rhod * this->_Ms;
    this->_QQ = (this->_q + this->_HUMIDITY)/(1.0-(this->_q +this->_HUMIDITY))/this->_rs10;

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
        monflux << this->_t_m << " " << this->_x_p_m << " " << this->_v_p_m << " " << this->_r_p_m*1e6 << " " << this->_m_p_m*1e9 << " " << this->_T_p_m - 273.15 <<" "<<this->_Ms<<" "<< this->_QQ<<std::endl;
        monflux.close();
    } else {
        std::cerr << "Erreur : impossible d'ouvrir le fichier " << n_file << std::endl;
    }
}


#define _SPRAY_CPP
#endif