#ifndef _DROP_CPP

#include "Drop.h"

Drop::Drop(DataFile* df, Function* fct) :
_df(df), _fct(fct)
{
   
} 

void Drop::Initialize()
{
    this->_t = 0.0;
    this->_x_p = 0.0;
    this->_v_p = 0.0;
    this->_r_p = _df->Get_r_p_0();
    this->_m_p = (4.0/3.0)*std::acos(-1.0)*std::pow(this->_r_p,3)*_df->Get_rho_p();
    this->_T_p = _df->Get_T_p_0();
}

void Drop::Update()
{
    double x_p_old = this->_x_p;
    double v_p_old = this->_v_p;
    double r_p_old = this->_r_p;
    double m_p_old = this->_m_p;
    double T_p_old = this->_T_p;

    double dt = _fct->tau_p(this->_r_p,this->_m_p);

    this->_x_p += dt * v_p_old;
    this->_v_p += (dt/m_p_old) * _fct->F(r_p_old,v_p_old,m_p_old);
    this->_r_p += dt * _fct->R(r_p_old,v_p_old,m_p_old,T_p_old);
    this->_m_p += dt * _fct->M(r_p_old,v_p_old,m_p_old,T_p_old);
    this->_T_p += dt * _fct->T(r_p_old,v_p_old,m_p_old,T_p_old);

    this->_t += dt;
}

void Drop::Display()
{
    std::cout << "t = " << this->_t << " [s] " << std::endl;
    std::cout << "x_p = " << this->_x_p << " [m] " << std::endl;
    std::cout << "v_p = " << this->_v_p << " [m/s] " << std::endl;
    std::cout << "r_p = " << this->_r_p << " [m] " << std::endl;
    std::cout << "m_p = " << this->_m_p << " [kg] " << std::endl;
    std::cout << "T_p = " << this->_T_p << " [K] " << std::endl;
}

void Drop::Save(std::string n_drop)
{
    std::string n_file = "../res/" + n_drop + ".dat";
    
    std::ofstream monflux;
    monflux.open(n_file, std::ios::app);  
    if (monflux.is_open()) {
        monflux << this->_t << " " << this->_x_p << " " << this->_v_p << " " << this->_r_p*1e6 << " " << this->_m_p << " " << this->_T_p - 273.15 << std::endl;
        monflux.close();
    } else {
        std::cerr << "Erreur : impossible d'ouvrir le fichier " << n_file << std::endl;
    }
}


#define _DATA_FILE_CPP
#endif