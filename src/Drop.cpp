#ifndef _DROP_CPP
#include "Drop.h"


Drop::Drop(DataFile* df, Function* fct) :
_df(df), _fct(fct)
{
   
} 

void Drop::Initialize()
{
    std::random_device rd;
    std :: default_random_engine seed(rd());
    std::uniform_real_distribution<double> u1{0, _df->Get_L()}, u2{1e-4, 1e-3};
    this->_t = 0.0;
    this->_x_p = u1(seed);
    this->_v_p = 0.0;
    this->_r_p = _fct->acceptation_rejet(seed);
    this->_m_p = (4.0/3.0)*std::acos(-1.0)*std::pow(this->_r_p,3)*_df->Get_rho_p();
    this->_m_s = (4.0/3.0)*std::acos(-1.0)*std::pow(this->_r_p,3)*_df->Get_rho_p()*_df->Get_Salinity_p()/1000.0;
    this->_T_p = _df->Get_T_p_0();
}

void Drop::Update()
{
    double x_p_old = this->_x_p;
    double v_p_old = this->_v_p;
    double r_p_old = this->_r_p;
    double m_p_old = this->_m_p;
    double T_p_old = this->_T_p;

    //double dt = _fct->tau_p(this->_r_p,this->_m_p);
    double dt;
    if(this->_t < 2e-4)
    {
        dt = 1e-5;
    }
    else if(this->_t < 2e-3)
    {
        dt = 1e-4;
    }
    else
    {
        dt = 2e-3;
    }

    switch (_df->Get_cas())
    {
    case 1:
        if (this->_x_p <= _df->Get_L()) this->_x_p += dt * v_p_old;
        else this->_x_p = 0.0;
        this->_v_p += (dt/m_p_old) * _fct->F(r_p_old,v_p_old,m_p_old);
        this->_r_p += dt * _fct->R(r_p_old,v_p_old,m_p_old,this->_m_s,T_p_old);
        this->_m_p += dt * _fct->M(r_p_old,v_p_old,m_p_old, this->_m_s,T_p_old);
        this->_T_p += dt * _fct->T(r_p_old,v_p_old,m_p_old,this->_m_s,T_p_old);

        break;
    default:
        if (this->_x_p <= _df->Get_L()) this->_x_p += dt * v_p_old;
        else this->_x_p = 0.0;
        this->_v_p = this->_v_p * exp(-dt/_fct->tau_p(this->_r_p, this->_m_p)) + _df->Get_U_air() 
                                                *(1 - exp(-dt/_fct->tau_p(this->_r_p, this->_m_p)) );    
        this->_r_p += dt * _fct->R(r_p_old,v_p_old,m_p_old,this->_m_s,T_p_old);
        this->_m_p += dt * _fct->M(r_p_old,v_p_old,m_p_old,this->_m_s,T_p_old);
        this->_T_p = this->_T_p*exp(-dt/_fct->tau_t(this->_r_p, this->_v_p, this->_m_p, this->_T_p))
                                            +_fct->b(this->_r_p, this-> _v_p, this->_m_p, this->_m_s, this->_T_p)
                            *(1 - exp(-dt/_fct->tau_t(this->_r_p, this->_v_p, this->_m_p, this->_T_p)));


        break;
    }

    if(this->_x_p > _df->Get_L())
    {
       this->_x_p -= _df->Get_L();
    }
    
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
        monflux << this->_t << " " << this->_x_p << " " << this->_v_p << " " << this->_r_p*1e6 << " " << this->_m_p*1e9 << " " << this->_T_p - 273.15 << std::endl;
        monflux.close();
    } else {
        std::cerr << "Erreur : impossible d'ouvrir le fichier " << n_file << std::endl;
    }
}


#define _DATA_FILE_CPP
#endif