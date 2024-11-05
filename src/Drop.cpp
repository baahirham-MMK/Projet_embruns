#ifndef _DROP_CPP
#include <cmath>
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
    this->_m_p = (4.0/3.0)*3.14159265358979323846264338327950288*std::pow(this->_r_p,3);
    this->_T_p = _df->Get_T_p_0();
}

void Drop::Update(int cas)
{
    double x_p_old = this->_x_p;
    double v_p_old = this->_v_p;
    double r_p_old = this->_r_p;
    double m_p_old = this->_m_p;
    double T_p_old = this->_T_p;

    double dt = 1e-15;

    switch (cas)
    {
    case 1:
        this->_x_p += dt * v_p_old;
        this->_v_p += (dt/m_p_old) * _fct->F(r_p_old,v_p_old,m_p_old);
        this->_r_p += dt * _fct->R(r_p_old,v_p_old,m_p_old,T_p_old);
        this->_m_p += dt * _fct->M(r_p_old,v_p_old,m_p_old,T_p_old);
        this->_T_p += dt * _fct->T(r_p_old,v_p_old,m_p_old,T_p_old);

        break;
    default:
        this->_x_p += dt * v_p_old;
        this->_v_p = this->_v_p * exp(-dt/_fct->tau_p(this->_r_p)) + (_df->Get_U_air() + _df->Get_g()*_fct->tau_p(this->_r_p))
                     *(1 - exp(-dt/_fct->tau_p(this->_r_p)) );    
        this->_r_p += dt * _fct->R(r_p_old,v_p_old,m_p_old,T_p_old);
        this->_m_p += dt * _fct->M(r_p_old,v_p_old,m_p_old,T_p_old);
        //this->_T_p = this->_T_p*exp(-dt/_fct->tau_T())+b(this->_T_p)*(1 - exp(-dt/_fct->tau_T()))
        this->_T_p += dt * _fct->T(r_p_old,v_p_old,m_p_old,T_p_old);


        break;
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

#define _DATA_FILE_CPP
#endif