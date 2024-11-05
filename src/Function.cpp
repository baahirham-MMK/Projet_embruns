#ifndef _FUNCTION_CPP

#include "Function.h"
#include <cmath>

Function::Function(DataFile* df) :
_df(df)
{
   
}

double Function::Re_p(const double r_p, const double v_p) const
{
    double v_s = v_p - _df->Get_U_air();
    return (2.0*r_p*std::fabs(v_s))/_df->Get_nu_air();
}

double Function::tau_p(const double r_p) const 
{
    return (2.0*std::pow(r_p,2)*_df->Get_rho_p())/(9.0*_df->Get_mu_air());
}

double Function::F(const double r_p, const double v_p, const double m_p) const
{
    double tau_D = tau_p(r_p);
    double f = (m_p/tau_D)*(_df->Get_U_air() - v_p);
    return m_p*_df->Get_g() + f;
}

double Function::f_v(const double r_p, const double v_p) const
{
    return 1.0 + std::sqrt(Re_p(r_p,v_p))/4.0;
}

double Function::f_h(const double r_p, const double v_p) const
{
    return f_v(r_p,v_p);
}

double Function::D_v_e(const double r_p, const double v_p) const
{
    double pi = 3.14159265358979323846264338327950288;
    double a = f_v(r_p,v_p)*_df->Get_D_v();
    double b = r_p/(r_p + _df->Get_Delta_v());
    double c = (_df->Get_D_v()/(r_p*_df->Get_alpha_c()))*std::sqrt((2.0*pi*_df->Get_M_w())/(_df->Get_R_g()*_df->Get_T_air()));
    return a/(b + c);
}

double Function::k_a_e(const double r_p, const double v_p) const
{
    double pi = 3.14159265358979323846264338327950288;
    double a = f_h(r_p,v_p)*_df->Get_k_a();
    double b = r_p/(r_p + _df->Get_Delta_T());
    double c = (_df->Get_k_a()/(r_p*_df->Get_alpha_T()*_df->Get_rho_air_sec()*_df->Get_c_p_air_sec()))*std::sqrt((2.0*pi*_df->Get_M_a())/(_df->Get_R_g()*_df->Get_T_air()));
    return a/(b + c);
}

double Function::M(const double r_p, const double v_p, const double m_p, const double T_p) const
{
    double pi = 3.14159265358979323846264338327950288;
    double a = (4.0*pi*r_p*D_v_e(r_p,v_p)*_df->Get_M_w()*_df->Get_p_v_sat_T_air())/(_df->Get_R_g()*_df->Get_T_air());
    double b = _df->Get_Q_RH() - (_df->Get_T_air()/T_p)*std::exp(((_df->Get_L_v()*_df->Get_M_w())/_df->Get_R_g())*((1.0/_df->Get_T_air()) - (1.0/T_p)) + (2.0*_df->Get_M_w()*_df->Get_Gamma_p())/(_df->Get_R_g()*_df->Get_rho_w()*r_p*T_p) - (_df->Get_I()*_df->Get_Phi_s()*_df->Get_m_s()*(_df->Get_M_w()/_df->Get_M_s()))/(m_p - _df->Get_m_s()));
    return a*b;
}

double Function::R(const double r_p, const double v_p, const double m_p, const double T_p) const
{
    double pi = 3.14159265358979323846264338327950288;
    return M(r_p,v_p,m_p,T_p)/(4.0*pi*std::pow(r_p,2)*_df->Get_rho_w());
}

double Function::T(const double r_p, const double v_p, const double m_p, const double T_p) const
{
    double pi = 3.14159265358979323846264338327950288;
    double a = (4.0*pi*r_p*k_a_e(r_p,v_p)*(_df->Get_T_air() - T_p))/(m_p*_df->Get_c_p_s());
    return a + (_df->Get_L_v()/(m_p*_df->Get_c_p_s()))*M(r_p,v_p,m_p,T_p);
}

double Function::b(const double r_p, const double v_p, const double m_p, const double T_p) const
{   
    double pi = 3.14159265358979323846264338327950288;
    double a = (_df->Get_T_air()/T_p)*std::exp((_df->Get_L_v()*_df->Get_M_w()/_df->Get_R_g())*((1.0/_df->Get_T_air())-(1.0/T_p)));
    double b = _df->Get_Q_RH() + (2.0*_df->Get_M_w()*_df->Get_Gamma_p())/(_df->Get_R_g()*_df->Get_rho_w()*r_p*T_p) - (_df->Get_I()*_df->Get_Phi_s()*_df->Get_m_s()*(_df->Get_M_w()/_df->Get_M_s()))/(m_p - _df->Get_m_s());
    double c = (4.0*pi*r_p*D_v_e(r_p,v_p)*_df->Get_M_w()*_df->Get_p_v_sat_T_air())/(_df->Get_R_g()*_df->Get_T_air());
    return _df->Get_T_air() + (_df->Get_L_v()*c/(4.0*pi*r_p*k_a_e(r_p,v_p)))*(b-a);
}

double Function::tau_t(const double r_p, const double v_p, const double m_p, const double T_p) const
{
    double pi = 3.14159265358979323846264338327950288;
    return m_p*_df->Get_c_p_s()/(4*pi*r_p*k_a_e(r_p,v_p));
}

#define _FUNCTION_CPP
#endif
