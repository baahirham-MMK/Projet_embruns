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

double Function::tau_p(const double r_p, const double m_p) const 
{
    return (2.0*std::pow(r_p,2)*(m_p/(4./3.*std::acos(-1.0)*std::pow(r_p,3))))/(9.0*_df->Get_mu_air());
}

double Function::F(const double r_p, const double v_p, const double m_p) const
{
    double tau_D = tau_p(r_p,m_p);
    double f = (m_p/tau_D)*(_df->Get_U_air() - v_p);
    return f;
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
    double pi = std::acos(-1.0);
    double a = f_v(r_p,v_p)*_df->Get_D_v();
    double b = r_p/(r_p + _df->Get_Delta_v());
    double c = (_df->Get_D_v()/(r_p*_df->Get_alpha_c()))*std::sqrt((2.0*pi*_df->Get_M_w())/(_df->Get_R_g()*_df->Get_T_air()));
    return a/(b + c);
}

double Function::k_a_e(const double r_p, const double v_p) const
{
    double pi = std::acos(-1.0);
    double a = f_h(r_p,v_p)*_df->Get_k_a();
    double b = r_p/(r_p + _df->Get_Delta_T());
    double c = (_df->Get_k_a()/(r_p*_df->Get_alpha_T()*_df->Get_rho_air_sec()*_df->Get_c_p_air_sec()))*std::sqrt((2.0*pi*_df->Get_M_a())/(_df->Get_R_g()*_df->Get_T_air()));
    return a/(b + c);
}

double Function::M(const double r_p, const double v_p, const double m_p, const double m_s, const double T_p) const
{
    double pi = std::acos(-1.0);
    double a = (4.0*pi*r_p*D_v_e(r_p,v_p)*_df->Get_M_w()*_df->Get_p_v_sat_T_air())/(_df->Get_R_g()*_df->Get_T_air());
    double b = _df->Get_Q_RH() - (_df->Get_T_air()/T_p)*std::exp(((_df->Get_L_v()*_df->Get_M_w())/_df->Get_R_g())*((1.0/_df->Get_T_air()) - (1.0/T_p)) + (2.0*_df->Get_M_w()*_df->Get_Gamma_p())/(_df->Get_R_g()*_df->Get_rho_w()*r_p*T_p) - (_df->Get_I()*_df->Get_Phi_s()*m_s*(_df->Get_M_w()/_df->Get_M_s()))/(m_p - m_s));
    return a*b;
}

double Function::R(const double r_p, const double v_p, const double m_p, const double m_s, const double T_p) const
{
    double pi = std::acos(-1.0);
    return M(r_p,v_p,m_p,m_s,T_p)/(4.0*pi*std::pow(r_p,2)*_df->Get_rho_w());
}

double Function::T(const double r_p, const double v_p, const double m_p, const double m_s, const double T_p) const
{
    double pi = std::acos(-1.0);
    double a = (4.0*pi*r_p*k_a_e(r_p,v_p)*(_df->Get_T_air() - T_p))/(m_p*_df->Get_c_p_s());
    return a + (_df->Get_L_v()/(m_p*_df->Get_c_p_s()))*M(r_p,v_p,m_p,m_s,T_p);
}

double Function::b(const double r_p, const double v_p, const double m_p, const double m_s, const double T_p) const
{   
    double pi =  std::acos(-1.0);
    return _df->Get_T_air() + (_df->Get_L_v()*M(r_p,v_p,m_p,T_p)/(4*pi*r_p*k_a_e(r_p,v_p)));
}

double Function::tau_t(const double r_p, const double v_p, const double m_p, const double T_p) const
{
    double pi =  std::acos(-1.0);
    return m_p*_df->Get_c_p_s()/(4*pi*r_p*k_a_e(r_p,v_p));
}

double Function::dFdr(const double R) const {

    double dFdr_F, r_80;
    double Uwind10(15.0);

    r_80 = 0.518 * std::pow(R * 1e6, 0.976);

    if (r_80 > 0.9 && r_80 < 15) {
        dFdr_F = 75.7672 / 2.0 * 3.8e-6 * std::pow(Uwind10, 3.4) * std::pow(r_80, -0.024) * 
                        std::pow(10, 4.405 - 2.646 * std::log10(r_80) - 3.156 * std::pow(std::log10(r_80), 2) +
                                    8.902 * std::pow(std::log10(r_80), 3) - 4.482 * std::pow(std::log10(r_80), 4));
    } else if (r_80 >= 15 && r_80 <= 37.5) {
        dFdr_F = 75.7672 / 2.0 * 3.8e-6 * std::pow(Uwind10, 3.4) * std::pow(r_80, -0.024) * 
                        1.02e4 * std::pow(r_80, -1);
    } else if (r_80 >= 37.5 && r_80 <= 100) {
        dFdr_F = 75.7672 / 2.0 * 3.8e-6 * std::pow(Uwind10, 3.4) * std::pow(r_80, -0.024) * 
                        6.95e6 * std::pow(r_80, -2.8);
    } else if (r_80 > 100) {
        dFdr_F = 75.7672 / 2.0 * 3.8e-6 * std::pow(Uwind10, 3.4) * std::pow(r_80, -0.024) * 
                        1.75e17 * std::pow(r_80, -8);
    } else {
        dFdr_F = NAN;
    }

    if (r_80 > 0.5 && r_80 < 0.9) {
        double a = 34623;
        dFdr_F = a * std::pow(r_80, -3);
    }

    dFdr_F *= 1e6;

    return dFdr_F;
}

double Function::vp4(const double R) const {
    const double g = 9.81;
    const double lambda = 68e-9;
    const double rhop = 1024.0;
    const double rhoa = 1.19448;
    const double nua = 1.51905e-05;
    double Vt_p, Cd_p, RE, CU;

    CU = 1 + lambda / R * (1.257 + 0.4 * std::exp(-0.55 * 2 * R / lambda));

    double vp2;
    vp2 = 0.01 * g * std::pow(R * 2, 2) * (rhop - rhoa) / (18 * nua * rhoa);

    for (int iter = 0; iter < 200; ++iter) {
        
        RE = 2 * R * vp2 / nua;

        if (RE < 500) {
            Cd_p = 24.0 / RE * (1 + 0.15 * std::pow(RE, 0.687) + 
                        0.175 * std::pow(1 + 4.25e4 * std::pow(RE, -1.16), -1));
        } else {
            Cd_p = 0.85 - 9.76e-4 * RE + 1.091e-6 * std::pow(RE, 2) - 
                        6.84e-10 * std::pow(RE, 3) + 2.72e-13 * std::pow(RE, 4) -
                        6.68e-17 * std::pow(RE, 5) + 9.83e-21 * std::pow(RE, 6) -
                        7.96e-25 * std::pow(RE, 7) + 2.73e-29 * std::pow(RE, 8);
        }

        Cd_p /= CU;

        double f = Cd_p * RE / 24.0;
        vp2 = 0.5 * (g * std::pow(R * 2, 2) * (rhop - rhoa) / (18 * nua * rhoa * f)) + 0.5 * vp2;

        if (2 * R > 1e-2) {
            vp2 = 9.36;
        }
    }

    Vt_p = vp2;

    return Vt_p;
}

double Function::Vdp(const double Vt_p) const {

    double Vd_p;
    double Uwind10(15.0);

    Vd_p = Vt_p/(1-std::exp(-Vt_p/1e-3/Uwind10));

    return Vd_p;
}

double Function::dCd_r(const double dFdr, const double Vdp) const {
    
    double dCdr;

    dCdr = dFdr/Vdp;

    return dCdr;
}

double Function::normalised_N_r(const double dCdr) const {
    double K(8.55422e+08);
    return dCdr/K;
}

double Function::acceptation_rejet(std::default_random_engine& seed) const {

    double max(19412.8);
    double a(1e-6), b(1e-3), c((b-a)*max), x, y;

    std::uniform_real_distribution<double> u1{0, 1}, u2{a, b};

    x = u1(seed);
    y = u2(seed);

    while (x > (b-a)*normalised_N_r(dCd_r(dFdr(y),Vdp(y)))/c){
        x = u1(seed);
        y = u2(seed);
    }

    return y;
}

#define _FUNCTION_CPP
#endif