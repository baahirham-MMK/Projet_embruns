#include "Function.h"
#include <cmath>
#include <limits> 
#include <cstdio> 

namespace Function {

double Re_p(const DataFile& df, const double r_p, const double v_p)
{
    double v_s = v_p - df.Get_U_air();
    return (2.0 * r_p * std::fabs(v_s)) / df.Get_nu_air();
}

double tau_p(const DataFile& df, const double r_p, const double m_p)
{
    double rho_p = m_p / ((4.0 / 3.0) * M_PI * std::pow(r_p, 3));
    return (2.0 * std::pow(r_p, 2) * rho_p) / (9.0 * df.Get_mu_air());
}

double F(const DataFile& df, const double r_p, const double v_p, const double m_p)
{
    double tau_D = tau_p(df, r_p, m_p);
    double f = (m_p / tau_D) * (df.Get_U_air() - v_p);
    return f;
}

double f_v(const DataFile& df, const double r_p, const double v_p)
{
    return 1.0 + std::sqrt(Re_p(df, r_p, v_p)) / 4.0;
}

double f_h(const DataFile& df, const double r_p, const double v_p)
{
    return f_v(df, r_p, v_p);
}

double D_v_e(const DataFile& df, const double r_p, const double v_p)
{
    double pi = M_PI;
    double a = f_v(df, r_p, v_p) * df.Get_D_v();
    double b = r_p / (r_p + df.Get_Delta_v());
    double c = (df.Get_D_v() / (r_p * df.Get_alpha_c())) * std::sqrt((2.0 * pi * df.Get_M_w()) / (df.Get_R_g() * df.Get_T_air()));
    return a / (b + c);
}

double k_a_e(const DataFile& df, const double r_p, const double v_p)
{
    double pi = M_PI;
    double a = f_h(df, r_p, v_p) * df.Get_k_a();
    double b = r_p / (r_p + df.Get_Delta_T());
    double c = (df.Get_k_a() / (r_p * df.Get_alpha_T() * df.Get_rho_air_sec() * df.Get_c_p_air_sec())) * std::sqrt((2.0 * pi * df.Get_M_a()) / (df.Get_R_g() * df.Get_T_air()));
    return a / (b + c);
}

double M(const DataFile& df, const double r_p, const double v_p, const double m_p, const double m_s, const double T_p)
{
    double pi = M_PI;
    double a = (4.0 * pi * r_p * D_v_e(df, r_p, v_p) * df.Get_M_w() * df.Get_p_v_sat_T_air()) / (df.Get_R_g() * df.Get_T_air());
    double exponent = ((df.Get_L_v() * df.Get_M_w()) / df.Get_R_g()) * ((1.0 / df.Get_T_air()) - (1.0 / T_p))
                      + (2.0 * df.Get_M_w() * df.Get_Gamma_p()) / (df.Get_R_g() * df.Get_rho_w() * r_p * T_p)
                      - (df.Get_I() * df.Get_Phi_s() * m_s * (df.Get_M_w() / df.Get_M_s())) / (m_p - m_s);
    double b = df.Get_Q_RH() - (df.Get_T_air() / T_p) * std::exp(exponent);
    return a * b;
}

double R(const DataFile& df, const double r_p, const double v_p, const double m_p, const double m_s, const double T_p)
{
    double pi = M_PI;
    return M(df, r_p, v_p, m_p, m_s, T_p) / (4.0 * pi * std::pow(r_p, 2) * df.Get_rho_w());
}

double T(const DataFile& df, const double r_p, const double v_p, const double m_p, const double m_s, const double T_p)
{
    double pi = M_PI;
    double a = (4.0 * pi * r_p * k_a_e(df, r_p, v_p) * (df.Get_T_air() - T_p)) / (m_p * df.Get_c_p_s());
    return a + (df.Get_L_v() / (m_p * df.Get_c_p_s())) * M(df, r_p, v_p, m_p, m_s, T_p);
}

double b(const DataFile& df, const double r_p, const double v_p, const double m_p, const double m_s, const double T_p)
{
    double pi = M_PI;
    return df.Get_T_air() + (df.Get_L_v() * M(df, r_p, v_p, m_p, m_s, T_p)) / (4.0 * pi * r_p * k_a_e(df, r_p, v_p));
}

double tau_t(const DataFile& df, const double r_p, const double v_p, const double m_p, const double T_p)
{
    double pi = M_PI;
    return (m_p * df.Get_c_p_s()) / (4.0 * pi * r_p * k_a_e(df, r_p, v_p));
}

double dFdr(const DataFile& df, const double R)
{
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
    } else if (r_80 > 37.5 && r_80 <= 100) {
        dFdr_F = 75.7672 / 2.0 * 3.8e-6 * std::pow(Uwind10, 3.4) * std::pow(r_80, -0.024) *
                 6.95e6 * std::pow(r_80, -2.8);
    } else if (r_80 > 100) {
        dFdr_F = 75.7672 / 2.0 * 3.8e-6 * std::pow(Uwind10, 3.4) * std::pow(r_80, -0.024) *
                 1.75e17 * std::pow(r_80, -8);
    } else if (r_80 > 0.5 && r_80 <= 0.9) {
        double a = 34623;
        dFdr_F = a * std::pow(r_80, -3);
    } else {
        dFdr_F = std::numeric_limits<double>::quiet_NaN();
    }

    dFdr_F *= 1e6;

    return dFdr_F;
}

double vp4(const DataFile& df, const double R)
{
    const double g = df.Get_g();
    const double lambda = 68e-9;
    const double rhop = df.Get_rho_p();
    const double rhoa = df.Get_rho_air();
    const double nua = df.Get_nu_air();
    double Vt_p, Cd_p, RE, CU;

    CU = 1 + lambda / R * (1.257 + 0.4 * std::exp(-0.55 * 2 * R / lambda));

    double vp2;
    vp2 = 0.01 * g * std::pow(2 * R, 2) * (rhop - rhoa) / (18 * nua * rhoa);

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
        vp2 = 0.5 * (g * std::pow(2 * R, 2) * (rhop - rhoa) / (18 * nua * rhoa * f)) + 0.5 * vp2;

        if (2 * R > 1e-2) {
            vp2 = 9.36;
        }
    }

    Vt_p = vp2;

    return Vt_p;
}

double Vdp(const DataFile& df, const double Vt_p)
{
    double Uwind10(15.0);
    double Vd_p;
    Vd_p = Vt_p / (1 - std::exp(-Vt_p / 1e-3 / Uwind10));

    return Vd_p;
}

double dCd_r(const DataFile& df, const double dFdr, const double Vdp)
{
    double dCdr;

    dCdr = dFdr / Vdp;

    return dCdr;
}

double normalised_N_r(const DataFile& df, const double dCdr)
{
    double K(8.55422e+08);
    return dCdr / K;
}

double acceptation_rejet(const DataFile& df, std::default_random_engine& seed)
{
    double max(19412.8);
    double a(1e-6), b(1e-3), c((b - a) * max), x, y;

    std::uniform_real_distribution<double> u1{0, 1}, u2{a, b};

    x = u1(seed);
    y = u2(seed);

    while (x > (b - a) * normalised_N_r(df, dCd_r(df, dFdr(df, y), Vdp(df, vp4(df, y)))) / c) {
        x = u1(seed);
        y = u2(seed);
    }

    return y;
}

double rho_0(const DataFile& df)
{
    double dlogr = 0.1;
    double ri = 1e-6;
    double rho_0 = 0.0;
    int i = 1;

    while (ri <= 1e-3)
    {
        double dri = ri * (std::pow(10, dlogr) - 1);
        double dFdr_val = dFdr(df, ri);
        double Vdp_val = Vdp(df, vp4(df, ri));
        double dCd_r_val = dCd_r(df, dFdr_val, Vdp_val);
        double m_p_i = (4.0 / 3.0) * M_PI * std::pow(ri, 3) * df.Get_rho_p();

        rho_0 += ri * dri * dCd_r_val * m_p_i;

        ri = std::pow(10, -6 + (i * dlogr));
        i += 1;
    }

    // printf("rho_0 = %lf\n", rho_0);

    return rho_0;
}

} // namespace Function

