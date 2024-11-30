#ifndef _FUNCTION_H
#define _FUNCTION_H

#include "DataFile.h"
#include <random>

namespace Function {

    // Les fonctions utilisées dans les régions parallèles OpenACC 
    #pragma acc routine seq
    double Re_p(const DataFile& df, const double r_p, const double v_p);

    #pragma acc routine seq
    double tau_p(const DataFile& df, const double r_p, const double m_p);

    #pragma acc routine seq
    double F(const DataFile& df, const double r_p, const double v_p, const double m_p);

    #pragma acc routine seq
    double f_v(const DataFile& df, const double r_p, const double v_p);

    #pragma acc routine seq
    double f_h(const DataFile& df, const double r_p, const double v_p);

    #pragma acc routine seq
    double D_v_e(const DataFile& df, const double r_p, const double v_p);

    #pragma acc routine seq
    double k_a_e(const DataFile& df, const double r_p, const double v_p);

    #pragma acc routine seq
    double M(const DataFile& df, const double r_p, const double v_p, const double m_p, const double m_s, const double T_p);

    #pragma acc routine seq
    double R(const DataFile& df, const double r_p, const double v_p, const double m_p, const double m_s, const double T_p);

    #pragma acc routine seq
    double T(const DataFile& df, const double r_p, const double v_p, const double m_p, const double m_s, const double T_p);

    #pragma acc routine seq
    double b(const DataFile& df, const double r_p, const double v_p, const double m_p, const double m_s, const double T_p);

    #pragma acc routine seq
    double tau_t(const DataFile& df, const double r_p, const double v_p, const double m_p, const double T_p);

    #pragma acc routine seq
    double dFdr(const DataFile& df, const double R);

    #pragma acc routine seq
    double vp4(const DataFile& df, const double R);

    #pragma acc routine seq
    double Vdp(const DataFile& df, const double Vt_p);

    #pragma acc routine seq
    double dCd_r(const DataFile& df, const double dFdr, const double Vdp);

    #pragma acc routine seq
    double normalised_N_r(const DataFile& df, const double dCdr);

    #pragma acc routine seq
    double rho_0(const DataFile& df); 

    double acceptation_rejet(const DataFile& df, std::default_random_engine& seed);

}

#endif // _FUNCTION_H

