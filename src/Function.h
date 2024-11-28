#ifndef _FUNCTION_H

#include "DataFile.h"

class Function {
    private:
   
    const DataFile* _df;

    public: // Méthodes et opérateurs de la classe
        Function(DataFile* df);
        double Re_p(const double r_p, const double v_p) const;
        double tau_p(const double r_p, const double m_p) const; //Masse volumique constante ???
        double F(const double r_p, const double v_p, const double m_p) const;
        double f_v(const double r_p, const double v_p) const;
        double f_h(const double r_p, const double v_p) const;
        double D_v_e(const double r_p, const double v_p) const;
        double k_a_e(const double r_p, const double v_p) const;
        double M(const double r_p, const double v_p, const double m_p, const double m_s, const double T_p) const;
        double R(const double r_p, const double v_p, const double m_p, const double m_s, const double T_p) const;
        double T(const double r_p, const double v_p, const double m_p, const double m_s, const double T_p) const;
        double b(const double r_p, const double v_p, const double m_p, const double m_s, const double T_p) const;
        double tau_t(const double r_p, const double v_p, const double m_p, const double T_p) const;
        double dFdr(const double R) const;
        double vp4(const double R) const;
        double Vdp(const double Vt_p) const;
        double dCd_r(const double dFdr, const double Vdp) const;
        double normalised_N_r(const double dCdr) const;
        double acceptation_rejet(std::default_random_engine& seed) const;
};

void charge(int Me, int N, int Np, int &iBeg, int &iEnd);

#define _FUNCTION_H
#endif