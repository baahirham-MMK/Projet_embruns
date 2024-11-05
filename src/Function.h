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
        double M(const double r_p, const double v_p, const double m_p, const double T_p) const;
        double R(const double r_p, const double v_p, const double m_p, const double T_p) const;
        double T(const double r_p, const double v_p, const double m_p, const double T_p) const;
        double b(const double r_p, const double v_p, const double m_p, const double T_p) const;
        double tau_t(const double r_p, const double v_p, const double m_p, const double T_p) const;
};

#define _FUNCTION_H
#endif