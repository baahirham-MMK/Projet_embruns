#ifndef _DROP_H

#include "Function.h"
// Définition de la classe

class Drop {

    private:
        const DataFile* _df;
        const Function* _fct;
        double _t, _x_p, _v_p, _m_p, _r_p, _T_p;

    public: // Méthodes et opérateurs de la classe
        Drop(DataFile* df, Function* fct);
        void Initialize();
        void Update(); 
        void Display();
        void Save(std::string n_sol);
        const double Get_t() const {return _t;};
        const double Get_x_p() const {return _x_p;};
        const double Get_v_p() const {return _v_p;};
        const double Get_m_p() const {return _m_p;};
        const double Get_r_p() const {return _r_p;};
        const double Get_T_p() const {return _T_p;};
};

#define _DROP_H
#endif