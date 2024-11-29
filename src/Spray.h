#ifndef _SPRAY_H

#include "Drop.h"
// Définition de la classe

class Spray {

    private:
        DataFile* _df;
        Function* _fct;
        std::vector<Drop*> _spray;
        double _t_m, _x_p_m, _v_p_m, _m_p_m, _r_p_m, _T_p_m;

    public: // Méthodes et opérateurs de la classe
        Spray(DataFile* df, Function* fct);
        ~Spray();
        void Initialize();
        void Update(); 
        void Display();
        void Save(std::string n_sol);
        const double Get_t_m() const {return _t_m;};
};

#define _SPRAY_H
#endif