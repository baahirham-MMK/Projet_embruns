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
        void Update(int cas);  // cas est la varaible de choix, Expo ou EE
        void Display();
};

#define _DROP_H
#endif