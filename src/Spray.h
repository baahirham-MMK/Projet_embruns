#ifndef _SPRAY_H

#include "Drop.h"
// Définition de la classe

class Spray {

    private:
        DataFile* _df;
        Function* _fct;
        std::vector<Drop*> _spray;
        double _t_m, _x_p_m, _v_p_m, _m_p_m, _r_p_m, _T_p_m;
        double _M0, _Ms;

        double _Abstemp10;      // Température absolue à 10 m (K)
        double _qs10;           // Humidité spécifique de saturation à 10 m
        double _rs10;           // Rapport de mélange de saturation à 10 m
        double _r10;            // Rapport de mélange réel à 10 m
        double _q10;            // Humidité spécifique réelle à 10 m
        double _q;              // Humidité spécifique réelle (initiale)
        double _rhod;
        double _HUMIDITY;       // Humidité supplémentaire due à l'évaporation
        double _QQ;             // Nouvelle humidité relative après évaporation

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