#ifndef SPRAY_H
#define SPRAY_H

#include "DataFile.h"
#include "Function.h"
#include <vector>

class Spray {
public:
    Spray(DataFile* df);
    ~Spray();

    void Initialize();
    void Update();
    void Display();
    void Save(std::string n_drop);

private:
    DataFile* _df;
    int _N; // Number of droplets

    // Particle properties
    std::vector<double> _t;
    std::vector<double> _x_p;
    std::vector<double> _v_p;
    std::vector<double> _r_p;
    std::vector<double> _m_p;
    std::vector<double> _m_s;
    std::vector<double> _T_p;

    // Mean values
    double _t_m;
    double _x_p_m;
    double _v_p_m;
    double _r_p_m;
    double _m_p_m;
    double _T_p_m;
};

#endif // SPRAY_H

