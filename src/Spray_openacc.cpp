#include "Spray.h"
#include "Function.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>

Spray::Spray(DataFile* df) :
    _df(df), _N(0)
{
}

Spray::~Spray(){
}

void Spray::Initialize()
{
    this->_t_m = 0.0;
    this->_x_p_m = 0.0;
    this->_v_p_m = 0.0;
    this->_r_p_m = 0.0;
    this->_m_p_m = 0.0;
    this->_T_p_m = 0.0;

    double Mtot = Function::rho_0(*_df) * std::pow(_df->Get_L(), 3);

    std::vector<double> t_temp;
    std::vector<double> x_p_temp;
    std::vector<double> v_p_temp;
    std::vector<double> r_p_temp;
    std::vector<double> m_p_temp;
    std::vector<double> m_s_temp;
    std::vector<double> T_p_temp;

    // Generate random numbers on the host
    std::random_device rd;
    std::default_random_engine seed(rd());
    std::uniform_real_distribution<double> u1{0, _df->Get_L()};

    double m_p_total = 0.0;

    // Initialize droplets until total mass exceeds Mtot
    while (m_p_total < Mtot)
    {
        double x_p_i = u1(seed);
        double r_p_i = Function::acceptation_rejet(*_df, seed);

        double volume = (4.0 / 3.0) * M_PI * std::pow(r_p_i, 3);
        double rho_p = _df->Get_rho_p();
        double m_p_i = volume * rho_p;
        double m_s_i = volume * rho_p * _df->Get_Salinity_p() / 1000.0;
        double T_p_i = _df->Get_T_p_0();

        t_temp.push_back(0.0);
        x_p_temp.push_back(x_p_i);
        v_p_temp.push_back(0.0);
        r_p_temp.push_back(r_p_i);
        m_p_temp.push_back(m_p_i);
        m_s_temp.push_back(m_s_i);
        T_p_temp.push_back(T_p_i);

        // Update total mass
        m_p_total += m_p_i;

        // Update mean values
        this->_t_m += 0.0;
        this->_x_p_m += x_p_i;
        this->_v_p_m += 0.0;
        this->_r_p_m += r_p_i;
        this->_m_p_m += m_p_i;
        this->_T_p_m += T_p_i;
    }

    // Number of droplets
    _N = t_temp.size();

    // Resize the arrays
    _t.resize(_N);
    _x_p.resize(_N);
    _v_p.resize(_N);
    _r_p.resize(_N);
    _m_p.resize(_N);
    _m_s.resize(_N);
    _T_p.resize(_N);

    // Copy data to class members
    for (int i = 0; i < _N; ++i) {
        _t[i] = t_temp[i];
        _x_p[i] = x_p_temp[i];
        _v_p[i] = v_p_temp[i];
        _r_p[i] = r_p_temp[i];
        _m_p[i] = m_p_temp[i];
        _m_s[i] = m_s_temp[i];
        _T_p[i] = T_p_temp[i];
    }

    // Calculate mean values
    this->_t_m /= double(_N);
    this->_x_p_m /= double(_N);
    this->_v_p_m /= double(_N);
    this->_r_p_m /= double(_N);
    this->_m_p_m /= double(_N);
    this->_T_p_m /= double(_N);

    std::cout << "Nombre de gouttes : " << _N << std::endl;
}

void Spray::Update()
{
    // Re-initialize mean values
    this->_t_m = 0.0;
    this->_x_p_m = 0.0;
    this->_v_p_m = 0.0;
    this->_r_p_m = 0.0;
    this->_m_p_m = 0.0;
    this->_T_p_m = 0.0;

    double dt = 1e-5;

    double U_air = _df->Get_U_air();
    double L = _df->Get_L();
    int cas = _df->Get_cas();
    double g = _df->Get_g();

    // Parallelize the loop with OpenACC
    #pragma acc parallel loop reduction(+:_t_m,_x_p_m,_v_p_m,_r_p_m,_m_p_m,_T_p_m)
    for (int i = 0; i < _N; ++i) {
        double x_p_old = this->_x_p[i];
        double v_p_old = this->_v_p[i];
        double r_p_old = this->_r_p[i];
        double m_p_old = this->_m_p[i];
        double m_s_old = this->_m_s[i];
        double T_p_old = this->_T_p[i];

        double tau_p = Function::tau_p(*_df, r_p_old, m_p_old);

        switch (cas)
        {
        case 1:
        {
            if (x_p_old <= L) this->_x_p[i] += dt * v_p_old;
            else this->_x_p[i] = 0.0;

            double F_val = Function::F(*_df, r_p_old, v_p_old, m_p_old);
            double R_val = Function::R(*_df, r_p_old, v_p_old, m_p_old, m_s_old, T_p_old);
            double M_val = Function::M(*_df, r_p_old, v_p_old, m_p_old, m_s_old, T_p_old);
            double T_val = Function::T(*_df, r_p_old, v_p_old, m_p_old, m_s_old, T_p_old);

            this->_v_p[i] += (dt / m_p_old) * F_val;
            this->_r_p[i] += dt * R_val;
            this->_m_p[i] += dt * M_val;
            this->_T_p[i] += dt * T_val;
            break;
        }
        default:
        {
            this->_x_p[i] += dt * v_p_old;
            this->_v_p[i] = v_p_old * exp(-dt / tau_p) + (U_air + g * tau_p) * (1 - exp(-dt / tau_p));

            double R_val = Function::R(*_df, r_p_old, v_p_old, m_p_old, m_s_old, T_p_old);
            double M_val = Function::M(*_df, r_p_old, v_p_old, m_p_old, m_s_old, T_p_old);
            double tau_t_val = Function::tau_t(*_df, r_p_old, this->_v_p[i], this->_m_p[i], T_p_old);
            double b_val = Function::b(*_df, r_p_old, this->_v_p[i], this->_m_p[i], m_s_old, T_p_old);

            this->_r_p[i] += dt * R_val;
            this->_m_p[i] += dt * M_val;
            this->_T_p[i] = T_p_old * exp(-dt / tau_t_val) + b_val * (1 - exp(-dt / tau_t_val));
            break;
        }
        }

        this->_t[i] += dt;

        // Update mean values
        this->_t_m += this->_t[i];
        this->_x_p_m += this->_x_p[i];
        this->_v_p_m += this->_v_p[i];
        this->_r_p_m += this->_r_p[i];
        this->_m_p_m += this->_m_p[i];
        this->_T_p_m += this->_T_p[i];
    }

    // Calculate mean values
    this->_t_m /= double(_N);
    this->_x_p_m /= double(_N);
    this->_v_p_m /= double(_N);
    this->_r_p_m /= double(_N);
    this->_m_p_m /= double(_N);
    this->_T_p_m /= double(_N);
}

void Spray::Display()
{
    std::cout << "t = " << this->_t_m << " [s] " << std::endl;
    std::cout << "x_p = " << this->_x_p_m << " [m] " << std::endl;
    std::cout << "v_p = " << this->_v_p_m << " [m/s] " << std::endl;
    std::cout << "r_p = " << this->_r_p_m << " [m] " << std::endl;
    std::cout << "m_p = " << this->_m_p_m << " [kg] " << std::endl;
    std::cout << "T_p = " << this->_T_p_m << " [K] " << std::endl;
}

void Spray::Save(std::string n_drop)
{
    std::string n_file = "../res/" + n_drop + ".dat";

    std::ofstream monflux;
    monflux.open(n_file, std::ios::app);
    if (monflux.is_open()) {
        monflux << this->_t_m << " " << this->_x_p_m << " " << this->_v_p_m << " " << this->_r_p_m * 1e6 << " " << this->_m_p_m * 1e9 << " " << this->_T_p_m - 273.15 << std::endl;
        monflux.close();
    } else {
        std::cerr << "Erreur : impossible d'ouvrir le fichier " << n_file << std::endl;
    }
}

