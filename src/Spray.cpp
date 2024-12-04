#ifndef _SPRAY_CPP
#include "Spray.h"


Spray::Spray(DataFile* df, Function* fct) :
_df(df), _fct(fct)
{

} 

Spray::~Spray(){
    int N = this->_spray.size();
    for (int i = 0; i < N; ++i){
        delete this->_spray[i];
    }
}

void Spray::Initialize()
{
    this->_t_m = 0.0;
    this->_x_p_m = 0.0;
    this->_v_p_m = 0.0;
    this->_r_p_m = 0.0;
    this->_m_p_m = 0.0;
    this->_T_p_m = 0.0;
    this->_M_n = 0.0;
    this->_M_s = 0.0;
    this->_Vol = 0.0;

    double Mtot = _fct->rho_0()*pow(_df->Get_L(),3);
    int N, Me, Np;
    MPI_Comm_rank(MPI_COMM_WORLD, &Me);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);

    while(this->_M_n < Mtot) 
    {
        this->_spray.push_back(new Drop(this->_df, this->_fct));
        N = this->_spray.size();
        this->_spray[N-1]->Initialize();
        this->_t_m += this->_spray[N-1]->Get_t();
        this->_x_p_m += this->_spray[N-1]->Get_x_p();
        this->_v_p_m += this->_spray[N-1]->Get_v_p();
        this->_r_p_m += this->_spray[N-1]->Get_r_p();
        this->_m_p_m += this->_spray[N-1]->Get_m_p();
        this->_T_p_m += this->_spray[N-1]->Get_T_p();
        MPI_Allreduce(&this->_m_p_m, &this->_M_n, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&this->_r_p_m, &this->_Vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    printf("Me = %d, Nombre de goutte : %d\n", Me, N);
    if(Me == 0){
        printf("Mtot = %e [kg], M_0 = %e [kg]\n", Mtot, this->_M_n);
    }

    this->_t_m *= (1./double(N));
    this->_x_p_m *= (1./double(N));
    this->_v_p_m *= (1./double(N));
    this->_r_p_m *= (1./double(N));
    this->_m_p_m *= (1./double(N));
    this->_T_p_m *= (1./double(N));
    this->_M_0 = this->_M_n;
    this->_Vol = 4.0/3.0*std::acos(-1.0)*std::pow(this->_Vol,3);
    this->_M_s = (this->_M_0 - this->_M_n);
}

void Spray::Update()
{
    this->_t_m = 0.0;
    this->_x_p_m = 0.0;
    this->_v_p_m = 0.0;
    this->_r_p_m = 0.0;
    this->_m_p_m = 0.0;
    this->_T_p_m = 0.0;
    this->_M_n = 0.0;
    this->_Vol = 0.0;

    int N = this->_spray.size();
    
    for (int i = 0; i < N; ++i) 
    {
        this->_spray[i]->Update();
        this->_t_m += this->_spray[i]->Get_t();
        this->_x_p_m += this->_spray[i]->Get_x_p();
        this->_v_p_m += this->_spray[i]->Get_v_p();
        this->_r_p_m += this->_spray[i]->Get_r_p();
        this->_m_p_m += this->_spray[i]->Get_m_p();
        this->_T_p_m += this->_spray[i]->Get_T_p();
    }

    this->_t_m *= (1./double(N));
    this->_x_p_m *= (1./double(N));
    this->_v_p_m *= (1./double(N));
    MPI_Allreduce(&this->_r_p_m, &this->_Vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    this->_r_p_m *= (1./double(N));
    MPI_Allreduce(&this->_m_p_m, &this->_M_n, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    this->_m_p_m *= (1./double(N));
    this->_T_p_m *= (1./double(N));
    this->_Vol = 4.0/3.0*std::acos(-1.0)*std::pow(this->_Vol,3);
    this->_M_s = (this->_M_0 - this->_M_n);
}

void Spray::Display()
{
    std::cout << "t = " << this->_t_m << " [s] " << std::endl;
    std::cout << "x_p = " << this->_x_p_m << " [m] " << std::endl;
    std::cout << "v_p = " << this->_v_p_m << " [m/s] " << std::endl;
    std::cout << "r_p = " << this->_r_p_m << " [m] " << std::endl;
    std::cout << "m_p = " << this->_m_p_m << " [kg] " << std::endl;
    std::cout << "M_n = " << this->_M_n << " [kg] " << std::endl;
    std::cout << "M_s = " << this->_M_s << " [kg] " << std::endl;
    std::cout << "T_p = " << this->_T_p_m << " [K] " << std::endl;
}

void Spray::Save(std::string n_drop)
{
    int Np, Me;
    MPI_Comm_rank(MPI_COMM_WORLD, &Me);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);
    double t_m_tot(0.0), x_p_m_tot(0.0), v_p_m_tot(0.0), r_p_m_tot(0.0), m_p_m_tot(0.0), T_p_m_tot(0.0), humidity, abs_hum, QQ;

    std::string n_file = "../res/" + n_drop + ".dat";

    MPI_Allreduce(&this->_t_m, &t_m_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&this->_x_p_m, &x_p_m_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&this->_v_p_m, &v_p_m_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&this->_r_p_m, &r_p_m_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&this->_m_p_m, &m_p_m_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&this->_T_p_m, &T_p_m_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    t_m_tot /= double(Np);
    x_p_m_tot /= double(Np);
    v_p_m_tot /= double(Np);
    r_p_m_tot /= double(Np);
    m_p_m_tot /= double(Np);
    T_p_m_tot /= double(Np);

    humidity = (1.0-_df->Get_q10())/_df->Get_rho_air()*(this->_M_s/(std::pow(_df->Get_L(),3)));
    QQ = (_df->Get_q10()+humidity)/(1-(_df->Get_q10()+humidity))/_df->Get_rs10();
    
    if (Me == 0){
        std::ofstream monflux;
        monflux.open(n_file, std::ios::app);  
        if (monflux.is_open()) {
            monflux << this->_t_m << " " << x_p_m_tot << " " << v_p_m_tot << " " << r_p_m_tot*1e6 << " " << m_p_m_tot*1e9 << " " << T_p_m_tot - 273.15 << " " << this->_M_s << " " << QQ*100 << std::endl;
            monflux.close();
        } else {
            std::cerr << "Erreur : impossible d'ouvrir le fichier " << n_file << std::endl;
        }
    }
}


#define _SPRAY_CPP
#endif