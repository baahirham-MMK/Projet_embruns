#ifndef _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <random>
#include <algorithm>
#include <mpi.h>

// Définition de la classe

class DataFile {

   private:
      const std::string _file_name;
      double _U_air, _T_air_celcius, _T_p_0_celcius, _Q_RH, _g, _M_w, _M_s, _M_a, _R_g, _Salinity_w, _Salinity_p, _p_0;
      double _Delta_v, _alpha_c, _Delta_T, _alpha_T, _T_f, _L;
      int _I, _cas;

   public: // Méthodes et opérateurs de la classe
   DataFile(std::string file_name);
   const double Get_U_air() const {return _U_air;};
   const double Get_T_air_celcius() const {return _T_air_celcius;};
   const double Get_T_p_0_celcius() const {return _T_p_0_celcius;};
   const double Get_Q_RH() const {return _Q_RH;};
   const double Get_g() const {return _g;};
   const double Get_k_a() const {return 0.02411*(1+3.309e-3*(Get_T_p_0()-273.15) -1.441e-6*std::pow(Get_T_p_0()-273.15,2));};
   const double Get_M_w() const {return _M_w;};
   const double Get_M_s() const {return _M_s;};
   const double Get_M_a() const {return _M_a;};
   const double Get_R_g() const {return _R_g;};
   const double Get_p_0() const {return _p_0;};
   const double Get_Salinity_w() const {return _Salinity_w;};
   const double Get_Salinity_p() const {return _Salinity_p;};
   const double Get_Delta_v() const {return _Delta_v;};
   const double Get_alpha_c() const {return _alpha_c;};
   const double Get_Delta_T() const {return _Delta_T;};
   const double Get_alpha_T() const {return _alpha_T;};
   const int Get_I() const {return _I;};
   const double Get_T_air() const {return _T_air_celcius + 273.15;};
   const double Get_T_p_0() const {return _T_p_0_celcius + 273.15;};
   const double Get_qs0() const {return 0.62197*611.2/Get_p_0()*std::exp((17.67*(Get_T_p_0()-273.15))/(Get_T_p_0()-29.66));};
   const double Get_rs0() const {return Get_qs0()/(1-Get_qs0());};
   const double Get_r0() const {return 98.0/100.0*Get_rs0();};
   const double Get_q0() const {return Get_r0()/(1+Get_r0());};
   const double Get_Tv0() const {return Get_T_p_0()*(1+0.6078*Get_q0());};
   const double Get_rho_air() const {return 1.2929*273.15/Get_Tv0();};
   const double Get_nu_air() const {return 1.33e-5 +(0.0084*(Get_Tv0()-273.15))*1e-5;};
   const double Get_mu_air() const {return Get_nu_air()*Get_rho_air();};
   const double Get_p_air() const {return Get_p_0()-Get_rho_air()*Get_g()*10.0;};
   const double Get_p_v_sat_T_air() const {return 2.33e3;}; // ????
   const double Get_L_v() const {return (25.00895-0.02274*(Get_T_air()-273.15))*1e5;};
   const double Get_D_v() const {return 2.11e-5*std::pow((Get_T_air()/273.15),1.94)*(1013.25/(Get_p_air()/100.0));};
   const double Get_Gamma_p() const {return (75.63 - 0.144*Get_T_p_0_celcius() + 0.221*Get_Salinity_w())*1e-3;};
   const double Get_Phi_s() const {return 0.91154614+1.7317496707e-4*Get_Salinity_p()+4.7616058412e-6*std::pow(Get_Salinity_p(),2)-9.2541509027e-9*std::pow(Get_Salinity_p(),3)+7.3475024678e-12*std::pow(Get_Salinity_p(),4);};
   const double Get_rho_w() const {return (999.842594+6.793952e-2*Get_T_p_0_celcius()-9.095290e-3*std::pow(Get_T_p_0_celcius(),2)+1.001685e-4*std::pow(Get_T_p_0_celcius(),3)-1.120083e-6*std::pow(Get_T_p_0_celcius(),4)+6.536332e-9*std::pow(Get_T_p_0_celcius(),5));};
   const double Get_rho_p() const {return Get_rho_w() + Get_Salinity_p()*(0.824493-4.0899e-3*Get_T_p_0_celcius() +7.6438e-5*std::pow(Get_T_p_0_celcius(),2)-8.2467e-7*std::pow(Get_T_p_0_celcius(),3)+5.3875e-9*std::pow(Get_T_p_0_celcius(),4))+std::pow(Get_Salinity_p(),(3./2.))*(-5.72466e-3+1.0227e-4*Get_T_p_0_celcius() -1.6546e-6*std::pow(Get_T_p_0_celcius(),2)) +4.8314e-4*std::pow(Get_Salinity_p(),2);};
   const double Get_c_p_s() const {return 4217.4 -3.720283*(Get_T_air()-273.15)+0.1412855*std::pow((Get_T_air()-273.15),2)-2.654387e-3*std::pow((Get_T_air()-273.15),3) +2.093236e-5*std::pow((Get_T_air()-273.15),4) + Get_Salinity_p()*(-7.6444+0.107276*(Get_T_air()-273.15)-1.3839e-3*std::pow(Get_T_p_0_celcius(),2))+std::pow(Get_Salinity_p(),(3/2))*(0.17709-4.0772e-3*(Get_T_air()-273.15)+5.3539e-5*std::pow(Get_T_p_0_celcius(),2));};
   const double Get_c_p_air_sec() const {return 1.9327e-10*std::pow(Get_T_p_0(),4)-7.9999e-7*std::pow(Get_T_p_0(),3)+1.1407e-3*std::pow(Get_T_p_0(),2)-4.4890e-1*Get_T_p_0()+1.0575e+3;};
   const double Get_rho_air_sec() const {return 1.2929*273.13/Get_T_air();};
   const int Get_cas() const {return _cas;};
   const double Get_T_f() const {return _T_f;};
   const double Get_L() const {return _L;};
};

#define _DATA_FILE_H
#endif