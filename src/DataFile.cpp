#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include "../data/toml.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

DataFile::DataFile(std::string file_name)
: _file_name(file_name)
{
   // Lecture du fichier de données
   auto config = toml::parse(file_name);

   // Other
   const auto& parameter = toml::find(config, "parameter");
   this->_U_air = toml::find<double>(parameter, "U_air");
   this->_T_air_celcius = toml::find<double>(parameter, "T_air_celcius");
   this->_T_p_0_celcius = toml::find<double>(parameter, "T_p_0_celcius");
   this->_Q_RH = toml::find<double>(parameter, "Q_RH");
   this->_g = toml::find<double>(parameter, "g");
   this->_M_w = toml::find<double>(parameter, "M_w");
   this->_M_s = toml::find<double>(parameter, "M_s");
   this->_M_a = toml::find<double>(parameter, "M_a");
   this->_R_g = toml::find<double>(parameter, "R_g");
   this->_p_0 = toml::find<double>(parameter, "p_0");
   this->_Salinity_w = toml::find<double>(parameter, "Salinity_w");
   this->_Salinity_p = toml::find<double>(parameter, "Salinity_p");
   this->_r_p_0 = toml::find<double>(parameter, "r_p_0");
   this->_Delta_v = toml::find<double>(parameter, "Delta_v");
   this->_alpha_c = toml::find<double>(parameter, "alpha_c");
   this->_Delta_T = toml::find<double>(parameter, "Delta_T");
   this->_alpha_T = toml::find<double>(parameter, "alpha_T");
   this->_I = toml::find<int>(parameter, "I");
} 

#define _DATA_FILE_CPP
#endif