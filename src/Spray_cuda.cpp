#ifndef _SPRAY_CPP
#define _SPRAY_CPP

#include "Spray.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <iostream>

// Déclaration des variables globales pour CUDA
__device__ double d_t_m;
__device__ double d_x_p_m;
__device__ double d_v_p_m;
__device__ double d_r_p_m;
__device__ double d_m_p_m;
__device__ double d_T_p_m;

// Kernel pour initialiser les gouttes
__global__ void initialize_kernel(Drop** d_spray, DataFile* d_df, Function* d_fct, int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) {
        // Allocation dynamique sur le GPU
        d_spray[i] = new Drop(d_df, d_fct);
        d_spray[i]->Initialize();
        // d_spray[i]->Display(); // Affichage sur le GPU n'est pas recommandé

        // Accumulation atomique des valeurs
        atomicAdd(&d_t_m, d_spray[i]->Get_t());
        atomicAdd(&d_x_p_m, d_spray[i]->Get_x_p());
        atomicAdd(&d_v_p_m, d_spray[i]->Get_v_p());
        atomicAdd(&d_r_p_m, d_spray[i]->Get_r_p());
        atomicAdd(&d_m_p_m, d_spray[i]->Get_m_p());
        atomicAdd(&d_T_p_m, d_spray[i]->Get_T_p());
    }
}

// Kernel pour mettre à jour les gouttes
__global__ void update_kernel(Drop** d_spray, int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) {
        d_spray[i]->Update();

        // Accumulation atomique des valeurs
        atomicAdd(&d_t_m, d_spray[i]->Get_t());
        atomicAdd(&d_x_p_m, d_spray[i]->Get_x_p());
        atomicAdd(&d_v_p_m, d_spray[i]->Get_v_p());
        atomicAdd(&d_r_p_m, d_spray[i]->Get_r_p());
        atomicAdd(&d_m_p_m, d_spray[i]->Get_m_p());
        atomicAdd(&d_T_p_m, d_spray[i]->Get_T_p());
    }
}

Spray::Spray(DataFile* df, Function* fct) :
    _df(df), _fct(fct)
{

}

Spray::~Spray(){
    // Les objets alloués sur le GPU doivent être libérés sur le GPU
    // On suppose que la méthode ~Drop() est marquée __device__
    int N = _df->Get_N();
    Drop** d_spray;

    cudaMalloc((void**)&d_spray, N * sizeof(Drop*));
    cudaMemcpy(d_spray, _spray.data(), N * sizeof(Drop*), cudaMemcpyHostToDevice);

    // Libération des gouttes sur le GPU
    int blockSize = 256;
    int gridSize = (N + blockSize - 1) / blockSize;
    
    // Libération du tableau de gouttes
    cudaFree(d_spray);
}

void Spray::Initialize()
{
    this->_t_m = 0.0;
    this->_x_p_m = 0.0;
    this->_v_p_m = 0.0;
    this->_r_p_m = 0.0;
    this->_m_p_m = 0.0;
    this->_T_p_m = 0.0;

    int N = _df->Get_N();
    this->_spray.resize(N);

    // Allocation sur le GPU
    Drop** d_spray;
    cudaMalloc((void**)&d_spray, N * sizeof(Drop*));

    DataFile* d_df;
    Function* d_fct;
    cudaMalloc((void**)&d_df, sizeof(DataFile));
    cudaMalloc((void**)&d_fct, sizeof(Function));

    cudaMemcpy(d_df, _df, sizeof(DataFile), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fct, _fct, sizeof(Function), cudaMemcpyHostToDevice);

    // Initialisation des variables globales sur le GPU
    cudaMemset(&d_t_m, 0, sizeof(double));
    cudaMemset(&d_x_p_m, 0, sizeof(double));
    cudaMemset(&d_v_p_m, 0, sizeof(double));
    cudaMemset(&d_r_p_m, 0, sizeof(double));
    cudaMemset(&d_m_p_m, 0, sizeof(double));
    cudaMemset(&d_T_p_m, 0, sizeof(double));

    // Lancement du kernel
    int blockSize = 256;
    int gridSize = (N + blockSize - 1) / blockSize;
    initialize_kernel<<<gridSize, blockSize>>>(d_spray, d_df, d_fct, N);
    cudaDeviceSynchronize();

    // Récupération des sommes
    double h_t_m, h_x_p_m, h_v_p_m, h_r_p_m, h_m_p_m, h_T_p_m;
    cudaMemcpyFromSymbol(&h_t_m, d_t_m, sizeof(double), 0, cudaMemcpyDeviceToHost);
    cudaMemcpyFromSymbol(&h_x_p_m, d_x_p_m, sizeof(double), 0, cudaMemcpyDeviceToHost);
    cudaMemcpyFromSymbol(&h_v_p_m, d_v_p_m, sizeof(double), 0, cudaMemcpyDeviceToHost);
    cudaMemcpyFromSymbol(&h_r_p_m, d_r_p_m, sizeof(double), 0, cudaMemcpyDeviceToHost);
    cudaMemcpyFromSymbol(&h_m_p_m, d_m_p_m, sizeof(double), 0, cudaMemcpyDeviceToHost);
    cudaMemcpyFromSymbol(&h_T_p_m, d_T_p_m, sizeof(double), 0, cudaMemcpyDeviceToHost);

    this->_t_m = h_t_m / double(N);
    this->_x_p_m = h_x_p_m / double(N);
    this->_v_p_m = h_v_p_m / double(N);
    this->_r_p_m = h_r_p_m / double(N);
    this->_m_p_m = h_m_p_m / double(N);
    this->_T_p_m = h_T_p_m / double(N);

    cudaMemcpy(_spray.data(), d_spray, N * sizeof(Drop*), cudaMemcpyDeviceToHost);

    // Libération de la mémoire sur le GPU
    cudaFree(d_spray);
    cudaFree(d_df);
    cudaFree(d_fct);
}

void Spray::Update()
{
    this->_t_m = 0.0;
    this->_x_p_m = 0.0;
    this->_v_p_m = 0.0;
    this->_r_p_m = 0.0;
    this->_m_p_m = 0.0;
    this->_T_p_m = 0.0;

    int N = _df->Get_N();

    // Allocation sur le GPU
    Drop** d_spray;
    cudaMalloc((void**)&d_spray, N * sizeof(Drop*));
    cudaMemcpy(d_spray, _spray.data(), N * sizeof(Drop*), cudaMemcpyHostToDevice);

    // Initialisation des variables globales sur le GPU
    cudaMemset(&d_t_m, 0, sizeof(double));
    cudaMemset(&d_x_p_m, 0, sizeof(double));
    cudaMemset(&d_v_p_m, 0, sizeof(double));
    cudaMemset(&d_r_p_m, 0, sizeof(double));
    cudaMemset(&d_m_p_m, 0, sizeof(double));
    cudaMemset(&d_T_p_m, 0, sizeof(double));

    // Lancement du kernel
    int blockSize = 256;
    int gridSize = (N + blockSize - 1) / blockSize;
    update_kernel<<<gridSize, blockSize>>>(d_spray, N);
    cudaDeviceSynchronize();

    // Récupération des sommes
    double h_t_m, h_x_p_m, h_v_p_m, h_r_p_m, h_m_p_m, h_T_p_m;
    cudaMemcpyFromSymbol(&h_t_m, d_t_m, sizeof(double), 0, cudaMemcpyDeviceToHost);
    cudaMemcpyFromSymbol(&h_x_p_m, d_x_p_m, sizeof(double), 0, cudaMemcpyDeviceToHost);
    cudaMemcpyFromSymbol(&h_v_p_m, d_v_p_m, sizeof(double), 0, cudaMemcpyDeviceToHost);
    cudaMemcpyFromSymbol(&h_r_p_m, d_r_p_m, sizeof(double), 0, cudaMemcpyDeviceToHost);
    cudaMemcpyFromSymbol(&h_m_p_m, d_m_p_m, sizeof(double), 0, cudaMemcpyDeviceToHost);
    cudaMemcpyFromSymbol(&h_T_p_m, d_T_p_m, sizeof(double), 0, cudaMemcpyDeviceToHost);

    this->_t_m = h_t_m / double(N);
    this->_x_p_m = h_x_p_m / double(N);
    this->_v_p_m = h_v_p_m / double(N);
    this->_r_p_m = h_r_p_m / double(N);
    this->_m_p_m = h_m_p_m / double(N);
    this->_T_p_m = h_T_p_m / double(N);

    cudaMemcpy(_spray.data(), d_spray, N * sizeof(Drop*), cudaMemcpyDeviceToHost);

    // Libération de la mémoire sur le GPU
    cudaFree(d_spray);
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
        monflux << this->_t_m << " " << this->_x_p_m << " " << this->_v_p_m << " " << this->_r_p_m*1e6 << " " << this->_m_p_m*1e9 << " " << this->_T_p_m - 273.15 << std::endl;
        monflux.close();
    } else {
        std::cerr << "Erreur : impossible d'ouvrir le fichier " << n_file << std::endl;
    }
}

#endif

