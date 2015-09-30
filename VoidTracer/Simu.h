/* 
 * File:   Simu.h
 * Author: Paul de Fromont
 *
 * Created on 22 octobre 2013, 15:58
 */

#ifndef SIMU_H
#define	SIMU_H

#include "Tools.h"
#include "WaitingBar.h"
#include "FOFReaderLib/FOFReaderLib.h"
#include "FOFReaderLib/FOFFiles/FOFParticles.h"

#define Pi 3.14159265359
struct SimuData
{
    float unit_l,unit_d,unit_t,unit_m,unit_v;
    float a_exp,H0,omega_m,omega_l,omega_k,omega_b;
};

using namespace std;

class Simu 
{
public:
    Simu();
    Simu(int const Boxlen, int const Npart,string const cosmo);
    Simu(int const Boxlen, int const Npart,string const cosmo,int const output);
    virtual ~Simu();
   
    bool load(int const Boxlen, int const Npart,string const cosmo,int output = -1);
    string saveHalosPositions(const int min_particles = 0, const int max_particles = -1)  const;
    void profileAnalysis(const string position_file,const string directory_name = "",int const NobjectsMax = -1);
    
    //setters et getters
    void setDataDirectory(const string data){
        _data_path = data;
    }
    string getDataDirectory()   const{return _data_path;}
    void setIsOverDensity(bool isOverDensity){ _isOverDensity = isOverDensity;}
    bool isOverDensity(){return _isOverDensity;}
    int getBoxlen() const{return _boxlen;}
    int getNpart()  const{return _npart;}
    string getCosmo()   const{return _cosmo;}
    int getOutput() const{return _output;}
    
private:    
    float getR1_ramses(vector<float> & f_profile, vector<float> & r_ramses,bool isHalo);
    float getR1_delta(vector<float> & f_profile, vector<float> & r_ramses,bool isHalo);
    void ProfileAroundPosition(FVector Position ,vector<float> & f,vector<float> & v,FOFMultiCube & multi,vector<float> const radius_ramses);

    //les propriétés de la simu
    string _data_path,_save_name,_cosmo;
    string _simu_name;
    int _boxlen,_npart,_output,_Numpy_cores;
    bool _isOverDensity;
    
    //les propriétés physiques de la simu
    SimuData _data;
    
    //les constantes physiques
    float Mpc,c0,SunMass; 
    
    //les constantes de tracage
    int _Npoints;
    float _dr_ramses,_r0_ramses;
};

#endif	/* SIMU_H */