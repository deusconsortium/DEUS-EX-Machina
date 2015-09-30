/* 
 * File:   Tools.h
 * Author: pdefromont
 *
 * Created on 9 septembre 2015, 10:25
 */

#ifndef TOOLS_H
#define	TOOLS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include "stdlib.h"  
#include "math.h"   
#include <sstream>
#include <sys/stat.h>

struct FVector{
    float x,y,z;
};

class Tools{
public:
    Tools();
    static std::string IntToString(int const number,bool formated=false,int nb_chiffres=5);
    static int StringToInt(std::string const Snumber);
    static float Random(float min, float max);
    
    //gestion des tableaux
    static int getMinimumIndex(const std::vector<float> & x, const std::vector<float> & y);
    static void getMeanAndSigma(float & mean,float & sigma,const std::vector<float> tab);
    static float getXvalue(const float y0, const std::vector<float> & x, const std::vector<float> & y);
    
    //gestion des fichiers
    static int getNbcolumn(std::string const file_path);
    static int getNbLines(std::string const file_path);
    static bool fileExists(std::string const file_path);
    static bool createFolder(std::string const folder_name);
    
    //fonctions particulières
    static int getLastOutput(const int boxlen, const int npart, const std::string cosmo);
};

#endif	/* TOOLS_H */

