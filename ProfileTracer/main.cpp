/* 
 * File:   main.cpp
 * Author: pdefromont
 *
 *
 * Created on 9 septembre 2015, 10:23
 */

#include <cstdlib>

#include "Simu.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    
    //interface utilisateur
    /*int boxlen, npart,output,Max;
    string cosmo,name;
    bool bo;
    string path = "none";
            
    Simu my_simu = Simu();*/
    
    int outputs[8] = {56,52,45,39,25,18,6,1};
    int boxlen = 2592;
    int npart = 2048;
    
    Simu my_simu = Simu();
    string name = "work_Z_200_210_"+Tools::IntToString(boxlen)+"_"+Tools::IntToString(npart)+"_output"+Tools::IntToString(outputs[0]);
    my_simu.load(boxlen,npart,"lcdmw5",outputs[0]);
    my_simu.setIsOverDensity(true);
    string path = my_simu.saveHalosPositions(200,210);
    my_simu.profileAnalysis(path,name,10000);
    
    for(unsigned int i=1; i<8 ; i++)
    {
        Simu my_simu = Simu();
        string name = "work_Z_200_210_"+Tools::IntToString(boxlen)+"_"+Tools::IntToString(npart)+"_output"+Tools::IntToString(outputs[i]);
        my_simu.load(boxlen,npart,"lcdmw5",outputs[i]);
        my_simu.setIsOverDensity(true);
        my_simu.profileAnalysis(path,name);
    }
    
    /*
    
    cout<<endl<<"##################################################"<<endl;
    cout<<"#                ProfileTracer v1.0              #"<<endl;
    cout<<"##################################################"<<endl<<endl;
    cout<<"Enter Boxlen : ";
    cin >> boxlen;
    cout<<"Enter Npart : ";
    cin >> npart;
    cout<<"Enter cosmo : ";
    cin >> cosmo;
    cout << "Enter Output : ";
    cin >> output;
    my_simu.load(boxlen,npart,cosmo,output);
    
    cout<<"Is it over-density ? (0/1) ";
    cin >> bo;
    my_simu.setIsOverDensity(bo);
    if(bo)
    {
        cout << "Do you want to compute halo position ? (0/1) ";
        cin >> bo;
        if(bo)
        {
            int min,max;
            cout << "Enter min halo mass : ";
            cin >> min;
            cout << "Enter max halo mass : ";
            cin >> max;
            path = my_simu.saveHalosPositions(min,max);
        }
    }
    else
    {
        cout << "Do you want to compute voids position ? (0/1) ";
        cin >> bo;
        if(bo)
        {
            float seuil_min,seuil_mean,density_max;
            cout << "Enter seuil min : ";
            cin >> seuil_min;
            cout << "Enter seuil mean : ";
            cin >> seuil_mean;
            cout << "Enter max density : ";
            cin >> density_max;
            path = my_simu.saveVoidPositions(seuil_min,seuil_mean,density_max);
        }
    }
    if(path == "none")
    {
        cout << "Enter file name of the position File : ";
        cin >> path;
    }
    cout << "Enter output name (can be '') : ";
    cin >> name;
    cout << "Enter max objects : ";
    cin >> Max;
    my_simu.profileAnalysis(path,name,Max);
    
    */
    return 0;
}

