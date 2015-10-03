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
    int boxlen, npart,output,Max;
    string cosmo,name;
    bool bo;
    string path = "none";
            
    Simu my_simu = Simu();
    
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
    return 0;
}

