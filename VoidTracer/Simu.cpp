/* 
 * File:   Simu.cpp
 * Author: Paul
 * 
 * Created on 22 octobre 2013, 15:58
 */

#include "Simu.h"

Simu::Simu():_data_path("/data_bingo/Babel/")
{
    _simu_name = "";
    _save_name = "";
    _cosmo = "";
    _Numpy_cores = 32;
    _boxlen = 0;
    _npart = 0;
    _output = 0;
    _DrCoarseGrid = 2.0;
    _R0CoarseGrid = 2.0,
    _isOverDensity = false;
    _RmaxCoarseGrid = 100.;
    _Numpy_cores = 32;
    
    Mpc = 3.08567758e+22;
    SunMass	= 1.989e+30;
    c0 = 299792458.0; 
}

Simu::Simu(int const Boxlen, int const Npart,string const cosmo,int const output)
{
    _simu_name = "";
    _save_name = "";
    _cosmo = cosmo;
    _Numpy_cores = 32;
    _boxlen = Boxlen;
    _npart = Npart;
    _output = output;
    _DrCoarseGrid = 2.0;
    _R0CoarseGrid = 2.0,
    _isOverDensity = false;
    _RmaxCoarseGrid = 100.;
    _Numpy_cores = 32;
    
    Mpc = 3.08567758e+22;
    SunMass	= 1.989e+30;
    c0 = 299792458.0; 
    
    load(Boxlen,Npart,cosmo,output);
}
 
Simu::~Simu() 
{
}

bool Simu::load(int const Boxlen, int const Npart,string cosmo,int output)
{    
    //conversion de cosmo :
    _boxlen = Boxlen;
    _npart = Npart;
    _cosmo = cosmo;
    _output = output;    
    
    _simu_name= _data_path + "boxlen"+Tools::IntToString(_boxlen)+"_n"+Tools::IntToString(_npart)+"_"+_cosmo;   
    _save_name = "boxlen"+Tools::IntToString(_boxlen)+"_n"+Tools::IntToString(_npart)+"_"+_cosmo;

    //chargement du fichier info    
    string info_path = _simu_name + "/cube"+ Tools::IntToString(_output,true) + "/info"+Tools::IntToString(_output,true)+".txt";

    ifstream info_file( info_path.c_str() ); 
    if ( info_file )
    {        
        string ligne,poubelle,utile; // variable contenant chaque ligne lue
        for(int i(0); i<9 ; i++)
                getline( info_file, ligne);
        
        info_file >> poubelle >> poubelle >> _data.a_exp;
        info_file >> poubelle >> poubelle >> _data.H0;
        info_file >> poubelle >> poubelle >> _data.omega_m;
        info_file >> poubelle >> poubelle >> _data.omega_l;
        info_file >> poubelle >> poubelle >> _data.omega_k;
        info_file >> poubelle >> poubelle >> _data.omega_b;
        info_file >> poubelle >> poubelle >> _data.unit_l;
        info_file >> poubelle >> poubelle >> _data.unit_d;
        info_file >> poubelle >> poubelle >> _data.unit_t;
        
        if(_data.unit_t == 0.0)
        {
            cout<<endl<<"\n\nERROR !!  no time unit"<<endl<<endl;
            return false;
        }            

        //transcription de la densite en masse :
        
        double d0 = (double) _npart;
        d0 = d0*d0*d0;
        
        //conversion en unité standard : unit_m en Mo ; unot_v en 'c'
        
        _data.unit_l = _data.unit_l/100.0;
        _data.unit_d = _data.unit_d*1000.0;        
        _data.unit_m = 1/SunMass*pow(_data.unit_l,3.0)*_data.unit_d/d0;        
        _data.unit_v = _data.unit_l/(c0*_data.unit_t);
        
        info_file.close();
    }
    else
    {
        cout<<endl<<"\n\nERROR !!  file "<<info_path<<" doesn't exists. leaving"<<endl<<endl;
        return false;
    }
    
    cout<<endl;
    cout<<"_____________________________________"<<endl<<endl;
    cout<<"        Simu well loaded "<<endl;
    cout<<"_____________________________________"<<endl<<endl;
    cout<<"     - Npart = "<<_npart <<endl;
    cout<<"     - Boxlen = "<<_boxlen <<endl;
    cout<<"     - cosmo = "<<_cosmo <<endl;
    cout<<"     - H0 = "<<_data.H0<<" km/s/Mpc"<<endl;
    cout<<"     - a0 = "<<_data.a_exp<<endl;
    cout<<"     - Mp = "<<_data.unit_m<<" [Mo]"<<endl;
    cout<<endl<<"_____________________________________"<<endl<<endl;
        
    return true;
}

void Simu::ProfileAroundPosition(FVector Position ,vector<float> & f,vector<float> & v,FOFMultiCube & multi,vector<float> const radius_ramses)
{
    int taille = radius_ramses.size();
    if(taille == 0 || f.size() != radius_ramses.size())
    {
        cout<<endl<<"ERROR : f.size() is null - or f has a different size of r - for ProfileAroundPosition. Leaving"<<endl<<endl;
        return;
    }
    
    float X = Position.x;
    float Y = Position.y;
    float Z = Position.z;
    
    float Mass[taille];
    float Speed[taille];
    float Nb_part[taille];

    for(int i(0) ; i < f.size() ; i++)
    {
        Mass[i] = 0.0;
        Speed[i] = 0.0;
        Nb_part[i] = 0.0;
    }
    
    float Rmax = radius_ramses[radius_ramses.size() - 1];
    
    DEUSArea Zone(X,Y,Z,2.*Rmax);

    for(int cube (0) ; cube < multi.nCubes() ; cube ++) 
    {
        if(Zone.cubeIntersect(multi.cubes(cube)->boundaries()))
        {
            for(int particle(0) ; particle < multi.cubes(cube)->npart(); particle ++)
            {
                float xp = multi.cubes(cube)->posX(particle);
                float yp = multi.cubes(cube)->posY(particle);
                float zp = multi.cubes(cube)->posZ(particle);
                
                float vx = multi.cubes(cube)->velX(particle);
                float vy = multi.cubes(cube)->velY(particle);
                float vz = multi.cubes(cube)->velZ(particle);
                
                if(Zone.particuleIsInside(xp,yp,zp))
                {
                    float D = sqrt((xp-X)*(xp-X) + (yp-Y)*(yp-Y) + (zp-Z)*(zp-Z));
                    float Dr = radius_ramses[3] - radius_ramses[2];
                    float Reff = 0.620350491/_npart;
                    
                    float vr = ((xp-X)*vx + (yp-Y)*vy + (zp-Z)*vz)/D;
                    //on cherche l'index correspondant
                    
                    for(unsigned int index(0) ; index < radius_ramses.size()  ;  index++)
                    {
                        unsigned int ri = radius_ramses.size()-1 - index;
                        float a = (radius_ramses[ri] - D)/Reff;
                        float b = (D - radius_ramses[ri])/Dr;
                        if( b>=-1.0 and b<=1.0)
                        {
                            Speed[ri] += vr;
                            Nb_part[ri] += 1.0;
                        }
                        
                        if(b<= 0)
                            Mass[ri] += 1.0;
                        else
                            break;
                        /*
                        if(a < -1)
                            break;
                        else if(a >= 1)
                        {
                            Mass[ri] += 1.0;
                        }
                        else if(a > -1.0 && a < 1.0)
                        {
                            Mass[ri] += (2.0 - a)*(1.0 + a)*(1.0 + a)/4.0;
                        }*/
                    }
                }
            }
        }
    }
    
    //normalisation :
    
    double d0 = (double) _npart;
    d0 = d0*d0*d0;  
    
    for(int i(0) ; i < f.size() ; i++)
    {
        f[i] = 3.0*Mass[i]/(4.0*Pi*d0*pow(radius_ramses[i],3.0));
        v[i] = Speed[i]*_data.unit_v/Nb_part[i];
    }
}

string Simu::saveHalosPositions(const int min_particles, const int max_particles)   const
{
    string directory_path = _simu_name + "/post/fof/output" + Tools::IntToString(_output,true) + "/";
    DEUSHalos* Halos = new DEUSHalos(directory_path.c_str());
    
    if(Halos == NULL)
        cout <<"\n\nERROR !!  No halos founded !"<<endl;
    else{
        string file_name = "halo_boxlen"+Tools::IntToString(_boxlen)+"_npart"+Tools::IntToString(_npart)+"_"+_cosmo+"_output"+Tools::IntToString(_output)+".txt";
        ofstream new_file(("data/PositionFiles/"+file_name).c_str(), ios::out | ios::trunc);
        if(!new_file)
            cout<<"\n\nERROR !!  Impossible to create file "<<file_name<<endl;
        else
        {
            int _nb_halos = Halos->nHalos();
            for(int h(0) ; h < _nb_halos ; h++)
            {
                DEUSHalo* temp = Halos->halos(h);

                if(temp->mass() >= min_particles && (max_particles == -1 || (max_particles > 0 && temp->mass() <= max_particles)))
                    new_file << temp->x() << "\t" << temp->y() << "\t" << temp->z()  << endl;  
                
            }
            new_file.close();
        }
        cout<<"Halo position well registered !" << endl;
        delete Halos;
        return file_name;
    }
    return "";
}

void Simu::profileAnalysis(const string position_file,const string directory_name,int const Nmax)
{
    const unsigned int NumProc = _Numpy_cores;
    _save_name = _save_name + directory_name;
    
    string file_path = "data/PositionFiles/" + position_file;
    string new_name = "data/PositionFiles/original.txt";
    
    int Nobjects = Tools::getNbLines(file_path);
    ifstream file(file_path.c_str());
    
    if(!file)
    {
        cout << "\nERROR !! No file "<<file_path<<" founded ... Leaving."<<endl;
        return;
    }
    
    //gestion du nombre de vides
    if(Nmax >= 0 and Nmax < Nobjects)
    {
        vector<int> used_voids;
        cout << "Selecting randomly " << Nmax << " positions on " << Nobjects << " ..."<< endl;
        //création du tableau des valeurs disponibles
        vector<unsigned int> original(Nobjects);
        for(unsigned int i(0) ; i < Nobjects ; i++)
            original[i] = i;
        
        //on remplit le tableau
        WaitingBar bar1(Nmax);
        for(unsigned int i(0) ; i < Nmax ; i++)
        {
            unsigned int index = floor(Tools::Random(0,original.size()));
            used_voids.push_back(index);
            original.erase(original.begin() + index);
            bar1.Up();
        }

        vector<float> X(Nobjects);
        vector<float> Y(Nobjects);
        vector<float> Z(Nobjects);
        
        cout<<endl<<"reading original file ..." << endl;
        
        string ligne;
        unsigned int i = 0;
        while(getline(file,ligne))
        {
            file >> X[i] >> Y[i] >> Z[i];
            i++;
        }
                
        //on ferme le fichier et on crée le nouveau
        file.close();
        
        cout << "renaming it as 'original.txt'" << endl;
        rename(file_path.c_str(),new_name.c_str());
        
        //on recrée le nouveau fichier ...
        cout<<"creating new_file named '"<<file_path<<"'"<<endl;
        ofstream new_file(file_path.c_str(), ios::out | ios::trunc);
        if (!new_file)
        {
            cout<<"ERROR : impossible to create file "<<file_path<<" ! Leaving"<<endl;
            rename(new_name.c_str(),file_path.c_str());
            return;
        }
        cout<<"filling the new file ..." << endl;
        for(unsigned int i(0) ; i < Nmax ; i++)
            new_file << X[used_voids[i]] << "\t" << Y[used_voids[i]] << "\t" << Z[used_voids[i]]  << endl;
        cout<<"done. Closing files"<<endl;
        new_file.close();
    } 
    else
        file.close();
    
    int _Npoints = floor((_RmaxCoarseGrid - _R0CoarseGrid)/_DrCoarseGrid) + 1;
    
    vector<float> r_ramses(_Npoints);
    float Rcell = 1.0/_npart;
    r_ramses[0] = _R0CoarseGrid*Rcell;    
    
    for(unsigned int i(1) ; i < _Npoints ; i++)
        r_ramses[i] = r_ramses[i-1] + _DrCoarseGrid*Rcell;
    
    vector<FVector> void_position;
    
    cout<<"Looking for file '"<<file_path<<"'"<<endl;

    ifstream part_file(file_path.c_str(), ios::out);
    if (!part_file)
    {
        cout<<"ERROR : no file '"<<file_path<<"' exists ! Leaving"<<endl;
        return;
    }
    
    if(Tools::fileExists(new_name))
    {
        rename(file_path.c_str(),(file_path.substr(0,file_path.size()-4)+"_"+Tools::IntToString(Nmax,0,false)+".txt").c_str());
        rename(new_name.c_str(),file_path.c_str());
    }
    
    string ligne;
    unsigned int i=0;
    while(getline(part_file,ligne))
    {
        float gx,gy,gz;   
        part_file >> gx >> gy >> gz;

        void_position.push_back(FVector());
        void_position[i].x = gx;
        void_position[i].y = gy;
        void_position[i].z = gz;
        i++;
    }

    part_file.close();

    cout << "Position file well readed, founded N = " << void_position.size() << " elements" << endl;
    
    
    string path = _simu_name + "/cube" + Tools::IntToString(_output,true) + "/";
    cout<<"Looking in '"<< path <<"'"<< endl;
    cout<<"Loading all the simulation particles ..." <<endl;
    FOFMultiCube multi(path,FOFParticles::DONT_READ_PARTICLES);
    
    if(multi.nCubes() == 0)
    {
        cout<< "\n\n!!! No Cube founded ... Leaving !!!\n"<<endl;
        cout<<"path : "<<path<<endl;
        
        return;
    }
    
    WaitingBar bar(multi.nCubes());
    #pragma omp parallel for num_threads(NumProc)
    for(unsigned int i=0; i < multi.nCubes() ; i++)
    {
        multi.cubes(i)->readParticles(FOFParticles::READ_POS | FOFParticles::READ_IDS | FOFParticles::READ_VEL);
        bar.Up();
    }
    
    cout<<endl<<"------------------------------------"<<endl<<endl;
    cout << "Begining the profile tracer with : \n\t- Npoints = " << _Npoints << "\n\t- dr = " << _DrCoarseGrid <<"/npart\n\t- N = " << void_position.size() << endl;

    vector<float> _f[void_position.size()];
    vector<float> _v[void_position.size()];
    for(unsigned int i(0) ; i < void_position.size() ; i++)
        for(unsigned int r(0) ; r < _Npoints ; r++)
        {
            _f[i].push_back(0.0);
            _v[i].push_back(0.0);
        }
    cout<<"------------------------------------"<<endl<<endl;
    cout<<"\ncomputing profiles ... "<<endl;
    
    bar.Reset(void_position.size());
    #pragma omp parallel for num_threads(NumProc)
    for(unsigned int i=0 ; i < void_position.size() ; i++)
    {
        ProfileAroundPosition(void_position[i],_f[i],_v[i],multi,r_ramses);
        bar.Up();
    }

    for(unsigned int i(0); i < multi.nCubes() ; i++)
        multi.cubes(i)->releaseParticles();

    cout<<endl<<"------------------------------------"<<endl<<endl;

    cout<< endl << "Profile computation done ! Saving profiles ..." <<endl;
    
    FILE* save_file = fopen((_save_name + ".bin").c_str(),"wb");
    if(save_file != NULL)
    {
        fwrite( &Nmax , sizeof(int) , 1 , save_file);
        int N = r_ramses.size();
        fwrite( &N , sizeof(int) , 1 , save_file);
        fwrite( &_R0CoarseGrid , sizeof(int) , 1 , save_file);
        fwrite( &_DrCoarseGrid , sizeof(int) , 1 , save_file);
        fwrite( &r_ramses[0] , sizeof(float) , r_ramses.size() , save_file);
        for(int i(0) ; i < Nmax ; i++){
            fwrite( &_f[i][0] , sizeof(float) , _f[i].size() , save_file);
            fwrite( &_v[i][0] , sizeof(float) , _v[i].size() , save_file);
        }
        fclose(save_file);
    }
}

float Simu::getR1DensityCoarseGrid(vector<float> & f, vector<float> & r,bool isOverDensity)
{
    vector<float> delta(f.size());
    for(unsigned int i(0) ; i < f.size()-1 ; i++)
    {
        delta[i] = (f[i+1]*pow(r[i+1],3.) - f[i]*pow(r[i],3.))/(3.*r[i]*r[i]*(r[i+1]-r[i]));
    }    
    delta[f.size() - 1] = delta[f.size() - 2];
    return getR1MassCoarseGrid(delta,r,isOverDensity);
}

float Simu::getR1MassCoarseGrid(vector<float> & f_profile, vector<float> & r,bool isOverDensity)
{
    if(f_profile[0] > 1.0)
    {
        if(isOverDensity){
            for(unsigned int i(0) ; i < f_profile.size() ; i++)
            {
                if(f_profile[i] <= 1.0 and i > 0)
                {
                    return r[i - 1] + (1.0 - f_profile[i-1])*(r[i] - r[i-1])/(f_profile[i] - f_profile[i-1]);
                }
            }
        }
        else    
            return -1.0;
    }
    else
    {
        if(!isOverDensity){
            for(unsigned int i(0) ; i < f_profile.size() ; i++)
            {
                if(f_profile[i] >= 1.0 and i > 0)
                {
                    return r[i - 1] + (1.0 - f_profile[i-1])*(r[i] - r[i-1])/(f_profile[i] - f_profile[i-1]);
                }
            }
        }
        else
            return -1.0;
    }
}
