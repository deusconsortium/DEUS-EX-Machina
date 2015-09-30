#include "Tools.h"

using namespace std;

Tools::Tools(){
    srand(time(NULL));
}

string Tools::IntToString(int const number,bool formated,int nb_chiffres){
    string retour,chiffre;
    ostringstream oss;
    oss<<number;

    chiffre=oss.str();
    
    if(formated){
        retour="_";
        for(int i(0) ; i<nb_chiffres-chiffre.size() ; i++ )
                retour+="0";
        retour+=chiffre;
        return retour;
    }
    
    return chiffre;
}

int Tools::StringToInt(const string Snumber){
    std::istringstream iss(Snumber);
    
    int nombre;
    iss >> nombre; 
    
    return nombre;
}

float Tools::Random(float min, float max){
    return ( rand()/(double)RAND_MAX ) * (max-min) + min;
}

int Tools::getMinimumIndex(const std::vector<float> & x, const std::vector<float> & y){
    if(x.size() != y .size())
        return -1;
    
    int index = 0;
    for(int i(0); i< x.size() ; i++)
    {
        if(y[i] < y[index])
            index = i;
    }
 
    return index;
}



float Tools::getXvalue(float const y0, const vector<float> & x, const vector<float> & y)
{
    if(x.size() != y .size())
        return -1.0;
    
    int i = 0;
    if(y[0] > y0)
    {
        while(i < x.size() && y[i] > y0)
        {
            i++;
        }
    }
    else
    {
        while(i < x.size() && y[i] < y0)
        {
            i++;
        }
    }
    
    if(i!= x.size() -1)
        return ((x[i+1]-x[i])*y0 - y[i]*x[i+1] + y[i+1]*x[i])/(y[i+1]-y[i]);
    
    return -1.0;   
}

void Tools::getMeanAndSigma(float & mean, float & sigma,const vector<float> tab)
{
    mean = 0.0;
    sigma = 0.0;
    
    for( int k(0) ; k< tab.size() ; k++)
        mean+= tab[k];
    
    if(tab.size() > 0)
        mean = mean/(tab.size());
    
    for(int k(0) ; k < tab.size() ; k++)
        sigma+= (tab[k] - mean)*(tab[k] - mean);
    
    if(tab.size() > 1)
        sigma = sqrt(sigma/(tab.size() - 1));
}

int Tools::getNbcolumn(string const file_path)
{
    ifstream file(file_path.c_str());
    if(file)
    {
        string ligne; 
        //on lit la deuxième ligne
        std::getline(file,ligne);
        std::getline(file,ligne);
        
        istringstream iss(ligne);
        int c= 0;
        
        do
        {
            string sub;
            iss >> sub;
            if(sub != "" && sub != "\n")
                c++;
        } while (iss);
        
        file.close();
        return c;
    }
    return 0;
}

int Tools::getNbLines(string const file_path)
{
    ifstream file(file_path.c_str());
    if(file)
    {
        string ligne; //Création d'une chaine de caractere
        int nbLignes = 0;
        while(std::getline(file, ligne)) //On lit chaque ligne du fichier que l'on stoke dans "ligne"
              nbLignes++;
        
        file.close(); //On ferme le fichier
        return nbLignes;
    }
    else
    {
        cout << "ERROR ! File doesn't exists  for getNbLines(" << file_path <<")" <<endl;
        return 0;
    }
}

bool Tools::fileExists(string const file_path)
{
    ifstream file(file_path.c_str());
    if(file)
    {
        file.close();
        return true;
    }
    return false;
}

bool Tools::createFolder(std::string const folder_name)
{
    const int dir_err = mkdir(folder_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (dir_err == -1)
        return false;
    return true;
}

int Tools::getLastOutput(const int boxlen, const int npart, const std::string cosmo)
{
    ifstream output_file( "data/lastOutput.txt" ); 
    if(output_file){
        string ligne;
        vector<int> fboxlen;
        vector<int> fnpart;
        vector<int> lout;
        vector<int> rout;
        int b,n,outl,outr;
        while(getline(output_file,ligne))
        {
            output_file >> b >> n >> outl >> outr;
            fboxlen.push_back(b);
            fnpart.push_back(n);
            lout.push_back(outl);
            rout.push_back(outr);
        }
        output_file.close();
        string c = cosmo.substr(0,1);
        
        if(c=="l"){
            for(unsigned int i(0) ; i <fboxlen.size() ; i++)
                if(fboxlen[i]==boxlen && fnpart[i]==npart)
                    return lout[i];
        }
        else if(c=="r"){
            for(unsigned int i(0) ; i <fboxlen.size() ; i++)
                if(fboxlen[i]==boxlen && fnpart[i]==npart)
                    return rout[i];
        }
    }
    else
        cout << "ERROR : No Output file !!" << endl;
    return -1;
}