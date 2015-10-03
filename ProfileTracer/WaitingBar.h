/* 
 * File:   WaitingBar.h
 * Author: pdefromont
 *
 * Created on 9 septembre 2015, 10:30
 */

#ifndef WAITINGBAR_H
#define	WAITINGBAR_H

using namespace std;

class WaitingBar
{
    public:
    WaitingBar()
    {
        _Lenght = 0;
        _N = 0;
        _k = 0;
        _n = 0;
        _dn = 0;
    }
    WaitingBar(int Nmax,int lenght = 50)
    {
        _Lenght = lenght;
        _N = Nmax;
        _k = 1;
        _n = 0;
        _dn = _N/_Lenght;
        for(unsigned int i(0) ; i < _Lenght ; i++)
            cout << "_";
        cout<<endl;
    }
    ~WaitingBar(){}
    
    void Reset(int Nmax,int lenght = 50)
    {
        _Lenght = lenght;
        _N = Nmax;
        _k = 1;
        _n = 0;
        _dn = _N/_Lenght;
        for(unsigned int i(0) ; i < _Lenght ; i++)
            cout << "_";
        cout<<endl;
    }
    void Up()
    {
        _n++;
        if(_n >= _k*_dn)
        {
            cout << "*";
            cout.flush();
            _k++;
        }
    }
    private:
        
    int _N,_Lenght,_k,_n;
    float _dn;
};

#endif	/* WAITINGBAR_H */

