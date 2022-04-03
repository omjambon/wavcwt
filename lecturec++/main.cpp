#include <iostream>
#include <fstream> // Fichier d'entete pour la gestion des fichiers
#include <vector>
#include <complex>
#include<cmath>
#include <algorithm>

using namespace std;

std::vector<double> getDefaultScales(int n, double ds)
{
    int s0 = 2; // Plus petit scale utile
    double max_scale = n /(sqrt(2)*s0); //Determine longest useful scale for wavelet
    if (max_scale <= 1)
    {
        max_scale = n / 2;
    }
    double nv = 1 / ds;
    max_scale = floor(nv * log2(max_scale));
    double a0 = pow(2,ds);
    std::vector<double> scales;
    for (int value = 0; value < max_scale + 1; value ++)
    {
        scales.push_back(s0*pow(a0,value));
    }
    return scales;
}

int main()
{
    std::ifstream ifs("pasnoisy.iq", std::ios::binary | std::ios::in);
    std::vector<std::complex<double>> v(320);
    ifs.read(reinterpret_cast<char*>(v.data()), 500*sizeof(double));
    ifs.close();
    for (int i=0; i<v.size();i++)
    {
       //std::cout << i << v[i] << std::endl;
    }

    std::string binFileName = "sortiec++.bin";
    std::ofstream rawFile(binFileName, std::ios::binary);
    rawFile.write(reinterpret_cast<const char*>(v.data()), 500*sizeof(double));

    int n_orig = v.size();
    int nv = 10;
    double ds = 1.0 / nv;
    int fs = 40e6;
    double dt = 1.0 / fs;
    int padvalue = n_orig/2; //TODO : PadValue mettre en double arrondi

    std::vector<std::complex<double>> x;
    std::vector<std::complex<double>> range1( &v[0], &v[0]+padvalue );
    std::vector<std::complex<double>> range2( &v[0]+padvalue , &v[0]+v.size() );
    std::reverse(range1.begin(), range1.end());
    std::reverse(range2.begin(), range2.end());
    range1.insert(range1.end(), v.begin(), v.end()); //Concatenate range et v
    range1.insert(range1.end(), range2.begin(), range2.end()); //Concatenate range2 et v
    x = range1;
    for (int i=0; i<x.size();i++)
    {
       std::cout << x[i] << std::endl;
    }



    int n = x.size();
    std::cout << range1.size() <<" ds "<< ds <<" padvalue "<<padvalue<< std::endl;
    std::vector<double> wavscales = getDefaultScales(n_orig,ds);

    for (int i=0; i<wavscales.size();i++)
    {
    std::cout<<" wavscales : "<<wavscales[i]<< std::endl;
    }
    int num_scales = wavscales.size();
    std::cout<<" num_scales : "<<num_scales<< std::endl;

    std::vector<double> omega;
    for (int value = 1; value < floor(n/2) + 1; value++)
    {
        omega.push_back(value*(2 * M_PI) / n);

    }

        for (int value = 0; value < omega.size(); value++)
    {
        std::cout<<" omega : "<<omega[value]<< std::endl;
    }


    cout << "Hello world!" << endl;
    return 0;
}


