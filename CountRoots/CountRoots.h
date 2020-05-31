#ifndef ABSTRACT_ALGEBRA2_COUNTROOTS_H
#define ABSTRACT_ALGEBRA2_COUNTROOTS_H

#include <iostream>
#include "GaloisField/Polynomial.h"
#include "utils.h"
#include "string.h"

unsigned int CountRoots(gf::Polynomial<unsigned int> p)
{
    std::string s = p.toString();
    std::map<size_t, unsigned int> pol = toPolynomial<unsigned int>(s,'x');
    auto itr = pol.begin();

    unsigned int arr[p.q()-1][p.q()-1];
    for(auto i = 0; i < p.q()-1; i++)
        if (itr->first == i && itr!=pol.end())
        {
            arr[0][i] = itr->second;
            itr++;
        }
        else
            arr[0][i] = 0;
    for(auto i = 1; i < p.q()-1; i++)
    {
        for(auto j = 0; j < p.q()-2; j++)
            arr[i][j] = arr[i-1][j+1];
        arr[i][p.q()-2] = arr[i-1][0];
    }


    /*for(int i = 0; i < p.q()-1; i++)
    {
        for(int j = 0; j < p.q()-1; j++)
            std::cout<<arr[i][j]<<" ";
        std::cout<<std::endl;
    }*/

    unsigned int q = p.q();
    unsigned int n = p.q()-1;
    unsigned int r = 0;

    for(auto i = 0; i < n; i++)
    {
        for(auto j = i; j < n; j++)
            if (arr[j][i]>0)
            {
                r++;
                for(auto temp = i; temp < n; temp++)
                    std::swap(arr[i][temp],arr[j][temp]);
                break;
            }

        for(auto j = i + 1; j < n; j++)
            if (arr[j][i])
            while(arr[j][i] != 0)
                for(auto temp = i; temp < n; temp++)
                    arr[j][temp] = (arr[j][temp] + arr[i][temp]) % q;
    }

    /*std::cout<<std::endl;
    for(int i = 0; i < p.q()-1; i++)
    {
        for(int j = 0; j < p.q()-1; j++)
            std::cout<<arr[i][j]<<" ";
        std::cout<<std::endl;
    }*/
    return q - 1 - r;
}
#endif //ABSTRACT_ALGEBRA2_COUNTROOTS_H
