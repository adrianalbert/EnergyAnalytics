/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cCyclicVector.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _CCYCLICVECTOR_H_
#define _CCYCLICVECTOR_H_
#pragma once
#include <vector>
#include "Hmm.h"

template<class T>
class cCyclicVector : public std::vector<T>
{
public:
    typename std::vector<T>::reference operator[](int __n)
    {
        return std::vector<T>::at(__n % std::vector<T>::size());
    }

    cCyclicVector & operator=(const cCyclicVector &rhs)
    {
        if (this == &rhs)      // Same object?
                return *this;

        this->clear();
        for (int i=0;i<(int)rhs.size();i++)
            this->push_back(rhs.at(i));

        return *this;
    }
};
#endif //_CCYCLICVECTOR_H_
