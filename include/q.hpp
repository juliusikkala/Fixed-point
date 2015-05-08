/*
The MIT License (MIT)

Copyright (c) 2015 Julius Ikkala

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include <cmath>
#include <cstdio>
#include "internal.hpp"
#ifndef FIXED_POINT_Q_H_
    #error Please do not include this header directly; include q.h
#endif

template<unsigned f, typename I>
fp::q<f, I>::q(double d)
{
    int exp=0;
    double tmp=frexp(d, &exp);
    i=lroundl(ldexp(tmp,exp+f));
}
template<unsigned f, typename I>
fp::q<f, I>::q(long double d)
{
    int exp=0;
    long double tmp=frexp(d, &exp);
    i=llroundl(ldexp(tmp,exp+f));
}
template<unsigned f, typename I>
fp::q<f, I>::q(int u)
{
    i=fixed_point_internal::signed_lsh<f>((I)u);
}

template<unsigned f, typename I>
fp::q<f, I>::operator long double() const
{
    return ldexp((long double)i, -(int)f);
}

template<unsigned f, typename I>
fp::q<f, I> fp::q<f, I>::operator + (q<f, I> b) const
{
    q<f, I> t;
    t.i=i+b.i;
    return t;
}
template<unsigned f, typename I>
fp::q<f, I> fp::q<f, I>::operator - (q<f, I> b) const
{
    q<f, I> t;
    t.i=i-b.i;
    return t;
}
template<unsigned f, typename I>
template<unsigned fb>
fp::q<f, I> fp::q<f, I>::operator * (q<fb, I> b) const
{
    q<f, I> t;
    t.i=fixed_point_internal::no_overflow_mul_rsh<fb>(i, b.i);
    return t;
}
template<unsigned f, typename I>
fp::q<f, I> fp::q<f, I>::operator / (q<f, I> b) const
{
    static const I msb=(I)1<<(sizeof(I)*8-1);
    I abs_b=b.i<0?-b.i:b.i;
    unsigned lz=fixed_point_internal::clz(abs_b);
    I d=(abs_b<<lz);/*0.5-1.0*/
    q<sizeof(I)*8+1, typename std::make_unsigned<I>::type> e;/*Must always be below 0.5*/
    e.i=(~d+1)<<1;
    q<sizeof(I)*8-1, typename std::make_unsigned<I>::type> r(1);
    for(unsigned i=0;i<sizeof(I)-1;++i)
    {
        r=r+r*e;
        e=e*e;
    }
    q<f, I> t;
    t.i=fixed_point_internal::no_overflow_mul_rsh<sizeof(I)*8-f-1>(
        r.i>>(sizeof(I)*8-lz-(d==msb)),
        (typename std::make_unsigned<I>::type)(this->i<0?-this->i:this->i)
    );
    t.i=(b.i^this->i)&msb?-t.i:t.i;
    return t;
}
template<unsigned f, typename I>
fp::q<f, I> fp::q<f, I>::operator - () const
{
    q<f, I> t;
    t.i=-i;
    return t;
}
template<unsigned f, typename I>
bool fp::q<f, I>::operator >= (q<f, I> b) const{return i>=b.i;}
template<unsigned f, typename I>
bool fp::q<f, I>::operator > (q<f, I> b) const{return i>b.i;}
template<unsigned f, typename I>
bool fp::q<f, I>::operator <= (q<f, I> b) const{return i<=b.i;}
template<unsigned f, typename I>
bool fp::q<f, I>::operator < (q<f, I> b) const{return i<b.i;}
template<unsigned f, typename I>
bool fp::q<f, I>::operator == (q<f, I> b) const{return i==b.i;}
template<unsigned f, typename I>
bool fp::q<f, I>::operator != (q<f, I> b) const{return i!=b.i;}

template<unsigned f, typename I>
fp::q<f, I> fp::abs(q<f, I> x)
{
    return (x.q>=0?x:-x);
}
