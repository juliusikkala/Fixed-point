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
fp::q<f, I>::q(int u)
{
    i=fp_internal::signed_lsh<f>((I)u);
}
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
template<unsigned fb>
fp::q<f, I>::q(q<fb, I> a)
{
    i=fp_internal::signed_rsh<fb-f>(a.i);
}
template<unsigned f, typename I>
fp::q<f, I>::operator long double() const
{
    return ldexp((long double)i, -(int)f);
}

template<unsigned f, typename I>
template<unsigned fb>
fp::q<f, I> fp::q<f, I>::operator + (q<fb, I> b) const
{
    q<f, I> t;
    t.i=i+fp_internal::signed_rsh<fb-f>(b.i);
    return t;
}
template<unsigned f, typename I>
template<unsigned fb>
fp::q<f, I> fp::q<f, I>::operator - (q<fb, I> b) const
{
    q<f, I> t;
    t.i=i-fp_internal::signed_rsh<fb-f>(b.i);
    return t;
}
template<unsigned f, typename I>
template<unsigned fb>
fp::q<f, I> fp::q<f, I>::operator * (q<fb, I> b) const
{
    q<f, I> t;
    t.i=fp_internal::no_overflow_mul_rsh<fb>(i, b.i);
    return t;
}
template<unsigned f, typename I>
template<unsigned fb>
fp::q<f, I> fp::q<f, I>::operator / (q<fb, I> b) const
{
    static const I msb=(I)1<<(sizeof(I)*8-1);//Most significant bit for the type
    //Make b positive so that leading zeroes can be properly computed
    I abs_b=b.i<0?-b.i:b.i;
    unsigned lz=fp_internal::clz(abs_b);//Amount of leading zeroes
    //normalize b to [0.5, 1.0[, where all digits are after radix
    I d=(abs_b<<lz);
    q<sizeof(I)*8+1, typename std::make_unsigned<I>::type> e;
    e.i=(~d+1)<<1;//[0, 0.5[
    //r is the reciprocal of d
    q<sizeof(I)*8-1, typename std::make_unsigned<I>::type> r(1);
    for(unsigned i=0;i<sizeof(I)-1;++i)
    {
        r=r+r*e;
        e=e*e;
    }
    q<f, I> t;
    t.i=fp_internal::no_overflow_mul_rsh(//adjust the radix point of (this*r)
        r.i,
        (typename std::make_unsigned<I>::type)(this->i<0?-this->i:this->i),
        sizeof(i)*16-fb-lz-(d==msb)-1
    );
    t.i=(b.i^this->i)&msb?-t.i:t.i;//set correct sign
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
template<unsigned fb>
fp::q<f, I> fp::q<f, I>::operator % (q<fb, I> b) const
{
    q<f, I> t;
    t.i=i%fp_internal::signed_rsh<fb-f>(b.i);
    return t;
}

template<unsigned f, typename I>
template<unsigned fb>
bool fp::q<f, I>::operator >= (q<fb, I> b) const
{
    static const unsigned bits=sizeof(I)*8;
    I ma=fp_internal::signed_rsh<f>(i);
    I mb=fp_internal::signed_rsh<fb>(b.i);
    return (ma>mb)||
           ((ma==mb)&&
           ((typename std::make_unsigned<I>::type)
            fp_internal::signed_lsh<bits-f>(i)
            >=(typename std::make_unsigned<I>::type)
            fp_internal::signed_lsh<bits-fb>(b.i)));
}
template<unsigned f, typename I>
template<unsigned fb>
bool fp::q<f, I>::operator > (q<fb, I> b) const
{
    static const unsigned bits=sizeof(I)*8;
    I ma=fp_internal::signed_rsh<f>(i);
    I mb=fp_internal::signed_rsh<fb>(b.i);
    return (ma>mb)||
           ((ma==mb)&&
           ((typename std::make_unsigned<I>::type)
            fp_internal::signed_lsh<bits-f>(i)
            >(typename std::make_unsigned<I>::type)
            fp_internal::signed_lsh<bits-fb>(b.i)));
}
template<unsigned f, typename I>
template<unsigned fb>
bool fp::q<f, I>::operator <= (q<fb, I> b) const
{
    static const unsigned bits=sizeof(I)*8;
    I ma=fp_internal::signed_rsh<f>(i);
    I mb=fp_internal::signed_rsh<fb>(b.i);
    return (ma<mb)||
           ((ma==mb)&&
           ((typename std::make_unsigned<I>::type)
            fp_internal::signed_lsh<bits-f>(i)
            <=(typename std::make_unsigned<I>::type)
            fp_internal::signed_lsh<bits-fb>(b.i)));
}
template<unsigned f, typename I>
template<unsigned fb>
bool fp::q<f, I>::operator < (q<fb, I> b) const
{
    static const unsigned bits=sizeof(I)*8;
    I ma=fp_internal::signed_rsh<f>(i);
    I mb=fp_internal::signed_rsh<fb>(b.i);
    return (ma<mb)||
           ((ma==mb)&&
           ((typename std::make_unsigned<I>::type)
            fp_internal::signed_lsh<bits-f>(i)
            <(typename std::make_unsigned<I>::type)
            fp_internal::signed_lsh<bits-fb>(b.i)));
}
template<unsigned f, typename I>
template<unsigned fb>
bool fp::q<f, I>::operator == (q<fb, I> b) const
{
    return fp_internal::signed_rsh<f-fb>(i)==b.i&&
           fp_internal::signed_rsh<fb-f>(b.i)==i;
}
template<unsigned f, typename I>
template<unsigned fb>
bool fp::q<f, I>::operator != (q<fb, I> b) const
{
    return fp_internal::signed_rsh<f-fb>(i)!=b.i||
           fp_internal::signed_rsh<fb-f>(b.i)!=i;
}

template<unsigned f, typename I>
fp::q<f, I> fp::abs(q<f, I> x)
{
    return (x.q>=0?x:-x);
}
