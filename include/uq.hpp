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
#ifndef _FIXED_POINT_UQ_H_
    #error Please do not include this header directly; include uq.h
#endif
/*==============================================================================
 *   Helper functions
 *=============================================================================/

/*Ridiculously long name to avoid namespace clashes*/
namespace fixed_point_internal
{
    template<typename T>
    unsigned clz(T x)
    {
        unsigned i=0;
        while(!(x&(((T)1)<<(sizeof(T)*8-1)))&&(x<<=1)) ++i;
        return i;
    }
    template<int shift, bool rsh, bool zero, typename T>
    struct shifter
    {
        static T op(T x);
    };
    template<int shift, typename T>
    struct shifter<shift, true, false, T>
    {
        static T op(T x){return x>>shift;};
    };
    template<int shift, typename T>
    struct shifter<shift, false, false, T>
    {
        static T op(T x){return x<<-shift;};
    };
    template<int shift, bool rsh, typename T>
    struct shifter<shift, rsh, true, T>
    {
        static T op(T x){return 0;};
    };
    template<int shift, typename T>
    T signed_rsh(T x)
    {
        return shifter<shift, (shift>=0), (std::abs(shift)>=sizeof(T)*8), T>::op(x);
    }
    template<int shift, typename T>
    T signed_lsh(T x)
    {
        return signed_rsh<-shift, T>(x);
    }
    template<int shift, typename T>
    T no_overflow_mul_rsh(T a, T b)
    {
        static const int bits=sizeof(T)*8;
        static const T lowmask=(((T)1)<<(bits/2))-1;
        T a1=a>>bits/2;
        T a2=a&lowmask;
        T b1=b>>bits/2;
        T b2=b&lowmask;
        return signed_rsh<shift-bits>(a1*b1)+
               signed_rsh<shift-bits/2>(a1*b2)+
               signed_rsh<shift-bits/2>(a2*b1)+
               signed_rsh<shift>(a2*b2);
    }
}
template<unsigned f, typename I>
fp::uq<f, I>::uq(double d)
{
    int exp=0;
    double tmp=frexp(d, &exp);
    q=lroundl(ldexp(tmp,exp+f));
}
template<unsigned f, typename I>
fp::uq<f, I>::uq(long double d)
{
    int exp=0;
    long double tmp=frexp(d, &exp);
    q=llroundl(ldexp(tmp,exp+f));
}
template<unsigned f, typename I>
fp::uq<f, I>::uq(long unsigned u)
{
    q=fixed_point_internal::signed_lsh<f>((I)u);
}
template<unsigned f, typename I>
fp::uq<f, I>::uq(int u)
{
    q=fixed_point_internal::signed_lsh<f>((I)u);
}

template<unsigned f, typename I>
fp::uq<f, I>::operator long double() const
{
    return ldexp((long double)q, -(int)f);
}

template<unsigned f, typename I>
fp::uq<f, I> fp::uq<f, I>::operator + (uq<f, I> b) const
{
    uq<f, I> t;
    t.q=q+b.q;
    return t;
}
template<unsigned f, typename I>
fp::uq<f, I> fp::uq<f, I>::operator - (uq<f, I> b) const
{
    uq<f, I> t;
    t.q=q-b.q;
    return t;
}
template<unsigned f, typename I>
template<unsigned fb>
fp::uq<f, I> fp::uq<f, I>::operator * (uq<fb, I> b) const
{
    uq<f, I> t;
    t.q=fixed_point_internal::no_overflow_mul_rsh<fb>(q, b.q);
    return t;
}
template<unsigned f, typename I>
fp::uq<f, I> fp::uq<f, I>::operator / (uq<f, I> b) const
{
    unsigned lz=fixed_point_internal::clz(b.q);
    I d=(b.q<<lz);//0x80000000-0xFFFFFFFF
    uq<sizeof(I)*8+1, I> e;/*Must always be below 0.5*/
    e.q=(~d+1)<<1;
    uq<sizeof(I)*8-1, I> q(1);
    for(unsigned i=0;i<7;++i)/*TODO: fix loop length, a formula exists*/
    {
        q=q+q*e;
        e=e*e;
    }
    uq<f, I> r;
    r.q=fixed_point_internal::no_overflow_mul_rsh<sizeof(I)*8-f-1>(q.q>>(sizeof(I)*8-lz-(d==((I)1<<(sizeof(I)*8-1)))), this->q);
    return r;
}
