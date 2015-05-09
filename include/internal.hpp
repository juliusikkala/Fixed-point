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
#ifndef FIXED_POINT_INTERNAL_H_
#define FIXED_POINT_INTERNAL_H_
    #include <type_traits>
    #include <cmath>
    /*Ridiculously long name to avoid namespace clashes*/
    namespace fp_internal
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
            return shifter<
                shift, (shift>=0), (std::abs(shift)>=sizeof(T)*8), T
            >::op(x);
        }
        template<typename T>
        T signed_rsh(T x, int shift)
        {
            return std::abs(shift)<sizeof(T)*8?(shift<0?x<<-shift:x>>shift):0;
        }
        template<int shift, typename T>
        T signed_lsh(T x)
        {
            return signed_rsh<-shift, T>(x);
        }
        template<typename T>
        T signed_lsh(T x, int shift)
        {
            return std::abs(shift)<sizeof(T)*8?(shift<0?x>>-shift:x<<shift):0;
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
        template<typename T>
        T no_overflow_mul_rsh(T a, T b, int shift)
        {
            static const int bits=sizeof(T)*8;
            static const T lowmask=(((T)1)<<(bits/2))-1;
            T a1=a>>bits/2;
            T a2=a&lowmask;
            T b1=b>>bits/2;
            T b2=b&lowmask;
            return signed_rsh(a1*b1, shift-bits)+
                   signed_rsh(a1*b2, shift-bits/2)+
                   signed_rsh(a2*b1, shift-bits/2)+
                   signed_rsh(a2*b2, shift);
        }
    }
#endif
