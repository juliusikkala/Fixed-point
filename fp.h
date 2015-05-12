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
#ifndef FIXED_POINT_H_
#define FIXED_POINT_H_
    #include <type_traits>
    #include <cmath>
    #include <cstdint>
    //==========================================================================
    //Fixed-point numbers
    //==========================================================================
    namespace fp
    {
        /*The fixed-point number class. f is the amount of fractional bits, I is
          the internal type. If I is signed, q will be signed.*/
        template<unsigned f, typename I>
        struct q
        {
            I i;
            
            q(int u=0);/*Integer literals are int*/
            q(double d);
            q(long double d);
            template<unsigned fb>
            q(q<fb, I> a);
            
            operator long double() const;
            template<unsigned fb>
            q<f, I> operator + (q<fb, I> b) const;
            template<unsigned fb>
            q<f, I> operator - (q<fb, I> b) const;
            template<unsigned fb>
            q<f, I> operator * (q<fb, I> b) const;
            template<unsigned fb>
            q<f, I> operator / (q<fb, I> b) const;
            
            q<f, I> operator - () const;
            template<unsigned fb>
            q<f, I> operator % (q<fb, I> b) const;
            
            template<unsigned fb>
            q<f, I> & operator = (q<fb, I> b);
            q<f, I> & operator = (int b);
            q<f, I> & operator = (double b);
            q<f, I> & operator = (long double b);
            
            template<unsigned fb>
            bool operator >= (q<fb, I> b) const;
            template<unsigned fb>
            bool operator > (q<fb, I> b) const;
            template<unsigned fb>
            bool operator <= (q<fb, I> b) const;
            template<unsigned fb>
            bool operator < (q<fb, I> b) const;
            template<unsigned fb>
            bool operator == (q<fb, I> b) const;
            template<unsigned fb>
            bool operator != (q<fb, I> b) const;
        };
        typedef q<16, uint32_t> UQ16_16;
        typedef q<32, uint64_t> UQ32_32;
        typedef q<16, uint64_t> UQ48_16;
        typedef q<16, int32_t> Q16_16;
        typedef q<32, int64_t> Q32_32;
        typedef q<16, int64_t> Q48_16;
        
        /*Returns the absolute value of x*/
        template<unsigned f, typename I>
        q<f, I> abs(q<f, I> x);
    }
    //==========================================================================
    //Internal helper functions and classes
    //==========================================================================
    namespace fp_internal
    {
        /*Count leading zeroes*/
        template<typename T>
        unsigned clz(T x)
        {
            unsigned i=0;
            while(!(x&(((T)1)<<(sizeof(T)*8-1)))&&(x<<=1)) ++i;
            return i;
        }
        /*Used to select the correct shift operator in compile-time*/
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
        /*Signed right shift. Accepts negative values of shift and never
          complains about too big shifts. Compile-time version.*/
        template<int shift, typename T>
        T signed_rsh(T x)
        {
            return shifter<
                shift, (shift>=0), (std::abs(shift)>=sizeof(T)*8), T
            >::op(x);
        }
        /*Signed right shift, run time version.*/
        template<typename T>
        T signed_rsh(T x, int shift)
        {
            return std::abs(shift)<sizeof(T)*8?(shift<0?x<<-shift:x>>shift):0;
        }
        /*Signed left shift, compile-time version*/
        template<int shift, typename T>
        T signed_lsh(T x)
        {
            return signed_rsh<-shift, T>(x);
        }
        /*Signed left shift, run time version*/
        template<typename T>
        T signed_lsh(T x, int shift)
        {
            return std::abs(shift)<sizeof(T)*8?(shift<0?x>>-shift:x<<shift):0;
        }
        /*Multiplies and simultaneously right-shifts the argument values,
          without losing the high bits. Compile-time version*/
        template<int shift, typename T>
        T mul_rsh(T a, T b)
        {
            static const T msb=(T)1<<(sizeof(T)*8-1);
            static const int bits=sizeof(T)*8;
            static const T lowmask=(((T)1)<<(bits/2))-1;
            typedef typename std::make_unsigned<T>::type U;
            U abs_a=a>=0?a:-a;
            U abs_b=b>=0?b:-b;
            U a1=abs_a>>bits/2;
            U a2=abs_a&lowmask;
            U b1=abs_b>>bits/2;
            U b2=abs_b&lowmask;
            U res=signed_rsh<shift-bits>(a1*b1)+
                  signed_rsh<shift-bits/2>(a1*b2)+
                  signed_rsh<shift-bits/2>(a2*b1)+
                  signed_rsh<shift>(a2*b2);
            return (a^b)&msb?-((T)res):(T)res;
        }
        /*Run time version*/
        template<typename T>
        T mul_rsh(T a, T b, int shift)
        {
            static const T msb=(T)1<<(sizeof(T)*8-1);
            static const int bits=sizeof(T)*8;
            static const T lowmask=(((T)1)<<(bits/2))-1;
            typedef typename std::make_unsigned<T>::type U;
            U abs_a=a>=0?a:-a;
            U abs_b=b>=0?b:-b;
            U a1=abs_a>>bits/2;
            U a2=abs_a&lowmask;
            U b1=abs_b>>bits/2;
            U b2=abs_b&lowmask;
            U res=signed_rsh(a1*b1, shift-bits)+
                  signed_rsh(a1*b2, shift-bits/2)+
                  signed_rsh(a2*b1, shift-bits/2)+
                  signed_rsh(a2*b2, shift);
            return (a^b)&msb?-((T)res):(T)res;
        }
    }
    //==========================================================================
    //Implementation
    //==========================================================================
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
        t.i=fp_internal::mul_rsh<fb>(i, b.i);
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
        t.i=fp_internal::mul_rsh(//adjust the radix point of (this*r)
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
    fp::q<f, I> & fp::q<f, I>::operator = (q<fb, I> b)
    {
        i=fp_internal::signed_rsh<fb-f>(b.i);
        return *this;
    }
    template<unsigned f, typename I>
    fp::q<f, I> & fp::q<f, I>::operator = (int b)
    {
        i=fp_internal::signed_lsh<f>((I)b);
        return *this;
    }
    template<unsigned f, typename I>
    fp::q<f, I> & fp::q<f, I>::operator = (double b)
    {
        int exp=0;
        double tmp=frexp(b, &exp);
        i=llroundl(ldexp(tmp,exp+f));
        return *this;
    }
    template<unsigned f, typename I>
    fp::q<f, I> & fp::q<f, I>::operator = (long double b)
    {
        int exp=0;
        long double tmp=frexp(b, &exp);
        i=llroundl(ldexp(tmp,exp+f));
        return *this;
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
#endif
