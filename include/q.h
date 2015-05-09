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

#ifndef FIXED_POINT_Q_H_
#define FIXED_POINT_Q_H_
    #include <cstdint>
    #include <type_traits>
    namespace fp
    {
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
        
        /*Returns the absolute value of x, returns the number unchanged if it's
          type is unsigned.*/
        template<unsigned f, typename I>
        q<f, I> abs(q<f, I> x);
    }
    #include "q.hpp"
#endif
