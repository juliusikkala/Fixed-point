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

#ifndef _FIXED_POINT_UQ_H_
#define _FIXED_POINT_UQ_H_
    #include <cstdint>
    namespace fp
    {
        template<unsigned f, typename I>
        struct uq
        {
            I q;
            
            uq(double d);
            uq(long double d);
            uq(long unsigned u);
            uq(int u=0);/*Integer literals are int*/
            operator long double() const;
            uq<f, I> operator + (uq<f, I> b) const;
            uq<f, I> operator - (uq<f, I> b) const;
            template<unsigned fb>
            uq<f, I> operator * (uq<fb, I> b) const;
            uq<f, I> operator / (uq<f, I> b) const;
        };
        typedef uq<16, uint32_t> UQ16_16;
        typedef uq<32, uint64_t> UQ32_32;
        typedef uq<16, uint64_t> UQ48_16;
    }
    #include "uq.hpp"
#endif
