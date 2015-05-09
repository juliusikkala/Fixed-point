# Fixed-point numbers
**This is a WIP, and doesn't even have a proper name yet.** It works already,
though.

A header only C++11 fixed-point number library.
* The number of fractional bits can be set and it does not have to be a power of
  two
* The underlying integer type can be freely chosen
* Supports both signed and unsigned numbers
* Allows for easy casting between floating-point numbers and fixed-point numbers
* Precise*

<sup>\* or at least should be. Please create an issue with the values you tested
if you notice imprecision.</sup>

Floating-point numbers are usually better than fixed-point numbers, but if you
need constant resolution, consider using fixed-point numbers.
**This library is not meant for embedded devices.** Instead, it is meant to be
used _only_ when you need constant resolution and high precision and performance
is not an issue.

##Example
```c++
#include "fp.h"
...
fp::q<17, int64_t> a(18.12);/*17 fractional bits, signed*/
fp::q<17, int64_t> b(1.2339);
/*This should print the value. Notice that it isn't exactly correct, the
  precision is limited by the number of fractional bits, and every division and
  multiplication may cause rounding errors, again, due to the number of
  fractional bits.*/
printf("%f\n", (double)((a+b)/(a-b)*b));
...
```
