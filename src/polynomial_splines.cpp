#include <polynomial_splines.h>


int cutoff_fact(int const& n,
                int const& k)
{
    int prod = int(1.0);
    for(int i = int(0); i < k; i+=int(1.0))
    {
        prod = prod*(n-i);
    }
    return prod;
};