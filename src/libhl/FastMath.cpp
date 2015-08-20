#ifdef HUGIN_LITE

#include <math.h>

#include "FastMath.h"

namespace libhl {

/* atan2 approximation by Hastings that has |err| < 0.0005 */
double hl_atan2(double y, double x)
{
    const static double M_PI_BY2 = M_PI / 2.0;

    if (x == 0.0)
    {
        if (y > 0.0) return M_PI_BY2;
        if (y == 0.0) return 0.0;
        return -M_PI_BY2;
    }

    double atan;
    double z = y / x;
    if (fabs(z) < 1.0)
    {
        atan = z / (1.0 + 0.28 * z * z);
        if (x < 0.0)
        {
            if (y < 0.0) return atan - M_PI;
            return atan + M_PI;
        }
    }
    else
    {
        atan = M_PI_BY2 - z / (z * z + 0.28);
        if (y < 0.0) return atan - M_PI;
    }

    return atan;
}

/* exp ~= (1 + x / n)^n for large n, we choose n = 1024 */
double hl_exp(double x)
{
    x = 1.0 + x / 1024.0;
    for(int i = 0;i < 10;i++)
    {
        x *= x;
    }

    return x;
}

} // namespace

#endif
