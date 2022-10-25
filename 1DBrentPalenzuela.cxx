#include "con2primFactory.hxx"
#include <assert.h>

bool tol(double a, double b)
{
    using std::abs;
    double rel_error;
    if (a == 0)
    {
        rel_error = abs(b);
    }
    else if (b == 0)
    {
        rel_error = abs(a);
    }
    else
    {
        rel_error = 2.0 * abs(a - b) / (abs(a) + abs(b));
    }
    if (rel_error < 1.E-10)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// See <https://en.wikipedia.org/wiki/Brent%27s_method>
double brent(double (*f)(double x, con2primFactory &plasma), con2primFactory &plasma)
{
    using std::abs;
    using std::max;
    using std::min;
    using std::swap;

    double qPalenzuela = plasma.ConservedVars[TAU] / plasma.ConservedVars[D];
    double rPalenzuela = plasma.Ssq / pow(plasma.ConservedVars[D], 2);
    double sPalenzuela = plasma.Bsq / plasma.ConservedVars[D];
    double tPalenzuela = plasma.BiSi / pow(plasma.ConservedVars[D], 3. / 2.);

    double xPalenzuela_lowerBound = 1.0 + qPalenzuela - sPalenzuela;
    double xPalenzuela_upperBound = 2.0 + 2.0 * qPalenzuela - sPalenzuela;

    double a = xPalenzuela_lowerBound;
    double b = xPalenzuela_upperBound;
    double c = xPalenzuela_upperBound;

    int iters = 0;
    double fa = f(a, plasma);
    double fb = f(b, plasma);

    if (abs(fa) < abs(fb))
    {
        swap(a, b);
        swap(fa, fb);
    }
    if (fb == 0)
        return b;
    if (fa * fb >= 0)
    {
        // Root is not bracketed
        iters = plasma.max_iterations;
        if (abs(fa) < abs(fb))
        {
            return a;
        }
        else
        {
            return b;
        }
    }
    c = a;
    double fc = fa;
    bool mflag = true;
    double d;

    while (fb != 0 && !tol(a, b) && iters < plasma.max_iterations)
    {
        double s;
        if (fa != fc && fb != fc)
            // inverse quadratic interpolation
            s = (a * fb * fc) / ((fa - fb) * (fa - fc)) +
                (b * fa * fc) / ((fb - fa) * (fb - fc)) +
                (c * fa * fb) / ((fc - fa) * (fc - fb));
        else
            // secant method
            s = (a + b) / 2 - (fa + fb) / 2 * (b - a) / (fb - fa);
        double u = (3 * a + b) / 4;
        double v = b;
        if (u > v)
            swap(u, v);
        bool cond1 = !(u <= s && s <= v);
        bool cond2 = mflag && abs(s - b) >= abs(b - c) / 2;
        bool cond3 = !mflag && abs(s - b) >= abs(c - d) / 2;
        bool cond4 = mflag && tol(c, b);
        bool cond5 = !mflag && tol(c, d);
        if (cond1 || cond2 || cond3 || cond4 || cond5)
        {
            // bisection
            s = (a + b) / 2;
            mflag = true;
        }
        else
        {
            mflag = false;
        }
        double fs = f(s, plasma);
        // `d` is assigned for the first time here; it won't be used above on the
        // first iteration because `mflag` is set
        d = c;
        c = b;
        fc = fb;
        if (fa * fs < 0)
        {
            b = s;
            fb = fs;
        }
        else
        {
            a = s;
            fa = fs;
        }
        // CCTK_VINFO("iters=%d mflag=%d   a=%.17g b=%.17g c=%.17g d=%.17g fa=%.17g"
        //            "fb=%.17g fc=%.17g",
        //            iters, int(mflag), double(a), double(b), double(c), double(d),
        //            double(fa), double(fb), double(fc));
        assert(fa * fb <= 0);
        if (abs(fa) < abs(fb))
        {
            swap(a, b);
            swap(fa, fb);
        }
        ++iters;
    }
    std::cout << iters << std::endl;

    if (abs(fa) < abs(fb))
    {
        return a;
    }
    else
    {
        return b;
    }
}

double funcRoot_1DPalenzuela(double x, con2primFactory &plasma)
{

    // computes f(x) from x and q,r,s,t
    const double qPalenzuela = plasma.ConservedVars[TAU] / plasma.ConservedVars[D];
    const double rPalenzuela = plasma.Ssq / pow(plasma.ConservedVars[D], 2);
    const double sPalenzuela = plasma.Bsq / plasma.ConservedVars[D];
    const double tPalenzuela = plasma.BiSi / pow(plasma.ConservedVars[D], 3. / 2.);

    // (i)
    double Wminus2 = 1.0 - (x * x * rPalenzuela + (2 * x + sPalenzuela) * tPalenzuela * tPalenzuela) / (x * x * (x + sPalenzuela) * (x + sPalenzuela));

    Wminus2 = fmin(fmax(Wminus2, 1e-10), 1 - 1e-10);
    const double W_loc = pow(Wminus2, -0.5);

    // (ii)
    double rho_loc = plasma.ConservedVars[D] / W_loc;

    // (iii)
    double eps_loc = W_loc - 1.0 + (1.0 - W_loc * W_loc) * x / W_loc + W_loc * (qPalenzuela - sPalenzuela + tPalenzuela * tPalenzuela / (2 * x * x) + sPalenzuela / (2 * W_loc * W_loc));

    // (iv)
    double P_loc = plasma.get_Press_funcRhoEps(rho_loc, eps_loc);

    return (x - (1.0 + eps_loc + P_loc / rho_loc) * W_loc);
}

/***************************************************************************
1DBrentPalenzuela C2P
------------------------------------
1D-Brent Palenzuela scheme for c2p.
Sources: Nielsen+2014, Palenzuela+2015, Section 3.4.1 of Siegel+2018,
https://en.wikipedia.org/wiki/Brent%27s_method
****************************************************************************/
void Con2Prim_1DBrentPalenzuela(con2primFactory &plasma)
{ // Send con2primFactory object as reference to modify it, and because we can not instantiate abstract class

    plasma.Failed_1DBrentPalenzuela = 1;

    // find x, this is the recovery process
    plasma.xPalenzuela_Sol = brent(*funcRoot_1DPalenzuela, plasma);
    std::cout << plasma.xPalenzuela_Sol << "\n";
    plasma.xPalenzuelaToPrim();
    return;
}