#ifndef LOGADD_H_INCLUDED
#define LOGADD_H_INCLUDED

#define LZERO  (-1.0E10) // log(0)
#define LSMALL (-0.5E10) // log values < LSMALL are set to LZERO
#define minLogExp -log(-LZERO) // ~=-23

double LogAdd(double x, double y) {
    double temp, diff, z;
    if (x < y)
    {
        temp = x; x = y; y = temp;
    }
    diff = y-x; // notice that diff <= 0
    if (diff < minLogExp)   // if y' is far smaller than x'
        return (x < LSMALL) ? LZERO : x;
    else
    {
        z = exp(diff);
        return x + log(1.0 + z);
    }
}

#endif // LOGADD_H_INCLUDED
