#ifndef HMM_H_INCLUDED
#define HMM_H_INCLUDED

#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>

using namespace std;

class HMM {

public:
    HMM(int a, int b, int c, int d) {
        N = a; M = b; T = c; L = d;

        initi_set.resize(N, 0);
        trans_set.resize(N, vector<double>(N, 0));
        state_set.resize(N, vector<double>(M, 0));

        forwa_set.resize(T, vector<double>(N, 0));
        backw_set.resize(T, vector<double>(N, 0));
        deco1_set.resize(T, vector<double>(N, 0));
        deco2_set.resize(T, vector<int>   (N, 0));
        lear1_set.resize(T, vector<double>(N, 0));
        lear2_set.resize(T, vector< vector<double> >(N, vector<double>(N, 0)));
    }
    double forward   (int* o, int T);
    double Logforward(int* o, int T);
    double bakward   (int* o, int T);
    double Logbakward(int* o, int T);
    double decode    (int* o, int T, int* q);
    void   learn     (int* o, int T);
    void   Loglearn  (int* o, int T);
    void   clear     ();

    void   set_initi(vector<double> &a);
    void   set_trans(vector< vector<double> > &a);
    void   set_state(vector< vector<double> > &a);
    void   print_HMM();
    void   print_HMM(string);
private:

    int N, M, T, L;

    vector<double>                     initi_set; /// N
    vector< vector<double> >           trans_set; /// N x N
    vector< vector<double> >           state_set; /// N x M

    vector< vector<double> >           forwa_set; /// T x N
    vector< vector<double> >           backw_set; /// T x N
    vector< vector<double> >           deco1_set; /// T x N
    vector< vector<int> >              deco2_set; /// T x N
    vector< vector<double> >           lear1_set; /// T x N
    vector< vector< vector<double> > > lear2_set; /// T x N x N
};

#define LZERO  (-1.0E10) // log(0)
#define LSMALL (-0.5E10) // log values < LSMALL are set to LZERO
#define minLogExp -log(-LZERO) // ~=-23

double LogAdd(double x, double y);

#endif // HMM_H_INCLUDED
