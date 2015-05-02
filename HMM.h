#ifndef HMM_H_INCLUDED
#define HMM_H_INCLUDED
#include <iostream>
#include <vector>
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
    double forward(int* o, int T);
    double bakward(int* o, int T);
/// double decode (int* o, int T, int* q);
    void   learn  (int* o, int T);

    void set_initi(vector<double> &a);
    void set_trans(vector< vector<double> > &a);
    void set_state(vector< vector<double> > &a);
    void print_HMM();
private:
/// Constant
    int N, M, T, L;
/// Probability
    vector<double>                     initi_set; /// N
    vector< vector<double> >           trans_set; /// N x N
    vector< vector<double> >           state_set; /// N x M
/// Memory Set
    vector< vector<double> >           forwa_set; /// T x N
    vector< vector<double> >           backw_set; /// T x N
    vector< vector<double> >           deco1_set; /// T x N
    vector< vector<int> >              deco2_set; /// T x N
    vector< vector<double> >           lear1_set; /// T x N
    vector< vector< vector<double> > > lear2_set; /// T x N x N
};

#endif // HMM_H_INCLUDED
