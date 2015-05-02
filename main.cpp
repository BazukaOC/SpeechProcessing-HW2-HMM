#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include "HMM.h"

using namespace std;

int N = 3, M = 3, T = 80, L = 5;

double initi[3]    =  {0.333333333333333, 0.333333333333333, 0.333333333333334};
double trans[3][3] = {{0.34, 0.33, 0.33}, {0.33, 0.34, 0.33}, {0.33, 0.33, 0.34}};
double state[3][3] = {{0.34, 0.33, 0.33}, {0.33, 0.34, 0.33}, {0.33, 0.33, 0.34}};

HMM hmm[2] = HMM(N, M, T, L);

void initialize(vector<double> &initi_set,
                vector< vector<double> > &trans_set,
                vector< vector<double> > &state_set) {
    for(int i = 0; i < N; ++i)
        initi_set[i] = initi[i];

    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            trans_set[i][j] = trans[i][j];

    for(int i = 0; i < N; ++i)
        for(int j = 0; j < M; ++j)
            state_set[i][j] = state[i][j];
}

void train(string dataset[2][9]) {
    for (int i = 0; i < 2; ++i) {
        for (int l = 0; l < L; ++l) {
            int o[80], length = 0;
            for (int j = 0; j < 9; ++j) {
                for (unsigned int k = 0; k < dataset[i][j].size(); ++k) {
                    o[length] = dataset[i][j][k] - 'A';
                    length++;
                }
            }
            hmm[i].learn(o, length);
        }
    }
}

int main()
{
    vector<double>           initi_set(N, 0);
    vector< vector<double> > trans_set(N, vector<double>(N, 0));
    vector< vector<double> > state_set(N, vector<double>(M, 0));

    initialize(initi_set, trans_set, state_set);

    string dataset[2][9] = {{"ABBCABCAABC", "ABCABC", "ABCAABC", "BBABCAB", "BCAABCCAB", "CACCABCA", "CABCABCA", "CABCA", "CABCA"},
                            {"BBBCCBC",     "CCBABB", "AACCBBB", "BBABBAC", "CCAABBAB",  "BBBCCBAA", "ABBBBABA", "CCCCC", "BBAAA"}};

    for(int i = 0; i < 2; ++i) {
        hmm[i].set_initi(initi_set);
        hmm[i].set_trans(trans_set);
        hmm[i].set_state(state_set);
    }
    train(dataset);
//    hmm[0].print_HMM();
//    hmm[1].print_HMM();
    return 0;
}
