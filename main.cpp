#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <string>
#include <stdio.h>
#include "HMM.h"

using namespace std;

int N = 3, M = 3, T = 80, L = 0;

double initi[3]    =  {0.34, 0.33, 0.33};
double trans[3][3] = {{0.34, 0.33, 0.33}, {0.33, 0.34, 0.33}, {0.33, 0.33, 0.34}};
double state[3][3] = {{0.34, 0.33, 0.33}, {0.33, 0.34, 0.33}, {0.33, 0.33, 0.34}};

vector<double>           initi_set(N, 0);
vector< vector<double> > trans_set(N, vector<double>(N, 0));
vector< vector<double> > state_set(N, vector<double>(M, 0));

HMM hmm[2] = HMM(N, M, T, L);

void HMMInitialize() {
    for(int i = 0; i < 2; ++i)
        hmm[i].clear();

    for(int i = 0; i < N; ++i)
        initi_set[i] = initi[i];

    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            trans_set[i][j] = trans[i][j];

    for(int i = 0; i < N; ++i)
        for(int j = 0; j < M; ++j)
            state_set[i][j] = state[i][j];

    for(int i = 0; i < 2; ++i) {
        hmm[i].set_initi(initi_set);
        hmm[i].set_trans(trans_set);
        hmm[i].set_state(state_set);
    }
}

void OMMInitialize() {
    for(int i = 0; i < 2; ++i)
        hmm[i].clear();

    for(int i = 0; i < N; ++i)
        initi_set[i] = initi[i];

    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            trans_set[i][j] = trans[i][j];

    state_set[0][0] = 1.0; state_set[0][1] = 0.0; state_set[0][2] = 0.0;
    state_set[1][0] = 0.0; state_set[1][1] = 1.0; state_set[1][2] = 0.0;
    state_set[2][0] = 0.0; state_set[2][1] = 0.0; state_set[2][2] = 1.0;

    for(int i = 0; i < 2; ++i) {
        hmm[i].set_initi(initi_set);
        hmm[i].set_trans(trans_set);
        hmm[i].set_state(state_set);
    }
}

void train(string dataset[2][9]) {
    for (int i = 0; i < 2; ++i) {
        for (int l = 0; l < L; ++l) {
            int o[T], length = 0;
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

void Logtrain(string dataset[2][9]) {
    for (int i = 0; i < 2; ++i) {
        for (int l = 0; l < L; ++l) {
            int o[T], length = 0;
            for (int j = 0; j < 9; ++j) {
                for (unsigned int k = 0; k < dataset[i][j].size(); ++k) {
                    o[length] = dataset[i][j][k] - 'A';
                    length++;
                }
            }
            hmm[i].Loglearn(o, length);
        }
    }
}

int  recognize(int* o, int T) {
    int ans = -1, q[12] = {0};
    double p = -1e9;
    for (int i = 0; i < 2; ++i) {
        double pp = hmm[i].decode(o, T, q);
        if (pp > p) p = pp, ans = i;
    }
    return ans;
}

void insideTest(string dataset[2][9]) {
    int o[15];
    for (int i = 0; i < 2; ++i) {
        for (int j = 0, k = 0; j < 9; ++j) {
            cout << "    ";
            for (k = 0; k < dataset[i][j].size(); ++k) {
                cout << dataset[i][j][k];
                o[k]  = dataset[i][j][k] - 'A';
            }
            cout << " model " << (char)('A' + recognize(o, k)) << " " << endl;
        }
        cout << endl;
    }
}

int main()
{
    string dataset[2][9] = {{"ABBCABCAABC", "ABCABC", "ABCAABC", "BBABCAB", "BCAABCCAB", "CACCABCA", "CABCABCA", "CABCA", "CABCA"},
                            {"BBBCCBC",     "CCBABB", "AACCBBB", "BBABBAC", "CCAABBAB",  "BBBCCBAA", "ABBBBABA", "CCCCC", "BBAAA"}};
    bool isHMM = true;
    int select = 0;
    do {
        isHMM ? HMMInitialize() : OMMInitialize();

        cout << "1. Show model parameters." << endl;
        cout << "2. Inside test." << endl;
        cout << "3. sequences belong to which one?" << endl;
        cout << "4. Switch to Observable Markov Models." << endl;
        cout << "0. exit. >>> ";
        cin  >> select;
        switch(select) {
            case 1: {
            /// P1. Please specify the model parameters after the first and 50th iterations of Baum-Welch training
                cout << "[ORIGIN MODEL PARAMETERS]" << endl;
                hmm[0].print_HMM("HMM_A--------------------------------");
                hmm[1].print_HMM("HMM_B--------------------------------");
                cout << "[TRAINED ONCE MODEL PARAMETERS]" << endl;
                L = 1;  isHMM ? Logtrain(dataset) : train(dataset);
                hmm[0].print_HMM("HMM_A--------------------------------");
                hmm[1].print_HMM("HMM_B--------------------------------");
                cout << "[TRAINED 50 TIMES MODEL PARAMETERS]" << endl;
                L = 49; isHMM ? Logtrain(dataset) : train(dataset);
                hmm[0].print_HMM("HMM_A--------------------------------");
                hmm[1].print_HMM("HMM_B--------------------------------");
                break;
            }
            case 2: {
            /// P2. Please show the recognition results by using the above training sequences as the testing data (The so-called inside testing).
                cout << "[TRAINED ONCE]" << endl;
                L = 1;  isHMM ? Logtrain(dataset) : train(dataset);
                insideTest(dataset);
                cout << "[TRAINED 50 TIMES]" << endl;
                L = 49; isHMM ? Logtrain(dataset) : train(dataset);
                insideTest(dataset);
                break;
            }
            case 3: {
            /// P3. Which class do the following testing sequences belong to?
                int o[15], k, s;
                string str1("ABCABCCAB");
                string str2("AABABCCCCBBB");
                for(k = 0; k < str1.size(); ++k)
                    o[k] = str1[k] - 'A';
                for(s = 0; s < str2.size(); ++s)
                    o[s] = str2[s] - 'A';
                cout << "[TRAINED ONCE]" << endl;
                L = 1;  isHMM ? Logtrain(dataset) : train(dataset);
                cout << "    " << str1 << " model " << (char)(recognize(o, k) + 'A') << endl;
                cout << "    " << str2 << " model " << (char)(recognize(o, s) + 'A') << endl;
                cout << "[TRAINED 50 TIMES]" << endl;
                L = 49; isHMM ? Logtrain(dataset) : train(dataset);
                cout << "    " << str1 << " model " << (char)(recognize(o, k) + 'A') << endl;
                cout << "    " << str2 << " model " << (char)(recognize(o, s) + 'A') << endl;
                break;
            }
            case 4: {
            /// P4. What are the results if Observable Markov Models were instead used in P1, P2 and P3?
                if (isHMM) {
                    isHMM = false;
                    cout << "Switch to Observable Markov Models." << endl;
                } else {
                    isHMM = true;
                    cout << "Switch to Hidden Markov Models." << endl;
                }
                break;
            }
        }
        system("pause");
        system("cls");
    } while (select);

    return 0;
}
