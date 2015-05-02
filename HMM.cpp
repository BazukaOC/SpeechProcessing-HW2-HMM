#include "HMM.h"
#include <stdio.h>

void   HMM::set_initi(vector<double> &a) {
    initi_set = a;
}

void   HMM::set_trans(vector< vector<double> > &a) {
    trans_set = a;
}

void   HMM::set_state(vector< vector<double> > &a) {
    state_set = a;
}

double HMM::forward(int* o, int T) {
    for (int t = 0; t < T; ++t) {
        for (int j = 0; j < N; ++j) {
            if (t == 0) {
                forwa_set[t][j] = initi_set[j] * state_set[j][o[t]];
            }
            else {

                double p = 0;
                for (int i = 0; i < N; ++i) {
                    p += forwa_set[t-1][i] * trans_set[i][j];
                }
                forwa_set[t][j] = p * state_set[j][o[t]];
            }
        }
    }
    double p = 0;
    for (int i = 0; i < N; ++i) {
        p += forwa_set[T-1][i];
    }
    return p;
}

double HMM::bakward(int* o, int T) {
    for (int t = T - 1; t >= 0; --t) {
        for (int i = 0; i < N; ++i) {
            if (t == T-1)
                backw_set[t][i] = 1.0;
            else
            {
                double p = 0;
                for (int j = 0; j < N; ++j) {
                    p += trans_set[i][j] * state_set[j][o[t+1]] * backw_set[t+1][j];
                }
                backw_set[t][i] = p;
            }
        }
    }
    double p = 0;
    for (int j = 0; j < N; ++j) {
        p += initi_set[j] * state_set[j][o[0]] * backw_set[0][j];
    }
    return p;
}

void   HMM::learn(int* o, int T) {
    forward(o, T);
    bakward(o, T);

    for (int t = 0; t < T; ++t) {
        double p = 0;
        for (int i = 0; i < N; ++i) {
            p += forwa_set[t][i] * backw_set[t][i];
        }
        //assert(p != 0);

        for (int i = 0; i < N; ++i) {
            lear1_set[t][i] = forwa_set[t][i] * backw_set[t][i] / p;
        }
    }

    for (int t = 0; t < T - 1; ++t) {
        double p = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                p += forwa_set[t][i] * trans_set[i][j] * state_set[j][o[t+1]] * backw_set[t+1][j];
            }
        }
        //assert(p != 0);

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                lear2_set[t][i][j] = forwa_set[t][i] * trans_set[i][j] * state_set[j][o[t+1]] * backw_set[t+1][j] / p;
            }
        }
    }

    // 更新Π
    for (int i = 0; i < N; ++i) {
        initi_set[i] = lear1_set[0][i];
    }

    // 更新A
    for (int i = 0; i < N; ++i) {
        double p2 = 0;
        for (int t = 0; t < T - 1; ++t) {
            p2 += lear1_set[t][i];
        }
        //assert(p2 != 0);

        for (int j = 0; j < N; ++j) {
            double p1 = 0;
            for (int t = 0; t < T - 1; ++t) {
                p1 += lear2_set[t][i][j];
            }
            trans_set[i][j] = p1 / p2;
        }
    }

    // 更新B
    for (int i = 0; i < N; ++i) {
        double p[M], p2 = 0;
        for(int j = 0; j < M; ++j) {
            p[j] = 0;
        }
        for (int t = 0; t < T; ++t) {
            p[o[t]] += lear1_set[t][i];
            p2 += lear1_set[t][i];
        }
        //assert(p2 != 0);

        for (int k = 0; k < M; ++k) {
            state_set[i][k] = p[k] / p2;
        }
    }
}

void   HMM::print_HMM() {
    printf("initi_set:\n");
    for(int j = 0; j < N; ++j) {
        printf("    %.16f, ", initi_set[j]);
    }
    printf("\n");
    printf("state_set\n");
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < M; ++j) {
            printf("    %.16f, ", state_set[i][j]);
        }
        printf("\n");
    }
    printf("trans_set:\n");
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            printf("    %.16f, ", trans_set[i][j]);
        }
        printf("\n");
    }
}
