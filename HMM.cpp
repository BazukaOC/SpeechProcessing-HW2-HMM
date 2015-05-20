#include "HMM.h"

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

double HMM::Logforward(int* o, int T) {
    for (int t = 0; t < T; ++t) {
        for (int j = 0; j < N; ++j) {
            if (t == 0) {
                forwa_set[t][j] = log(initi_set[j]) + log(state_set[j][o[t]]);
            }
            else {

                double p = forwa_set[t-1][0] + log(trans_set[0][j]);
                for (int i = 1; i < N; ++i) {
                    double temp = forwa_set[t-1][i] + log(trans_set[i][j]);
                    p = LogAdd(p, temp);
                }
                forwa_set[t][j] = p + log(state_set[j][o[t]]);
            }
        }
    }
    double p = forwa_set[T-1][0];
    for (int i = 1; i < N; ++i) {
        p = LogAdd(p, forwa_set[T-1][i]);
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

double HMM::Logbakward(int* o, int T) {
    for (int t = T - 1; t >= 0; --t) {
        for (int i = 0; i < N; ++i) {
            if (t == T-1)
                backw_set[t][i] = 0.0;
            else
            {
                double p = log(trans_set[i][0]) + log(state_set[0][o[t+1]]) + backw_set[t+1][0];
                for (int j = 1; j < N; ++j) {
                    double temp = log(trans_set[i][j]) + log(state_set[j][o[t+1]]) + backw_set[t+1][j];
                    p = LogAdd(p, temp);
                }
                backw_set[t][i] = p;
            }
        }
    }
    double p = log(initi_set[0]) + log(state_set[0][o[0]]) + backw_set[0][0];
    for (int j = 1; j < N; ++j) {
        double temp = log(initi_set[j]) + log(state_set[j][o[0]]) + backw_set[0][j];
        p = LogAdd(p, temp);
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

void   HMM::Loglearn(int* o, int T) {
    Logforward(o, T);
    Logbakward(o, T);

    for (int t = 0; t < T; ++t) {
        double p = LZERO;
        for (int i = 0; i < N; ++i) {
            double temp = forwa_set[t][i] + backw_set[t][i];
            p = LogAdd(p, temp);
        }
        //assert(p != 0);

        for (int i = 0; i < N; ++i) {
            lear1_set[t][i] = forwa_set[t][i] + backw_set[t][i] - p;
        }
    }

    for (int t = 0; t < T - 1; ++t) {
        double p = LZERO;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double temp = forwa_set[t][i] + log(trans_set[i][j]) + log(state_set[j][o[t+1]]) + backw_set[t+1][j];
                p = LogAdd(p, temp);
            }
        }
        //assert(p != 0);

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                lear2_set[t][i][j] = forwa_set[t][i] + log(trans_set[i][j]) + log(state_set[j][o[t+1]]) + backw_set[t+1][j] - p;
            }
        }
    }

    // 更新Π
    for (int i = 0; i < N; ++i) {
        initi_set[i] = exp(lear1_set[0][i]);
    }

    // 更新A
    for (int i = 0; i < N; ++i) {
        double p2 = LZERO;
        for (int t = 0; t < T - 1; ++t) {
            p2 = LogAdd(p2, lear1_set[t][i]);
        }
        //assert(p2 != 0);

        for (int j = 0; j < N; ++j) {
            double p1 = LZERO;
            for (int t = 0; t < T - 1; ++t) {
                p1 = LogAdd(p1, lear2_set[t][i][j]);
            }
            trans_set[i][j] = exp(p1 - p2);
        }
    }

    // 更新B
    for (int i = 0; i < N; ++i) {
        double p[M], p2 = LZERO;
        for(int j = 0; j < M; ++j) {
            p[j] = LZERO;
        }
        for (int t = 0; t < T; ++t) {
            p[o[t]] = LogAdd(p[o[t]], lear1_set[t][i]);
            p2 = LogAdd(p2, lear1_set[t][i]);
        }
        //assert(p2 != 0);

        for (int k = 0; k < M; ++k) {
            state_set[i][k] = exp(p[k] - p2);
        }
    }
}

void   HMM::print_HMM() {
    printf("initi_set:\n");
    for(int j = 0; j < N; ++j) {
        printf("    %.8f, ", initi_set[j]);
    }
    printf("\n");
    printf("state_set\n");
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < M; ++j) {
            printf("    %.8f, ", state_set[i][j]);
        }
        printf("\n");
    }
    printf("trans_set:\n");
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            printf("    %.8f, ", trans_set[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void   HMM::print_HMM(string title) {
    printf("initi_set: ");
    cout << title << endl;
    for(int j = 0; j < N; ++j) {
        printf("    %.8f, ", initi_set[j]);
    }
    printf("\n");
    printf("state_set\n");
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < M; ++j) {
            printf("    %.8f, ", state_set[i][j]);
        }
        printf("\n");
    }
    printf("trans_set:\n");
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            printf("    %.8f, ", trans_set[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

double HMM::decode(int* o, int T, int* q) {
    for (int t = 0; t < T; ++t) {
        for (int j = 0; j < N; ++j) {
            if (t == 0)
                deco1_set[t][j] = initi_set[j] * state_set[j][o[t]];
            else
            {
                double p = -1e9;
                for (int i = 0; i < N; ++i)
                {
                    double w = deco1_set[t-1][i] * trans_set[i][j];
                    if (w > p) {
                        p = w;
                        deco2_set[t][j] = i;
                    }
                }
                deco1_set[t][j] = p * state_set[j][o[t]];
            }
        }
    }

    double p = -1e9;
    for (int j = 0; j < N; ++j) {
        if (deco1_set[T-1][j] > p) {
            p = deco1_set[T-1][j];
            q[T-1] = j;
        }
    }

    for (int t = T - 1; t > 0; --t) {
        q[t-1] = deco2_set[t][q[t]];
    }

    return p;
}

void   HMM::clear() {
    for(int i = 0; i < N; ++i)
        initi_set[i] = 0;

    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            trans_set[i][j] = 0;

    for(int i = 0; i < N; ++i)
        for(int j = 0; j < M; ++j)
            state_set[i][j] = 0;

    for(int i = 0; i < T; ++i)
        for(int j = 0; j < N; ++j) {
            forwa_set[i][j] = 0;
            backw_set[i][j] = 0;
            deco1_set[i][j] = 0;
            deco2_set[i][j] = 0;
            lear1_set[i][j] = 0;
        }
    for(int i = 0; i < T; ++i)
        for(int j = 0; j < N; ++j)
            for(int k = 0; k < N; ++k)
                lear2_set[i][j][k] = 0;
}
