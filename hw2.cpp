const int N = 3, M = 3, T = 15;
double init_set[N], a[N][N], b[N][M];  // HMM
double α[T][N], β[T][N];        // evaluation problem
double δ[T][N]; int ψ[T][N];    // decoding problem
double γ[T][N], ξ[T][N][N];     // learning problem
 
void learn(int* o, int T)
{
    forward(o, T);
    backward(o, T);
 
    for (int t=0; t<T; ++t)
    {
        double p = 0;
        for (int i=0; i<N; ++i)
            p += α[t][i] * β[t][i];
        assert(p != 0);
 
        for (int i=0; i<N; ++i)
            γ[t][i] = α[t][i] * β[t][i] / p;
    }
 
    for (int t=0; t<T-1; ++t)
    {
        double p = 0;
        for (int i=0; i<N; ++i)
            for (int j=0; j<N; ++j)
                p += α[t][i] * a[i][j] * b[j][o[t+1]] * β[t+1][j];
        assert(p != 0);
 
        for (int i=0; i<N; ++i)
            for (int j=0; j<N; ++j)
                ξ[t][i][j] = α[t][i] * a[i][j] * b[j][o[t+1]] * β[t+1][j] / p;
    }
 
    // 更新init_set
    for (int i=0; i<N; ++i)
        init_set[i] = γ[0][i];
 
    // 更新A
    for (int i=0; i<N; ++i)
    {
        double p2 = 0;
        for (int t=0; t<T-1; ++t)
            p2 += γ[t][i];
        assert(p2 != 0);
 
        for (int j=0; j<N; ++j)
        {
            double p1 = 0;
            for (int t=0; t<T-1; ++t)
                p1 += ξ[t][i][j];
            a[i][j] = p1 / p2;
        }
    }
 
    // 更新B
    for (int i=0; i<N; ++i)
    {
        double p[M] = {0}, p2 = 0;
        for (int t=0; t<T; ++t)
        {
            p[o[t]] += γ[t][i];
            p2 += γ[t][i];
        }
        assert(p2 != 0);
 
        for (int k=0; k<M; ++k)
            b[i][k] = p[k] / p2;
    }
}