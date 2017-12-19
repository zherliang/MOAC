#include <iostream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <queue>
#include <cmath>
#include <map>
#include <stdlib.h>
#include <time.h>
#define eps 0.0000000001
#define maxn 3009
#define inf 0x3f3f3f3f
#define mod 1000000009
#define ll long long
using namespace std;
int n;//表示维度
int max_itr = 1000, itr = 1;//表示最大迭代次数
pair<int, double> dis[1009];//存储每个点的最小距离的点
vector<int> v[1009];//存储每种类别的点的标号
double diss;//需设定的邻居的距离参数
struct ee
{
	int wei;//表示维数
	double x[10];
} P[1009];
//目标的数组
double G0 = 0.02, G = 0.0;
int K;//表示要分成的种类数
int sui[109], ha[1009];
//表示随机数标号
int max_itt = 1;
double ddis[109], F[109], D[109], TF[109];
//F表示的是Force数组,TF是Total_Force数组
ee V[1009], A[1009];
double M[1009];
double AA = 1, BB = 0.5;
//需设定的参数（自己调试）
double R[1009][11], L[1009][11];
//上下限界
int main()
{
	//freopen("f:\\iris1.txt","r",stdin);
	int wei;
	srand((unsigned)time(NULL));
	scanf("%d%d%d", &n, &K, &wei);
	//n表示点的个数
	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= wei; j++)
		{
			scanf("%lf", &P[i].x[j]);
			double mia = 9999999, mib = 9999999;
			//这里需要设置一个最大值
			for (int op = 1; op <= n; op++)
			{
				if (P[op].x[j]<P[i].x[j] && mia>P[i].x[j] - P[op].x[j])
					mia = P[op].x[j];
				if (P[op].x[j] >= P[i].x[j] && mib>P[op].x[j] - P[i].x[j])
					mib = P[op].x[j];
			}
			L[i][j] = mia, R[i][j] = mib;
		}
	}
	for (int i = 1;i <= K;i++)
	{
		int tty = (rand() % (n)) + 1;
		while (ha[tty])
			tty = (rand() % (n)) + 1;
		sui[i] = tty;
		// printf("%d  ",sui[i]);
		ha[sui[i]] = i;
	}
	for (int i = 1;i <= n;i++)
		dis[i].second = 10000000009;
	while (itr <= max_itr)
	{
		for (int i = 1; i <= K; i++)
		{
			for (int j = 1; j <= n; j++)
			{
				double tmp = 0;
				for (int ij = 1; ij <= wei; ij++)
					tmp += (P[j].x[ij] - P[sui[i]].x[ij])*(P[j].x[ij] - P[sui[i]].x[ij]);
				if (tmp<dis[j].second)
					dis[j].first = sui[i], dis[j].second = tmp;
			}
		}
		for (int i = 1; i <= n; i++)
		{
			// printf("%d %d\n",ha[dis[i].first],i);
			v[ha[dis[i].first]].push_back(i);
		}
		if (itr == max_itt)
		{
			break;
		}
		//ddis中保存的是Bi
		for (int i = 1; i <= K; i++)
		{
			ddis[i] = 0;//记录内部簇数组
			int yu = sui[i];
			double vp = 0;
			for (unsigned int j = 0; j<v[yu].size(); j++)
			{
				int y = v[yu][j];
				double vr = 0;
				for (int o = 1; o <= wei; o++)
					vr += (P[yu].x[o] - P[y].x[o])*(P[yu].x[o] - P[y].x[o]);
				vr = sqrt(vr);
				vp += vr;
			}
			ddis[i] = vp;
		}
		//前二步调试部分
		/*
		*/
		for (int i = 1;i <= K;i++)
			M[i] = AA + BB*ddis[i];
		//...............................................F的计算
		for (int i = 1; i <= K; i++)
		{
			int yu = sui[i];
			F[i] = 0.0;
			for (unsigned int j = 0; j<v[i].size(); j++)
			{
				int y = v[yu][j];
				double t1 = 0.0, t2 = 0.0, mp = 0.0;
				for (int op = 1;op <= wei;op++)
					mp += abs((P[y].x[op] - P[i].x[op]) / (R[y][op] - L[y][op]));
				mp /= 1.0*wei;
				for (int op = 1; op <= wei; op++)
					t1 += (P[y].x[op] - P[i].x[op])*ddis[i];
				F[i] += t1 / mp;
			}
		}
		//前四步调试部分
		//...............................................TF,M的计算
		for (int i = 1; i <= n; i++)
			TF[i] = G / (1.0*v[i].size())*F[i];
		G = 1.0*G0*(1 - itr) / max_itr;
		//最后一步计算
		for (int i = 1; i <= K; i++)
		{
			int RAND = 1;
			int yu = sui[i];
			for (int j = 1;j <= wei;j++)
			{
				A[yu].x[j] = TF[i] / M[i] * RAND;
				V[yu].x[j] += A[yu].x[j];
			}
			//for()
			// A[yu]=TF[i]/M[i]*RAND;
			// V[yu]+=A[yu];
			//P[i].
		}
		for (int i = 1;i <= K;i++)
			for (int j = 1;j <= wei;j++)
				P[i].x[j] += V[i].x[j];
		itr++;
	}
	for (int i = 1;i <= K;i++)
	{
		printf("第%d类的点有这些:", i);
		for (unsigned int j = 0;j<v[i].size();j++)
			printf("%d ", v[i][j]);
		printf("\n");
	}
	system("pause");
	return 0;
}
