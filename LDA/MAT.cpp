#include"MAT.h"
#include<iostream>
using namespace std;


MAT::MAT(){
		row = 1000;
		int n = 1001;
		num = new double*[n];
		if (num == NULL)
		{
			cout << "NO ENOUGH SPACE!" << endl;
			exit(0);
		}
		for (int i = 0; i <= n; i++) {
			num[i] = new double[21];
			if (num[i] == NULL)
			{
				cout << "NO ENOUGH SPACE!" << endl;
				exit(0);
			}
		}
}

MAT::MAT(const MAT &a)
{
	row = a.row;
	col = a.col;
	num = new double*[(row + 1)];
	for (int i = 0; i <= row; i++)
		num[i] = new double[(col + 1)];
	for (int i = 1; i <= row; i++)
		for (int j = 1; j <= col; j++)
			num[i][j] = a.num[i][j];

};

//建立行为r列为c的矩阵
MAT::MAT(int r, int c){
		row = r;
		col = c;
		num = new double*[(r + 1)];
		if (num == NULL)
		{
			cout << "NO ENOUGH SPACE!" << endl;
			exit(0);
		}
		for (int i = 0; i <= r; i++) {
			num[i] = new double[(c + 1)];
			if (num[i] == NULL)
			{
				cout << "NO ENOUGH SPACE!" << endl;
				exit(0);
			}
		}
}

MAT::~MAT()
{
	for (int i = 0; i <= row; i++)
	{
		if (num[i] != NULL)
		{
			delete[] num[i];
			num[i] = NULL;
		}
	}
	if (num != NULL)
	{
		delete[] num;
		num = NULL;
	}
};




//返回矩阵第i行的向量
MAT MAT::getrow(int i)
	{
		MAT c(1, col);
		for (int j = 1; j <= col; j++)
			c.num[1][j] = num[i][j];
		return c;

};

//求矩阵的同一列均值
MAT MAT::getmid()
	{
		MAT c(1, col);
		for (int j = 1; j <= col; j++)
		{
			double sum = 0;
			for (int i = 1; i <= row; i++)
				sum += num[i][j];
			c.num[1][j] = sum / double(row);
		}
		return c;
};

//转置矩阵
MAT MAT::trans()
{
		MAT c(col, row);

		for (int i = 1; i <= row; i++)
			for (int j = 1; j <= col; j++)

				c.num[j][i] = num[i][j];
		return c;

};

//求矩阵的逆矩阵
MAT MAT::getNi()/* 参考引用网络代码  http://blog.163.com/beckhero@126/blog/static/365417832009111424347517/ */
{
	int n = row;
	MAT res(col, col);
	int i, j, k, m = 2 * n;
	double mik, temp;
	double **a = new double*[n];
	double **B = new double*[n];

	for (i = 0; i<n; i++)
	{
		a[i] = new double[2 * n];
		B[i] = new double[n];
	}

	for (i = 0; i<n; i++)
	{
		for (j = 0; j<n; j++)
		{
			if (i == j)
				B[i][j] = 1.0;
			else
				B[i][j] = 0.0;
		}
	}        //初始化B=E

	for (i = 0; i<n; i++)
		for (j = 0; j<n; j++)
			a[i][j] = num[i + 1][j + 1];  //复制A到a，避免改变A的值
	for (i = 0; i<n; i++)
		for (j = n; j<m; j++)
			a[i][j] = B[i][j - n];  //复制B到a，增广矩阵

	for (k = 1; k <= n - 1; k++)
	{
		for (i = k + 1; i <= n; i++)
		{
			mik = a[i - 1][k - 1] / a[k - 1][k - 1];
			for (j = k + 1; j <= m; j++)
			{
				a[i - 1][j - 1] -= mik*a[k - 1][j - 1];
			}
		}
	}        //顺序高斯消去法化左下角为零

	for (i = 1; i <= n; i++)
	{
		temp = a[i - 1][i - 1];
		for (j = 1; j <= m; j++)
		{
			a[i - 1][j - 1] /= temp;
		}
	}        //归一化

	for (k = n - 1; k >= 1; k--)
	{
		for (i = k; i >= 1; i--)
		{
			mik = a[i - 1][k];
			for (j = k + 1; j <= m; j++)
			{
				a[i - 1][j - 1] -= mik*a[k][j - 1];
			}
		}
	}        //逆序高斯消去法化增广矩阵左边为单位矩阵

	for (i = 0; i<n; i++)
		for (j = 0; j<n; j++)
			B[i][j] = a[i][j + n];  //取出求逆结果

	delete[]a;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			res.num[i + 1][j + 1] = B[i][j];
		}
	}
	delete[]B;

	return res;

};



//矩阵相加
MAT MAT::add(const MAT &b)
	
{
		MAT c(row, col);
		for (int i = 1; i <= row; i++)
			for (int j = 1; j <= col; j++)

				c.num[i][j] = num[i][j] + b.num[i][j];
		return c;
	
};

//矩阵相减
MAT MAT::sub(const MAT &b)
	
{
		MAT c(row, col);

		for (int i = 1; i <= row; i++)
			for (int j = 1; j <= col; j++)

				c.num[i][j] = num[i][j] - b.num[i][j];
		return c;
	
};

//矩阵相乘
MAT MAT::multi(const MAT &b)
	
{
		MAT c(row, b.col);
		for (int i = 1; i <= c.row; i++)
			for (int j = 1; j <= c.col; j++)
			{
				double m = 0;
				for (int n = 1; n <= col; n++)
					m += num[i][n] * b.num[n][j];
				c.num[i][j] = m;
			}

		return c;

};




//求类内散度矩阵sw
MAT MAT::SW()
{

		MAT d = getrow(1).sub(getmid()).trans();
		MAT c = d.multi(d.trans());

		for (int i = 2; i <= row; i++)
		{
			MAT e = getrow(i).sub(getmid()).trans();
			c = c.add(e.multi(e.trans()));

		}

		return c;
};


//矩阵信息输出
void MAT::print()
{
		cout << "----------------------------" << endl;
		for (int i = 1; i <= row; i++)
		{
			cout << "    ";
			for (int j = 1; j <= col; j++)
			{

				cout << num[i][j] << " ";
			}
			cout << endl;
		}
		cout << "----------------------------" << endl;
};

//把本矩阵中行数除以chushu余数为yushu的行赋给a,其余赋给b
void MAT::spilt(MAT *a, MAT *b, int chushu, int yushu)
{
		(*a).row = (*b).row = 0;
		(*a).col = (*b).col = col;
		for (int i = 1; i <= row; i++)
		{
			if (i%chushu == yushu)
			{
				++(*a).row;
				for (int j = 0; j <= col; j++)
					(*a).num[(*a).row][j] = num[i][j];
			}
			else
			{
				++(*b).row;
				for (int j = 0; j <= col; j++)
					(*b).num[(*b).row][j] = num[i][j];
			}
		}
};



//等号重载
MAT & MAT::operator =(const MAT &a)
{

		if (this == &a)
		{
			return *this;
		}
		for (int i = 0; i <= row; i++)
		{
			if (num[i] != NULL)
			{
				delete[]  num[i];
				num[i] = NULL;
			}
		}
		if (num != NULL)
		{
			delete[] num;
			num = NULL;
		}

		row = a.row;
		col = a.col;
		num = new double*[(row + 1)];
		for (int i = 0; i <= row; i++)
			num[i] = new double[(col + 1)];
		for (int i = 1; i <= row; i++)
			for (int j = 1; j <= col; j++)
				num[i][j] = a.num[i][j];
		return *this;
};


