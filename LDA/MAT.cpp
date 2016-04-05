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

//������Ϊr��Ϊc�ľ���
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




//���ؾ����i�е�����
MAT MAT::getrow(int i)
	{
		MAT c(1, col);
		for (int j = 1; j <= col; j++)
			c.num[1][j] = num[i][j];
		return c;

};

//������ͬһ�о�ֵ
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

//ת�þ���
MAT MAT::trans()
{
		MAT c(col, row);

		for (int i = 1; i <= row; i++)
			for (int j = 1; j <= col; j++)

				c.num[j][i] = num[i][j];
		return c;

};

//�����������
MAT MAT::getNi()/* �ο������������  http://blog.163.com/beckhero@126/blog/static/365417832009111424347517/ */
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
	}        //��ʼ��B=E

	for (i = 0; i<n; i++)
		for (j = 0; j<n; j++)
			a[i][j] = num[i + 1][j + 1];  //����A��a������ı�A��ֵ
	for (i = 0; i<n; i++)
		for (j = n; j<m; j++)
			a[i][j] = B[i][j - n];  //����B��a���������

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
	}        //˳���˹��ȥ�������½�Ϊ��

	for (i = 1; i <= n; i++)
	{
		temp = a[i - 1][i - 1];
		for (j = 1; j <= m; j++)
		{
			a[i - 1][j - 1] /= temp;
		}
	}        //��һ��

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
	}        //�����˹��ȥ��������������Ϊ��λ����

	for (i = 0; i<n; i++)
		for (j = 0; j<n; j++)
			B[i][j] = a[i][j + n];  //ȡ��������

	delete[]a;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			res.num[i + 1][j + 1] = B[i][j];
		}
	}
	delete[]B;

	return res;

};



//�������
MAT MAT::add(const MAT &b)
	
{
		MAT c(row, col);
		for (int i = 1; i <= row; i++)
			for (int j = 1; j <= col; j++)

				c.num[i][j] = num[i][j] + b.num[i][j];
		return c;
	
};

//�������
MAT MAT::sub(const MAT &b)
	
{
		MAT c(row, col);

		for (int i = 1; i <= row; i++)
			for (int j = 1; j <= col; j++)

				c.num[i][j] = num[i][j] - b.num[i][j];
		return c;
	
};

//�������
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




//������ɢ�Ⱦ���sw
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


//������Ϣ���
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

//�ѱ���������������chushu����Ϊyushu���и���a,���ำ��b
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



//�Ⱥ�����
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


