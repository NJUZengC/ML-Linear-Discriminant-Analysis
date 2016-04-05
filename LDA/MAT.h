#ifndef __MAT_H_
#define __MAT_H


//README  二维矩阵数组下标从1开始！！


class MAT{

public:
	int row;//行数
	int col;//列数
	double **num;//二维数组

	MAT();
	MAT(const MAT &a);
	MAT(int r, int c);//建立行为r列为c的矩阵
	~MAT();

	MAT getrow(int i);//返回矩阵第i行的向量
	MAT getmid();//求矩阵的同一列均值
	MAT trans();//转置矩阵
	MAT getNi();//求矩阵的逆矩阵
	/* 参考引用网络代码  http://blog.163.com/beckhero@126/blog/static/365417832009111424347517/ */
	

	MAT add(const MAT &b);//矩阵相加
	MAT sub(const MAT &b);//矩阵相减
	MAT multi(const MAT &b);//矩阵相乘

	MAT SW();//求类内散度矩阵sw

	void print();//矩阵信息输出
	void spilt(MAT *a, MAT *b, int chushu, int yushu);//把本矩阵中行数除以chushu余数为yushu的行赋给a,其余赋给b
	
	MAT & operator =(const MAT &a);//等号重载


};

#endif