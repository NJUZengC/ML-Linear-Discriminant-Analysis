#ifndef __MAT_H_
#define __MAT_H


//README  ��ά���������±��1��ʼ����


class MAT{

public:
	int row;//����
	int col;//����
	double **num;//��ά����

	MAT();
	MAT(const MAT &a);
	MAT(int r, int c);//������Ϊr��Ϊc�ľ���
	~MAT();

	MAT getrow(int i);//���ؾ����i�е�����
	MAT getmid();//������ͬһ�о�ֵ
	MAT trans();//ת�þ���
	MAT getNi();//�����������
	/* �ο������������  http://blog.163.com/beckhero@126/blog/static/365417832009111424347517/ */
	

	MAT add(const MAT &b);//�������
	MAT sub(const MAT &b);//�������
	MAT multi(const MAT &b);//�������

	MAT SW();//������ɢ�Ⱦ���sw

	void print();//������Ϣ���
	void spilt(MAT *a, MAT *b, int chushu, int yushu);//�ѱ���������������chushu����Ϊyushu���и���a,���ำ��b
	
	MAT & operator =(const MAT &a);//�Ⱥ�����


};

#endif