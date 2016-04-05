#include"MAT.h"
#include<string>
#include<fstream>
#include<cstdio>
#include<cstdlib>
#include<math.h>
#include<iostream>
using namespace std;

//��������������
MAT melonpos(20,2);
MAT melonneg(20,2);

//fourclass����������
MAT fourpos(600,2);
MAT fourneg(600,2);

//heart����������
MAT heartpos(200,13);
MAT heartneg(200,13);


void initmelon()//��ʼ�������ļ�����
{
	melonpos.row = melonneg.row = 0;
	melonpos.col = 2;
	melonneg.col = 2;
	ifstream input;
	input.open("��������3.0a.txt");
	if (!input )
	{
		cout << "�ļ�������";
		exit(-1);
	}
	char s[100],s1[100];
	input.getline(s,100);
	input.getline(s1, 100);
	while (!input.fail())
	{
		double ma[3];
		int a; input >> a;
		input >> ma[0];
		input >> ma[1];
	
		char *p=new char[10];
		input >> p;
		if (strcmp(p, "") == 0)
			break;
		if (a<=8)
		{
			melonpos.row++;
			melonpos.num[melonpos.row][1] = ma[0];
			melonpos.num[melonpos.row][2] = ma[1];
		}
		else 
		{
			melonneg.row++;
			
			melonneg.num[melonneg.row][1] = ma[0];
			melonneg.num[melonneg.row][2] = ma[1];
		}
	}
	input.close();
}

void moduleMelon()//����ģ�����ģ��
{

	initmelon();
	MAT SW = melonneg.SW().add(melonpos.SW());
	MAT w = SW.getNi().multi(melonneg.getmid().sub(melonpos.getmid()).trans());
	cout <<endl<< "         w is��" << endl;
	w.print();
	cout << endl<< "   ��ѵ������" ;
	double mida = melonneg.getmid().multi(w).num[1][1];

	double midb = melonpos.getmid().multi(w).num[1][1];

	int falsenum = 0;
	for (int i = 1; i <= melonneg.row; i++)
	{
		double n = melonneg.getrow(i).multi(w).num[1][1];
		if (fabs(n - mida) > fabs(n - midb))
		{
			falsenum++;
		}
	}

	for (int i = 1; i <= melonpos.row; i++)
	{
		double n = melonpos.getrow(i).multi(w).num[1][1];
		if (fabs(n - mida) < fabs(n - midb))
		{
			falsenum++;
		}
	}
	double errorrate = double(falsenum) / (melonpos.row + melonneg.row);
	cout << errorrate << endl;
}


void initfour()//��ʼ��four�ļ�����
{
	fourpos.row = fourneg.row = 0;
	fourpos.col = 2;
	fourneg.col = 2;
	ifstream input;
	input.open("fourclass.csv");
	if (!input)
	{
		cout << "�ļ�������";
		exit(-1);
	}
	int line = 0;
	while (!input.fail())
	{
		line++;
		if (line > 862)
			break;
		double ma[3];
		char s;
		input >> ma[0];
		input >> s;

		input >> ma[1];
		input >> s;
		int a;
		input >> a;

		if (a == 1)
		{
			fourpos.row++;

			fourpos.num[fourpos.row][1] = ma[0];
			fourpos.num[fourpos.row][2] = ma[1];
		}
		else
		{
			fourneg.row++;

			fourneg.num[fourneg.row][1] = ma[0];
			fourneg.num[fourneg.row][2] = ma[1];
		}

	}
	input.close();

}

void moduleFour(int flag)//Fourģ�����ģ��
{
	initfour();
	cout << endl;
	if (flag == 0)
		cout << "               ʮ�۽��淨      " << endl << endl;
	else
		cout << "                ��һ��      " << endl << endl;
	cout << "    ������Ϊ: " << fourneg.row << "    ������Ϊ:" << fourpos.row << endl;
	MAT pracpos(fourpos.row, fourneg.col), testpos(fourpos.row, fourneg.col), pracneg(fourneg.row, fourneg.col), testneg(fourneg.row, fourneg.col), SW(fourneg.col, fourneg.col), w(fourneg.col, fourneg.col);
	int testclass = 10;
	if (flag == 1) testclass = fourneg.row + fourpos.row;
	double *error = new double[testclass];

	for (int k = 0; k < testclass; k++){
		if (flag == 0){
			fourpos.spilt(&testpos, &pracpos, 10, k);
			fourneg.spilt(&testneg, &pracneg, 10, k);
		}
		else{
			fourpos.spilt(&testpos, &pracpos, fourpos.row, (k<fourpos.row) ? k : fourpos.row);
			fourneg.spilt(&testneg, &pracneg, fourneg.row, (k<fourpos.row) ? fourneg.row : (k - fourpos.row));
		}
		SW = pracneg.SW().add(pracpos.SW());
		w = SW.getNi().multi(pracneg.getmid().sub(pracpos.getmid()).trans());
		/*cout << "w ����Ϊ��" << endl;
		w.print();*/

		double mida = pracneg.getmid().multi(w).num[1][1];

		double midb = pracpos.getmid().multi(w).num[1][1];

		int falsenum = 0;
		for (int i = 1; i <= testneg.row; i++)
		{
			double n = testneg.getrow(i).multi(w).num[1][1];
			if (fabs(n - mida) > fabs(n - midb))
			{
				falsenum++;
			}
		}

		for (int i = 1; i <= testpos.row; i++)
		{
			double n = testpos.getrow(i).multi(w).num[1][1];
			if (fabs(n - mida) < fabs(n - midb))
			{
				falsenum++;
			}
		}

		
		double errorrate = double(falsenum) / (testpos.row + testneg.row);
		error[k] = errorrate;
	}
	double miderro = 0;
	for (int j = 0; j < testclass; j++)
	{
		miderro += error[j];
	}
	miderro /= double(testclass);
	cout << "    �����ʾ�ֵΪ: " << miderro << endl;

	double biaozuncha = 0;
	for (int j = 0; j < testclass; j++)
	{
		biaozuncha += pow(error[j]-miderro,2);
	}
	biaozuncha /= double(testclass);
	biaozuncha = sqrt(biaozuncha);
	cout << "    �����ʱ�׼��ֵΪ: " << biaozuncha << endl;

	delete[] error;

}


void initheart()//��ʼ��heart�ļ�����
{
	heartpos.row = heartneg.row = 0;
	heartpos.col = 13;
	heartneg.col = 13;
	ifstream input;
	input.open("heart.csv");
	if (!input)
	{
		cout << "�ļ�������";
		exit(-1);
	}
	int line = 0;
	while (!input.fail())
	{
		line++;
		if (line > 270)
			break;
		double ma[14];
		char s;
		int time = 1;
		while (time < 14)
		{
			input >> ma[time];
			time++;
			input >> s;
		}

		int a;
		input >> a;

		if (a == 1)
		{
			heartpos.row++;
			for (int j = 1; j <= 13; j++){
				heartpos.num[heartpos.row][j] = ma[j];

			}
		}
		else
		{
			heartneg.row++;
			for (int j = 1; j <= 13; j++){
				heartneg.num[heartneg.row][j] = ma[j];

			}
		}

	}
	input.close();

}

void moduleHeart(int flag)//heartģ�����ģ��

{
	initheart();
	cout << endl;
	if (flag == 0)
		cout << "               ʮ�۽��淨      " << endl<<endl;
	else
		cout << "                ��һ��      " << endl<<endl;
	cout << "    ������Ϊ: "<< heartneg.row << "    ������Ϊ:" << heartpos.row << endl;
	MAT pracpos(heartpos.row, heartneg.col), testpos(heartpos.row, heartneg.col), pracneg(heartneg.row, heartneg.col), testneg(heartneg.row, heartneg.col), SW(heartneg.col, heartneg.col), w(heartneg.col, heartneg.col);
	int testclass = 10;
	if (flag == 1) testclass = heartneg.row + heartpos.row;
	double *error = new double[testclass];
	
	for (int k = 0; k < testclass; k++){
		if (flag == 0){
			heartpos.spilt(&testpos, &pracpos, 10, k);
			heartneg.spilt(&testneg, &pracneg, 10, k);
		}
		else{
			heartpos.spilt(&testpos, &pracpos, heartpos.row, (k<heartpos.row)?k:heartpos.row);
			heartneg.spilt(&testneg, &pracneg, heartneg.row, (k<heartpos.row)?heartneg.row:(k-heartpos.row));
		}
		SW = pracneg.SW().add(pracpos.SW());
		w = SW.getNi().multi(pracneg.getmid().sub(pracpos.getmid()).trans());
		/*cout << "w ����Ϊ��" << endl;
		w.print();*/

		double mida = pracneg.getmid().multi(w).num[1][1];

		double midb = pracpos.getmid().multi(w).num[1][1];

		int falsenum = 0;
		for (int i = 1; i <= testneg.row; i++)
		{
			double n = testneg.getrow(i).multi(w).num[1][1];
			if (fabs(n - mida) > fabs(n - midb))
			{
				falsenum++;
			}
		}

		for (int i = 1; i <= testpos.row; i++)
		{
			double n = testpos.getrow(i).multi(w).num[1][1];
			if (fabs(n - mida) < fabs(n - midb))
			{
				falsenum++;
			}
		}

		
		double errorrate = double(falsenum) / (testpos.row + testneg.row);
		error[k] = errorrate;
	}
	double miderro = 0;
	for (int j = 0; j < testclass; j++)
	{
		miderro += error[j];
	}
	miderro /= double(testclass);
	cout << "    �����ʾ�ֵΪ: " << miderro << endl;

	double biaozuncha = 0;
	for (int j = 0; j < testclass; j++)
	{
		biaozuncha += pow(error[j] - miderro, 2);
	}
	biaozuncha /= double(testclass);
	biaozuncha = sqrt(biaozuncha);
	cout << "    �����ʱ�׼��ֵΪ: " << biaozuncha << endl;

	delete[] error;
}

int main()
{
	while (1){
		int m = 0;
	A:	do{
		    cout << endl << "----------------------------" << endl;
			cout << "  ";

			cout<< "ѡ��ѡ��: " << endl;
			cout << "       0: �������� " << endl;
			cout << "       1: Fourclass " << endl;
			cout << "       2: Heart " << endl;
			cout << "       3: �˳� " << endl;
			cin >> m;
		} while (m < 0 || m>3);
		switch (m){
		case 0:  moduleMelon(); break;
		case 1:case 2:{
			while (1){
				int n = 0;
				do{
					cout << endl << "----------------------------" << endl;
					cout << "  ";
					if (m == 1)
						cout << "Fourclass";
					else
						cout << "Heart";
					cout<<"ѡ��ѡ��: " << endl;
					cout << "       0: ʮ�۽��淨 " << endl;
					cout << "       1:  ��һ�� " << endl;
					cout << "       2:  ������һ�� " << endl;
					cin >> n;
				} while (n < 0 || n>2);
				switch (n)
				{
				case 0: if (m == 1) moduleFour(0); else moduleHeart(0); break;
				case 1: if (m == 1) moduleFour(1); else moduleHeart(1); break;
				case 2: goto A; break;
				}
			}
			break;
		}
		case 3: return 0; break;

		}
	}

	return 0;
}