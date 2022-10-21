#ifndef DATAPROCESS_H
#define DATAPROCESS_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <QString>
#include <opencv2/opencv.hpp>
#include <QFile>
#include <QFileDialog>
#include <QDir>
#include <QTextStream>
#include <QTextCodec>
#include <cmath>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <random>
#include <QtGlobal>
#include <QTime>
#include <QRandomGenerator>
#include <QList>
#include <QStringList>
#include <armadillo>

using namespace arma;

using namespace cv;
using namespace std;

namespace tool
{
  class DataProcess;
}

class DataProcess
{

public:
  int scale;
  DataProcess(int);
  ~DataProcess();

  mat vvToMatrixXd(const vector<vector<double>> &);
  vector<vector<double>> MatrixXdToVV(mat);

  void randomSeperate(const vector<double> &, vector<vector<double>> &, int);
  void getVectorFromFile(QString, vector<double> &);
  unsigned long long strToULongLong(QString);

  double strTodouble(QString); // Prevent loss of accuracy
  void getInvertibleMatrix(vector<vector<double>> &, int);
  vector<vector<double>> getInvertibleMatrix(int);

  void getIdentityMatrix(vector<vector<double>> &, int);
  vector<int> getRandomOrder(int);
  vector<vector<double>> getDiagonalMatrix(const vector<double> &, const vector<int> &);
  vector<vector<double>> getDiagonalMatrix(const vector<double> &);
  vector<vector<double>> getRandomLowerTriangleMatrix(int);

  vector<int> getRandomOrderVector(const vector<int> &);
  vector<double> getRandomOrderVector(const vector<double> &);
  vector<double> fixVectorOrder(const vector<double> &, const vector<int> &);

  vector<int> getNonzeroRandomVector(int);
  vector<int> getFixedSumRandomIntVector(int, int);
  vector<double> getFixedSumRandomDoubleVector(double, int);
  vector<vector<double>> getEmbedRandomVV(const vector<vector<double>> &, double);
  vector<vector<double>> getEmbedRandomVV(const vector<int> &, double);

  vector<int> getSign(const vector<double> &, const vector<double> &);
  vector<int> getSign(const vector<double> &, int);

  vector<int> getNiVector(const vector<vector<double>> &);

  unsigned long long getMatrixHash(const vector<vector<double>> &, int, int);
  unsigned long long getMatrixHash(const vector<vector<int>> &, int, int);

  unsigned long long getVectorHash(const vector<double> &, int, int);
  unsigned long long getVectorHash(const vector<int> &, int, int);

  bool echoHashVectortoFile(const vector<unsigned long long> &, QString, QString);
  bool echoVectortoFile(const vector<double> &, QString, QString);
  bool echoVectortoFile(const vector<int> &, QString, QString);
  bool echoHashVtoFile(const vector<unsigned long long> &, QString, QString);
  bool echoVVtoFile(const vector<vector<double>> &, QString, QString);
  bool echoVVtoFile(const vector<vector<int>> &, QString, QString);
  bool echoQStrtoFile(const QString, QString, QString);
  vector<vector<double>> getDoubleVVfromFile(vector<int> &, QString, QString);
  vector<vector<double>> getDoubleVVfromFile(int, QString, QString);

  vector<vector<int>> getIntegerVVfromFile(vector<int> &, QString, QString);
  vector<vector<int>> getIntegerVVfromFile(int, QString, QString);

  vector<int> getIntegerVectorfromFile(int, QString, QString);
  vector<unsigned long long> getHashVectorfromFile(int, QString, QString);
  vector<double> getDoubleVectorfromFile(int, QString, QString);
  QString getQstrFromFile(QString, QString);
  QString getFileFolderPath(const QString);

  vector<vector<double>> matrixMultiplication(const vector<vector<double>> &, const vector<vector<double>> &);
  double getMatrixValue(const vector<vector<double>> &);

private:
};

#endif