#ifndef WORKTHREAD_H
#define WORKTHREAD_H

#define PRIMENUMBER_P 17
#define PRIMENUMBER_MOD 6151

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
#include <QThread>
#include <QObject>
#include <armadillo>
#include <QFileInfo>
#include <QCoreApplication>

using namespace arma;

using namespace cv;
using namespace std;

class DataProcess;

namespace mythread
{
  class workThread;
}

struct non_sigularMatrix
{
  vector<vector<double>> matri;
};

struct embed_res
{
  vector<vector<double>> alter_y;
  unsigned long long hash_w;
};

struct query_Template
{
  vector<vector<double>> Tyi;
};

struct sub_Template
{
  vector<vector<double>> Cxi;
};

struct param
{
  int seperateCount;
  double xita;
  int dimension;
};
struct privateKey
{
  vector<vector<double>> M1;
  vector<vector<double>> M2;
  int scale;
  vector<int> randomEx;
};
class workThread : public QObject
{
  // int id;
  Q_OBJECT
public:
  workThread(QObject *parent = nullptr);
  ~workThread();

  bool isWorking;
  bool isSent;
public slots:
  void startWork();
  void doWork();
  void getXita(double);
  void getDimension(int);
  void getRegisterVector(vector<double>);
  void getWorkMode(int);
  void getChallengeVector(vector<double>);
  void getUserID(QString);

signals:
  void workFinished();
  void workStart();
  void sendMessage(QString);

private:
  bool registerStage(const vector<double> &);
  bool loginStage(const vector<double> &);

  vector<privateKey> keyGenerationMode(vector<int> &);
  vector<sub_Template> encodingMode(const vector<privateKey> &, const vector<vector<double>> &);
  embed_res embeddingMode(const vector<privateKey> &, const vector<double> &, vector<int>);
  vector<query_Template> tokenGenerationMode(const vector<privateKey> &, const vector<vector<double>> &);
  vector<int> extractionMod(const vector<double> &, unsigned long long);
  vector<double> decodingMode(const vector<sub_Template> &, const vector<query_Template> &);
  void PrintMatrix(const vector<vector<double>> &);
  bool registerDataSaving(QString, const vector<privateKey> &, const vector<non_sigularMatrix> &, const vector<unsigned long long> &, const vector<sub_Template> &);
  vector<non_sigularMatrix> getRegisterRandomMatrixP(QString, const vector<int> &);
  vector<sub_Template> getFeatureTemplateCx(QString, const vector<int> &);
  vector<privateKey> getRegisterPrivateKey(QString);
  void PrintMatrix(const vector<vector<double>> &, QString);

  void PrintVector(const vector<double> &);
  void PrintVector(const vector<double> &, QString);
  void PrintVector(const vector<int> &);
  void PrintVector(const vector<int> &, QString);

  int workMode;
  vector<double> registerVector;
  vector<double> challengeVector;
  DataProcess *tool;
  param systemParam;
  vector<sub_Template> featureMatrix_Cx;
  vector<non_sigularMatrix> randomMatrix_P;
  vector<unsigned long long> hash_P; // h1 = hash(p)
  unsigned long long hash_W;         // h2 = hash(w)

  vector<non_sigularMatrix> randomMatrix_Q;
  vector<non_sigularMatrix> challengeMatrix_U;
  vector<non_sigularMatrix> queryMatrix_V;
  QString userID;
};

#endif // WORKTHREAD_H