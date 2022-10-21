#ifndef MAINWINDOW_H
#define MAINWINDOW_H
  
#define FEATURE_SCALE 512
#define DEVIATION_VALUE 0.0000000001
#define TEST_XITA 0.5

#include <QMainWindow>
#include <QDebug.h>
#include <QDesktopWidget.h>
#include <QDialog.h>
#include <QFileDialog.h>
#include <QFile>
#include <QTextStream>
#include <QImage.h>
#include <QMainWindow>
#include <QPushButton.h>
#include <QString.h>
#include <QThread>
#include <QTimer.h>
#include <Qlabel.h>
#include <QtCore/qmath.h>
#include <cmath>
#include <ctime>

#include <cstring>
#include <opencv2/opencv.hpp>
#include <qprocess.h>
#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include <vector>
#include "./src/dataProcess.h"
#include <armadillo>
#include <QMessageBox>

using namespace arma;

using namespace cv;
using namespace std;

class DataProcess;

QT_BEGIN_NAMESPACE
namespace Ui
{
  class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
  Q_OBJECT

public slots:
  void getsendMessage(QString);
  void registerButton_clicked();
  void loginButton_clicked();
  void cleanup_clicked();

public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

signals:
  void sendXita(double);
  void sendDimension(int);
  void sendRegisterVector(vector<double>);
  void sendWorkmode(int);
  void sendChallengeVector(vector<double>);
  void sendUserID(QString);

private:
  Ui::MainWindow *ui;
  void PrintLog(const QString, QString, QString);
  void setQssFile();
  void settleRegisterThread(QString, QString);

  void PrintVector(const vector<double> &);
  void PrintVector(const vector<int> &);

  void PrintVV(const vector<vector<double>> &);
  void PrintVV(const vector<vector<int>> &);
  void PrintMatrix(const vector<vector<double>> &);
  vector<int> getNiVector(const vector<vector<double>> &);
  bool embedVVTest(const vector<vector<double>> &, const vector<int> &, double);
  bool invertibleMTest();
  vector<double> getFeatureChara(QString);
  bool fileSystemInitial();

  DataProcess *tool;

  double xita;
  int featureDimension;
  int busyFlag;
};

#endif // MAINWINDOW_H
