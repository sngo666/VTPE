#include "mainwindow.h"
#include "./ui_mainwindow.h"

#include "./src/dataProcess.h"
#include "./src/workThread.h"

using namespace cv;
using namespace std;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow)
{
  ui->setupUi(this);
  this->xita = TEST_XITA;
  // fiz mainwindow size
  this->setFixedSize(this->width(), this->height());
  // load template qss
  this->setQssFile();

  this->fileSystemInitial();

  // TODO: create connection
  connect(ui->registerButton, SIGNAL(clicked()), this, SLOT(registerButton_clicked()));
  connect(ui->loginButton, SIGNAL(clicked()), this, SLOT(loginButton_clicked()));
  connect(ui->cleanup, SIGNAL(clicked()), this, SLOT(cleanup_clicked()));

  //
  busyFlag = 0;

  // this->PrintVV(vv);
  // vector<vector<double>> a = tool->getEmbedRandomVV(vv, 0.5);
  // this->PrintVector(a[0]);
  // this->PrintVector(a[1]);

  // test embed VV model
  // vector<double> featureVector = this->getFeatureChara("./xjc2-512.txt");
  // int flag = 0;
  // int testTime = 10;
  // for (int j = 0; j < testTime; j++)
  // {
  //   vector<vector<double>> vv;
  //   tool->randomSeperate(featureVector, vv, 100);
  //   vector<int> ni = this->getNiVector(vv);
  //   vector<vector<double>> m_vv = tool->getEmbedRandomVV(vv, 0.5);

  //   if (!this->embedVVTest(m_vv, ni, 0.5))
  //     flag++;
  // }
  // this->PrintLog("[log] embedVV test " + QString::number(testTime) + " times tried,  result: " + QString::number(flag) + " times failed!", "blue", "4");
}
void MainWindow::settleRegisterThread(QString workMode, QString filePath)
{
  vector<double> featureVector = this->getFeatureChara(filePath);

  QThread *m_workerThread = new QThread();
  workThread *worker = new workThread();

  connect(m_workerThread, &QThread::started, worker, &workThread::startWork);
  connect(worker, &workThread::workFinished, m_workerThread, &QThread::quit);
  connect(m_workerThread, &QThread::finished, m_workerThread,
          &QThread::deleteLater);

  // establish connection
  connect(this, SIGNAL(sendDimension(int)), worker, SLOT(getDimension(int)));
  connect(this, SIGNAL(sendXita(double)), worker, SLOT(getXita(double)));
  connect(this, SIGNAL(sendRegisterVector(vector<double>)), worker, SLOT(getRegisterVector(vector<double>)));
  connect(this, SIGNAL(sendWorkmode(int)), worker, SLOT(getWorkMode(int)));
  connect(this, SIGNAL(sendChallengeVector(vector<double>)), worker, SLOT(getChallengeVector(vector<double>)));
  connect(this, SIGNAL(sendUserID(QString)), worker, SLOT(getUserID(QString)));

  connect(worker, SIGNAL(sendMessage(QString)), this, SLOT(getsendMessage(QString)));

  QTimer *tr = new QTimer(this);

  tr->start(10);

  connect(tr, &QTimer::timeout, [=]()
          {
    if (worker->isSent == true)
    {
      m_workerThread->quit();
      m_workerThread->wait();

      worker->isSent = false;
    } });

  if (workMode == "register")
    emit sendWorkmode(1);
  else if (workMode == "login")
    emit sendWorkmode(2);
  else
  {
    this->PrintLog("[warn] stopped! undefined workmode.", "purple", "4");
    return;
  }
  emit sendRegisterVector(featureVector);
  if (ui->userID->toPlainText().isEmpty())
  {
    this->PrintLog("[warn] stopped! put in ID plz.", "purple", "4");

    return;
  }
  emit sendUserID(ui->userID->toPlainText());
  // this->PrintLog("[log] user:" + ui->userID->toPlainText(), "green", "4");

  emit sendDimension(featureVector.size());
  emit sendXita(xita);

  worker->moveToThread(m_workerThread);

  m_workerThread->start();
}
void MainWindow::cleanup_clicked()
{
  if (QMessageBox::Yes == QMessageBox::question(nullptr, "confirm", "sure to clean register DB?"))
  {
    QString filePath = QCoreApplication::applicationDirPath() + "/usr";
    // this->PrintLog(filePath, "green", "4");

    QDir dir(filePath);

    dir.setFilter(QDir::AllEntries | QDir::NoDotAndDotDot);
    QFileInfoList fileList = dir.entryInfoList();
    foreach (QFileInfo file, fileList)
    {
      if (file.isFile())
      {
        QFile::remove(file.absoluteFilePath());
      }
    }

    return;
  }
}
void MainWindow::registerButton_clicked()
{
  QString filePath = QFileDialog::getOpenFileName(this, tr("Put in data"), ".", tr("data Files(*.txt)"));
  if (filePath.isEmpty())
  {
    this->PrintLog("[warn] none file!", "purple", "4");
    return;
  }
  this->settleRegisterThread("register", filePath);
}

// void MainWindow::registerButton_clicked()
// {
//   QString filePath1 = QFileDialog::getOpenFileName(this, tr("Put in data"), ".", tr("data Files(*.txt)"));
//   if (filePath1.isEmpty())
//   {
//     this->PrintLog("[warn] none file!", "purple", "4");
//     return;
//   }
//   // this->settleRegisterThread("register", filePath);

//   vector<double> regV = this->getFeatureChara(filePath1);

//   double sum, squareSum;
//   for (int i = 0; i < regV.size(); i++)
//   {
//     sum += regV[i];
//     squareSum += regV[i] * regV[i];
//   }
//   this->PrintLog("[log] sum = " + QString::number(sum), "purple", "4");
//   this->PrintLog("[log] squareSum = " + QString::number(squareSum), "purple", "4");
// }

void MainWindow::loginButton_clicked()
{
  QString filePath = QFileDialog::getOpenFileName(this, tr("Put in data"), ".", tr("data Files(*.txt)"));
  if (filePath.isEmpty())
  {
    this->PrintLog("[warn] none file!", "purple", "4");
    return;
  }

  this->settleRegisterThread("login", filePath);
}

void MainWindow::setQssFile()
{
  QFile file(QString("%1/qss/myqss.qss").arg(QDir::currentPath()));
  file.open(QFile::ReadOnly);
  this->PrintLog("[log] load qss...", "blue", "4");

  this->PrintLog(QStringLiteral("[log] qss loaded."), "blue", "4");

  qApp->setStyleSheet(file.readAll());
}

vector<double> MainWindow::getFeatureChara(QString filePath) // an interface to get face feature vector
{

  QString file(filePath);
  vector<double> m_registerVector;
  tool->getVectorFromFile(file, m_registerVector);
  tool = new DataProcess(m_registerVector.size());

  return m_registerVector;
}

void MainWindow::PrintLog(const QString qstr, QString color, QString size)
{
  ui->MsgBox->moveCursor(QTextCursor::End, QTextCursor::MoveAnchor);
  ui->MsgBox->append("<font size=\"" + size + "\"" + " color=\"" + color +
                     "\"" + ">" + qstr + "</font>");
}

void MainWindow::PrintVector(const vector<double> &m_v)
{
  this->PrintLog("[log] print vector--------------------", "red", "4");

  for (int i = 0; i < m_v.size(); i++)
  {
    this->PrintLog(QString::number(m_v[i], 'g', 20), "blue", "3");
  }
  this->PrintLog("[log] " + QString::number(m_v.size()) + " elements printed end---------", "red", "4");
}

void MainWindow::PrintVector(const vector<int> &m_v)
{
  this->PrintLog("[log] print vector--------------------", "red", "4");

  for (int i = 0; i < m_v.size(); i++)
  {
    this->PrintLog(QString::number(m_v[i], 'g', 15), "blue", "3");
  }
  this->PrintLog("[log] " + QString::number(m_v.size()) + " elements printed end----------", "red", "4");
}

void MainWindow::PrintVV(const vector<vector<double>> &m_v)
{
  int count = 0;
  this->PrintLog("[log] print VV--------------------", "red", "4");

  for (int i = 0; i < m_v.size(); i++)
  {
    for (int j = 0; j < m_v[i].size(); j++)
    {
      this->PrintLog(QString::number(m_v[i][j], 'g', 15), "blue", "3");
      count++;
    }

    this->PrintLog(QString::number(m_v[i].size()) + " elements in this v printed", "blue", "3");
  }
  this->PrintLog("[log] " + QString::number(m_v.size()) + " elements printed end", "red", "4");
  this->PrintLog("[log] " + QString::number(count) + "  elements in total", "red", "4");
}

void MainWindow::PrintVV(const vector<vector<int>> &m_v)
{
  int count = 0;
  this->PrintLog("[log] print VV--------------------", "red", "4");

  for (int i = 0; i < m_v.size(); i++)
  {
    for (int j = 0; j < m_v[i].size(); j++)
    {
      this->PrintLog(QString::number(m_v[i][j], 'g', 15), "blue", "3");
      count++;
    }

    this->PrintLog(QString::number(m_v[i].size()) + " elements in this v printed", "blue", "3");
  }
  this->PrintLog("[log] " + QString::number(m_v.size()) + " elements printed end", "red", "4");
  this->PrintLog("[log] " + QString::number(count) + "  elements in total", "red", "4");
}

void MainWindow::PrintMatrix(const vector<vector<double>> &m_v)
{
  this->PrintLog("[log] print VV--------------------", "red", "4");

  for (int i = 0; i < m_v.size(); i++)
  {
    QString str;
    for (int j = 0; j < m_v.size(); j++)
    {
      str += QString::number(m_v[i][j]) + " ";
    }
    this->PrintLog(str, "blue", "4");

    // this->PrintLog(QString::number(i + 1) + "th row printed---------", "blue", "3");
  }
  this->PrintLog("[log] " + QString::number(m_v.size()) + " * " + QString::number(m_v.size()) + "  elements in total", "red", "4");
}

vector<int> MainWindow::getNiVector(const vector<vector<double>> &m_vv)
{
  vector<int> ni;
  for (int i = 0; i < m_vv.size(); i++)
  {
    ni.push_back(m_vv[i].size());
  }
  return ni;
}
bool MainWindow::embedVVTest(const vector<vector<double>> &m_vv, const vector<int> &ni, double xita)
{
  double mSum = 0, xitaSum = 0;
  bool flag = true;
  for (int i = 0; i < ni.size(); i++)
  {
    mSum += m_vv[1][i];
    xitaSum += m_vv[0][i];
    if (abs(m_vv[1][i] + m_vv[0][i]) <= ni[i])
    {
      this->PrintLog("[log] " + QString::number(i) + "th elements out of range!", "red", "4");

      flag = false;
    }
    this->PrintLog("[log] xita_i + m_i =" + QString::number(m_vv[1][i] + m_vv[0][i]) + ", while n_i =" + QString::number(ni[i]), "red", "4");

    // else
    // this->PrintLog("[log] " + QString::number(i) + "th elements in range!", "red", "4");
  }

  this->PrintLog("[log] mSum = " + QString::number(mSum), "red", "4");
  this->PrintLog("[log] xitaSum = " + QString::number(xitaSum, (char)103, 8), "red", "4");

  if (mSum < (0 - DEVIATION_VALUE) || mSum > (0 + DEVIATION_VALUE))
  {
    flag = false;
  }
  if (xitaSum < (xita - DEVIATION_VALUE) || xitaSum > (xita + DEVIATION_VALUE))
  {
    flag = false;
  }

  return flag;
}
bool MainWindow::fileSystemInitial()
{

  QString IDLibraryPath = QCoreApplication::applicationDirPath() + "/usr"; // remember to get absolute dir
  QDir dir(IDLibraryPath);
  if (dir.exists())
  {
    return false;
  }
  else
  {
    bool flag = dir.mkdir(IDLibraryPath);
    return flag;
  }
}
void MainWindow::getsendMessage(QString msg)
{
  this->PrintLog(msg, "green", "3");
}

bool MainWindow::invertibleMTest()
{
  int count = 0;
  for (int i = 0; i < 500; i++)
  {
    vector<vector<double>> m_vv;
    tool->getInvertibleMatrix(m_vv, 15);
    mat m_m(m_vv.size(), m_vv.size());
    m_m = tool->vvToMatrixXd(m_vv);
    this->PrintLog(QString::number(m_m(0, 0)), "blue", "3");

    if (det(m_m) == 0)
      count++;
  }
  this->PrintLog(QString::number(count) + " random matrix's determinant is zero !", "blue", "3");

  return (count == 0);
}

MainWindow::~MainWindow()
{
  delete ui;
}
