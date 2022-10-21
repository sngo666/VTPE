#include "dataProcess.h"

DataProcess::DataProcess(int n)
{
  scale = n;
}

void DataProcess::randomSeperate(const vector<double> &inputV, vector<vector<double>> &outputVV, int num)
{
  int scale = inputV.size();

  vector<int> m_segNum = this->getRandomOrder(inputV.size() - 1);
  vector<int> segNum;

  for (int i = 0; i < num - 1; i++)
  {
    segNum.push_back(m_segNum[i]);
  }
  sort(segNum.begin(), segNum.end());

  int count = 0;
  for (int i = 0; i < num - 1; i++)
  {
    vector<double> m_v;
    for (; count <= segNum[i]; count++)
    {
      m_v.push_back(inputV[count]);
    }
    outputVV.push_back(m_v);
  }

  vector<double> m_v;
  while (count <= inputV.size() - 1)
  {
    m_v.push_back(inputV[count]);
    count++;
  }
  if (m_v.size() != 0)
    outputVV.push_back(m_v);
}

// void DataProcess::randomSeperate(const vector<double> &inputV, vector<vector<double>> &outputVV, int num)
// {
//   int scale = inputV.size();
//   int count = 0;
//   vector<int> testCount;

//   vector<int> segNum;
//   for (int i = 0; i < num - 1; i++)
//   {

//     int ran = (int)(QRandomGenerator::global()->bounded(1, scale - num + i));
//     segNum.push_back(ran);
//     scale -= ran;
//   }
//   segNum.push_back(scale);

//   vector<int> order = this->getRandomOrder(num);
//   for (int i = 0; i < num; i++)
//   {
//     vector<double> m_v;
//     // m_segNum[order[i]] = segNum[i];
//     for (int j = 0; j < segNum[order[i]]; j++)
//     {
//       m_v.push_back(inputV[count]);
//       // testCount.push_back(count);
//       count++;
//     }
//     outputVV.push_back(m_v);
//   }
// }

unsigned long long DataProcess::getMatrixHash(const vector<vector<double>> &m_vv, int p, int mod)
{
  QString qstr = "";

  for (int i = 0; i < m_vv.size(); i++)
  {
    if (m_vv.size() != m_vv[i].size())
      return 0;
    for (int j = 0; j < m_vv[i].size(); j++)
    {
      qstr += QString::number(m_vv[i][j], (char)103, 15);
    }
  }
  string str = qstr.toStdString();
  vector<unsigned long long> hashArray;
  hashArray.push_back((int)str[0]);
  for (int i = 1; i < str.size(); i++)
  {
    hashArray.push_back(hashArray[i - 1] * p + (int)str[i] % mod);
  }
  return hashArray[hashArray.size() - 1];
}
unsigned long long DataProcess::getMatrixHash(const vector<vector<int>> &m_vv, int p, int mod) // p & mod must be Prime,p < mod,and make p & mod large as you can.
{
  QString qstr = "";

  for (int i = 0; i < m_vv.size(); i++)
  {
    if (m_vv.size() != m_vv[i].size())
      return 0;
    for (int j = 0; j < m_vv[i].size(); j++)
    {

      qstr += QString::number(m_vv[i][j]);
    }
  }
  string str = qstr.toStdString();
  vector<unsigned long long> hashArray;
  hashArray.push_back((int)str[0]);
  for (int i = 1; i < str.size(); i++)
  {
    hashArray.push_back(hashArray[i - 1] * p + (int)str[i] % mod);
  }
  return hashArray[hashArray.size() - 1];
}

unsigned long long DataProcess::getVectorHash(const vector<double> &m_v, int p, int mod)
{
  QString qstr = "";
  for (int i = 0; i < m_v.size(); i++)
  {
    qstr += QString::number(m_v[i]);
  }
  vector<double> hashArray;
  string str = qstr.toStdString();

  hashArray.push_back((double)str[0]);
  for (int i = 1; i < str.size(); i++)
  {
    hashArray.push_back(hashArray[i - 1] * p + (int)str[i] % mod);
  }
  return hashArray[hashArray.size() - 1];
}

unsigned long long DataProcess::getVectorHash(const vector<int> &m_v, int p, int mod) // p & mod must be Prime,p < mod,and make p & mod large as you can.
{
  QString qstr = "";
  for (int i = 0; i < m_v.size(); i++)
  {
    qstr += QString::number(m_v[i]);
  }
  vector<unsigned long long> hashArray;
  string str = qstr.toStdString();

  hashArray.push_back((int)str[0]);
  for (int i = 1; i < str.size(); i++)
  {
    hashArray.push_back(hashArray[i - 1] * p + (int)str[i] % mod);
  }
  return hashArray[hashArray.size() - 1];
}

vector<int> DataProcess::getSign(const vector<double> &m_v1, const vector<double> &m_v2)
{
  vector<int> res;
  if (m_v1.size() != m_v2.size())
    return res;
  for (int i = 0; i < m_v1.size(); i++)
  {
    if (m_v1[i] + m_v2[i] > 0)
      res.push_back(1);
    else
      res.push_back(-1);
  }
  return res;
}

vector<int> DataProcess::getSign(const vector<double> &m_v, int num)
{
  vector<int> res;
  for (int i = 0; i < m_v.size(); i++)
  {
    if (m_v[i] > 0)
      res.push_back(1 * num);
    else
      res.push_back(-1 * num);
  }
  return res;
}

void DataProcess::getVectorFromFile(QString filePath, vector<double> &m_v)
{
  QFileInfo fileinfo(filePath);
  if (!fileinfo.isFile())
    return;

  QFile file(filePath);
  QString strAll;
  QStringList strList;
  if (file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QTextCodec *codec = QTextCodec::codecForName("utf8");
    QTextStream stream(&file);
    stream.setCodec(codec);
    strAll = stream.readAll();
  };

  strList = strAll.split('\n');
  foreach (QString str, strList)
  {
    if (!str.isEmpty())
      m_v.push_back(str.toDouble());
    // m_v.push_back(this->strTodouble(str));
  }
  file.close();
}

unsigned long long DataProcess::strToULongLong(QString qstr)
{
  string str = qstr.toStdString();
  unsigned long long res = str[0] - '0';
  for (int i = 1; i < str.size(); i++)
  {
    int temp = (int)(str[i] - '0');
    res = res * 10 + temp;
  }

  return res;
}

// Prevent loss of accuracy
double DataProcess::strTodouble(QString str)
{
  if (str.isEmpty())
    return 0.0;
  string s;
  double num;
  double tempSign;
  if (str[0] == '-')
  {
    s = str.mid(1).toStdString();
    num = (double)(-1 * ((char)s[0] - '0'));
    tempSign = -1;
  }
  else
  {
    s = str.toStdString();
    num = (double)((char)s[0] - '0');
    tempSign = 1;
  }

  int flag = 0;
  double sign, magn = 0;

  for (int i = 1; i < s.size(); i++)
  {
    if (s[i] == 'e')
    {
      sign = (s[i + 1] == '+') ? 10 : 0.1;
      flag = i + 2;
      break;
    }

    if (s[i] == '.')
    {
      flag = 1;
      continue;
    }
    else if (flag == 0)
    {
      num = num * 10 + tempSign * (double)((char)s[i] - '0');
    }
    else
    {

      num += tempSign * (double)((char)s[i] - '0') * pow(0.1, flag);
      flag++;
    }
  }

  for (int i = flag; i < s.size(); i++)
  {
    if (s[i] == '0')
      continue;
    magn = magn * pow(10, (i - flag)) + ((char)s[i] - '0');
  }
  return num * pow(sign, magn);
}

void DataProcess::getInvertibleMatrix(vector<vector<double>> &m_vv, int scale) // put in a null vv plz,not an identity M
{
  this->getIdentityMatrix(m_vv, scale);
  int n = m_vv.size();
  int i = 0;
  int j = 0;
  int pattern = 0;
  int randMagn, randRow, randTempRow, flag = 0;

  vector<double> vTemp(n);
  vector<int> vState(n, 0);

  int transformTime = (int)(QRandomGenerator::global()->bounded(0, 501) + 500); // time of operate: 500~1000 times
  for (i = 0; i < transformTime; ++i)
  {
    randRow = QRandomGenerator::global()->bounded(0, n);
    // operate row : 0~n-1
    pattern = QRandomGenerator::global()->bounded(1, 4); // 1~3
    if (vState[randRow] == 0)
      switch (pattern)
      {
      case 1: // + -
      {
        randTempRow = (int)(QRandomGenerator::global()->bounded(0, n)); // 0~n-1
        for (j = 0; j < n; j++)
        {
          vTemp[j] = m_vv[randRow][j] + m_vv[randTempRow][j];
          if (vTemp[j] > scale || vTemp[j] < (-1) * scale)
            vState[randRow] = 1;
        }
        m_vv[randRow] = vTemp;

        break;
      }
      case 2: // swap
      {
        randTempRow = (int)(QRandomGenerator::global()->bounded(0, n)); // 0~n-1
        if (randTempRow != randTempRow)
        {
          vector<double> swapV = m_vv[randRow];
          m_vv[randRow] = m_vv[randTempRow];
          m_vv[randTempRow] = swapV;
        }
        break;
      }
      case 3: // */=
      {
        randMagn = (int)(QRandomGenerator::global()->bounded(1, 10)) * (QRandomGenerator::global()->bounded(0, 2) ? (1) : (-1)); // 1~3 / -1~-3
        for (j = 0; j < n; j++)
        {
          m_vv[randRow][j] *= randMagn;
          if (m_vv[randRow][j] > scale || m_vv[randRow][j] < (-1) * scale)
            vState[randRow] = 1;
        }
        break;
      }
      }
  }
}

vector<vector<double>> DataProcess::getInvertibleMatrix(int scale)
{
  vector<vector<double>> m_vv;
  this->getIdentityMatrix(m_vv, scale);
  int n = m_vv.size();
  int i = 0;
  int j = 0;
  int pattern = 0;
  int randMagn, randRow, randTempRow, flag = 0;

  vector<double> vTemp(n);
  vector<int> vState(n, 0);

  int transformTime = (int)(QRandomGenerator::global()->bounded(100, 301) + 200); // time of operate: 300~500 times
  for (i = 0; i < transformTime; ++i)
  {
    randRow = QRandomGenerator::global()->bounded(0, n);
    // operate row : 0~n-1
    pattern = QRandomGenerator::global()->bounded(1, 4); // 1~3
    if (vState[randRow] == 0)
      switch (pattern)
      {
      case 1: // + -
      {
        randTempRow = (int)(QRandomGenerator::global()->bounded(0, n)); // 0~n-1
        for (j = 0; j < n; j++)
        {
          vTemp[j] = m_vv[randRow][j] + m_vv[randTempRow][j];
          if (vTemp[j] > scale || vTemp[j] < (-1) * scale)
            vState[randRow] = 1;
        }
        m_vv[randRow] = vTemp;

        break;
      }
      case 2: // swap
      {
        randTempRow = (int)(QRandomGenerator::global()->bounded(0, n)); // 0~n-1
        if (randTempRow != randTempRow)
        {
          vector<double> swapV = m_vv[randRow];
          m_vv[randRow] = m_vv[randTempRow];
          m_vv[randTempRow] = swapV;
        }
        break;
      }
      case 3: // */=
      {
        randMagn = (int)(QRandomGenerator::global()->bounded(1, 4)) * (QRandomGenerator::global()->bounded(0, 2) ? (1) : (-1)); // 1~3 / -1~-3
        for (j = 0; j < n; j++)
        {
          m_vv[randRow][j] *= randMagn;
          if (m_vv[randRow][j] > scale || m_vv[randRow][j] < (-1) * scale)
            vState[randRow] = 1;
        }
        break;
      }
      }
  }
  return m_vv;
}

void DataProcess::getIdentityMatrix(vector<vector<double>> &m_vv, int n)
{
  int r = 0;
  int c = 0;

  for (r = 0; r < n; ++r)
  {
    vector<double> m_v;
    for (c = 0; c < n; ++c)
    {
      if (r == c)
        m_v.push_back(1);
      else
        m_v.push_back(0);
    }
    m_vv.push_back(m_v);
  }
}

vector<int> DataProcess::getNiVector(const vector<vector<double>> &m_vv)
{
  vector<int> resV;
  for (int i = 0; i < m_vv.size(); i++)
    resV.push_back(m_vv[i].size());
  return resV;
}

vector<vector<double>> DataProcess::getDiagonalMatrix(const vector<double> &m_v, const vector<int> &randomFun)
{
  vector<vector<double>> resVV;
  vector<double> operateV = m_v;
  int num = m_v.size();

  vector<double> ranV = this->fixVectorOrder(m_v, randomFun);
  for (int i = 0; i < num; i++)
  {
    vector<double> row_v;
    for (int j = 0; j < num; j++)
    {
      if (i == j)
        row_v.push_back(ranV[i]);
      else
        row_v.push_back(0);
    }
    resVV.push_back(row_v);
  }
  return resVV;
}

vector<vector<double>> DataProcess::getDiagonalMatrix(const vector<double> &m_v)
{
  vector<vector<double>> resVV;
  vector<double> operateV = m_v;
  int num = m_v.size();

  for (int i = 0; i < num; i++)
  {
    vector<double> row_v;
    for (int j = 0; j < num; j++)
    {
      if (i == j)
        row_v.push_back(operateV[i]);
      else
        row_v.push_back(0);
    }
    resVV.push_back(row_v);
  }
  return resVV;
}

vector<vector<double>> DataProcess::getRandomLowerTriangleMatrix(int size)
{
  vector<vector<double>> resVV;

  for (int i = 0; i < size; i++)
  {
    vector<double> row_v;
    for (int j = 0; j < size; j++)
    {
      if (i == j)
        row_v.push_back(1);
      else if (i > j)
        row_v.push_back((double)(QRandomGenerator::global()->bounded(0, 20)));
      else if (i < j)
        row_v.push_back(0);
    }
    resVV.push_back(row_v);
  }
  return resVV;
}

vector<int> DataProcess::getRandomOrder(int n) // start with 0
{

  vector<int> order(n);
  for (int i = 0; i < n; i++)
  {
    order[i] = i;
  }

  for (int i = 0; i < n; i++)
  {
    int j = QRandomGenerator::global()->bounded(0, n); // 0~n-1
    int temp = order[i];
    order[i] = order[j];
    order[j] = temp;
  }
  return order;
}

vector<int> DataProcess::getRandomOrderVector(const vector<int> &m_v)
{
  vector<int> order = this->getRandomOrder(m_v.size());

  vector<int> m_m_v;

  for (int i = 0; i < m_v.size(); i++)
  {
    m_m_v.push_back(m_v[order[i]]);
  }
  return m_m_v;
}

vector<double> DataProcess::getRandomOrderVector(const vector<double> &m_v)
{

  vector<int> order = this->getRandomOrder(m_v.size());

  vector<double> m_m_v;

  for (int i = 0; i < m_v.size(); i++)
  {
    m_m_v.push_back(m_v[order[i]]);
  }
  return m_m_v;
}

vector<double> DataProcess::fixVectorOrder(const vector<double> &m_v, const vector<int> &order)
{

  vector<double> m_m_v;

  for (int i = 0; i < m_v.size(); i++)
  {
    m_m_v.push_back(m_v[order[i]]);
  }
  return m_m_v;
}

vector<int> DataProcess::getFixedSumRandomIntVector(int sum, int num)
{

  vector<int> m_v;
  int m_sum;
  for (int i = 0; i < num - 1; i++)
  {
    int ran = (QRandomGenerator::global()->bounded(0, sum / num * 3 + 1)) * (QRandomGenerator::global()->bounded(0, 2) ? 1 : (-1));
    m_v.push_back(ran);
    m_sum += ran;
  }
  m_v.push_back((-1) * m_sum);
  this->getRandomOrderVector(m_v);
  return m_v;
}

vector<int> DataProcess::getNonzeroRandomVector(int num)
{

  vector<int> m_v;
  for (int i = 0; i < num; i++)
  {
    int ran = (QRandomGenerator::global()->bounded(1, 20)) * (QRandomGenerator::global()->bounded(2) ? 1 : (-1));
    m_v.push_back(ran);
  }
  return m_v;
}

vector<double> DataProcess::getFixedSumRandomDoubleVector(double sum, int num)
{

  vector<double> m_v;
  double m_sum;

  if (num <= 1)
  {
    m_v.push_back(sum);
    return m_v;
  }

  for (int i = 0; i < num - 1; i++)
  {
    double ran = (double)(((int)QRandomGenerator::global()->bounded(1001)) / (1000.0)) + QRandomGenerator::global()->bounded(10) * ((QRandomGenerator::global()->bounded(2) == 1) ? 1.0 : (-1.0));
    // 0~0.999 + a integer num * random sign
    m_v.push_back(ran);
    m_sum += ran;
  }
  m_v.push_back(sum - m_sum);
  // this->getRandomOrderVector(m_v);
  return m_v;
}

vector<vector<double>> DataProcess::getEmbedRandomVV(const vector<int> &Ni, double xita)
{
  int num = Ni.size();
  vector<double> vecXita = this->getFixedSumRandomDoubleVector(xita, num);
  vector<vector<double>> embedVV, mRange;
  vector<double> mRang_more, mRang_less, mVec;
  double offset = 0;

  embedVV.push_back(vecXita);

  vector<double> absSum;
  for (int i = 0; i < num; i++)
  {
    mRang_more.push_back((-1) * vecXita[i] + Ni[i]);
    mRang_less.push_back((-1) * vecXita[i] - Ni[i]);
  }

  for (int i = 0; i < num - 1; i++)
  {
    int pattern = QRandomGenerator::global()->bounded(0, 2); // 0~1
    if (i == 0)
      pattern = 0;
    if (i == 1)
      pattern = 1;
    switch (pattern)
    {
    case 0: // get radom double fromm Rang_more
    {
      double m = mRang_more[i] + (double)(QRandomGenerator::global()->bounded(1, 1001) / ((double)1000) + (double)(QRandomGenerator::global()->bounded(0, 3)));
      offset += m;
      mVec.push_back(m);
      break;
    }
    case 1: // get radom double fromm Rang_less
    {
      double m = mRang_less[i] - (double)(QRandomGenerator::global()->bounded(1, 1001) / ((double)1000) + (double)(QRandomGenerator::global()->bounded(0, 3)));
      offset += m;
      mVec.push_back(m);
      break;
    }
    }
  }
  mVec.push_back((-1) * offset);

  if (mVec[num - 1] > mRang_less[num - 1] && mVec[num - 1] < mRang_more[num - 1])
  {
    mVec[0] += mRang_more[num - 1] + 1 - mVec[num - 1];
    mVec[num - 1] = mRang_more[num - 1] + 1;
  }
  double tempSum = 0.0;
  for (int i = 0; i < num; i++)
  {
    tempSum += mVec[i];
  }
  if (tempSum > 0)
    mVec[1] -= tempSum;
  else if (tempSum < 0)
    mVec[0] += tempSum;

  embedVV.push_back(mVec);
  return embedVV;
}

vector<vector<double>> DataProcess::getEmbedRandomVV(const vector<vector<double>> &m_vv, double xita)
{
  int num = m_vv.size();
  vector<double> vecXita = this->getFixedSumRandomDoubleVector(xita, num);
  vector<vector<double>> embedVV, mRange;
  vector<double> mRang_more, mRang_less, mVec;
  double offset = 0;

  embedVV.push_back(vecXita);

  vector<double> absSum;
  for (int i = 0; i < num; i++)
  {
    mRang_more.push_back((-1) * vecXita[i] + m_vv[i].size());
    mRang_less.push_back((-1) * vecXita[i] - m_vv[i].size());
  }

  for (int i = 0; i < num - 1; i++)
  {
    int pattern = QRandomGenerator::global()->bounded(0, 2); // 0~1
    if (i == 0)
      pattern = 0;
    if (i == 1)
      pattern = 1;
    switch (pattern)
    {
    case 0: // get radom double fromm Rang_more
    {
      double m = mRang_more[i] + (double)(QRandomGenerator::global()->bounded(1, 1001) / ((double)1000) + (double)(QRandomGenerator::global()->bounded(0, 3)));
      offset += m;
      mVec.push_back(m);
      break;
    }
    case 1: // get radom double fromm Rang_less
    {
      double m = mRang_less[i] - (double)(QRandomGenerator::global()->bounded(1, 1001) / ((double)1000) + (double)(QRandomGenerator::global()->bounded(0, 3)));
      offset += m;
      mVec.push_back(m);
      break;
    }
    }
  }
  mVec.push_back((-1) * offset);

  if (mVec[num - 1] > mRang_less[num - 1] && mVec[num - 1] < mRang_more[num - 1])
  {
    mVec[0] += mRang_more[num - 1] + 1 - mVec[num - 1];
    mVec[num - 1] = mRang_more[num - 1] + 1;
  }
  double tempSum = 0.0;
  for (int i = 0; i < num; i++)
  {
    tempSum += mVec[i];
  }
  if (tempSum > 0)
    mVec[1] -= tempSum;
  else if (tempSum < 0)
    mVec[0] += tempSum;

  embedVV.push_back(mVec);
  return embedVV;
}

mat DataProcess::vvToMatrixXd(const vector<vector<double>> &m_vv)
{
  mat m_m(m_vv.size(), m_vv.size(), fill::zeros);

  for (int i = 0; i < m_vv.size(); i++)
  {
    for (int j = 0; j < m_vv.size(); j++)
      m_m(i, j) = m_vv[i][j];
  }
  return m_m;
}

vector<vector<double>> DataProcess::MatrixXdToVV(mat m_m)
{
  vector<vector<double>> m_vv;

  for (int i = 0; i < m_m.n_rows; i++)
  {
    vector<double> m_v;
    for (int j = 0; j < m_m.n_cols; j++)
    {
      m_v.push_back(m_m(i, j));
    }
    m_vv.push_back(m_v);
  }
  return m_vv;
}

bool DataProcess::echoHashVectortoFile(const vector<unsigned long long> &m_v, QString filePath, QString descriptors)
{
  QFile file(filePath);
  if (!file.exists())
    return false;
  file.open(QIODevice::Append | QIODevice::Text);
  if (file.isOpen())
  {
    QTextStream in(&file);
    in << descriptors + QString(":") << Qt::endl;
    for (int i = 0; i < m_v.size(); i++)
    {
      in << QString::number(m_v[i]) << Qt::endl;
    }

    file.close();
    return true;
  }
  else
    return false;
}

bool DataProcess::echoVectortoFile(const vector<double> &m_v, QString filePath, QString descriptors)
{
  QFile file(filePath);
  if (!file.exists())
    return false;
  file.open(QIODevice::Append | QIODevice::Text);
  if (file.isOpen())
  {
    QTextStream in(&file);
    in << descriptors + QString(":") << Qt::endl;
    for (int i = 0; i < m_v.size(); i++)
    {
      in << QString::number(m_v[i]) << Qt::endl;
    }

    file.close();
    return true;
  }
  else
    return false;
}

bool DataProcess::echoVectortoFile(const vector<int> &m_v, QString filePath, QString descriptors)
{
  QFile file(filePath);
  if (!file.exists())
    return false;
  file.open(QIODevice::Append | QIODevice::Text);
  if (file.isOpen())
  {
    QTextStream in(&file);
    in << descriptors + QString(":") << Qt::endl;
    for (int i = 0; i < m_v.size(); i++)
    {
      in << QString::number(m_v[i]) << Qt::endl;
    }

    file.close();
    return true;
  }
  else
    return false;
}
bool DataProcess::echoHashVtoFile(const vector<unsigned long long> &m_v, QString filePath, QString descriptors)
{
  QFile file(filePath);
  if (!file.exists())
    return false;
  file.open(QIODevice::Append | QIODevice::Text);
  if (file.isOpen())
  {
    QTextStream in(&file);
    in << descriptors + QString(":") << Qt::endl;
    for (int i = 0; i < m_v.size(); i++)
    {
      in << QString::number(m_v[i]) << Qt::endl;
    }

    file.close();
    return true;
  }
  else
    return false;
}
bool DataProcess::echoVVtoFile(const vector<vector<double>> &m_vv, QString filePath, QString descriptors)
{
  QFile file(filePath);
  if (!file.exists())
    return false;
  file.open(QIODevice::Append | QIODevice::Text);
  if (file.isOpen())
  {
    QTextStream in(&file);
    in << descriptors + QString(":") << Qt::endl;
    for (int i = 0; i < m_vv.size(); i++)
    {
      for (int j = 0; j < m_vv[i].size(); j++)
        in << QString::number(m_vv[i][j]) << Qt::endl;
    }

    file.close();
    return true;
  }
  else
    return false;
}

bool DataProcess::echoVVtoFile(const vector<vector<int>> &m_vv, QString filePath, QString descriptors)
{
  QFile file(filePath);
  if (!file.exists())
    return false;
  file.open(QIODevice::Append | QIODevice::Text);
  if (file.isOpen())
  {
    QTextStream in(&file);
    in << descriptors + QString(":") << Qt::endl;
    for (int i = 0; i < m_vv.size(); i++)
    {
      for (int j = 0; j < m_vv[i].size(); j++)
        in << QString::number(m_vv[i][j]) << Qt::endl;
    }

    file.close();
    return true;
  }
  else
    return false;
}

bool DataProcess::echoQStrtoFile(const QString m_qstr, QString filePath, QString descriptors)
{
  QFile file(filePath);
  if (!file.exists())
    return false;
  file.open(QIODevice::Append | QIODevice::Text);
  if (file.isOpen())
  {
    QTextStream in(&file);
    in << descriptors + QString(":") << Qt::endl;
    in << m_qstr << Qt::endl;

    file.close();
    return true;
  }
  else
    return false;
}

vector<vector<double>> DataProcess::getDoubleVVfromFile(vector<int> &Ni, QString filePath, QString descriptors)
{
  vector<vector<double>> m_vv;

  QFileInfo fileinfo(filePath);
  if (!fileinfo.isFile())
    return m_vv;

  QFile file(filePath);
  QString strAll;
  QStringList strList;
  if (file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QTextCodec *codec = QTextCodec::codecForName("utf8");
    QTextStream stream(&file);
    stream.setCodec(codec);
    strAll = stream.readAll();
  };

  strList = strAll.split('\n');
  int flag = 0;
  int count = -1;
  for (int i = 0; i < strList.size(); i++)
  {
    QString str = strList[i];
    if (str.size() - 1 == descriptors.size())
    {
      for (int j = 0; j < descriptors.size(); j++)
      {
        if (str[j] != descriptors[j])
          count = -1;
      }
      count = i + 1;
    }
  }

  if (count != -1)
  {
    for (int i = 0; i < Ni.size(); i++)
    {
      vector<double> m_v;
      for (int j = 0; j < Ni[i]; j++)
      {
        m_v.push_back(strList[count].toDouble());
        count++;
      }
      m_vv.push_back(m_v);
    }
  }
  file.close();
  return m_vv;
}

vector<vector<double>> DataProcess::getDoubleVVfromFile(int m_n, QString filePath, QString descriptors)
{
  vector<vector<double>> m_vv;

  QFileInfo fileinfo(filePath);
  if (!fileinfo.isFile())
    return m_vv;

  QFile file(filePath);
  QString strAll;
  QStringList strList;
  if (file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QTextCodec *codec = QTextCodec::codecForName("utf8");
    QTextStream stream(&file);
    stream.setCodec(codec);
    strAll = stream.readAll();
  };

  strList = strAll.split('\n');
  int flag = 0;
  int count = -1;
  for (int i = 0; i < strList.size(); i++)
  {
    QString str = strList[i];
    if (str.size() - 1 == descriptors.size())
    {
      for (int j = 0; j < descriptors.size(); j++)
      {
        if (str[j] != descriptors[j])
          count = -1;
      }
      count = i + 1;
    }
  }

  if (count != -1)
  {
    for (int i = 0; i < m_n; i++)
    {
      vector<double> m_v;
      for (int j = 0; j < m_n; j++)
      {
        m_v.push_back(strList[count].toDouble());
        count++;
      }
      m_vv.push_back(m_v);
    }
  }
  file.close();
  return m_vv;
}

vector<vector<int>> DataProcess::getIntegerVVfromFile(vector<int> &Ni, QString filePath, QString descriptors)
{
  vector<vector<int>> m_vv;

  QFileInfo fileinfo(filePath);
  if (!fileinfo.isFile())
    return m_vv;

  QFile file(filePath);
  QString strAll;
  QStringList strList;
  if (file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QTextCodec *codec = QTextCodec::codecForName("utf8");
    QTextStream stream(&file);
    stream.setCodec(codec);
    strAll = stream.readAll();
  };

  strList = strAll.split('\n');
  int flag = 0;
  int count = -1;
  for (int i = 0; i < strList.size(); i++)
  {
    QString str = strList[i];
    if (str.size() - 1 == descriptors.size())
    {
      for (int j = 0; j < descriptors.size(); j++)
      {
        if (str[j] != descriptors[j])
          count = -1;
      }
      count = i + 1;
    }
  }

  if (count != -1)
  {
    for (int i = 0; i < Ni.size(); i++)
    {
      vector<int> m_v;
      for (int j = 0; j < Ni[i]; j++)
      {
        m_v.push_back(strList[count].toInt());
        count++;
      }
      m_vv.push_back(m_v);
    }
  }
  file.close();
  return m_vv;
}

vector<vector<int>> DataProcess::getIntegerVVfromFile(int m_n, QString filePath, QString descriptors)
{
  vector<vector<int>> m_vv;

  QFileInfo fileinfo(filePath);
  if (!fileinfo.isFile())
    return m_vv;

  QFile file(filePath);
  QString strAll;
  QStringList strList;
  if (file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QTextCodec *codec = QTextCodec::codecForName("utf8");
    QTextStream stream(&file);
    stream.setCodec(codec);
    strAll = stream.readAll();
  };

  strList = strAll.split('\n');
  int flag = 0;
  int count = -1;
  for (int i = 0; i < strList.size(); i++)
  {
    QString str = strList[i];
    if (str.size() - 1 == descriptors.size())
    {
      for (int j = 0; j < descriptors.size(); j++)
      {
        if (str[j] != descriptors[j])
          count = -1;
      }
      count = i + 1;
    }
  }

  if (count != -1)
  {
    for (int i = 0; i < m_n; i++)
    {
      vector<int> m_v;
      for (int j = 0; j < m_n; j++)
      {
        m_v.push_back(strList[count].toInt());
        count++;
      }
      m_vv.push_back(m_v);
    }
  }
  file.close();
  return m_vv;
}

vector<int> DataProcess::getIntegerVectorfromFile(int m_n, QString filePath, QString descriptors)
{
  vector<int> m_v;

  QFileInfo fileinfo(filePath);
  if (!fileinfo.isFile())
    return m_v;

  QFile file(filePath);
  QString strAll;
  QStringList strList;
  if (file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QTextCodec *codec = QTextCodec::codecForName("utf8");
    QTextStream stream(&file);
    stream.setCodec(codec);
    strAll = stream.readAll();
  };

  strList = strAll.split('\n');
  int flag = 0;
  int count = -1;
  for (int i = 0; i < strList.size(); i++)
  {
    QString str = strList[i];
    if (str.size() - 1 == descriptors.size())
    {
      for (int j = 0; j < descriptors.size(); j++)
      {
        if (str[j] != descriptors[j])
          count = -1;
      }
      count = i + 1;
    }
  }

  if (count != -1)
  {
    for (int i = 0; i < m_n; i++)
    {
      m_v.push_back(strList[count].toInt());
      count++;
    }
  }
  file.close();
  return m_v;
}

vector<double> DataProcess::getDoubleVectorfromFile(int m_n, QString filePath, QString descriptors)
{
  vector<double> m_v;

  QFileInfo fileinfo(filePath);
  if (!fileinfo.isFile())
    return m_v;

  QFile file(filePath);
  QString strAll;
  QStringList strList;
  if (file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QTextCodec *codec = QTextCodec::codecForName("utf8");
    QTextStream stream(&file);
    stream.setCodec(codec);
    strAll = stream.readAll();
  };

  strList = strAll.split('\n');
  int flag = 0;
  int count = -1;
  for (int i = 0; i < strList.size(); i++)
  {
    QString str = strList[i];
    if (str.size() - 1 == descriptors.size())
    {
      for (int j = 0; j < descriptors.size(); j++)
      {
        if (str[j] != descriptors[j])
          count = -1;
      }
      count = i + 1;
    }
  }

  if (count != -1)
  {
    for (int i = 0; i < m_n; i++)
    {
      m_v.push_back(strList[count].toDouble());
      count++;
    }
  }
  file.close();
  return m_v;
}
vector<unsigned long long> DataProcess::getHashVectorfromFile(int m_n, QString filePath, QString descriptors)
{
  vector<unsigned long long> m_v;

  QFileInfo fileinfo(filePath);
  if (!fileinfo.isFile())
    return m_v;

  QFile file(filePath);
  QString strAll;
  QStringList strList;
  if (file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QTextCodec *codec = QTextCodec::codecForName("utf8");
    QTextStream stream(&file);
    stream.setCodec(codec);
    strAll = stream.readAll();
  };

  strList = strAll.split('\n');
  int flag = 0;
  int count = -1;
  for (int i = 0; i < strList.size(); i++)
  {
    QString str = strList[i];
    if (str.size() - 1 == descriptors.size())
    {
      for (int j = 0; j < descriptors.size(); j++)
      {
        if (str[j] != descriptors[j])
          count = -1;
      }
      count = i + 1;
    }
  }

  if (count != -1)
  {
    for (int i = 0; i < m_n; i++)
    {
      m_v.push_back(this->strToULongLong(strList[count]));
      count++;
    }
  }
  file.close();
  return m_v;
}

QString DataProcess::getQstrFromFile(QString filePath, QString descriptors)
{
  QFileInfo fileinfo(filePath);
  if (!fileinfo.isFile())
    return "";

  QFile file(filePath);
  QString strAll;
  QStringList strList;
  if (file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QTextCodec *codec = QTextCodec::codecForName("utf8");
    QTextStream stream(&file);
    stream.setCodec(codec);
    strAll = stream.readAll();
  };

  strList = strAll.split('\n');
  int flag = 0;
  int count = -1;
  for (int i = 0; i < strList.size(); i++)
  {
    QString str = strList[i];
    if (str.size() - 1 == descriptors.size())
    {
      for (int j = 0; j < descriptors.size(); j++)
      {
        if (str[j] != descriptors[j])
          count = -1;
      }
      count = i;
    }
  }
  if (count != -1)
  {
    file.close();
    return strList[++count];
  }
  file.close();
  return "";
}

QString DataProcess::getFileFolderPath(const QString qstr)
{
  QString Folder;
  if (qstr.isEmpty())
    return "";
  int ptr = qstr.size();
  for (int i = 0; i < qstr.size(); i++)
  {
    if (qstr[i] == '/')
      ptr = i;
  }
  for (int i = 0; i < ptr; i++)
  {
    Folder.append(qstr[i]);
  }
  return Folder;
}

vector<vector<double>> DataProcess::matrixMultiplication(const vector<vector<double>> &vv1, const vector<vector<double>> &vv2)
{
  vector<vector<double>> res_VV;
  if (vv1.size() != vv2[0].size() || vv1.size() == 0 || vv2.size() == 0)
    return res_VV;

  for (int i = 0; i < vv1.size(); i++)
  {

    vector<double> m_v;
    for (int j = 0; j < vv2[0].size(); j++)
    {
      double temp = 0;
      for (int k = 0; k < vv1[0].size(); k++)
      {
        temp += vv1[i][k] * vv2[k][j];
      }
      m_v.push_back(temp);
    }
    res_VV.push_back(m_v);
  }
  return res_VV;
}

double DataProcess::getMatrixValue(const vector<vector<double>> &m_vv)
{
  if (m_vv.size() == 1)
  {
    return m_vv[0][0];
  }
  double ans = 0;
  int i, j, k;
  for (i = 0; i < m_vv.size(); i++)
  {
    vector<vector<double>> temp;

    for (j = 0; j < m_vv.size() - 1; j++)
    {
      vector<double> m_v;
      for (k = 0; k < m_vv.size() - 1; k++)
      {
        m_v.push_back(m_vv[j + 1][(k >= i) ? k + 1 : k]);
      }
      temp.push_back(m_v);
    }
    double t = this->getMatrixValue(temp);

    if (i % 2 == 0)
    {
      ans += m_vv[0][i] * t;
    }
    else
    {
      ans -= m_vv[0][i] * t;
    }
  }
  return ans;
}

// void getAStart(double arcs[N][N], int n, double ans[N][N])
// {
//   if (n == 1)
//   {
//     ans[0][0] = 1;
//     return;
//   }
//   int i, j, k, t;
//   double temp[N][N];
//   for (i = 0; i < n; i++)
//   {
//     for (j = 0; j < n; j++)
//     {
//       for (k = 0; k < n - 1; k++)
//       {
//         for (t = 0; t < n - 1; t++)
//         {
//           temp[k][t] = arcs[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
//         }
//       }

//       ans[j][i] = getA(temp, n - 1);
//       if ((i + j) % 2 == 1)
//       {
//         ans[j][i] = -ans[j][i];
//       }
//     }
//   }
// }

// bool GetMatrixInverse(double src[N][N], int n, double des[N][N])
// {
//   double flag = getA(src, n);
//   double t[N][N];
//   if (flag == 0)
//   {
//     return false;
//   }
//   else
//   {
//     getAStart(src, n, t);
//     for (int i = 0; i < n; i++)
//     {
//       for (int j = 0; j < n; j++)
//       {
//         des[i][j] = t[i][j] / flag;
//       }
//     }
//   }
//   return true;
// }

DataProcess::~DataProcess()
{
}
