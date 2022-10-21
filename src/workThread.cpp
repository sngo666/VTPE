#include "dataProcess.h"
#include "workThread.h"

workThread::workThread(QObject *parent) : QObject(parent)
{
  isWorking = false;
  isSent = false;
}
void workThread::startWork()
{
  emit workStart();

  doWork();
}

void workThread::doWork()
{
  this->isWorking = true;
  emit sendMessage("workMode:" + QString::number(this->workMode) + " Start!");
  emit sendMessage("Xita: " + QString::number(this->systemParam.xita));
  emit sendMessage("Victor Dimension: " + QString::number(this->systemParam.dimension));
  emit sendMessage("userID: " + this->userID);

  switch (this->workMode)
  {
  // register stage
  case (1):
  {

    if (!this->registerStage(registerVector))
      emit sendMessage("register err!");
    else
      emit sendMessage("register completed.");

    break;
  }
  // challenge stage
  case (2):
  {
    int flag = 0;
    // for (int i = 0; i < 10; i++)

    // {
    if (!this->loginStage(registerVector))
    {
      emit sendMessage("login err!");
      flag++;
    }
    else
      emit sendMessage("successful authentication !!!");
    // }
    // emit sendMessage("[log] 10 times try, fail " + QString::number(flag) + " times");

    break;
  }
  }
  this->isSent = false;
  emit workFinished();
}

bool workThread::registerStage(const vector<double> &registerV)
{
  int num = registerV.size(), b;
  if (num > 500)
    b = QRandomGenerator::global()->bounded(128, num);
  else
    b = QRandomGenerator::global()->bounded(20, num); // for test

  vector<vector<double>> sepFeatureVV;
  tool->randomSeperate(registerV, sepFeatureVV, b);
  vector<int> niV = tool->getNiVector(sepFeatureVV);

  vector<privateKey> Sk = this->keyGenerationMode(niV);
  // this->PrintMatrix(Sk[0].M1, "SK_M1");
  // this->PrintMatrix(Sk[0].M2, "SK_M2");
  // this->PrintVector(Sk[0].randomEx, "Sk_randomorder");

  this->featureMatrix_Cx = this->encodingMode(Sk, sepFeatureVV);
  // this->PrintMatrix(featureMatrix_Cx[0].Cxi, "Cxi");

  for (int i = 0; i < b; i++)
  {
    vector<vector<double>> m_pi;
    tool->getInvertibleMatrix(m_pi, niV[i] + 3);
    non_sigularMatrix temp;
    temp.matri = m_pi;

    this->randomMatrix_P.push_back(temp);
  }

  for (int i = 0; i < b; i++)
  {
    this->hash_P.push_back(tool->getMatrixHash(this->randomMatrix_P[i].matri, PRIMENUMBER_P, PRIMENUMBER_MOD));
  }
  if (!registerDataSaving(this->userID, Sk, this->randomMatrix_P, this->hash_P, this->featureMatrix_Cx))
    return false;

  bool extractFlag = 1;
  vector<non_sigularMatrix> test_ns = this->getRegisterRandomMatrixP(this->userID, niV);
  for (int i = 0; i < test_ns.size(); i++)
  {
    for (int j = 0; j < test_ns[i].matri.size(); j++)
    {
      for (int k = 0; k < test_ns[i].matri.size(); k++)
      {
        if (test_ns[i].matri[j][k] != this->randomMatrix_P[i].matri[j][k])
          extractFlag = 0;
      }
    }
  }
  if (extractFlag && test_ns.size() != 0)
    emit sendMessage("[err] Matrix P extraction success!");
  else
  {
    emit sendMessage("[err] Matrix P extraction fail!");
    return false;
  }

  extractFlag = 1;
  vector<sub_Template> test_st = this->getFeatureTemplateCx(this->userID, niV);
  for (int i = 0; i < test_st.size(); i++)
  {
    for (int j = 0; j < test_st[i].Cxi.size(); j++)
    {
      for (int k = 0; k < test_st[i].Cxi.size(); k++)
      {
        if (test_st[i].Cxi[j][k] != this->featureMatrix_Cx[i].Cxi[j][k])
          extractFlag = 0;
      }
    }
  }
  if (extractFlag && test_st.size() != 0)
    emit sendMessage("[err] Template Cx extraction success!");
  else
  {
    emit sendMessage("[err] Template Cx extraction fail!");
    return false;
  }

  extractFlag = 1;
  vector<privateKey> m_Sk = this->getRegisterPrivateKey(this->userID);
  emit sendMessage("[log] " + QString::number(m_Sk.size()) + " privateKey eles has been extracted!");

  for (int i = 0; i < m_Sk.size(); i++)
  {
    if (m_Sk[i].scale != Sk[i].scale)
    {
      emit sendMessage("[err] Ni extracting error!");
      extractFlag = 0;
    }
    for (int j = 0; j < Sk[i].scale + 3; j++)
    {
      if (m_Sk[i].randomEx[j] != Sk[i].randomEx[j])
      {
        emit sendMessage("[err] randomOrder extracting error!");
        extractFlag = 0;
      }
      for (int k = 0; k < Sk[i].scale + 3; k++)
      {
        if (m_Sk[i].M1[j][k] != Sk[i].M1[j][k] || m_Sk[i].M2[j][k] != Sk[i].M2[j][k])
        {
          emit sendMessage("[err] M1/M2 extracting error!");
          extractFlag = 0;
        }
      }
    }
  }
  if (extractFlag == 0)
    return false;
  emit sendMessage("[log] Sk extracting successfully!");

  // this->PrintMatrix(Sk[0].M1);
  // this->PrintMatrix(Sk[0].M2);
  // this->PrintVector(Sk[0].randomEx);

  return true;
}

bool workThread::loginStage(const vector<double> &loginV)
{
  QString m_id = this->userID;
  vector<privateKey> m_sk = this->getRegisterPrivateKey(m_id);
  if (m_sk.size() == 0)
  {
    emit sendMessage("[err] user info query wrong!");
    return false;
  }
  int b = m_sk.size();

  vector<int> m_ni;
  for (int i = 0; i < b; i++)
  {
    m_ni.push_back(m_sk[i].scale);
  }

  this->featureMatrix_Cx = this->getFeatureTemplateCx(m_id, m_ni);
  this->randomMatrix_P = this->getRegisterRandomMatrixP(m_id, m_ni);

  for (int i = 0; i < b; i++)
  {
    vector<vector<double>> m_qi;
    tool->getInvertibleMatrix(m_qi, m_sk[i].scale + 3);
    non_sigularMatrix temp;
    temp.matri = m_qi;
    this->randomMatrix_Q.push_back(temp);
  }

  for (int i = 0; i < b; i++)
  {
    int ni = m_sk[i].scale + 3;

    vector<vector<double>> m_ui;
    mat m_qi(ni, ni);
    m_qi = tool->vvToMatrixXd(this->randomMatrix_Q[i].matri);
    mat m_cxi(ni, ni);
    m_cxi = tool->vvToMatrixXd(this->featureMatrix_Cx[i].Cxi);

    m_ui = tool->MatrixXdToVV(m_qi * m_cxi);
    non_sigularMatrix temp;
    temp.matri = m_ui;
    this->challengeMatrix_U.push_back(temp);
  }

  embed_res loginEmbed = this->embeddingMode(m_sk, loginV, m_ni);
  this->hash_W = loginEmbed.hash_w;
  for (int i = 0; i < b; i++)
  {
    if (loginEmbed.alter_y[i].size() != m_ni[i] + 3)
      emit sendMessage("[err] y' calculate error!");
  }

  vector<query_Template> Ty = this->tokenGenerationMode(m_sk, loginEmbed.alter_y);

  for (int i = 0; i < m_sk.size(); i++)
  {
    int ni = m_sk[i].scale + 3;
    vector<vector<double>> m_vi;

    mat m_pi(ni, ni);
    m_pi = tool->vvToMatrixXd(this->randomMatrix_P[i].matri);
    mat m_ui(ni, ni);
    m_ui = tool->vvToMatrixXd(this->challengeMatrix_U[i].matri);
    mat m_tyi(ni, ni);
    m_tyi = tool->vvToMatrixXd(Ty[i].Tyi);

    m_vi = tool->MatrixXdToVV(m_pi * m_ui * m_tyi);
    non_sigularMatrix temp;
    temp.matri = m_vi;
    this->queryMatrix_V.push_back(temp);
  }

  // vector<double> resVector_v;
  vector<double> resVector_v = this->decodingMode(this->featureMatrix_Cx, Ty);

  vector<int> signVector_w = this->extractionMod(resVector_v, this->hash_W);
  double sum_v = 0.0;
  for (int i = 0; i < resVector_v.size(); i++)
  {
    sum_v += resVector_v[i];
  }
  unsigned long long m_hashw = tool->getVectorHash(signVector_w, PRIMENUMBER_P, PRIMENUMBER_MOD);
  emit sendMessage("[log] sum_v = " + QString::number(sum_v));
  bool flag = true;
  if (sum_v < 0)
  {
    flag = false;
  }
  if (signVector_w.size() == 0)
  {
    emit sendMessage("[err] signV error ");
    flag = false;
  }

  if (m_hashw != this->hash_W)
  {
    emit sendMessage("[err] hash error");
    flag = false;
  }

  return flag;
}

vector<privateKey> workThread::keyGenerationMode(vector<int> &subVectorSize)
{
  vector<privateKey> m_sk;
  for (int i = 0; i < subVectorSize.size(); i++)
  {
    privateKey m_p;
    m_p.M1 = tool->getInvertibleMatrix(subVectorSize[i] + 3);
    m_p.M2 = tool->getInvertibleMatrix(subVectorSize[i] + 3);
    m_p.randomEx = tool->getRandomOrder(subVectorSize[i] + 3);
    m_p.scale = subVectorSize[i];
    m_sk.push_back(m_p);
  }

  // just for test
  int flag = 0, sum = 0;
  for (int i = 0; i < subVectorSize.size(); i++)
  {
    // emit sendMessage("scale: " + QString::number(m_sk[i].scale));

    flag++;
    sum += m_sk[i].scale;
  }

  emit sendMessage("[log] scale sum: " + QString::number(sum));
  emit sendMessage("[log] " + QString::number(flag) + " groups in sk.");

  return m_sk;
}

vector<sub_Template> workThread::encodingMode(const vector<privateKey> &m_sk, const vector<vector<double>> &seperatedVV)
{
  vector<sub_Template> featureTemp;
  double beta = ((double)QRandomGenerator::global()->bounded(1, 10));
  for (int i = 0; i < m_sk.size(); i++)
  {
    sub_Template m_sT;

    vector<double> operateV;
    for (int j = 0; j < seperatedVV[i].size(); j++)
    {
      operateV.push_back(beta * seperatedVV[i][j]);
    }
    operateV.push_back(beta * (-1.0));
    operateV.push_back(beta * (double)QRandomGenerator::global()->bounded(-10, 10));
    operateV.push_back(0);

    if (operateV.size() != (m_sk[i].scale + 3))
    {
      emit sendMessage("[err] scale wrong!");
      return featureTemp;
    }
    vector<vector<double>> Xi = tool->getDiagonalMatrix(operateV, m_sk[i].randomEx);
    vector<vector<double>> Sxi = tool->getRandomLowerTriangleMatrix(operateV.size()); // size: n+3

    // if (i == 0)
    // {
    //   this->PrintMatrix(Xi, "Xi");
    //   this->PrintMatrix(Sxi, "Sxi");
    // }

    mat m_M1(operateV.size(), operateV.size());
    m_M1 = tool->vvToMatrixXd(m_sk[i].M1);
    mat m_M2(operateV.size(), operateV.size());
    m_M2 = tool->vvToMatrixXd(m_sk[i].M2);
    mat m_Xi(operateV.size(), operateV.size());
    m_Xi = tool->vvToMatrixXd(Xi);
    mat m_Sxi(operateV.size(), operateV.size());
    m_Sxi = tool->vvToMatrixXd(Sxi);

    vector<vector<double>> m_Cxi = tool->MatrixXdToVV(m_M1 * m_Sxi * m_Xi * m_M2);

    // emit sendMessage(QString::number(m_M1(0, 0)));

    // this->PrintMatrix(Xi);
    // this->PrintMatrix(Sxi);
    // this->PrintMatrix(m_Cxi);
    m_sT.Cxi = m_Cxi;
    featureTemp.push_back(m_sT);
  }
  return featureTemp;
}

embed_res workThread::embeddingMode(const vector<privateKey> &m_sk, const vector<double> &y, vector<int> Ni)
{
  embed_res res;
  int num = this->systemParam.dimension;
  double xita = this->systemParam.xita;

  vector<vector<double>> ranMXita = tool->getEmbedRandomVV(Ni, xita);
  vector<int> sign_W = tool->getSign(ranMXita[0], ranMXita[1]);

  unsigned long long h = tool->getVectorHash(sign_W, PRIMENUMBER_P, PRIMENUMBER_MOD);
  res.hash_w = h;
  int count = 0;
  double aerfa = ((double)QRandomGenerator::global()->bounded(1, 10));
  vector<vector<double>> alter_yVV;

  for (int i = 0; i < m_sk.size(); i++)
  {
    // emit sendMessage(QString::number(i) + "th group");

    vector<double> yi;
    for (int j = 0; j < Ni[i]; j++)
    {
      // emit sendMessage(QString::number(y[count]));

      yi.push_back(y[count] * aerfa);

      count++;
    }
    yi.push_back((ranMXita[0][i] + ranMXita[1][i]) * aerfa);
    yi.push_back(0);
    yi.push_back(((double)QRandomGenerator::global()->bounded(-10, 10)) * aerfa);

    alter_yVV.push_back(yi);
  }

  res.alter_y = alter_yVV;
  return res;
}

vector<query_Template> workThread::tokenGenerationMode(const vector<privateKey> &m_sk, const vector<vector<double>> &m_alter_y)
{
  vector<query_Template> res_V;
  for (int i = 0; i < m_sk.size(); i++)
  {
    query_Template q_temp;
    int ni = m_sk[i].scale + 3;
    vector<vector<double>> Yi = tool->getDiagonalMatrix(m_alter_y[i], m_sk[i].randomEx);
    vector<vector<double>> Syi = tool->getRandomLowerTriangleMatrix(ni); // size: ni+3
    // if (i == 0)
    // {
    //   this->PrintVector(m_sk[i].randomEx);
    //   this->PrintVector(m_alter_y[i]);
    //   this->PrintMatrix(Yi);
    //   this->PrintMatrix(Syi);
    //   this->PrintMatrix(m_sk[i].M1);
    //   this->PrintMatrix(m_sk[i].M2);
    // }

    mat m_M1(ni, ni);
    m_M1 = tool->vvToMatrixXd(m_sk[i].M1);
    mat m_M2(ni, ni);
    m_M2 = tool->vvToMatrixXd(m_sk[i].M2);
    mat m_Yi(ni, ni);
    m_Yi = tool->vvToMatrixXd(Yi);
    mat m_Syi(ni, ni);
    m_Syi = tool->vvToMatrixXd(Syi);

    q_temp.Tyi = tool->MatrixXdToVV(inv(m_M2) * m_Yi * m_Syi * inv(m_M1));
    res_V.push_back(q_temp);
  }
  return res_V;
}

vector<double> workThread::decodingMode(const vector<sub_Template> &m_Cx, const vector<query_Template> &m_Ty)
{
  vector<double> res_V;
  for (int i = 0; i < m_Cx.size(); i++)
  {
    int ni = m_Cx[i].Cxi.size();
    if (m_Cx[i].Cxi.size() != m_Ty[i].Tyi.size())
      emit sendMessage("[err] mode:decoding scale wrong!");

    mat vi_matrix;
    mat m_cxi(ni, ni);
    m_cxi = tool->vvToMatrixXd(m_Cx[i].Cxi);
    mat m_tyi(ni, ni);
    m_tyi = tool->vvToMatrixXd(m_Ty[i].Tyi);

    res_V.push_back(trace(m_cxi * m_tyi));
    // emit sendMessage("[log] trace: " + QString::number(trace(vi_matrix)));
  }
  return res_V;
}

vector<int> workThread::extractionMod(const vector<double> &m_v, unsigned long long m_h)
{
  vector<int> res_V;
  vector<int> sign_V = tool->getSign(m_v, -1);
  unsigned long long hash_w = tool->getVectorHash(sign_V, PRIMENUMBER_P, PRIMENUMBER_MOD);
  if (hash_w != m_h)
    return res_V;
  else
    return sign_V;
}

bool workThread::registerDataSaving(QString id, const vector<privateKey> &m_sk, const vector<non_sigularMatrix> &m_P, const vector<unsigned long long> &m_hash_p, const vector<sub_Template> &m_Cx)
{
  QString filePath = QCoreApplication::applicationDirPath() + "/usr/." + id;
  QFileInfo userInfo(filePath);
  if (userInfo.isFile())
  {
    QFile::remove(filePath);
  }
  QFile userfile(filePath);
  bool flag = userfile.open(QIODevice::WriteOnly | QIODevice::Truncate);
  if (!flag)
  {
    emit sendMessage(filePath);

    emit sendMessage("[err] data establish wrong!");
    return flag;
  }
  if (userfile.isOpen())
  {
    QTextStream in(&userfile);
    in << QString("userID:") << Qt::endl;
    in << id << Qt::endl;

    in << QString("Ni:") << Qt::endl;
    for (int i = 0; i < m_sk.size(); i++)
    {
      in << QString::number(m_sk[i].scale) << Qt::endl;
    }

    in << QString("M1:") << Qt::endl;
    for (int i = 0; i < m_sk.size(); i++)
    {
      for (int j = 0; j < m_sk[i].scale + 3; j++)
      {
        for (int k = 0; k < m_sk[i].scale + 3; k++)
          in << QString::number(m_sk[i].M1[j][k]) << Qt::endl;
      }
    }

    in << QString("M2:") << Qt::endl;
    for (int i = 0; i < m_sk.size(); i++)
    {
      for (int j = 0; j < m_sk[i].scale + 3; j++)
      {
        for (int k = 0; k < m_sk[i].scale + 3; k++)
          in << QString::number(m_sk[i].M2[j][k]) << Qt::endl;
      }
    }

    in << QString("RandomOrder:") << Qt::endl;
    for (int i = 0; i < m_sk.size(); i++)
    {
      for (int j = 0; j < m_sk[i].scale + 3; j++)
      {
        in << QString::number(m_sk[i].randomEx[j]) << Qt::endl;
      }
    }

    in << QString("RandomMatrixVectorP:") << Qt::endl;
    for (int i = 0; i < m_sk.size(); i++)
    {
      for (int j = 0; j < m_sk[i].scale + 3; j++)
      {
        for (int k = 0; k < m_sk[i].scale + 3; k++)
          in << QString::number(m_P[i].matri[j][k]) << Qt::endl;
      }
    }

    in << QString("HashP:") << Qt::endl;
    for (int i = 0; i < m_sk.size(); i++)
    {
      in << QString::number(m_hash_p[i]) << Qt::endl;
    }

    in << QString("fearureTemplatCx:") << Qt::endl;
    for (int i = 0; i < m_sk.size(); i++)
    {
      for (int j = 0; j < m_sk[i].scale + 3; j++)
      {
        for (int k = 0; k < m_sk[i].scale + 3; k++)
          in << QString::number(m_Cx[i].Cxi[j][k], 'g', 20) << Qt::endl;
      }
    }

    userfile.close();
    return true;
  }

  return false;
}

vector<privateKey> workThread::getRegisterPrivateKey(QString userID)
{
  vector<privateKey> m_Sk;
  vector<int> m_Ni;
  QString qstr = "Ni:";

  QString filePath = QCoreApplication::applicationDirPath() + "/usr/." + userID;
  QFileInfo fileinfo(filePath);
  if (!fileinfo.isFile())
    return m_Sk;

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
    if (str.size() == qstr.size())
    {
      bool flag = 1;
      for (int j = 0; j < qstr.size(); j++)
      {
        if (str[j] != qstr[j])
          flag = 0;
      }
      if (flag)
        count = i + 1;
    }
  }

  if (count == -1)
    return m_Sk;

  emit sendMessage("[log] Ni rows in database: " + QString::number(count));

  while (strList[count] != "M1:")
  {
    m_Ni.push_back(strList[count].toInt());
    count++;
  }

  emit sendMessage("[log] M1 rows in database: " + QString::number(count));

  count++;
  vector<non_sigularMatrix> m_M1;
  for (int i = 0; i < m_Ni.size(); i++)
  {
    non_sigularMatrix m_vv;
    for (int j = 0; j < m_Ni[i] + 3; j++)
    {
      vector<double> m_v;
      for (int k = 0; k < m_Ni[i] + 3; k++)
      {
        m_v.push_back(strList[count].toDouble());
        count++;
      }
      m_vv.matri.push_back(m_v);
    }
    m_M1.push_back(m_vv);
  }

  emit sendMessage("[log] M2 rows in database: " + QString::number(count));

  count++;

  vector<non_sigularMatrix> m_M2;
  for (int i = 0; i < m_Ni.size(); i++)
  {
    non_sigularMatrix m_vv;
    for (int j = 0; j < m_Ni[i] + 3; j++)
    {
      vector<double> m_v;
      for (int k = 0; k < m_Ni[i] + 3; k++)
      {
        m_v.push_back(strList[count].toDouble());
        count++;
      }
      m_vv.matri.push_back(m_v);
    }
    m_M2.push_back(m_vv);
  }

  emit sendMessage("[log] randomOrder rows in database: " + QString::number(count));

  count++;

  vector<vector<int>> m_RO;
  for (int i = 0; i < m_Ni.size(); i++)
  {
    vector<int> m_v;
    for (int j = 0; j < m_Ni[i] + 3; j++)
    {
      m_v.push_back(strList[count].toInt());
      count++;
    }
    m_RO.push_back(m_v);
  }

  for (int i = 0; i < m_Ni.size(); i++)
  {
    privateKey m_key;
    m_key.M1 = m_M1[i].matri;
    m_key.M2 = m_M2[i].matri;
    m_key.randomEx = m_RO[i];
    m_key.scale = m_Ni[i];
    m_Sk.push_back(m_key);
  }

  return m_Sk;
}

vector<sub_Template> workThread::getFeatureTemplateCx(QString userID, const vector<int> &Ni)
{

  vector<sub_Template> m_Cx;
  QString qstr = "fearureTemplatCx:";

  QString filePath = QCoreApplication::applicationDirPath() + "/usr/." + userID;
  QFileInfo fileinfo(filePath);
  if (!fileinfo.isFile())
    return m_Cx;

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
    if (str.size() == qstr.size())
    {
      bool flag = 1;
      for (int j = 0; j < qstr.size(); j++)
      {
        if (str[j] != qstr[j])
          flag = 0;
      }
      if (flag)
        count = i + 1;
    }
  }

  if (count == -1)
    return m_Cx;
  emit sendMessage("[log] Cx rows in database: " + QString::number(count));

  for (int i = 0; i < Ni.size(); i++)
  {
    sub_Template m_st;
    for (int j = 0; j < Ni[i] + 3; j++)
    {
      vector<double> m_v;
      for (int k = 0; k < Ni[i] + 3; k++)
      {
        m_v.push_back(strList[count].toDouble());
        count++;
      }
      m_st.Cxi.push_back(m_v);
    }
    m_Cx.push_back(m_st);
  }

  file.close();
  return m_Cx;
}

vector<non_sigularMatrix> workThread::getRegisterRandomMatrixP(QString userID, const vector<int> &Ni)
{

  vector<non_sigularMatrix> m_P;
  QString qstr = "RandomMatrixVectorP:";

  QString filePath = QCoreApplication::applicationDirPath() + "/usr/." + userID;
  QFileInfo fileinfo(filePath);
  if (!fileinfo.isFile())
    return m_P;

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
    if (str.size() == qstr.size())
    {
      bool flag = 1;
      for (int j = 0; j < qstr.size(); j++)
      {
        if (str[j] != qstr[j])
          flag = 0;
      }
      if (flag)
        count = i + 1;
    }
  }

  if (count == -1)
    return m_P;
  emit sendMessage("[log] P rows in database: " + QString::number(count));

  for (int i = 0; i < Ni.size(); i++)
  {
    non_sigularMatrix m_ns;
    for (int j = 0; j < Ni[i] + 3; j++)
    {
      vector<double> m_v;
      for (int k = 0; k < Ni[i] + 3; k++)
      {
        m_v.push_back(strList[count].toDouble());
        count++;
      }
      m_ns.matri.push_back(m_v);
    }
    m_P.push_back(m_ns);
  }

  file.close();
  return m_P;
}

void workThread::PrintMatrix(const vector<vector<double>> &m_v)
{
  emit sendMessage("[log] print VV--------------------");

  for (int i = 0; i < m_v.size(); i++)
  {
    QString str;
    for (int j = 0; j < m_v.size(); j++)
    {
      str += QString::number(m_v[i][j]) + " ";
    }
    emit sendMessage(str);
  }
  emit sendMessage("[log] " + QString::number(m_v.size()) + " * " + QString::number(m_v.size()) + "  elements in total");
}

void workThread::PrintMatrix(const vector<vector<double>> &m_v, QString descri)
{
  emit sendMessage("[log] ----print " + descri + "---");

  for (int i = 0; i < m_v.size(); i++)
  {
    QString str;
    for (int j = 0; j < m_v.size(); j++)
    {
      str += QString::number(m_v[i][j]) + " ";
    }
    emit sendMessage(str);
  }
  emit sendMessage("[log] " + QString::number(m_v.size()) + " * " + QString::number(m_v.size()) + "  elements in " + descri);
}

void workThread::PrintVector(const vector<double> &m_v)
{
  emit sendMessage("[log] print Vector-----");

  QString str;
  for (int j = 0; j < m_v.size(); j++)
  {
    str += QString::number(m_v[j]) + " ";
  }
  emit sendMessage(str);
  emit sendMessage("[log] " + QString::number(m_v.size()) + " * " + QString::number(m_v.size()) + "  elements in total");
}

void workThread::PrintVector(const vector<double> &m_v, QString descri)
{
  emit sendMessage("[log] ----print " + descri + "---");

  QString str;
  for (int j = 0; j < m_v.size(); j++)
  {
    str += QString::number(m_v[j]) + " ";
  }
  emit sendMessage(str);
  emit sendMessage("[log] " + QString::number(m_v.size()) + " * " + QString::number(m_v.size()) + "  elements in " + descri);
}

void workThread::PrintVector(const vector<int> &m_v)
{
  emit sendMessage("[log] print Vector-----");

  QString str;
  for (int j = 0; j < m_v.size(); j++)
  {
    str += QString::number(m_v[j]) + " ";
  }
  emit sendMessage(str);
  emit sendMessage("[log] " + QString::number(m_v.size()) + " * " + QString::number(m_v.size()) + "  elements in total");
}

void workThread::PrintVector(const vector<int> &m_v, QString descri)
{
  emit sendMessage("[log] ----print " + descri + "---");

  QString str;
  for (int j = 0; j < m_v.size(); j++)
  {
    str += QString::number(m_v[j]) + " ";
  }
  emit sendMessage(str);
  emit sendMessage("[log] " + QString::number(m_v.size()) + " * " + QString::number(m_v.size()) + "  elements in " + descri);
}

void workThread::getXita(double m_xita)
{
  this->systemParam.xita = m_xita;
}

void workThread::getDimension(int m_dimension)
{
  this->systemParam.dimension = m_dimension;
}
void workThread::getRegisterVector(vector<double> m_v)
{
  this->registerVector = m_v;
  tool = new DataProcess(m_v.size());
}

void workThread::getChallengeVector(vector<double> m_v)
{
  this->challengeVector = m_v;
  tool = new DataProcess(m_v.size());
}

void workThread::getWorkMode(int m_mode)
{
  this->workMode = m_mode;
}

void workThread::getUserID(QString m_str)
{
  emit sendMessage("[log] get userID!");

  this->userID = m_str;
}
workThread::~workThread()
{
}