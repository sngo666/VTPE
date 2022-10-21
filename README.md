# 可验证阈值认证加密系统

使用任意生物特征识别模型生成一组长度不低于128位的特征向量模板，需要对特征向量进行z-score均值归一化以及矩阵单位化处理，并同时需要保证每个元素的范围[-1, 1]之间，利用PKI技术和数字水印以及质询机制防范随机攻击和重播攻击。
为了使用项目，你需要安装QT5版本(5.15.2+)以及armadillo数学库。
使用cmake进行build。
输入特征文件为n(n>128)行,每个元素占据独立一行，不包含任意特殊符号的txt文件（正负号除外）。
参考文献详见《Attacks and Countermeasures on Privacy-Preserving Biometric Authentication Schemes》
