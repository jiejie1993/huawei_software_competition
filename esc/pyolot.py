# coding=utf-8

import matplotlib.pyplot as plt
import sys
from pylab import * 
mpl.rcParams['font.sans-serif'] = ['SimHei']


filepath = sys.argv[1]
datas = []
with open(filepath, 'r') as f:
    lines = f.readlines()[0:]

for line in lines:
    datas.append(line.strip().split(','))
names, y1, y2, y3, y4 = [], [], [], [], []
for data in datas:
    names.append(int(data[0]))
    y1.append(float(data[1]))
    y2.append(float(data[2]))
#    y3.append(float(data[3]))
#    y4.append(float(data[4]))

x = range(len(names))
# plt.plot(x, y, 'ro-')
# plt.plot(x, y1, 'bo-')
# pl.xlim(-1, 11)  # 限定横轴的范围
# pl.ylim(-1, 110)  # 限定纵轴的范围
plt.plot(x, y1, marker='', mec='r', mfc='w')
plt.plot(x, y2, marker='', ms=10)
#plt.plot(x, y1, marker='', mec='r', mfc='w', label=u'y1')
#plt.plot(x, y2, marker='', ms=10, label=u'y2')
#plt.plot(x, y3, marker='.', ms=10, label=u'y3')
#plt.plot(x, y4, marker='.', ms=10, label=u'y4')
plt.legend()  # 让图例生效
plt.xticks(x, names, rotation=45)
plt.margins(0)
plt.subplots_adjust(bottom=0.15)
#plt.xlabel("x")  # X轴标签
#plt.ylabel("y")  # Y轴标签
plt.title("Result")  # 标题

plt.show()
