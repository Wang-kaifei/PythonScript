'''
Descripttion: 
version: 
Author: Kaifei
Date: 2024-01-23 21:43:35
LastEditors: Kaifei
LastEditTime: 2024-01-25 00:14:41
'''
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# 假设的数据，每个块有5个数据点

# data = np.random.rand(4, 5) * 100  # 创建一个4x5的数组用于示例
engine = ['MaxQuant', 'Open-pFind', 'MS-GF+', 'MSFragger', 'Comet', 'pAnno2']
data = [[0, 0, 0, 0, 0, 0], [0, 0.993, 0, 0, 0, 0.998], [0, 0.988, 0, 0, 0, 0.998], [0, 0, 0, 0, 0, 0]] #precision
# data = [[0, 0, 0, 0, 0, 0], [0, 0.868, 0, 0, 0, 0.973], [0, 0.741, 0, 0, 0, 0.955], [0, 0, 0, 0, 0, 0]] #recall
data = np.array(data)
name = ['Dong-Ecoli', 'Xu-Yeast', 'Gygi-Human', 'Arabidopsis thaliana']
# Seaborn的色卡
palette = sns.color_palette("Set2", 6)

# 组的数量和组内柱子的数量
num_groups = len(data)
num_bars_in_group = len(data[0])

# 设置每个组的柱形图和组间的间距
bar_width = 0.15
group_width = (bar_width + 0.04) * num_bars_in_group + 0.6  # 每组的宽度包括所有柱子和一个更大间隔

# 设置初始位置
ind = np.arange(num_groups) * group_width

# 生成图表
fig, ax = plt.subplots()
fig.set_size_inches(15, 5)
# 绘制每个块中的柱形图
for i in range(num_bars_in_group):
    bars = ax.bar(ind + i * (bar_width + 0.04), data[:, i], width=bar_width, color=palette[i], label=engine[i])
    ax.bar_label(bars, labels=[f'{val*100:.1f}%' for val in data[:, i]], label_type='edge', fontsize=12)

# 添加一些文本用于标签、标题和自定义x轴刻度等
# ax.set_xlabel('Groups')
ax.set_ylabel('Precision', fontsize=16)
# ax.set_title('Multiple bars in each group')
ax.set_xticks(ind + group_width / 2 - 2*(bar_width + 0.04))
ax.set_xticklabels([name[i] for i in range(num_groups)], fontsize=16)

# 添加图例
ax.legend()
plt.savefig('pre_pep.png', format='png', dpi=300, bbox_inches='tight')  # 保存为高分辨率的图片
# 展示图表
plt.show()