import matplotlib.pyplot as plt
import seaborn as sns

def set_matplotlib_rcparams():
    # 设置全局字体为Arial
    plt.rcParams['font.family'] = 'Arial'
    # 设置全局字体为加粗
    plt.rcParams['font.weight'] = 'bold'
    # 设置标题字号
    plt.rcParams['axes.titlesize'] = 16

    # 设置轴标签字号
    plt.rcParams['axes.labelsize'] = 18

    # 设置图例字号
    plt.rcParams['legend.fontsize'] = 12

    # 设置坐标轴刻度标签字号
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 15
    
def plot_stacked_bar_with_labels():
    # 示例数据
    categories = ['Part 1', 'Part 2', 'Part 3', 'Part 4', 'Part 5']
    series_labels = ['Series A', 'Series B']

    # data = {
    #     'Series A': [4756, 4662, 0, 10890],
    #     'Series B': [0, 0, 0, 0]
    # }

    data = {
        'Series A': [7220, 5964, 6434, 0, 10209],
        'Series B': [0, 0, 0, 0, 2383]
    }
    
    # 使用seaborn的配色样式
    sns.set_palette("deep")
    total_bars = len(categories)
    bar_width = 0.23  # 设置每个柱子的宽度为0.25
    bar_padding = 0.08  # 设置柱子之间的间距为0.1

    # 计算柱子的位置
    positions = []
    for i in range(total_bars):
        position = i * (bar_width + bar_padding)
        positions.append(position)

    fig, ax = plt.subplots()

    for i, series_label in enumerate(series_labels):
        if i == 0:
            ax.bar(positions, data[series_label], label=series_label, width=bar_width)
        else:
            bottom = [sum(x) for x in zip(*[data[sl] for sl in series_labels[:i]])]
            ax.bar(positions, data[series_label], bottom=bottom, label=series_label, width=bar_width)

        # 在每个柱子上添加数字标签
        for j, value in enumerate(data[series_label]):
            height = bottom[j] + value if i > 0 else value
            if value != 0 or i == 0:
                ax.text(positions[j], height, str(value), ha='center', va='bottom', fontsize=18)

    plt.subplots_adjust(left=0.1, right=0.5)  # 可根据需要调整left和right参数，增加或减少间距
    ax.set_ylabel('# PSM', fontweight='bold')
    plt.ylim(0, 14500)  # 范围根据您的数据进行调整
    # 设置x轴标签和刻度
    ax.set_xticks([pos for pos in positions])
    ax.set_xticklabels(categories, rotation=0)  # 将x轴刻度设置为categories列表中的字符串
    # 展示图表
    plt.tight_layout()
    plt.savefig(r"D:\OneDrive - mails.ucas.ac.cn\年中技术报告\1.jpg", dpi=300)
    plt.show()
    
if __name__ == "__main__":
    set_matplotlib_rcparams()  # 设置matplotlib的全局字体样式
    plot_stacked_bar_with_labels()

