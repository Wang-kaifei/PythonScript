# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""脚本功能：查看decoy鉴定数量随鉴定总数的变化情况"""
def get_node(filename):
    """函数功能：得到要画的点，step为步长，每一步画一个点"""
    with open(filename) as f:
        f = f.readlines()
    step = len(f) / 30
    x_node = [0]  #横坐标点，鉴定到的PSM总数
    y_node = [0]  #纵坐标点，鉴定到的decoy PSM总数
    decoy = 0
    target = 0
    for i in range(1, len(f)):
       # print type(f[i].strip().rsplit('\t', 1)[1])
        if f[i].strip().rsplit('\t', 1)[1] == "True" or f[i].strip().rsplit('\t', 1)[1] == "TRUE":
            target += 1
        else:
            decoy += 1
        if (target + decoy) % step == 0 or i == (len(f) - 1):
            x_node.append(target + decoy)
            y_node.append(decoy)
    print x_node
    print y_node
    return x_node, y_node

def decoy(filename):
    """函数功能：查看文件中有多少decoy符合FDR限制"""
    with open(filename) as f:
        f = f.readlines()
    res = 0
    t = 0
    for i in range(1, len(f)):
        segs = f[i].strip().split('\t')
        t += 1
        if float(segs[4]) > 0.01:
            break
       # print segs[-1]
        if segs[-1] != "True":
            res += 1
            
    print res, t, (res / float(t))
    
'''if __name__ == "__main__":
    filename = "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\complete1\\groupt3.spectra"
    output = "C:\\Users\\kfwang\\Desktop\\plot.txt"
    x, y = get_node(filename)
    res = ""
    for node in x:
        res += str(node) + ','
    res += '\n'
    for node in y:
        res += str(node) + ','
    res += '\n'
    with open(output, 'w') as f:
        f.write(res)
        f.close()'''
if __name__ == "__main__":
    filename = ["C:\\test_wkf\\two_step_test\\Ecoli\\database1\\complete1\\groupt0.spectra", "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\complete1\\groupt1.spectra", "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\complete1\\groupt2.spectra", "C:\\test_wkf\\two_step_test\\Ecoli\\database1\\complete1\\groupt3.spectra"]
    for f in filename:
        decoy(f)
