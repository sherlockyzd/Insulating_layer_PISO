'''
将此文件夹下的所有matlab文件合并成一个文件，然后保存到当前目录下的whole.m文件中；
之后再将whole.m文件转换成wholenew.py文件！
'''
import os
import re

def get_file_list():
    file_list = []
    for root, dirs, files in os.walk('.'):
        for file in files:
            if file.endswith('.m'):
                file_list.append(os.path.join(root, file))
    return file_list

def get_whole_file(file_list):
    whole_file = ''
    for file in file_list:
        try:
            with open(file, 'r', encoding='utf-8') as f:
                whole_file += f.read()
        except (IOError, OSError) as e:
            print(f"Error reading file {file}: {e}")
    return whole_file

def save_whole_file(whole_file):
    with open('whole.m', 'w') as f:
        f.write(whole_file)

def m2py():
    with open('whole.m', 'r') as f:
        whole_file = f.read()
    whole_file = re.sub(r'%.*\n', '', whole_file)
    whole_file = re.sub(r'\.\.\.', '', whole_file)
    whole_file = re.sub(r'\n', '', whole_file)
    whole_file = re.sub(r'function', 'def', whole_file)
    whole_file = re.sub(r'end', '', whole_file)
    whole_file = re.sub(r'=', '#', whole_file)
    whole_file = re.sub(r';', '', whole_file)
    whole_file = re.sub(r'\[', '(', whole_file)
    whole_file = re.sub(r'\]', ')', whole_file)
    whole_file = re.sub(r'\(', '(', whole_file)
    whole_file = re.sub(r'\)', ')', whole_file)
    whole_file = re.sub(r'\{', '(', whole_file)
    whole_file = re.sub(r'\}', ')', whole_file)
    whole_file = re.sub(r'\.\*', '*', whole_file)
    whole_file = re.sub(r'\.\^', '**', whole_file)
    whole_file = re.sub(r'\.', '*', whole_file)
    whole_file = re.sub(r'/', '//', whole_file)
    whole_file = re.sub(r'//', '/', whole_file)
    whole_file = re.sub(r'for', 'for', whole_file)
    whole_file = re.sub(r'while', 'while', whole_file)
    whole_file = re.sub(r'if', 'if', whole_file)
    whole_file = re.sub(r'else', 'else', whole_file)
    whole_file = re.sub(r'elseif', 'elif', whole_file)
    whole_file = re.sub(r'break', 'break', whole_file)
    whole_file = re.sub(r'continue', 'continue', whole_file)
    whole_file = re.sub(r'return', 'return', whole_file)
    whole_file = re.sub(r'zeros', 'np.zeros', whole_file)
    whole_file = re.sub(r'ones', 'np.ones', whole_file)
    whole_file = re.sub(r'eye', 'np.eye', whole_file)
    whole_file = re.sub(r'rand', 'np.random.rand', whole_file)
    whole_file = re.sub(r'randn', 'np.random.randn', whole_file)
    whole_file = re.sub(r'length', 'len', whole_file)
    whole_file = re.sub(r'plot', 'plt.plot', whole_file)
    whole_file = re.sub(r'hold', 'plt.hold', whole_file)
    whole_file = re.sub(r'figure', 'plt.figure', whole_file)
    whole_file = re.sub(r'xlabel', 'plt.xlabel', whole_file)
    whole_file = re.sub(r'ylabel', 'plt.ylabel', whole_file)
    whole_file = re.sub(r'title', 'plt.title', whole_file)
    whole_file = re.sub(r'legend', 'plt.legend', whole_file)
    whole_file = re.sub(r'subplot', 'plt.subplot', whole_file)
    whole_file = re.sub(r'semilogy', 'plt.semilogy', whole_file)
    whole_file = re.sub(r'semilogx', 'plt.semilogx', whole_file)
    whole_file = re.sub(r'loglog', 'plt.loglog', whole_file)
    whole_file = re.sub(r'grid', 'plt.grid', whole_file)
    whole_file = re.sub(r'axis', 'plt.axis', whole_file)
    whole_file = re.sub(r'xlim', 'plt.xlim', whole_file)
    whole_file = re.sub(r'ylim', 'plt.ylim', whole_file)
    whole_file = re.sub(r'close', 'plt.close', whole_file)
    whole_file = re.sub(r'subplot', 'plt.subplot', whole_file)
    whole_file = re.sub(r'pause', 'plt.pause', whole_file)
    whole_file = re.sub(r'print', 'print', whole_file)
    whole_file = re.sub(r'input', 'input', whole_file)
    whole_file = re.sub(r'fprintf', 'print', whole_file)
    whole_file = re.sub(r'fopen', 'open', whole_file)
    whole_file = re.sub(r'fclose', 'close', whole_file)

    with open('wholenew.py', 'w') as f:
        f.write(whole_file)




whole_file_list=get_file_list()
whole_file=get_whole_file(whole_file_list)
save_whole_file(whole_file)
