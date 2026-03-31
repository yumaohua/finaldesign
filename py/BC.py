import pandas as pd
import matplotlib.pyplot as plt

# 读取数据（没有表头，因此指定列名）
df_B = pd.read_csv('./data/B_row0.csv', header=None, names=['dist', 'Bmag'])
df_C = pd.read_csv('./data/C_row0.csv', header=None, names=['dist', 'Cmag'])

# 绘制 B_row0
plt.figure()
plt.plot(df_B['dist'], df_B['Bmag'], marker='.', linestyle='None')
plt.xlabel('Distance from segment 0 midpoint')
plt.ylabel('|Re(B₀ₗ)| + |Im(B₀ₗ)|')
plt.title('B matrix row 0 vs distance')
plt.savefig('./pic/B_row0_plot.svg')
plt.show()

# 绘制 C_row0
plt.figure()
plt.plot(df_C['dist'], df_C['Cmag'], marker='.', linestyle='None')
plt.xlabel('Distance from segment 0 midpoint')
plt.ylabel('|Re(C₀ₗ)| + |Im(C₀ₗ)|')
plt.title('C matrix row 0 vs distance')
plt.savefig('./pic/C_row0_plot.svg')
plt.show()
