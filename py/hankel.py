import pandas as pd
import matplotlib.pyplot as plt

# 读取数据
data = pd.read_csv('./dataf/error_data.txt')

# 提取数据
z = data['z']
error_nu_0 = data['error_nu_0']
error_nu_1 = data['error_nu_1']

# 绘制误差曲线
plt.figure(figsize=(10, 6))
plt.plot(z, error_nu_0, label='Error for nu = 0', linewidth=1.5)
plt.plot(z, error_nu_1, label='Error for nu = 1', linewidth=1.5, linestyle='--')

# 配置坐标轴和标题
plt.xlabel('z', fontsize=12)
plt.ylabel('Error (Absolute Difference)', fontsize=12)
plt.title('Error Between Custom Hankel Function and Boost Implementation', fontsize=14, pad=20)
plt.legend(fontsize=11, loc='upper right')
plt.grid(True, linestyle='--', alpha=0.7)

# 保存为 SVG 矢量图
save_path = './picf/hankel_error_comparison.png'
plt.savefig(
    save_path,
    dpi=300,
    bbox_inches='tight',
    format='png'
)

# 关闭图形
plt.close()

print(f"图片已保存至：{save_path}")
    