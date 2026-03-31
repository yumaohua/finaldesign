import numpy as np
import matplotlib.pyplot as plt

# 读取已过滤的中心点数据（假设文件中仅包含目标区域内的点）
data = np.loadtxt('./dataf/center_points.txt')
x = data[:, 0]
y = data[:, 1]

# 定义绘图区域参数（与 C++ 中定义的区域一致）
R = 10.0
center_x = R * np.cos(np.radians(60)) - 0.9 * R  # 60 度转换为弧度
center_y = R * np.sin(np.radians(60)) + 0.25 * R
width = height = R / 250  # 宽度和高度

# 创建图形
plt.figure(figsize=(8, 8))

# 绘制中心点（红色圆圈）
plt.scatter(x, y, s=15, edgecolor='red', facecolor='none', label='Filtered Center Points')

# 按顺序连接中点（单条蓝色曲线，假设点顺序正确）
plt.plot(x, y, 'b-', linewidth=0.8, alpha=0.6, label='Connected Path')

# 设置绘图范围（严格匹配 C++ 过滤区域）
plt.xlim(center_x - width/2, center_x + width/2)
plt.ylim(center_y - height/2, center_y + height/2)

# 添加标签和标题
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('Filtered Center Points and Connecting Path\n(Rounded Hexagon Corner Region)')

# 显示图例和网格
plt.legend()
plt.grid(True, linestyle='--', alpha=0.3)

# 保存图片（高质量，适合打印）
plt.savefig('./picf/filtered_center_points.png', dpi=600, bbox_inches='tight')

# 关闭图形（非必需，但释放内存）
plt.close()

print("图片已保存，包含过滤后的点及单条连接曲线。")