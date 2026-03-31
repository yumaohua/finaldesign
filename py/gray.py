# py/gray_fixed.py
import numpy as np
import matplotlib.pyplot as plt

# 1) 读入数据（跳过第一行表头）
data = np.loadtxt("dataf/compute_psi2_grid.txt", delimiter=",")
# data 按 x 先后、y 紧跟的顺序排列

# 2) 重建 x,y 坐标向量
xs = np.unique(data[:,0])   # 长度 400
ys = np.unique(data[:,1])   # 长度 300
nx, ny = xs.size, ys.size

# 3) 按 (nx,ny) 重塑后再转置，得到 (ny,nx) 矩阵
Z = data[:,2].reshape((nx, ny)).T
vmin = Z.min()
vmax = Z.max()
# 4) 画图
fig, ax = plt.subplots(figsize=(8,6), dpi=800)
im = ax.imshow(
    Z,
    extent=(xs[0], xs[-1], ys[0], ys[-1]),
    origin="lower",
    cmap="gray",
    interpolation="none",
    aspect="equal",
    vmin=0,
    vmax=0.3,
)

# 5) 叠加离散元中点
mid = np.loadtxt("dataf/elements_midpoints.txt", delimiter=",")
ax.scatter(mid[:,0], mid[:,1], c="k", s=1)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Near-field $|\\psi|^2$")
plt.colorbar(im, ax=ax, label="$|\\psi|^2$")
plt.tight_layout()
plt.savefig("pic/psi2_grid_fixed.png")
print("Saved psi2_grid_fixed.png")
