# # plot_midpoints.py
# import numpy as np
# import matplotlib.pyplot as plt

# # 读文件，把空行当分组标志
# pts1 = []
# pts2 = []
# with open("./data/cavities_midpoints.txt") as f:
#     group = pts1
#     for line in f:
#         s = line.strip()
#         if not s:
#             group = pts2
#             continue
#         x,y = map(float, s.split(","))
#         group.append((x,y))

# pts1 = np.array(pts1)
# pts2 = np.array(pts2)

# fig, ax = plt.subplots(figsize=(6,6))
# # 连线画出边界近似
# ax.plot(pts1[:,0], pts1[:,1], "-", lw=1, label="Cavity 1")
# ax.plot(pts2[:,0], pts2[:,1], "-", lw=1, label="Cavity 2")
# # 中心点标记
# ax.plot(-9,  2.5, "ro", label="Center 1")
# ax.plot( 9, -2.5, "bo", label="Center 2")

# ax.set_aspect("equal", "box")
# ax.set_xlabel("x / R")
# ax.set_ylabel("y / R")
# ax.legend(loc="upper right")
# plt.tight_layout()

# # 保存矢量图
# plt.savefig("./pic/cavities_midpoints.pdf")
# plt.savefig("./pic/cavities_midpoints.svg")
# print("Saved cavities_midpoints.pdf / .svg")
# plot_midpoints.py

import numpy as np
import matplotlib.pyplot as plt

# 从文件读入两组中点
pts1 = []
pts2 = []
with open("./data/cavities_midpoints.txt") as f:
    group = pts1
    for line in f:
        s = line.strip()
        if not s:
            group = pts2
            continue
        x, y = map(float, s.split(","))
        group.append((x, y))

pts1 = np.array(pts1)
pts2 = np.array(pts2)

fig, ax = plt.subplots(figsize=(6, 6))

# 只绘制散点，不连线
ax.scatter(pts1[:, 0], pts1[:, 1],
           s=5,          # 点的大小，可根据需要调整
           c='C0',       # 颜色
           marker='o',
           label="Cavity 1 midpoints",
           alpha=0.7)

ax.scatter(pts2[:, 0], pts2[:, 1],
           s=5,
           c='C1',
           marker='o',
           label="Cavity 2 midpoints",
           alpha=0.7)

# 绘制中心
ax.scatter([-0.9 * 10], [0.25 * 10], s=50, c='red', marker='x', label="Center 1")
ax.scatter([ 0.9 * 10], [-0.25 * 10], s=50, c='blue', marker='x', label="Center 2")

ax.set_aspect("equal", "box")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.legend(loc="best")
plt.tight_layout()

# 保存为矢量图
plt.savefig("./pic/cavities_midpoints.pdf")
plt.savefig("./pic/cavities_midpoints.svg")
print("Saved cavities_midpoints.pdf and .svg")
