# plot_psi2_circle.py
import numpy as np
import matplotlib.pyplot as plt
import os

# 1) 读取 (theta,psi2)，跳过表头
data = np.loadtxt("dataf/compute_psi2_circle.txt",
                  delimiter=",",
                  skiprows=1)  # <-- 跳过表头

theta = data[:,0]
psi2  = data[:,1]
deg   = theta * 180.0 / np.pi

# 2) 绘图
fig, ax = plt.subplots(figsize=(8,4))
ax.plot(deg, psi2, "-", lw=1)
ax.set_xlabel(r"$\theta$ (deg)")
ax.set_ylabel(r"$|\psi|^2$")
ax.set_title("Far‐field $|\psi|^2$ vs. θ")
ax.grid(True)
plt.tight_layout()

os.makedirs("pic", exist_ok=True)
plt.savefig("pic/psi2_circle.png", dpi=300)
plt.close()
print("Saved pic/psi2_circle.png")
