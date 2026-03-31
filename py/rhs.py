import numpy as np
import matplotlib.pyplot as plt
import os

# 创建保存图片的文件夹
# os.makedirs(save_dir, exist_ok=True)

# 1) 读数据（跳过 # 开头的注释行）
data = np.loadtxt('./data/rhs_compare.txt', comments='#')
idx        = data[:, 0]
rhs_re     = data[:, 1]
rhs_ori_re = data[:, 2]
rhs_im     = data[:, 3]
rhs_ori_im = data[:, 4]

# 2) 作图
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6))

ax1.plot(idx, rhs_re,     label='rhs real')
ax1.plot(idx, rhs_ori_re, label='rhs_ori real', linestyle='--')
ax1.set_ylabel('Real part')
ax1.legend()
ax1.grid(True)

ax2.plot(idx, rhs_im,     label='rhs imag')
ax2.plot(idx, rhs_ori_im, label='rhs_ori imag', linestyle='--')
ax2.set_xlabel('Element index')
ax2.set_ylabel('Imag part')
ax2.legend()
ax2.grid(True)

plt.suptitle('Comparison of RHS real & imag parts')
plt.tight_layout(rect=[0, 0, 1, 0.96])

# 保存图片到 figures/rhs_compare.png
save_path = os.path.join('./pic/rhs_compare.png')
plt.savefig(save_path, dpi=300)
