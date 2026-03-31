import os
import numpy as np
import matplotlib.pyplot as plt

# 1) 读数据
#    每行：l real(diff) imag(diff)
data = np.loadtxt('./data/dsum_diff.txt')
idx   = data[:, 0]
diff_re = data[:, 1]
diff_im = data[:, 2]

# 2) 作图
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

ax1.plot(idx, diff_re, label='Re(dsum - dsum2)')
ax1.set_ylabel('Real part')
ax1.legend()
ax1.grid(True)

ax2.plot(idx, diff_im, label='Im(dsum - dsum2)')
ax2.set_xlabel('Element index l')
ax2.set_ylabel('Imag part')
ax2.legend()
ax2.grid(True)

plt.suptitle('Difference between numerical dsum and analytical dsum2')
plt.tight_layout(rect=[0, 0, 1, 0.95])

# 3) 保存到文件夹
save_dir = 'pic'
os.makedirs(save_dir, exist_ok=True)
out_path = os.path.join(save_dir, 'dsum_diff.png')
plt.savefig(out_path)
plt.close(fig)
# 计算总和
total_re = np.sum(diff_re)
total_im = np.sum(diff_im)
total = total_re + 1j * total_im
print(f"Plot saved to {out_path}")
print(f"Sum of all (dsum - dsum2) over elements: {total.real:.6e} + {total.imag:.6e}j")