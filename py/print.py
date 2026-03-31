import numpy as np
import matplotlib.pyplot as plt

def plot_cross_section(filename="./data/cross_section.txt"):
    data = np.loadtxt(filename)
    k_vals, sigma_vals = data[:, 0], data[:, 1]

    plt.figure(figsize=(8, 4))
    plt.plot(k_vals, sigma_vals, label=r"$\sigma(k)$", color="C0")
    plt.xlabel("Wavenumber $k$")
    plt.ylabel("Total cross section $\sigma(k)$")
    plt.title("Scattering Cross Section")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("./pic/cross_section.png", dpi=300)
    plt.show()


# def plot_far_field(filename="far_field.txt"):
#     data = np.loadtxt(filename)
#     thetas = data[:, 0]
#     f_vals = data[:, 1] + 1j * data[:, 2]
#     intensity = np.abs(f_vals)**2

#     fig = plt.figure(figsize=(6,6))
#     ax = fig.add_subplot(111, polar=True)
#     ax.plot(thetas, intensity, label=r"$|f(\theta)|^2$", color="C1")
#     ax.set_title("Far-field Scattering Pattern", va='bottom')
#     ax.grid(True)
#     ax.legend()
#     plt.tight_layout()
#     plt.savefig("far_field.png", dpi=300)
#     plt.show()


if __name__ == "__main__":
    plot_cross_section()
    # plot_far_field()
