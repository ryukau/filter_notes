import numpy as np
import matplotlib.pyplot as plt


def compute_stability_region(Q_val, resolution=200, limit=3):
    re = np.linspace(-limit, limit, resolution)
    im = np.linspace(-limit, limit, resolution)
    Re, Im = np.meshgrid(re, im)
    F = Re + 1j * Im

    Trace = 2 - F**2 - F / Q_val
    Det = 1 - F / Q_val

    delta = np.sqrt(Trace**2 - 4 * Det)
    lam1 = (Trace + delta) / 2
    lam2 = (Trace - delta) / 2

    max_mag = np.maximum(np.abs(lam1), np.abs(lam2))
    stable = max_mag < 1.0
    return Re, Im, stable


def plot_for_Q():
    # Plotting for a few Q values
    plt.figure(figsize=(12, 5))

    # Case 1: Moderate Q
    Q1 = 0.707
    Re1, Im1, Stable1 = compute_stability_region(Q1)
    plt.subplot(1, 2, 1)
    plt.pcolormesh(Re1, Im1, Stable1, cmap="RdYlGn", shading="auto", alpha=0.6)
    plt.colorbar(label="Stable")
    plt.title(f"Stability Region for Q={Q1}")
    plt.xlabel("Re(f)")
    plt.ylabel("Im(f)")
    plt.grid(True, alpha=0.3)
    plt.axhline(0, color="black", linewidth=0.5)
    plt.axvline(0, color="black", linewidth=0.5)

    # Case 2: Higher Q (Resonant)
    Q2 = 2.0
    Re2, Im2, Stable2 = compute_stability_region(Q2)
    plt.subplot(1, 2, 2)
    plt.pcolormesh(Re2, Im2, Stable2, cmap="RdYlGn", shading="auto", alpha=0.6)
    plt.colorbar(label="Stable")
    plt.title(f"Stability Region for Q={Q2}")
    plt.xlabel("Re(f)")
    plt.ylabel("Im(f)")
    plt.grid(True, alpha=0.3)
    plt.axhline(0, color="black", linewidth=0.5)
    plt.axvline(0, color="black", linewidth=0.5)

    plt.tight_layout()
    # plt.savefig("chamberlin_stability.png")


def plot_with_circle():
    plt.figure(figsize=(10, 8))

    # Use a representative Q
    Q_val = 2.0
    Re, Im, Stable = compute_stability_region(Q_val, limit=3.0)

    # Plot stability region
    plt.pcolormesh(Re, Im, Stable, cmap="Greys", shading="auto", alpha=0.3)
    plt.contour(Re, Im, Stable, levels=[0.5], colors="green", linewidths=2)

    # Overlay "Twist" trajectories
    # Twist preserves magnitude, changes angle.
    # f = f_mag * exp(j * twist)
    angles = np.linspace(-np.pi, np.pi, 200)

    # Frequencies to check (corresponding to different Cutoffs)
    # f_mag approx 2 * pi * cutoff_norm
    # Normalized cutoffs: 0.05, 0.1, 0.2, 0.3
    f_mags = [0.2, 0.5, 1.0, 1.5, 1.9]

    for f_mag in f_mags:
        traj_real = f_mag * np.cos(angles)
        traj_imag = f_mag * np.sin(angles)

        # Check intersection with stability to color code?
        # Simple plot for now
        plt.plot(traj_real, traj_imag, label=f"|f| = {f_mag}", linestyle="--")

    plt.title(
        f'Stability Region of Complex Chamberlin Filter (Q={Q_val})\nDashed Lines = "Twist" Trajectories at const freq'
    )
    plt.xlabel("Re(f)")
    plt.ylabel("Im(f)")
    plt.axhline(0, color="black", alpha=0.5)
    plt.axvline(0, color="black", alpha=0.5)
    plt.legend(loc="upper right")
    plt.grid(True, alpha=0.3)
    plt.axis("equal")

    # plt.savefig('chamberlin_complex_stability.png')


plot_for_Q()
plot_with_circle()
plt.show()
