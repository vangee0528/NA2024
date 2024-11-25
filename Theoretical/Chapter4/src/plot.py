import numpy as np
import matplotlib.pyplot as plt

# Define functions
def cond_f(x):
    """Condition number for f(x) = 1 - e^(-x)."""
    return np.abs((x * np.exp(-x)) / (1 - np.exp(-x)))

def cond_A_upper_bound(x):
    """Estimated upper bound for cond_A(x)."""
    return np.exp(x) / x

# Generate x values avoiding division by zero
x = np.linspace(0.01, 1, 500)

# Compute cond_f(x) and the upper bound of cond_A(x)
cond_f_values = cond_f(x)
cond_A_values = cond_A_upper_bound(x)

# Plot results
plt.figure(figsize=(8, 6))
plt.plot(x, cond_f_values, label=r'$\text{cond}_f(x)$', color='blue', linewidth=2)
plt.plot(x, cond_A_values, label=r'Estimated upper bound of $\text{cond}_A(x)$', color='red', linestyle='--', linewidth=2)
plt.xlabel(r'$x$', fontsize=14)
plt.ylabel(r'Condition Number', fontsize=14)
plt.title(r'$\text{cond}_f(x)$ and Upper Bound of $\text{cond}_A(x)$ on $[0, 1]$', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True)
plt.show()
