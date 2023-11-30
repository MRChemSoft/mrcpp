#Plots 3 graphs in the following order: error, solution and initial condition
#in ONE WINDOW

import numpy as np
import matplotlib.pyplot as plt

# Load error data from the file
real_data = np.loadtxt('../build/Re_error.line')
imag_data = np.loadtxt('../build/Im_error.line')

# Plot error
plt.subplot(3, 1, 1)  # 3 rows, 1 column, first plot
plt.plot(real_data[:, 0], real_data[:, 1], label='Re_error(x)')
plt.plot(imag_data[:, 0], imag_data[:, 1], label='Im_error(x)')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Error')
plt.legend()
plt.grid(True)

# Load solution from the file
real_data = np.loadtxt('../build/Re_g_tree.line')
imag_data = np.loadtxt('../build/Im_g_tree.line')

# Plot solution
plt.subplot(3, 1, 2)  # 3 rows, 1 column, second plot
plt.plot(real_data[:, 0], real_data[:, 1], label='Reg(x)')
plt.plot(imag_data[:, 0], imag_data[:, 1], label='Img(x)')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Solution')
plt.legend()
plt.grid(True)

# Load initial data from the file
real_data = np.loadtxt('../build/Re_f_tree.line')
imag_data = np.loadtxt('../build/Im_f_tree.line')

# Plot initial condition
plt.subplot(3, 1, 3)  # 3 rows, 1 column, third plot
plt.plot(real_data[:, 0], real_data[:, 1], label='Ref(x)')
plt.plot(imag_data[:, 0], imag_data[:, 1], label='Imf(x)')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Initial condition')
plt.legend()
plt.grid(True)

plt.tight_layout()  # Adjust layout for better spacing
plt.show()
