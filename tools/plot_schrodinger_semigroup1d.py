#Plots 3 graphs in the following order: error, solution and initial condition
# Use `q` to skip to the next figure

import numpy as np
import matplotlib.pyplot as plt

# Load error data from the file
real_data = np.loadtxt('../build/Re_error.line')
imag_data = np.loadtxt('../build/Im_error.line')

# Extract and plot values
x = real_data[:, 0]
y = real_data[:, 1]
print("Left real boundary value:")
print(y[0])
plt.plot(x, y, label='Re_error(x)')
# Extract and plot values
x = imag_data[:, 0]
y = imag_data[:, 1]
print("Left imaginary boundary value:")
print(y[0])
plt.plot(x, y, label='Im_error(x)')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Error')
plt.legend()
plt.grid(True)
plt.show()


# Load solution from the file
real_data = np.loadtxt('../build/Re_g_tree.line')
imag_data = np.loadtxt('../build/Im_g_tree.line')

# Extract and plot values
x = real_data[:, 0]
y = real_data[:, 1]
print("Left real boundary value:")
print(y[0])
plt.plot(x, y, label='Reg(x)')
# Extract and plot values
x = imag_data[:, 0]
y = imag_data[:, 1]
print("Left imaginary boundary value:")
print(y[0])
plt.plot(x, y, label='Img(x)')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Solution')
plt.legend()
plt.grid(True)
plt.show()


# Load initial data from the file
real_data = np.loadtxt('../build/Re_f_tree.line')
imag_data = np.loadtxt('../build/Im_f_tree.line')

# Extract and plot values
x = real_data[:, 0]
y = real_data[:, 1]
print("Left real boundary value:")
print(y[0])
plt.plot(x, y, label='Ref(x)')
# Extract and plot values
x = imag_data[:, 0]
y = imag_data[:, 1]
print("Left imaginary boundary value:")
print(y[0])
plt.plot(x, y, label='Imf(x)')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Initial condition')
plt.legend()
plt.grid(True)
plt.show()
