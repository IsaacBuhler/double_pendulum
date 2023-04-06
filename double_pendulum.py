import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import animation

# Set the working directory where the map will be saved
os.chdir("C:/Users/Isaac/Documents/PROJECTS/PYTHON_PROGRAMMS/double_pendulum")

# Define the parameters of the double pendulum
m1 = 1.0  # mass of the first pendulum
m2 = 4.2  # mass of the second pendulum
L1 = 1.0  # length of the first pendulum
L2 = 3.0  # length of the second pendulum
g = 9.81  # gravitational constant

# Define the initial conditions
theta1_0 = np.pi/2   # initial angle of the first pendulum (in radians)
theta2_0 = np.pi/2   # initial angle of the second pendulum (in radians)
omega1_0 = 0.0       # initial angular velocity of the first pendulum (in radians/s)
omega2_0 = 0.0       # initial angular velocity of the second pendulum (in radians/s)

# Define the time vector
t0 = 0      # initial time
tf = 30     # final time
dt = 0.01   # time step
t = np.arange(t0, tf, dt)

# Define the function to calculate the derivatives
def f(theta1, theta2, omega1, omega2):
    dtheta1_dt = omega1
    dtheta2_dt = omega2
    domega1_dt = (-g*(2*m1 + m2)*np.sin(theta1) - m2*g*np.sin(theta1 - 2*theta2) - 2*np.sin(theta1 - theta2)*m2*(omega2**2*L2 + omega1**2*L1*np.cos(theta1 - theta2))) / (L1*(2*m1 + m2 - m2*np.cos(2*theta1 - 2*theta2)))
    domega2_dt = (2*np.sin(theta1 - theta2)*((omega1**2)*L1*(m1 + m2) + g*(m1 + m2)*np.cos(theta1) + (omega2**2)*L2*m2*np.cos(theta1 - theta2))) / (L2*(2*m1 + m2 - m2*np.cos(2*theta1 - 2*theta2)))
    return dtheta1_dt, dtheta2_dt, domega1_dt, domega2_dt

# Use the Runge-Kutta method to solve the differential equations
theta1 = np.zeros_like(t)
theta2 = np.zeros_like(t)
omega1 = np.zeros_like(t)
omega2 = np.zeros_like(t)
theta1[0] = theta1_0
theta2[0] = theta2_0
omega1[0] = omega1_0
omega2[0] = omega2_0

for i in range(1, len(t)):
    k1_theta1, k1_theta2, k1_omega1, k1_omega2 = f(theta1[i-1], theta2[i-1], omega1[i-1], omega2[i-1])
    k2_theta1, k2_theta2, k2_omega1, k2_omega2 = f(theta1[i-1] + 0.5*dt*k1_theta1, theta2[i-1] + 0.5*dt*k1_theta2, omega1[i-1] + 0.5*dt*k1_omega1, omega2[i-1] + 0.5*dt*k1_omega2)
    k3_theta1,k3_theta2, k3_omega1, k3_omega2 = f(theta1[i-1] + 0.5*dt*k2_theta1, theta2[i-1] + 0.5*dt*k2_theta2, omega1[i-1] + 0.5*dt*k2_omega1, omega2[i-1] + 0.5*dt*k2_omega2)
    k4_theta1, k4_theta2, k4_omega1, k4_omega2 = f(theta1[i-1] + dt*k3_theta1, theta2[i-1] + dt*k3_theta2, omega1[i-1] + dt*k3_omega1, omega2[i-1] + dt*k3_omega2)
    theta1[i] = theta1[i-1] + (dt/6)*(k1_theta1 + 2*k2_theta1 + 2*k3_theta1 + k4_theta1)
    theta2[i] = theta2[i-1] + (dt/6)*(k1_theta2 + 2*k2_theta2 + 2*k3_theta2 + k4_theta2)
    omega1[i] = omega1[i-1] + (dt/6)*(k1_omega1 + 2*k2_omega1 + 2*k3_omega1 + k4_omega1)
    omega2[i] = omega2[i-1] + (dt/6)*(k1_omega2 + 2*k2_omega2 + 2*k3_omega2 + k4_omega2)

x1 = L1*np.sin(theta1)
y1 = -L1*np.cos(theta1)
x2 = x1 + L2*np.sin(theta2)
y2 = y1 - L2*np.cos(theta2)

fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)


line, = ax.plot([], [], 'o-', lw=2)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    line.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
    return line,

ani = animation.FuncAnimation(fig, animate, frames=len(t), interval=dt*1000, blit=True, init_func=init)
plt.show()
'''

# Plot the trajectory of the double pendulum
ax.plot(x2, y2, '-r', lw=2)

# Set the axis labels
ax.set_xlabel('x')
ax.set_ylabel('y')

# Set the title of the plot
ax.set_title('Double Pendulum Motion, m1= '+str(m1)+', m2 = '+str(m2))

# Show the plot
plt.savefig('double_pendulum_motion.png')
plt.show()

'''