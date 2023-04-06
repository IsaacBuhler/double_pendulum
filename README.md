# double_pendulum
Modelling semi-chaotic double pendulum system using Runge-Kutta method

The motion of the two pendula can be described using the following equations:
$$\ddot{\theta}_1 = \frac{-g(2m_1+m_2)\sin{\theta_1} - m_2g\sin{(\theta_1 - 2\theta_2)} - 2\sin{(\theta_1 - \theta_2)}m_2(\dot{\theta}_2^2 L_2 + \dot{\theta}_1^2 L_1\cos{(\theta_1 - \theta_2)})}{L_1(2m_1 + m_2 - m_2\cos{(2\theta_1 - 2\theta_2)})} $$

$$\ddot{\theta}_2 = \frac{2\sin{(\theta_1 - \theta_2)}(\dot{\theta}_1^2 L_1(m_1+m_2) + g(m_1+m_2)\cos{\theta_1} + \dot{\theta}_2^2 L_2 m_2 \cos{(\theta_1 - \theta_2)})}{L_2(2m_1 + m_2 - m_2\cos{(2\theta_1 - 2\theta_2)})} $$

where $\theta_1$ and $\theta_2$ are the angular displacements of the first and second pendulums from the vertical, $\dot{\theta}_1$ and $\dot{\theta}_2$ are the corresponding angular velocities, $m_1$ and $m_2$ are the masses of the first and second pendulums, $L_1$ and $L_2$ are their respective lengths, and $g$ is the acceleration due to gravity.

To numerically solve the equations of motion, can use the fourth-order Runge-Kutta method:

$$
\begin{align}
k_1 &= hf(t_n, y_n) \\
k_2 &= hf(t_n + \frac{h}{2}, y_n + \frac{k_1}{2}) \\
k_3 &= hf(t_n + \frac{h}{2}, y_n + \frac{k_2}{2}) \\
k_4 &= hf(t_n + h, y_n + k_3) \\
y_{n+1} &= y_n + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)
\end{align}
$$

where $h$ is the time step, $t_n$ and $y_n$ are the time and state (i.e., the angular displacements and velocities) at step $n$, and $k_1$, $k_2$, $k_3$, and $k_4$ are intermediate values calculated at each step.

Using these equations and the Runge-Kutta method, we can simulate the motion of the double pendulum system and plot its trajectory over time!.
