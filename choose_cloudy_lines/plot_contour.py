import numpy as np
import matplotlib.pyplot as plt

#线：
Tc_Mo = 1.05
k_B = 8.617333262145 * 1e-5
Lamda_F = 0.462
n_Mo = 0.29 * 1e2
n_Cu = 0.125 * 1e2
d0 = 1./((np.pi / 2.) * k_B * Tc_Mo * Lamda_F**2. * n_Mo)
t = 0.135
t = 0.13
def Tc(d_Mo, d_Cu):
    power_index = (d_Cu * n_Cu) / (d_Mo * n_Mo)
    return Tc_Mo * (d_Mo / (d0 * 1.13 * (1. + 1. / power_index) * t))**power_index * 1000.

mesh_number = 1000
d_Mo = np.linspace(10., 80., mesh_number)
d_Cu = np.linspace(10., 400., mesh_number)
D_Cu, D_Mo = np.meshgrid(d_Cu, d_Mo)

fig = plt.figure(figsize=(6,5), dpi=144)
select = [50, 60, 70, 80, 90, 100, 110, 120, 130]
contourPlot = plt.contour(D_Cu, D_Mo, Tc(D_Mo, D_Cu), select)
levels = contourPlot.levels[1:-2:]
plt.clabel(contourPlot, levels, inline=True, fmt = 'Tc = %2.0f mK', fontsize=7)
plt.legend(loc='best')
plt.xlabel(r'$t_{\rm{Cu}}$ (nm)')
plt.ylabel(r'$t_{\rm{Mo}}$ (nm)')
x_major_locator = MultipleLocator(30)
y_major_locator = MultipleLocator(10)
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
plt.grid(True)
plt.show(fig)
plt.close(fig)


#实心：
def PSF(alpha, beta, long_angle, short_angle):
    alpha = np.deg2rad(alpha)
    beta = np.deg2rad(beta)
    long_angle = np.deg2rad(long_angle)
    short_angle = np.deg2rad(short_angle)
    psf = (1. - abs(np.tan(alpha))/np.tan(long_angle)) \
          *(1. - abs(np.tan(beta))/np.tan(short_angle)) \
          / np.sqrt(np.tan(alpha)**2 + np.tan(beta)**2 + 1)
    psf[(np.abs(alpha) > long_angle) | (np.abs(beta)> short_angle)] = 0
    return psf

a = np.arange(-6.7, 6.7, 0.05)
b = np.arange(-6.7, 6.7, 0.05)
a, b = np.meshgrid(a, b)
psf = PSF(a, b, 5.7, 5.7)
plt.contourf(a, b, psf)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()