import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages 
#*************************************************************************************
#Open files for graphocal output
pp = PdfPages('GraphOut.pdf')
#Fuctions for creating different types of plots and analysis
#To show th plot on the console during the calculations, activarecommand plt.show() below
def MA(Density,Counter,num_iterations,Window):
    Density1 = np.convolve(Counter, np.ones([Window]), mode='valid')/Window
    for i in range(0,num_iterations-Window):
        Density[i] = Density1[i]
    for i in range(num_iterations-Window,num_iterations):
        Density[i] = 0.0
    return;
def UserPlot2DScatter(x_pos, y_pos):
    plt.scatter(x_pos,y_pos, marker='o')
    pp.savefig()
    plt.show()
    return;
def UserPlot2Dline(nodenumber, z_pos,title,xlabel,ylabel):
    plt.plot(nodenumber,z_pos)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()
    pp.savefig()
    plt.show()
    return;
def QuiverPlot2D(x_pos, y_pos, x_direct, y_direct, color):
    fig, ax = plt.subplots()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(x_pos, y_pos, rstride=10, cstride=10)
    ax.quiver(x_pos, y_pos, x_direct, y_direct,color)
    ax.set_title('Quiver plot with Force Vactors at Each of the Nodes')
    pp.savefig()
    plt.show()
    return;
def QuiverPlot3D(x_pos, y_pos, z_pos, x_direct, y_direct, z_direct, color):
    fig, ax = plt.subplots()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(x_pos, y_pos, z_pos, rstride=10, cstride=10)
    ax.quiver(x_pos, y_pos, x_direct, y_direct,color)
    ax.quiver(x_pos, y_pos, z_pos, x_direct, y_direct, z_direct,color)
    ax.set_title('Quiver plot with Force Vactors at Each of the Nodes')
    pp.savefig()
    plt.show()
    fig = plt.figure()
    return;
def UserPlot3D(x_pos, y_pos, z_pos):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x_pos, y_pos, z_pos, marker='o')
    pp.savefig()
    plt.show()
    return;
def UserPlot3DSurf(x_pos, y_pos, z_pos):
    fig = plt.figure()
    ax = Axes3D(fig)
    surf = ax.plot_trisurf(x_pos, y_pos, z_pos, cmap=cm.jet, linewidth=0.1)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    pp.savefig()
    plt.show()
    return;
#*******************************************************************************************
