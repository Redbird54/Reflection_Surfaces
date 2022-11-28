from Reflections import *

##Settings for how to show plot(s)
raysOnly = True
indivPlots = False

nextPoint,nextDir = np.array([0,0]),np.array([1,4])

curv = Ellipse(1,2,0,2,nextPoint, nextDir,math.pi/3) 
nextPoint,nextDir = curv.reflect(raysOnly, indivPlots)

curv = Hyperbola(3,2,1,3,nextPoint,nextDir,math.pi/3)
nextPoint,nextDir = curv.reflect(raysOnly, indivPlots)

curv = Linear(nextPoint,nextDir, a=-3, b=5)
nextPoint,nextDir = curv.reflect(raysOnly, indivPlots)

curv = Parabola(1,4,-3,nextPoint,nextDir,math.pi/3)
nextPoint,nextDir = curv.reflect(raysOnly, indivPlots)


if not(indivPlots):
    t = np.linspace(0, 10, 500)
    plt.plot(nextPoint[0] + t*nextDir[0], nextPoint[1] + t*nextDir[1],'green')

    ##Print plot
    plt.grid(color='lightgray',linestyle='--')
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()