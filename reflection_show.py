from Reflections import *

##Settings for how to show plot(s)
raysOnly = True
indivPlots = False

nextPoint,nextDir = np.array([0,0]),np.array([1,4])


curv = Rotated_Ellipse(1,2,-1,3,math.pi/3,nextPoint, nextDir) 
nextPoint,nextDir = curv.reflect(raysOnly, indivPlots)

curv = Linear(-3,5,nextPoint,nextDir)
nextPoint,nextDir = curv.reflect(raysOnly, indivPlots)

curv = Rotated_Hyperbola(3,2,1,3,math.pi/3,nextPoint,nextDir)
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