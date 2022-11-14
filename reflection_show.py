from Reflections import *

##Settings for how to show plot(s)
raysOnly = True
simpleView = False
indivPlots = False

# curv1 = Parabola(1,-2,-5,np.array([-2,-2]),np.array([1,2])) #a,b,c,initPoint,initDir
# nextPoint,nextDir = curv1.reflect_vert(raysOnly, simpleView, indivPlots)	

# curv2 = Parabola(1,2,4,np.array([-6,-2]),np.array([9,1])) 
# nextPoint,nextDir = curv2.reflect_horiz(raysOnly, simpleView, indivPlots)

# curv3 = Rotated_Parabola(1,-2,-5,np.array([-2,-2]),np.array([1,2])) #a,b,c,initPoint,initDir
# nextPoint,nextDir = curv3.reflect(raysOnly, simpleView, indivPlots)	

# curv3 = Linear(1,0,np.array([-1,-3]),np.array([1,2])) #a,b,initPoint,initDir
# nextPoint,nextDir = curv3.reflect(raysOnly, simpleView, indivPlots)

# curv4 = Hyperbola(1,2,1,3,np.array([1,-2]),np.array([1,4])) #a,b,h,k,initPoint,initDir
# nextPoint,nextDir = curv4.reflect(raysOnly, simpleView, indivPlots)

# curv4 = Rotated_Hyperbola(1,2,1,3,np.array([1,-2]),np.array([1,4]))
# nextPoint,nextDir = curv4.reflect(raysOnly, simpleView, indivPlots)

# curv5 = Ellipse(1,2,1,3, np.array([0,0]),np.array([1,4])) 
# nextPoint,nextDir = curv5.reflect(raysOnly, simpleView, indivPlots)

curv5 = Rotated_Ellipse(1,2,-1,3, np.array([0,0]),np.array([1,4])) 
nextPoint,nextDir = curv5.reflect(raysOnly, simpleView, indivPlots)

# curv3 = Linear(-3,5,nextPoint,nextDir) #a,b,initPoint,initDir
# nextPoint,nextDir = curv3.reflect(raysOnly, simpleView, indivPlots)

# curv4 = Rotated_Hyperbola(3,2,1,3,nextPoint,nextDir)
# nextPoint,nextDir = curv4.reflect(raysOnly, simpleView, indivPlots)


if not(indivPlots):
    t = np.linspace(0, 10, 500)
    plt.plot(nextPoint[0] + t*nextDir[0], nextPoint[1] + t*nextDir[1],'green')

    ##Print plot
    plt.grid(color='lightgray',linestyle='--')
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.gca().set_aspect('equal', adjustable='box')
    # plt.savefig('/Users/epeairs/Desktop/chain3.pdf', bbox_inches='tight')
    plt.show()