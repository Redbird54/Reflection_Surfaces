from Reflections import *

##Settings for how to show plot(s)
raysOnly = True
indivPlots = False
interactions = 4


nextPoint,nextDir = np.array([0,0]),np.array([1,4])
objs = []

objs.append(Ellipse(1,2,-1,3,math.pi/3))
objs.append(Hyperbola(3,2,1,3,math.pi/3))
objs.append(Linear(a=-3, b=5))
# objs.append(Linear(point=np.array([0,5]), dir=np.array([1,-3]))
objs.append(Parabola(1,4,-3,math.pi/3))

for x in range(interactions):
    distSmall = -1
    currObj = Object()
    for obj in objs:
        dist = obj.get_distance(nextPoint,nextDir)
        if distSmall == -1 and dist > 0:
            distSmall = dist
            currObj = obj
        elif dist == -1:
            distSmall = distSmall
        else:
            if dist < distSmall:
                distSmall = dist
                currObj = obj

    nextPoint,nextDir = currObj.reflect(distSmall, nextPoint, nextDir, raysOnly, indivPlots)


if not(indivPlots):
    t = np.linspace(0, 10, 500)
    plt.plot(nextPoint[0] + t*nextDir[0], nextPoint[1] + t*nextDir[1],'green')

    ##Print plot
    plt.grid(color='lightgray',linestyle='--')
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()