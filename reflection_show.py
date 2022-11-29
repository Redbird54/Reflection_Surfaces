from Reflections import *

##Settings for how to show plot(s)
indivPlots = False
interactions = 4


nextPoint,nextDir = np.array([0,0]),np.array([1,4])
initialObjs = []
objs = []

initialObjs.append(Ellipse(1,2,-1,3,math.pi/3))
initialObjs.append(Hyperbola(3,2,1,3,math.pi/3))
initialObjs.append(Linear(0, 5, a=-3, b=5))
# initialObjs.append(Linear(0,5,point=np.array([0,5]), dir=np.array([1,-3]))
initialObjs.append(Parabola(1,5,-3,math.pi/3))


for obj in initialObjs:
    if not(objs):
        objs.append(obj)
    else:
        notOverlap = True
        for x in range(len(objs)):
            currCent = [abs(a - b) >= 5 for a, b in zip(obj.get_center(), objs[x].get_center())]
            if not(all(currCent)):
                notOverlap = False
        if notOverlap:
            objs.append(obj)

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

    nextPoint,nextDir = currObj.reflect(distSmall, nextPoint, nextDir, indivPlots)


if not(indivPlots):
    t = np.linspace(0, 10, 500)
    plt.plot(nextPoint[0] + t*nextDir[0], nextPoint[1] + t*nextDir[1],'green')

    ##Print plot
    plt.grid(color='lightgray',linestyle='--')
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()