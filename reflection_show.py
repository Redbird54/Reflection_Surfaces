from Reflections import *

##Settings for how to show plot(s)
indivPlots = False
interactions = 6
boxsize = 5


nextPoint,nextDir = np.array([-10,-40]),np.array([1,4])
initialObjs = []
objs = []

initialObjs.append(Ellipse(1,2,1,8,boxsize,math.pi/3))
initialObjs.append(Hyperbola(3,2,-15,-3,boxsize,5*math.pi/6))
initialObjs.append(Linear(4,-9,4,3,boxsize))
initialObjs.append(Parabola(1,17,13,boxsize,5*math.pi/3))


for obj in initialObjs:
    if not(objs):
        if any(abs(nextPoint-obj.get_center()) > (boxsize*1.5)):
            objs.append(obj)
            obj.show_curve(indivPlots)
    else:
        if any(abs(nextPoint-obj.get_center()) > (boxsize*1.5)):
            notOverlap = True
            for x in range(len(objs)):
                centerDiff = [abs(a - b) >= (boxsize*3) for a, b in zip(obj.get_center(), objs[x].get_center())]
                if not(any(centerDiff)):
                    notOverlap = False
            if notOverlap:
                objs.append(obj)
                obj.show_curve(indivPlots)

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
    t = np.linspace(0, 20, 500)
    plt.plot(nextPoint[0] + t*nextDir[0], nextPoint[1] + t*nextDir[1],'green')

    ##Print plot
    plt.grid(color='lightgray',linestyle='--')
    plt.xlim(-30, 30)
    plt.ylim(-30, 30)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()