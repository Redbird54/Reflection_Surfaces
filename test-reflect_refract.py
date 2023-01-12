from reflect_refract import *
import queue

##Hyperparameters for model
indivPlots = False
interactions = 8
boxsize = 5

##Setting the medium 1 to air, medium 2 to glass
n1 = 1.0003
n2 = 1.52

outputs = queue.Queue()
initialObjs = []
objs = []

nextPoint, nextRefl = np.array([-10,-40]),np.array([1,4])
initialObjs.append(Ellipse(1,2,1,8,boxsize,"reflection",math.pi/3))
initialObjs.append(Hyperbola(3,2,-15,-3,boxsize,"reflection",5*math.pi/6))
initialObjs.append(Linear(4,-9,4,3,boxsize,"reflection"))
initialObjs.append(Parabola(1,-15,13,boxsize,"refraction",5*math.pi/3))
# initialObjs.append(Lens(-15,13,7,10,5,boxsize,9*math.pi/20))
# initialObjs.append(Linear_Lens(-15,13,4,boxsize,18*math.pi/20))

for obj in initialObjs:
    if not(objs):
        if any(abs(nextPoint-obj.get_center()) > (boxsize*1.5)):
            objs.append(obj)
            if not(indivPlots):
                obj.show_curve()
    else:
        if any(abs(nextPoint-obj.get_center()) > (boxsize*1.5)):
            notOverlap = True
            for x in range(len(objs)):
                centerDiff = [abs(a - b) >= (boxsize*3) for a, b in zip(obj.get_center(), objs[x].get_center())]
                if not(any(centerDiff)):
                    notOverlap = False
            if notOverlap:
                objs.append(obj)
                if not(indivPlots):
                    obj.show_curve()


outputs.put([nextPoint,nextRefl,n1,n2])

for x in range(interactions):
    if outputs.empty():
        continue
    currInfo = outputs.get()
    distSmall = -1
    currObj = Object()
    for obj in objs:
        dist = obj.get_distance(currInfo[0],currInfo[1])
        if distSmall == -1 and dist > 0:
            distSmall = dist
            currObj = obj
        elif dist == -1:
            distSmall = distSmall
        else:
            if dist < distSmall:
                distSmall = dist
                currObj = obj

    nextRays = currObj.output(distSmall, currInfo[0], currInfo[1], currInfo[2], currInfo[3], indivPlots)
    for nextRay in nextRays: 
        if not(np.array_equal(nextRay[1],currInfo[1])):
            outputs.put([nextRay[0],nextRay[1],currInfo[2],currInfo[3]])
        elif (np.array_equal(nextRay[2],currInfo[1])):
            ##Show rays not interacting with any curves here
            t = np.linspace(0, 30, 500)
            plt.plot(nextRay[0][0] + t*nextRay[1][0], nextRay[0][1] + t*nextRay[1][1],'orange')
        if not(np.array_equal(nextRay[2],currInfo[1])):
            if not(any(np.isnan(nextRay[2]))):
                outputs.put([nextRay[0],nextRay[2],currInfo[3],currInfo[2]])



if not(indivPlots):
    ##Show rays currently in the queue here
    t = np.linspace(0, 30, 500)
    for x in range(outputs.qsize()):
        output = outputs.get()
        plt.plot(output[0][0] + t*output[1][0], output[0][1] + t*output[1][1],'green')

    ##Print plot
    plt.grid(color='lightgray',linestyle='--')
    plt.xlim(-30, 30)
    plt.ylim(-30, 30)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()