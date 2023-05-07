from reflect_refract3d import *
from encrypt_ecc import *
import queue

# Plotting shape
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

##Hyperparameters for model
indivPlots = False
interactions = 12
boxsize = 5
boxedge = False
intensity = 1

n1 = 1.0003
n2 = 1.52


outputs = queue.Queue()
reverse = queue.Queue()
initialObjs = []
objs = []


startPoint, startDir = np.array([0,-8,0.5]),np.array([0, 1, 0])

initialObjs.append(Ellipsoid(2,1,1,0,0,0,boxsize,"reflection",math.pi/3,3,-math.pi/4))
initialObjs.append(Ellipsoid(2,1,1,-16.365,-3.73,14.97,boxsize,"reflection",math.pi/3,math.pi/2,-math.pi/4))
# initialObjs.append(Polynomial2(2,1,0,3,5,-2,0,-1,3,2, 0,0,0.5,boxsize,"reflection",0,0,0))
# initialObjs.append(Polynomial3(2,1,0,3,5,-2,0,-1,3,2,4,-6,1,0,1,-5,3,3,-2, 0,0,0,boxsize,"reflection",0,0,0))

def get_valid_objs(initialObjs, objs, startPoint):
    for obj in initialObjs:
        if not(objs):
            if any(abs(startPoint-obj.get_center()) > (boxsize*1.5)):
                objs.append(obj)
                if not(indivPlots):
                    obj.show_curve(ax)
        else:
            if any(abs(startPoint-obj.get_center()) > (boxsize*1.5)):
                notOverlap = True
                for x in range(len(objs)):
                    centerDiff = [abs(a - b) >= (boxsize*3) for a, b in zip(obj.get_center(), objs[x].get_center())]
                    if not(any(centerDiff)):
                        notOverlap = False
                if notOverlap:
                    objs.append(obj)
                    if not(indivPlots):
                        obj.show_curve(ax)
    return objs


def find_rays(outputs, objs, numbInteractions):
    outPoint = []
    outRay = []
    outRay2 = []
    isBoxEdge = 0
    for x in range(numbInteractions):
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

        actionCount = x
        if x == 0:
            initMag = distSmall
        if currObj.get_type() == "none":
            ##Show rays not interacting with any curves here
            t = np.linspace(0, 5, 500)
            ax.plot(currInfo[0][0] + t*currInfo[1][0], currInfo[0][1] + t*currInfo[1][1], currInfo[0][2] + t*currInfo[1][2], 'orange')
            break
        else:
            nextRays = currObj.output(distSmall, currInfo[0], currInfo[1], currInfo[2], currInfo[3], currInfo[4], boxedge, ax, indivPlots)
            for nextRay in nextRays:
                outPoint = nextRay[0]
                outRay = currInfo[1]
                isBoxEdge = nextRay[4]
                if currObj.get_type() == "reflection":  #Assume always reflecting for now
                    outRay2 = nextRay[1]
                elif currObj.get_type() == "refraction":
                    outRay2 = nextRay[2]
                    
                if np.array_equal(nextRay[1], currInfo[1]) and np.array_equal(nextRay[2], currInfo[1]):
                    outputs.put([nextRay[0], nextRay[2], currInfo[3], currInfo[2], nextRay[2]])
                else:
                    if not(np.array_equal(nextRay[1], currInfo[1])):
                        outputs.put([nextRay[0], nextRay[1], currInfo[2], currInfo[3], nextRay[3]])
                    if not(np.array_equal(nextRay[2], currInfo[1])):
                        if not(any(np.isnan(nextRay[2]))):
                            outputs.put([nextRay[0], nextRay[2], currInfo[3], currInfo[2], nextRay[3]])
    return outputs, outPoint, outRay, outRay2, initMag, isBoxEdge, actionCount

objs = get_valid_objs(initialObjs, objs, startPoint)
outputs.put([startPoint, startDir, n1, n2, intensity])
outputs, outPoint, outRay, outRay2, mag, isBoxEdge, actionCount = find_rays(outputs, objs, interactions)
if actionCount > 0:
    if boxedge and isBoxEdge:
        inRay = -outRay2
    else:
        inRay = -outRay
    test1 = np.array([float(encrypt(str(outPoint[0]))), float(encrypt(str(outPoint[1]))), float(encrypt(str(outPoint[2])))])
    test2 = np.array([float(encrypt(str(inRay[0]))), float(encrypt(str(inRay[1]))), float(encrypt(str(inRay[2])))])
    test3 = float(encrypt(str(mag)))

    reverse.put([test1, test2, n1, n2, intensity])
    if actionCount > 1:
        tests, testPoint, testRay, testRay2, testMag, testFlag, testCount = find_rays(reverse, objs, actionCount)
    else: 
        testPoint = test1
        testRay2 = test2
    t = np.linspace(0, test3, 500)
    print(np.array([testPoint[0] + test3*testRay2[0], testPoint[1] + test3*testRay2[1], testPoint[2] + test3*testRay2[2]]))
    # plt.plot(testPoint[0] + t*testRay2[0], testPoint[1] + t*testRay2[1], 'orange')
    ax.plot(testPoint[0] + t*testRay2[0], testPoint[1] + t*testRay2[1], testPoint[2] + t*testRay2[2], 'orange')
else:
    raise Exception("No interactions occurred, please adjust input parameters")

if not(indivPlots):
    ##Show rays currently in the queue here
    t = np.linspace(0, 5, 500)
    for x in range(outputs.qsize()):
        output = outputs.get()
        plt.plot(output[0][0] + t*output[1][0], output[0][1] + t*output[1][1],output[0][2] + t*output[1][2],'green')
        # ax.plot(output[0][0] + t*output[1][0], output[0][1] + t*output[1][1],output[0][2] + t*output[1][2],'green')
    ax.set_xlim(-10,10)
    ax.set_ylim(-10,10)
    ax.set_zlim(-10,10)
    ax.set_aspect('equal')
    plt.show()
