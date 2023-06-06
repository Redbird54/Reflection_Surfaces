from reflect_refract import *
from general_curves import *
from encrypt_ecc import *
import queue

##Hyperparameters for model
indivPlots = True
interactions = 12
boxsize = 5
boxedge = False
globalbox = 100
intensity = 1

##Setting the medium 1 to air, medium 2 to glass
n1 = 1.0003
n2 = 1.52

outputs = queue.Queue()
reverse = queue.Queue()
initialObjs = []
objs = []

if interactions < 1:
    raise Exception("interactions must be positive")

startPoint, startDir = np.array([-5,-20]),np.array([1,4])
# initialObjs.append(Polynomial3(5, 2, 3, -1, 7, 2, 0, -4, 1, -2, 0, 0,boxsize,"reflection"))
# initialObjs.append(Polynomial2(5, 3, 7, 5, -1, 0, 0, 0,boxsize,"reflection"))
# initialObjs.append(Ellipse(1,2,1,8,boxsize,"reflection",math.pi/3))
# initialObjs.append(Hyperbola(3,2,-15,-3,boxsize,"reflection",5*math.pi/6))
# initialObjs.append(Linear(4,-9,4,3,boxsize,"reflection"))
# initialObjs.append(Parabola(1,-15,13,boxsize,"reflection",5*math.pi/3))
# initialObjs.append(Polynomial3(5, 2, 3, -1, 7, 2, 0, -4, 1, -2, 4, -9,boxsize,"reflection", math.pi/2))
initialObjs.append(Parabola(1,2,8,boxsize,"reflection",-math.pi/6))
# initialObjs.append(Lens(-15,13,7,10,5,boxsize,9*math.pi/20))
# initialObjs.append(Linear_Lens(-15,13,4,boxsize,18*math.pi/20))

def get_valid_objs(initialObjs, objs, startPoint):
    for obj in initialObjs:
        if not(objs):
            if any(abs(startPoint-obj.get_center()) > (boxsize*1.5)) and all(np.abs(obj.get_center())+(boxsize*1.5) <= globalbox):
                objs.append(obj)
                if not(indivPlots):
                    obj.show_curve()
        else:
            if any(abs(startPoint-obj.get_center()) > (boxsize*1.5)) and all(np.abs(obj.get_center())+(boxsize*1.5) <= globalbox):
                notOverlap = True
                for x in range(len(objs)):
                    centerDiff = [abs(a - b) >= (boxsize*3) for a, b in zip(obj.get_center(), objs[x].get_center())]
                    if not(any(centerDiff)):
                        notOverlap = False
                if notOverlap:
                    objs.append(obj)
                    if not(indivPlots):
                        obj.show_curve()
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
            ##Show rays not interacting with any curves here. Only used when both reflect and refract in single object
            # t = np.linspace(0, 30, 500)
            # plt.plot(currInfo[0][0] + t*currInfo[1][0], currInfo[0][1] + t*currInfo[1][1], 'green')

            ##When reflect and refract are exclusive, put current info back to print final ray at end
            outputs.put(currInfo)
            break
        else:
            nextRays = currObj.output(distSmall, currInfo[0], currInfo[1], currInfo[2], currInfo[3], currInfo[4], boxedge, indivPlots)
            for nextRay in nextRays:
                outPoint = nextRay[0]
                outRay = currInfo[1]
                isBoxEdge = nextRay[4]
                if currObj.get_type() == "reflection":  #Assume always reflecting for now
                    outRay2 = nextRay[1]
                elif currObj.get_type() == "refraction":
                    outRay2 = nextRay[2]
                    
                if np.array_equal(nextRay[1], currInfo[1]) and np.array_equal(nextRay[2], currInfo[1]):
                    outputs.put([nextRay[0], nextRay[2], currInfo[3], currInfo[2], nextRay[3]])
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
    encryptedMsgs = [encrypt(str(outPoint[0]), "x Position", globalbox), encrypt(str(outPoint[1]), "y Position", globalbox),
        encrypt(str(inRay[0]), "x Direction", globalbox), encrypt(str(inRay[1]), "y Direction", globalbox), 
        encrypt(str(mag), "Ray Magnitude", globalbox)]

    test1 = np.array([decrypt(encryptedMsgs[0]), decrypt(encryptedMsgs[1])])
    test2 = np.array([decrypt(encryptedMsgs[2]), decrypt(encryptedMsgs[3])])
    test3 = decrypt(encryptedMsgs[4])

    reverse.put([test1, test2, n1, n2, intensity])
    if actionCount > 1:
        tests, testPoint, testRay, testRay2, testMag, testFlag, testCount = find_rays(reverse, objs, actionCount)
    else: 
        testPoint = test1
        testRay2 = test2
    t = np.linspace(0, test3, 500)
    print(np.array([testPoint[0] + test3*testRay2[0], testPoint[1] + test3*testRay2[1]]))
    plt.plot(testPoint[0] + t*testRay2[0], testPoint[1] + t*testRay2[1], 'orange')
else:
    raise Exception("No interactions occurred, please adjust input parameters")

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
    plt.xlabel('x')
    plt.ylabel('y')
    lines = [Line2D([0], [0], color=c, linewidth=3) for c in ['red', 'black', 'green', 'orange']]
    labels = ['Curve Objects', 'Rays', 'Output Ray Forward Direction', 'Final Ray Reverse Direction']
    plt.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=2)
    plt.show()