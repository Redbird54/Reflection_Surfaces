from reflect_refract import *
from general_curves import *
import queue

##Hyperparameters for model
indivPlots = False
interactions = 8
boxsize = 5
intensity = 1
shape_count = 6

##Setting the medium 1 to air, medium 2 to glass
n1 = 1.0003
n2 = 1.52

outputs = queue.Queue()
initialObjs = []
objs = []

def get_valid_curves(initialObjs, objs, startPoint):
    for obj in initialObjs:
        if not(objs):
            if any(abs(startPoint-obj.get_center()) > (boxsize*1.5)):
                objs.append(obj)
                if not(indivPlots):
                    obj.show_curve()
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
                        obj.show_curve()
    return objs

def find_rays(outputs, objs):
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

        if currObj.get_type() == "none":
            ##Show rays not interacting with any curves here
            t = np.linspace(0, 30, 500)
            plt.plot(currInfo[0][0] + t*currInfo[1][0], currInfo[0][1] + t*currInfo[1][1], 'orange')
        else:
            nextRays = currObj.output(distSmall, currInfo[0], currInfo[1], currInfo[2], currInfo[3], currInfo[4], indivPlots)
            for nextRay in nextRays: 
                if np.array_equal(nextRay[1], currInfo[1]) and np.array_equal(nextRay[2], currInfo[1]):
                    outputs.put([nextRay[0], nextRay[2], currInfo[3], currInfo[2], nextRay[2]])
                else:
                    if not(np.array_equal(nextRay[1], currInfo[1])):
                        outputs.put([nextRay[0], nextRay[1], currInfo[2], currInfo[3], nextRay[2]])
                    if not(np.array_equal(nextRay[2], currInfo[1])):
                        if not(any(np.isnan(nextRay[2]))):
                            outputs.put([nextRay[0], nextRay[2], currInfo[3], currInfo[2], nextRay[2]])
    return outputs


def create_objs(h, k, numb):
    # Parabola: a, h, k, boxsize, outputType, theta=0
    # Linear: h, k, dx, dy, boxsize, outputType, notLens=True
    # Hyperbola: a, b, h, k, boxsize, outputType, theta=0
    # Ellipse: a, b, h, k, boxsize, outputType, theta=0, notLens=True
    # Lens: h, k, r1, r2, height, boxsize, theta=0
    # Linear Lens: h, k, height, boxsize, theta=0
    # Polynomial 3: a, b, c, d, e, f, g, h1, i, j, h, k, boxsize, outputType, theta=0
    # Pollynomial 2: a, b, c, d, e, f, h, k, boxsize, outputType, theta=0

    match numb:
        case 0:
            initialObjs.append(Polynomial3(5, 2, 3, -1, 7, 2, 0, -4, 1, -2, h, k, boxsize, "reflection"))
        case 1:
            initialObjs.append(Polynomial2(5, 3, 7, 5, -1, 0, h, k, boxsize, "reflection"))
        case 2:
            initialObjs.append(Ellipse(1, 2, h, k, boxsize, "reflection", math.pi/3))
        case 3:
            initialObjs.append(Hyperbola(3, 2, h, k, boxsize, "reflection", 5*math.pi/6))
        case 4:
            initialObjs.append(Linear(h, k, 4, 3, boxsize, "reflection"))
        case 5:
            initialObjs.append(Parabola(1, h, k, boxsize, "refraction", 5*math.pi/3))
        case 6:
            initialObjs.append(Lens(h, k, 7, 10, 5, boxsize, 9*math.pi/20))
        case 7:
            initialObjs.append(Linear_Lens(h, k, 4, boxsize, 18*math.pi/20))
        case _:
            pass

def generate_lattice():
    G = np.array([[15, 0], [0, 15]]).transpose()
    size = np.linspace(-10, 10, 21)
    numb = 0
    for i in size:
        for j in size:
            ans = np.matmul(G, np.array([i, j]))
            create_objs(ans[0], ans[1], numb)
            numb = (numb+1)%shape_count


startPoint, startRefl = np.array([0,0]),np.array([1,4])
generate_lattice()
objs = get_valid_curves(initialObjs, objs, startPoint)
outputs.put([startPoint,startRefl,n1,n2,intensity])
outputs = find_rays(outputs, objs)

if not(indivPlots):
    ##Show rays currently in the queue here
    t = np.linspace(0, 30, 500)
    for x in range(outputs.qsize()):
        output = outputs.get()
        plt.plot(output[0][0] + t*output[1][0], output[0][1] + t*output[1][1],'green')

    ##Print plot
    plt.grid(color='lightgray',linestyle='--')
    # plt.xlim(-157.5, 157.5)
    # plt.ylim(-157.5, 157,5)
    plt.xlim(-45, 45)
    plt.ylim(-45, 45)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()