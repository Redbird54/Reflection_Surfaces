import numpy as np 
from mpl_toolkits.mplot3d import axes3d
import pandas as pd
import math
import matplotlib.pyplot as plt
import numdifftools as nd
import queue

class Object:
    def __init__(self):
        pass

    def show_curve(self, ax):
        # Angles for polar coordinates:
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        # Cartesian coordinates from polar coordinates for ellipsoids:
        x0 = self.a * np.outer(np.cos(u), np.sin(v))
        y0 = self.b * np.outer(np.sin(u), np.sin(v))
        z0 = self.c * np.outer(np.ones_like(u), np.cos(v))

        xz, yz, zz = self.rotate(x0, y0, z0, self.thetaz)
        zy, xy, yy = self.rotate(zz, xz, yz, self.thetay)
        y, z, x = self.rotate(yy, zy, xy, self.thetax)

        x = x + self.h
        y = y + self.k
        z = z + self.l

        ax.plot_surface(x, y, z, color='red', alpha=0.75)

    def show_box(self, h, k, l): ##NOT YET UPDATED FOR 3D
        # PLOT CRYPTO RECTANGLE 
        plt.plot([h + (self.boxsize / 2), h + (self.boxsize / 2)], [k - (self.boxsize / 2), k + (self.boxsize / 2)], color = 'blue')
        plt.plot([h + (self.boxsize / 2), h - (self.boxsize / 2)], [k + (self.boxsize / 2), k + (self.boxsize / 2)], color = 'blue')
        plt.plot([h + (self.boxsize / 2), h - (self.boxsize / 2)], [k - (self.boxsize / 2), k - (self.boxsize / 2)], color = 'blue')
        plt.plot([h - (self.boxsize / 2), h - (self.boxsize / 2)], [k - (self.boxsize / 2), k + (self.boxsize / 2)], color = 'blue')

        plt.plot([h + (self.boxsize * 1.5), h + (self.boxsize * 1.5)], [k - (self.boxsize * 1.5), k + (self.boxsize * 1.5)], color = 'brown')
        plt.plot([h + (self.boxsize * 1.5), h - (self.boxsize * 1.5)], [k + (self.boxsize * 1.5), k + (self.boxsize * 1.5)], color = 'brown')
        plt.plot([h + (self.boxsize * 1.5), h - (self.boxsize * 1.5)], [k - (self.boxsize * 1.5), k - (self.boxsize * 1.5)], color = 'brown')
        plt.plot([h - (self.boxsize * 1.5), h - (self.boxsize * 1.5)], [k - (self.boxsize * 1.5), k + (self.boxsize * 1.5)], color = 'brown')

    def get_distance(self, initPoint, initDir):
        return -1

    def get_box_distance(self, initPoint, dir): ##NOT YET UPDATED FOR 3D
        edge = self.boxsize / 2
        if not(self.h + edge > initPoint[0] and self.h - edge < initPoint[0] and
            self.k + edge > initPoint[1] and self.k - edge < initPoint[1]): #Intercept not in box
            return -1
        
        if (self.h + edge - initPoint[0]) / dir[0] >= 0: #Crosses right side of box
            if initPoint[1] + (self.h + edge - initPoint[0])/dir[0] * dir[1] >= self.k - edge and initPoint[1] + (self.h + edge - initPoint[0])/dir[0] * dir[1] <= self.k + edge:
                outDist = (self.h + edge - initPoint[0]) / dir[0]
            elif (self.k + edge - initPoint[1])/dir[1] >= 0: #Crosses top of box before right side
                outDist = (self.k + edge - initPoint[1]) / dir[1]
            else: #Crosses bottom of box before right side
                outDist = (self.k - edge - initPoint[1]) / dir[1]
        elif (self.h - edge - initPoint[0]) / dir[0] >= 0: #Crosses left side of box
            if initPoint[1] + (self.h - edge - initPoint[0])/dir[0] * dir[1] >= self.k - edge and initPoint[1] + (self.h - edge - initPoint[0])/dir[0] * dir[1] <= self.k + edge:
                outDist = (self.h - edge - initPoint[0]) / dir[0]
            elif (self.k + edge - initPoint[1]) / dir[1] >= 0:
                outDist = (self.k + edge - initPoint[1]) / dir[1]
            else:
                outDist = (self.k - edge - initPoint[1]) / dir[1]
        else:
            outDist = -1
        return outDist

    def box_edge(self, initPoint, dir, ax):
        nextDist = self.get_distance(initPoint, dir)
        boxDist =  self.get_box_distance(initPoint, dir)
        if boxDist != -1 and (nextDist == -1 or boxDist <= nextDist):
            outPoint, outDist = initPoint + boxDist*dir, boxDist
        else:
            outPoint, outDist = initPoint, 0
        t = np.linspace(0, outDist, 100)
        ax.plot(initPoint[0] + t*dir[0], initPoint[1] + t*dir[1], initPoint[2] + t*dir[2], 'black')
        return outPoint

    def rotate(self, x1, x2, x3, theta):
        sin = np.sin(theta)
        cos = np.cos(theta)
        return (x1*cos + x2*sin), (x2*cos - x1*sin), x3

    def output(self, dist, initPoint, initDir, n1, n2, ax, isPlot):
        return [[initPoint, initDir, initDir]]

    def output_procedure(self, dist, initPoint, initDir, n1, n2, ax, isPlot):
        mu = n1 / n2

        ##Variable for plotting & incident ray 
        t = np.linspace(0, dist, 100)

        intercept = np.array([initPoint[0] + dist*initDir[0], initPoint[1] + dist*initDir[1], initPoint[2] + dist*initDir[2]])

        ##Plot where the ray actually goes (initial point until intercept)
        ax.plot(initPoint[0] + t*initDir[0], initPoint[1] + t*initDir[1], initPoint[2] + t*initDir[2], 'black')

        ##Normal vector is gradient of function
        normDir = nd.Gradient(self.func)(intercept)
        if np.dot(initDir, normDir) < 0:
            normDir = -normDir

        ##Direction of normal
        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        normNorm = normDir / np.linalg.norm(normDir)
        initNorm = initDir / np.linalg.norm(initDir)
        #Citation 1 
        outRefl = initDir - 2*(np.dot(initDir, normNorm) * normNorm)
        if np.sqrt(1 - ((mu**2) * (1 - ((np.dot(normNorm, initNorm))**2)))) >= 0:
            #Note: np.dot(normNorm, initNorm) is equal to cos(theta-i) in paper. Also normal vector defined to be in opposite direction of paper
            outRefr = mu*(initNorm - (normNorm * np.dot(normNorm, initNorm))) + (normNorm * np.sqrt(1 - ((mu**2) * (1 - ((np.dot(normNorm, initNorm))**2)))))
        else:
            outRefr = initDir ##THIS LINE PROBABLY CAUSES ISSUES WITH DECRYPTION##

        # ##Find output intercept with box
        # if self.notLens:
        #     if self.objType == "reflection":
        #         outPoint = self.box_edge(intercept, outRefl, ax)
        #     elif self.objType == "refraction":
        #         outPoint2 = self.box_edge(intercept, outRefr, ax)
        #     else:
        #         outPoint = self.box_edge(intercept, outRefl, ax)
        #         outPoint2 = self.box_edge(intercept, outRefr, ax)
        # else:
        #     outPoint = intercept
        #     outPoint2 = intercept

        outPoint = intercept
        outPoint2 = intercept
            
        ##Show plot
        if isPlot:
            self.show_curve(ax)
            ax.set_aspect('equal')
            t3 = np.linspace(0, self.boxsize, 500)
            if self.notLens:
                ax.plot(intercept[0] + t3*normDir[0], intercept[1] + t3*normDir[1], intercept[2] + t3*normDir[2], 'orange')
                if self.objType == "reflection" or self.objType == "both":
                    ax.plot(outPoint[0] + t3*outRefl[0], outPoint[1] + t3*outRefl[1], outPoint[2] + t3*outRefl[2], 'green')
                elif self.objType == "refraction" or self.objType == "both":
                    ax.plot(outPoint2[0] + t3*outRefr[0], outPoint2[1] + t3*outRefr[1], outPoint[2] + t3*outRefl[2], 'green')
            plt.show()

        if self.objType == "reflection":
            return [[outPoint, outRefl, initDir]]
        elif self.objType == "refraction":
            return [[outPoint2, initDir, outRefr]]
        else:
            return [[outPoint, outRefl, initDir],[outPoint2, initDir, outRefr]]


class Ellipsoid(Object):
    def __init__(self, a, b, c, h, k, l, boxsize, outputType, thetax=0, thetay=0, thetaz=0, notLens=True):
        if outputType == "reflection" or outputType == "refraction": # or outputType == "both":
            self.objType = outputType
        else:
            raise Exception("Type must be \"reflection\", \"refraction\"") #, or \"both\"")
        if a == 0 and b == 0 and c == 0:
            raise Exception("a, b, and c cannot all be 0")
        self.a = a
        self.b = b
        self.c = c
        self.h = h
        self.k = k
        self.l = l
        self.boxsize = boxsize
        self.thetax = thetax
        self.thetay = thetay
        self.thetaz = thetaz
        self.notLens = notLens

    def get_center(self):
        return self.h, self.k, self.l

    def get_distance(self, initPoint, initDir):
        a = self.a
        b = self.b
        c = self.c
        h = self.h
        k = self.k
        l = self.l
        xz, yz, zz = super().rotate((initPoint[0] - self.h), (initPoint[1] - self.k), (initPoint[2]- self.l), self.thetaz)
        zy, xy, yy = super().rotate(zz, xz, yz, self.thetay)
        y, z, x = super().rotate(yy, zy, xy, self.thetax)
        dxz, dyz, dzz = super().rotate(initDir[0], initDir[1], initDir[2], self.thetaz)
        dzy, dxy, dyy = super().rotate(dzz, dxz, dyz, self.thetay)
        dy, dz, dx = super().rotate(dyy, dzy, dxy, self.thetax)


        ##Find where ray and curve intersect
        a_term = (b**2)*(c**2)*(dx**2) + (a**2)*(c**2)*(dy**2) + (a**2)*(b**2)*(dz**2)
        b_term = (b**2)*(c**2)*(2*x*dx) + (a**2)*(c**2)*(2*y*dy) + (a**2)*(b**2)*(2*z*dz)
        c_term = (b**2)*(c**2)*(x**2) + (a**2)*(c**2)*(y**2) + (a**2)*(b**2)*(z**2) - (a**2)*(b**2)*(c**2)
        if b_term**2 - 4 * a_term * c_term < 0:
            return -1
        dists = np.array([(-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term), (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)])

        dists = dists[dists > 1e-12]
        goodDists = np.array([])
        for myDist in dists:
            intercept = np.array([initPoint[0] + myDist*initDir[0], initPoint[1] + myDist*initDir[1]])
            if not(any(abs(intercept - np.array([h, k])) > (self.boxsize * 1.5))):
                goodDists = np.append(goodDists, myDist)
        if len(goodDists) > 0:
            return min(goodDists)
        else:
            return -1

    def func(self, x):
        sinx = np.sin(self.thetax)
        cosx = np.cos(self.thetax)
        siny = np.sin(self.thetay)
        cosy = np.cos(self.thetay)
        sinz = np.sin(self.thetaz)
        cosz = np.cos(self.thetaz)

        return ((self.b**2)*(self.c**2) * ((((x[0]-self.h)*cosz + (x[1]-self.k)*sinz)*cosy - (x[2]-self.l)*siny)**2) 
        + (self.a**2)*(self.c**2) * ((((x[1]-self.k)*cosz - (x[0]-self.h)*sinz)*cosx + (((x[0]-self.h)*cosz + (x[1]-self.k)*sinz)*siny + (x[2]-self.l)*cosy)*sinx)**2) 
        + (self.a**2)*(self.b**2) * (((((x[0]-self.h)*cosz + (x[1]-self.k)*sinz)*siny + (x[2]-self.l)*cosy)*cosx - ((x[1]-self.k)*cosz - (x[0]-self.h)*sinz)*sinx)**2)
        - (self.a**2)*(self.b**2)*(self.c**2))

    def output(self, dist, initPoint, initDir, n1, n2, ax, isPlot):
        print('ELLIPSOID')
        return super().output_procedure(dist, initPoint, initDir, n1, n2, ax, isPlot)


def test3D():
    indivPlots = False
    n1 = 1.0003
    n2 = 1.52
    outputs = queue.Queue()
    objs = []
    nextPoint, nextRefl = np.array([0,-3,0.5]),np.array([0, 1, 0])

    # Plotting shape
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    objs.append(Ellipsoid(2,1,1,0,0,0,5,"reflection",math.pi/3,3,-math.pi/4))

    for obj in objs:
        if not(indivPlots):
            obj.show_curve(ax)

    outputs.put([nextPoint,nextRefl,n1,n2])
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

    nextRays = currObj.output(distSmall, currInfo[0], currInfo[1], currInfo[2], currInfo[3], ax, indivPlots)
    for nextRay in nextRays: 
        if not(np.array_equal(nextRay[1],currInfo[1])):
            outputs.put([nextRay[0],nextRay[1],currInfo[2],currInfo[3]])
        elif (np.array_equal(nextRay[2],currInfo[1])):
            ##Show rays not interacting with any curves here
            t = np.linspace(0, 5, 500)
            ax.plot(nextRay[0][0] + t*nextRay[1][0], nextRay[0][1] + t*nextRay[1][1], nextRay[0][2] + t*nextRay[1][2],'orange')
        if not(np.array_equal(nextRay[2],currInfo[1])):
            if not(any(np.isnan(nextRay[2]))):
                outputs.put([nextRay[0],nextRay[2],currInfo[3],currInfo[2]])
    if not(indivPlots):
        ax.set_aspect('equal')
        ##Show rays currently in the queue here
        t = np.linspace(0, 5, 500)
        for x in range(outputs.qsize()):
            output = outputs.get()
            plt.plot(output[0][0] + t*output[1][0], output[0][1] + t*output[1][1],output[0][2] + t*output[1][2],'green')
        plt.show()

test3D()
