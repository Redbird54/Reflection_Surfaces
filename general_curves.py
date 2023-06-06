import numpy as np 
import pandas as pd
import math
import matplotlib.pyplot as plt
import numdifftools as nd
import decimal
from matplotlib.lines import Line2D

class Object:
    def __init__(self):
        pass

    def get_type(self):
        return "none"

    def show_curve(self, boxshow):
        x = np.linspace(self.h - self.boxsize * 1.5, self.h + self.boxsize * 1.5, 1000)
        y = np.linspace(self.k - self.boxsize * 1.5, self.k + self.boxsize * 1.5, 1000)
        hypx, hypy = np.meshgrid(x, y)

        ##Equation of curve
        plt.contour(hypx, hypy, self.func([hypx, hypy]), [0], colors='red')

        if boxshow:
            self.show_box(self.h, self.k)

    def show_box(self, h, k):
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

    def get_box_distance(self, initPoint, dir):
        edge = self.boxsize / 2
        if not(self.h + edge > initPoint[0] and self.h - edge < initPoint[0] and
            self.k + edge > initPoint[1] and self.k - edge < initPoint[1]): #Intercept not in box
            return -1

        dists = np.array([(self.h + edge - initPoint[0]) / dir[0], (self.h - edge - initPoint[0]) / dir[0],
            (self.k + edge - initPoint[1]) / dir[1], (self.k - edge - initPoint[1]) / dir[1]])

        dists = dists[dists >= 0]
        return min(dists)

    def box_edge(self, initPoint, dir):
        nextDist = self.get_distance(initPoint, dir)
        boxDist =  self.get_box_distance(initPoint, dir)
        if boxDist != -1 and (nextDist == -1 or boxDist <= nextDist):
            outPoint, outDist = initPoint + boxDist*dir, boxDist
            isBoxEdge = 1
        else:
            outPoint, outDist = initPoint, 0
            isBoxEdge = 0
        t2 = np.linspace(0, outDist, 100)
        plt.plot(initPoint[0] + t2*dir[0], initPoint[1] + t2*dir[1], 'black')
        return outPoint, isBoxEdge

    def rotate(self, x1, x2, theta):
        sin = np.sin(theta)
        cos = np.cos(theta)
        return (x1*cos + x2*sin, -(x1*sin) + x2*cos)

    def output(self, dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot):
        return [[initPoint, initDir, initDir, intensity]]

    def output_procedure(self, dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow):
        mu = n1 / n2

        ##Variable for plotting & incident ray 
        t = np.linspace(0, dist, 100)

        intercept = np.array([initPoint[0] + dist*initDir[0], initPoint[1] + dist*initDir[1]])

        ##Plot where the ray actually goes (initial point until intercept)
        plt.plot(initPoint[0] + t*initDir[0], initPoint[1] + t*initDir[1], 'black')

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
            outRefr = (mu*initNorm) + (normNorm * np.sqrt(1 - ((mu**2) * (1 - ((np.dot(normNorm, initNorm))**2))))) - (mu * np.dot(normNorm, np.dot(normNorm, initNorm)))
        else:
            outRefr = initDir ##THIS LINE CAN CAUSE ISSUES WITH DECRYPTION, TOTAL INTERNAL REFLECTION##

        ##Find output intercept with box
        if self.notLens and boxedge:
            if self.objType == "reflection":
                #Used when we have transformation occurring on edge of a box
                outPoint, isBoxEdge = self.box_edge(intercept, outRefl)
            elif self.objType == "refraction":
                outPoint2, isBoxEdge = self.box_edge(intercept, outRefr)
            else:
                outPoint, isBoxEdge = self.box_edge(intercept, outRefl)
                outPoint2, isBoxEdge = self.box_edge(intercept, outRefr)
        else:
            outPoint = intercept
            outPoint2 = intercept
            isBoxEdge = 0
            
        ##Show plot
        if isPlot:
            plt.grid(color='lightgray', linestyle='--')
            plt.xlim(self.h-10, self.h+10)
            plt.ylim(self.k-10, self.k+10)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.xlabel('x')
            plt.ylabel('y')
            self.show_curve(boxshow)
            t3 = np.linspace(0, self.boxsize, 500)
            if self.notLens:
                plt.plot(intercept[0] + t3*normDir[0], intercept[1] + t3*normDir[1], 'orange')
                if self.objType == "reflection" or self.objType == "both":
                    plt.plot(outPoint[0] + t3*outRefl[0], outPoint[1] + t3*outRefl[1], 'green')
                elif self.objType == "refraction" or self.objType == "both":
                    plt.plot(outPoint2[0] + t3*outRefr[0], outPoint2[1] + t3*outRefr[1], 'green')
            lines = [Line2D([0], [0], color=c, linewidth=3) for c in ['red', 'black', 'green', 'orange']]
            labels = ['Curve Object', 'Input Ray', 'Output', 'Normal']
            if boxshow:
                lines.append(Line2D([0], [0], color='blue', linewidth=3))
                lines.append(Line2D([0], [0], color='brown', linewidth=3))
                labels.append('Inner Box')
                labels.append('Outer Box')
            plt.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=2)
            plt.show()

        if self.objType == "reflection":
            return [[outPoint, outRefl, initDir, intensity, isBoxEdge]]
        elif self.objType == "refraction":
            return [[outPoint2, initDir, outRefr, intensity, isBoxEdge]]
        else:
            Rs = np.linalg.norm((n1*np.cos(self.theta) - n2*sqrt(1-((n1 / n2)*np.sin(self.theta))**2)) / (n1*np.cos(self.theta) + n2*sqrt(1-((n1 / n2)*np.sin(self.theta))**2)))**2
            Rp = np.linalg.norm((n1*sqrt(1-((n1 / n2)*np.sin(self.theta))**2) - n2*np.cos(self.theta)) / (n1*sqrt(1-((n1 / n2)*np.sin(self.theta))**2) + n2*np.cos(self.theta)))**2
            R = 0.5*(Rs + Rp) #Unpolarized
            T = 1-R
            return [[outPoint, outRefl, initDir, R*intensity, isBoxEdge],[outPoint2, initDir, outRefr, T*intensity, isBoxEdge]]

class Polynomial3(Object):
    def __init__(self, a, b, c, d, e, f, g, h1, i, j, h, k, boxsize, outputType, theta=0):
        if outputType == "reflection" or outputType == "refraction": # or outputType == "both":
            self.objType = outputType
        else:
            raise Exception("Type must be \"reflection\", \"refraction\"") #, or \"both\"")
        if a == 0 and b == 0 and c == 0 and d == 0:
            raise Exception("a, b, c, and d cannot all have value of 0")
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f
        self.g = g
        self.h1 = h1
        self.i = i
        self.j = j
        self.h = h
        self.k = k       
        self.boxsize = boxsize
        self.theta = theta
        self.notLens = True

    def get_type(self):
        return self.objType

    def get_center(self):
        return self.h, self.k

    def cuberoot(self, x, y):
        # Paper for De Moivre's Theorem: https://web.pdx.edu/~caughman/Cindy%20501%20Final.pdf
        # Application of theorem: https://stackoverflow.com/questions/1361740/cubic-root-of-the-negative-number-on-python/1362288#1362288
        z = complex(x, y)
        mag = abs(z)
        arg = math.atan2(y,x)
        return mag**(1/3) * np.exp( 1j*arg/3 )

    def get_distance(self, initPoint, initDir):
        x, y = super().rotate((initPoint[0] - self.h), (initPoint[1] - self.k), self.theta)
        dx, dy = super().rotate(initDir[0], initDir[1], self.theta)

        A = (self.a * (dx**3)) + (self.b * (dy**3)) + (self.c * (dx**2) * dy) + (self.d * dx * (dy**2))
        B = ((self.a * 3 * x * (dx**2)) + (self.b * 3 * y * (dy**2)) + (self.c * (y * (dx**2) + 2 * x * dx * dy)) 
            + (self.d * (x * (dy**2) + 2 * y * dx * dy)) + (self.e * (dx**2)) + (self.f * (dy**2)) + (self.g * dx * dy))
        C = ((self.a * 3 * (x**2) * dx) + (self.b * 3 * (y**2) * dy) + (self.c * ((x**2) * dy + 2 * x * y * dx)) 
            + (self.d * ((y**2) * dx + 2 * x * y * dy)) + (self.e * 2 * x * dx) + (self.f * 2 * y * dy) 
            + (self.g * (x * dy + y * dx)) + (self.h1 * dx) + (self.i * dy))
        D = ((self.a * (x**3)) + (self.b * (y**3)) + (self.c * (x**2) * y) + (self.d * x * (y**2)) + (self.e * (x**2)) 
            + (self.f * (y**2)) + (self.g * x * y) + (self.h1 * x) + (self.i * y) + (self.j))

        dists = np.array([])

        Q = (3 * A * C - (B**2))/(9 * (A**2))
        R = (9 * A * B * C - 27 * (A**2) * D - 2 * (B**3))/(54 * (A**3))
        rootTest = (Q**3) + (R**2)
        if rootTest >= 0:
            S = np.cbrt(R + math.sqrt(rootTest))
            T = np.cbrt(R - math.sqrt(rootTest))
        else:
            S = self.cuberoot(R, math.sqrt(abs(rootTest)))
            T = self.cuberoot(R, -math.sqrt(abs(rootTest)))

        lambda1 = complex(S + T - (B / (3 * A)))
        lambda2 = complex(-(S + T)/2 + (S - T) * (1j * math.sqrt(3)) / 2 - (B / (3 * A)))
        lambda3 = complex(-(S + T)/2 - (S - T) * (1j * math.sqrt(3)) / 2 - (B / (3 * A)))

        if abs(lambda1.imag) <= 1e-11:
            lambda1 = lambda1.real
        if abs(lambda2.imag) <= 1e-11:
            lambda2 = lambda2.real
        if abs(lambda3.imag) <= 1e-11:
            lambda3 = lambda3.real

        if not(isinstance(lambda1, complex)):
            dists = np.append(dists, lambda1)

        if not(isinstance(lambda2, complex)):
            dists = np.append(dists, lambda2)

        if not(isinstance(lambda3, complex)):
            dists = np.append(dists, lambda3)

        dists = dists[dists > 1e-11]
        goodDists = np.array([])
        for myDist in dists:
            intercept = np.array([initPoint[0] + myDist*initDir[0], initPoint[1] + myDist*initDir[1]])
            if not(any(abs(intercept - np.array([self.h, self.k])) > (self.boxsize * 1.5))):
                goodDists = np.append(goodDists, myDist)
        if len(goodDists) > 0:
            return min(goodDists)
        else:
            return -1

    def func(self, input):
        x, y = super().rotate((input[0] - self.h), (input[1] - self.k), self.theta)
        return (self.a * (x**3) + self.b * (y**3) + self.c * (x**2) * y + self.d * x * (y**2) 
            + self.e * (x**2) + self.f * (y**2) + self.g * x * y + self.h1 * x + self.i * y + self.j)

    def output(self, dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow):
        print('3RD DEGREE POLYNOMIAL')
        return super().output_procedure(dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow)

class Polynomial2(Object):
    def __init__(self, a, b, c, d, e, f, h, k, boxsize, outputType, theta=0):
        if outputType == "reflection" or outputType == "refraction": # or outputType == "both":
            self.objType = outputType
        else:
            raise Exception("Type must be \"reflection\", \"refraction\"") #, or \"both\"")
        if a == 0 and b == 0 and c == 0:
            raise Exception("a, b, and c cannot all have a value of 0")
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f
        self.h = h
        self.k = k       
        self.boxsize = boxsize
        self.theta = theta
        self.notLens = True

    def get_type(self):
        return self.objType

    def get_center(self):
        return self.h, self.k

    def get_distance(self, initPoint, initDir):
        x, y = super().rotate((initPoint[0] - self.h), (initPoint[1] - self.k), self.theta)
        dx, dy = super().rotate(initDir[0], initDir[1], self.theta)

        A = (self.a * (dx**2)) + (self.b * (dy**2)) + (self.c * dx * dy)
        B = (self.a * 2 * x * dx) + (self.b * 2 * y * dy) + (self.c * (y * dx + x * dy)) + (self.d * dx) + (self.e * dy)       
        C = (self.a * (x**2)) + (self.b * (y**2)) + (self.c * x * y) + (self.d * x) + (self.e * y) + self.f

        if B**2 - 4 * A * C < 0:
            return -1
        if A == 0:
            return (-C)/B
        dists = np.array([(-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A), (-B - np.sqrt(B**2 - 4 * A * C)) / (2 * A)])

        dists = dists[dists > 1e-11]
        goodDists = np.array([])
        for myDist in dists:
            intercept = np.array([initPoint[0] + myDist*initDir[0], initPoint[1] + myDist*initDir[1]])
            if not(any(abs(intercept - np.array([self.h, self.k])) > (self.boxsize * 1.5))):
                goodDists = np.append(goodDists, myDist)
        if len(goodDists) > 0:
            return min(goodDists)
        else:
            return -1

    def func(self, input):
        x, y = super().rotate((input[0] - self.h), (input[1] - self.k), self.theta)
        return (self.a * (x**2) + self.b * (y**2) + self.c * x * y + self.d * x 
            + self.e * y + self.f)

    def output(self, dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow):
        print('2ND DEGREE POLYNOMIAL')
        return super().output_procedure(dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow)
