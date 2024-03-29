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
        t = np.linspace(0, outDist, 100)
        plt.plot(initPoint[0] + t*dir[0], initPoint[1] + t*dir[1], 'black')
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
            #Note: np.dot(normNorm, initNorm) is equal to cos(theta-i) in paper. Also normal vector defined to be in opposite direction of paper
            outRefr = mu*(initNorm - (normNorm * np.dot(normNorm, initNorm))) + (normNorm * np.sqrt(1 - ((mu**2) * (1 - ((np.dot(normNorm, initNorm))**2)))))
        else:
            outRefr = initDir ##THIS LINE PROBABLY CAUSES ISSUES WITH DECRYPTION##
            raise Exception ("TOTAL INTERNAL REFLECTION")

        ##Find output intercept with box
        if self.notLens and boxedge:
            if self.objType == "reflection":
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
            labels = ['Object', 'Input Ray', 'Output', 'Normal']
            if boxshow:
                lines.append(Line2D([0], [0], color='blue', linewidth=3))
                lines.append(Line2D([0], [0], color='brown', linewidth=3))
                labels.append('Local Bounding Box')
                labels.append('Object Limit Boundary')
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


class Parabola(Object):
    def __init__(self, a, h, k, boxsize, outputType, theta=0):
        if outputType == "reflection" or outputType == "refraction": # or outputType == "both":
            self.objType = outputType
        else:
            raise Exception("Type must be \"reflection\", \"refraction\"") #, or \"both\"")
        if a == 0:
            raise Exception("a cannot be 0")
        self.a = a
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

        ##Find where ray and curve intersect
        a_term = self.a * (dx**2)
        b_term = (2 * self.a * dx * x) - dy
        c_term = self.a * (x**2) - y
        if b_term**2 - 4 * a_term * c_term < 0:
            return -1
        if a_term == 0:
            return (-c_term)/b_term
        dists = np.array([(-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term), (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)])

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
        return (self.a * (x**2) - y)

    def output(self, dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow):
        return super().output_procedure(dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow)


class Linear(Object):
    def __init__(self, h, k, dx, dy, boxsize, outputType, notLens=True):
        if outputType == "reflection" or outputType == "refraction": # or outputType == "both":
            self.objType = outputType
        else:
            raise Exception("Type must be \"reflection\", \"refraction\"") #, or \"both\"")
        if dx == 0 and dy == 0:
            raise Exception("dx and dy cannot both be 0")
        self.h = h
        self.k = k
        self.surfPoint = np.array([h, k])
        self.surfDir = np.array([dx, dy])
        self.boxsize = boxsize
        self.notLens = notLens

    def get_type(self):
        return self.objType

    def get_center(self):
        return self.h, self.k

    def get_distance(self, initPoint, initDir):
        surfPoint = self.surfPoint
        surfDir = self.surfDir

        ##Find where ray and curve intersect
        if np.array_equal(surfDir / np.linalg.norm(surfDir), initDir / np.linalg.norm(initDir)):
            return -1
        distance = (surfDir[0] * (initPoint[1] - surfPoint[1]) + surfDir[1] * (surfPoint[0] - initPoint[0])) / (surfDir[1]*initDir[0] - surfDir[0]*initDir[1])

        if distance <= 1e-11:
            return -1

        intercept = np.array([initPoint[0] + distance*initDir[0], initPoint[1] + distance*initDir[1]])
        if any(abs(intercept - surfPoint) > (self.boxsize * 1.5)):
            return -1

        return distance

    def func(self, input):
        return (self.surfDir[0] * (self.surfPoint[1] - input[1]) + self.surfDir[1] * (input[0] - self.surfPoint[0]))

    def output(self, dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow):
        return super().output_procedure(dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow)


class Hyperbola(Object):
    def __init__(self, a, b, h, k, boxsize, outputType, theta=0):
        if outputType == "reflection" or outputType == "refraction": # or outputType == "both":
            self.objType = outputType
        else:
            raise Exception("Type must be \"reflection\", \"refraction\"") #, or \"both\"")
        if a == 0 and b == 0:
            raise Exception("a and b cannot both be 0")
        self.a = a
        self.b = b
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
        a = self.a
        b = self.b

        x, y = super().rotate((initPoint[0] - self.h), (initPoint[1] - self.k), self.theta)
        dx, dy = super().rotate(initDir[0], initDir[1], self.theta)

        ##Find where ray and curve intersect
        a_term = (b**2) * (dx**2) - (a**2) * (dy**2)
        b_term = 2 * ((b**2) * (x * dx) - (a**2) * (y * dy))
        c_term = (b**2) * (x**2) - (a**2) * (y**2) - (a**2)*(b**2)
        if b_term**2 - 4 * a_term * c_term < 0:
            return -1
        dists = np.array([(-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term), (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)])

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
        return ((self.b**2) * (x**2) 
            - (self.a**2) * (y**2) 
            - ((self.a**2) * (self.b**2)))

    def output(self, dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow):
        return super().output_procedure(dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow)


class Ellipse(Object):
    def __init__(self, a, b, h, k, boxsize, outputType, theta=0, notLens=True):
        if outputType == "reflection" or outputType == "refraction": # or outputType == "both":
            self.objType = outputType
        else:
            raise Exception("Type must be \"reflection\", \"refraction\"") #, or \"both\"")
        if a == 0 and b == 0:
            raise Exception("a and b cannot both be 0")
        self.a = a
        self.b = b
        self.h = h
        self.k = k
        self.boxsize = boxsize
        self.theta = theta
        self.notLens = notLens

    def get_type(self):
        return self.objType

    def get_center(self):
        return self.h, self.k

    def get_distance(self, initPoint, initDir):
        a = self.a
        b = self.b

        x, y = super().rotate((initPoint[0] - self.h), (initPoint[1] - self.k), self.theta)
        dx, dy = super().rotate(initDir[0], initDir[1], self.theta)

        ##Find where ray and curve intersect
        a_term = (b**2) * (dx**2) + (a**2) * (dy**2)
        b_term = 2 * ((b**2) * (x * dx) + (a**2) * (y * dy))
        c_term = (b**2) * (x**2) + (a**2) * (y**2) - (a**2)*(b**2)      
        if b_term**2 - 4 * a_term * c_term < 0:
            return -1
        dists = np.array([(-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term), (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)])

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
        return ((self.b**2) * (x**2) 
            + (self.a**2) * (y**2) 
            - ((self.a**2) * (self.b**2)))

    def output(self, dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow):
        return super().output_procedure(dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow)


#Convex/concave lens using two curved surfaces
class ParentLens:
    def __init__(self):
        pass
    
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
        t = np.linspace(0, outDist, 100)
        plt.plot(initPoint[0] + t*dir[0], initPoint[1] + t*dir[1], 'black')
        return outPoint, isBoxEdge

    def rotate(self, x1, x2, theta):
        sin = np.sin(theta)
        cos = np.cos(theta)
        return (x1*cos + x2*sin, -(x1*sin) + x2*cos)


class Lens(ParentLens): 
    def __init__(self, h, k, r1, r2, height, boxsize, theta=0):
        self.h = h
        self.k = k
        self.r1 = r1   # radius of first curvature
        self.r2 = r2
        if height == 0:
            raise Exception("Lens thickness cannot be zero")
        else:
            self.height = height # thickness of lens
        self.boxsize = boxsize
        if theta >= 0:
            self.theta = theta%math.pi
        else:
            while theta < 0:
                theta += math.pi
            self.theta = theta
        self.center1 = np.array([h + np.cos(self.theta) * (r1 - (height / 2)), k + np.sin(self.theta) * (r1 - (height / 2))])
        self.center2 = np.array([h + np.cos(self.theta) * ((height / 2) - r2), k + np.sin(self.theta) * ((height / 2) - r2)])
        self.lens1 = Ellipse(r1, r1, self.center1[0], self.center1[1], 10*self.boxsize, "refraction", notLens = False)
        self.lens2 = Ellipse(r2, r2, self.center2[0], self.center2[1], 10*self.boxsize, "refraction", notLens = False)
        
        if (height > 0):
            d = np.linalg.norm(self.center1 - self.center2)
            self.theta1 = math.acos((self.r1**2 + d**2 - self.r2**2)/(2 * self.r1 * d))
            self.theta2 = math.acos((self.r2**2 + d**2 - self.r1**2)/(2 * self.r2 * d))

            self.slope = np.array([-math.sin(theta), math.cos(theta)])
            self.point1 = np.array([self.center1[0]+self.r1*np.cos(self.theta + math.pi - self.theta1), self.center1[1]+self.r1*np.sin(self.theta + math.pi - self.theta1)])
            self.point2 = self.point1
        else:
            check1 = False
            check2 = False
            
            for test in [[0, -1], [1, -1], [0, 1], [1, 1]] :
                out1 = self.findBoxPoint(r1, self.center1, np.array([h, k]), test[0], test[1])
                if out1.any():
                    self.point1 = out1
                    check1 = True
                out2 = self.findBoxPoint(r2, self.center2, np.array([h, k]), test[0], test[1])
                if out2.any():
                    self.point2 = out2
                    check2 = True
                if check1 and check2:
                    break

            self.slope = np.array([-math.sin(theta), math.cos(theta)])

    def get_type(self):
        return "refraction"

    def findBoxPoint(self, r, center, middle, var, side):
        ##Var is 0  when testing for y, 1 when testing for x
        ##Side is -1 when testing left of bottom of box, 1 for top or right
        b = self.boxsize*1.5
        if np.sqrt(r**2 - (middle[var]+side*b-center[var])**2):
            test = center[(var+1)%2] + np.sqrt(r**2 - (middle[var]+side*b-center[var])**2)
            test2 = center[(var+1)%2] - np.sqrt(r**2 - (middle[var]+side*b-center[var])**2)
            if (middle[(var+1)%2]-b) <= test <= (middle[(var+1)%2]+b):
                if var == 0:
                    return np.array([middle[0]+side*b, test])
                elif var == 1:
                    return np.array([test, middle[1]+side*b])
            elif (middle[(var+1)%2]-b) <= test2 <= (middle[(var+1)%2]+b):
                if var == 0:
                    return np.array([middle[0]+side*b, test2])
                elif var == 1:
                    return np.array([test2, middle[1]+side*b])
        return np.array([])

    def get_center(self):
        return self.h, self.k

    def show_curve(self, boxshow):
        # super().show_box(self.h, self.k)
        if self.height > 0:
            t = np.linspace(self.theta + math.pi - self.theta1, self.theta + math.pi + self.theta1, 100)
            t2 = np.linspace(self.theta + self.theta2, self.theta - self.theta2, 100)
            
            ##Equation of curves
            plt.plot(self.center1[0]+self.r1*np.cos(t), self.center1[1]+self.r1*np.sin(t), color='red')
            plt.plot(self.center2[0]+self.r2*np.cos(t2), self.center2[1]+self.r2*np.sin(t2), color='red')
        else:
            x = np.linspace(self.h - self.boxsize * 1.5, self.h + self.boxsize * 1.5, 1000)
            y = np.linspace(self.k - self.boxsize * 1.5, self.k + self.boxsize * 1.5, 1000)
            hypx, hypy = np.meshgrid(x, y)

            plt.contour(hypx, hypy, (self.func(self.center1, [hypx, hypy], self.r1)), [0], colors='red')
            plt.contour(hypx, hypy, (self.func(self.center2, [hypx, hypy], self.r2)), [0], colors='red')

        if boxshow:
            super().show_box(self.h, self.k)
        
    def get_distance(self, initPoint, initDir):
        dist1 = self.lens1.get_distance(initPoint, initDir)
        if dist1 != -1:
            intercept1 = initPoint + dist1*initDir
        dist2 = self.lens2.get_distance(initPoint, initDir)
        if dist2 != -1:
            intercept2 = initPoint + dist2*initDir

        
        #Testing if intercepts are valid by testing if they are above or below line connecting where lenses intercept
        #If not finding 2nd intercept and testing it as well. Distance not -1 only if valid intercept found 
        if self.theta > 0:
            if dist1 != -1:
                if intercept1[1] > self.point1[1] + (intercept1[0]-self.point1[0]) * self.slope[1]/self.slope[0]:
                    nextDist = self.lens1.get_distance(intercept1, initDir)
                    intercept1 = intercept1 + nextDist*initDir
                    if (nextDist != -1 and intercept1[1] > self.point1[1] + (intercept1[0]-self.point1[0]) * self.slope[1]/self.slope[0]) or nextDist == -1:
                        dist1 = -1
                    else:
                        dist1 = dist1 + nextDist
            if dist2 != -1:
                if intercept2[1] < self.point2[1] + (intercept2[0]-self.point2[0]) * self.slope[1]/self.slope[0]:
                    nextDist = self.lens2.get_distance(intercept2, initDir)
                    intercept2 = intercept2 + nextDist*initDir
                    if (nextDist != -1 and intercept2[1] < self.point2[1] + (intercept2[0]-self.point2[0]) * self.slope[1]/self.slope[0]) or nextDist == -1:
                        dist2 = -1
                    else:
                        dist2 = dist2 + nextDist
        else: #Testing for left or right side when theta = 0
            if dist1 != -1 and intercept1[0] > self.point1[0]:
                nextDist = self.lens1.get_distance(intercept1, initDir)
                intercept1 = intercept1 + nextDist*initDir
                if (nextDist != -1 and intercept1[0] > self.point1[0]) or nextDist == -1:
                    dist1 = -1
                else:
                    dist1 = dist1 + nextDist
            if dist2 != -1 and intercept2[0] < self.point2[0]:
                nextDist = self.lens2.get_distance(intercept2, initDir)
                intercept2 = intercept2 + nextDist*initDir
                if (nextDist != -1 and intercept2[0] < self.point2[0]) or nextDist == -1:
                    dist2 = -1
                else:
                    dist2 = dist2 + nextDist

        if dist1 != -1 and dist2 == -1:
            return dist1
        elif dist1 == -1 and dist2 != -1:
            return dist2
        else:
            return min(dist1, dist2)

    def func(self, center, input, r):
        x, y = super().rotate((input[0] - self.center[0]), (input[1] - self.center[1]), self.theta)
        return (x**2 + y**2 - r**2)

    def output(self, dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow):
        intercept = initPoint + dist*initDir

        if -1e-11 <= (intercept[0] - self.center2[0])**2 + (intercept[1] - self.center2[1])**2 - (self.r2)**2 <= 1e-11:
            (nextPoint, nextRefl, nextRefr, intensity, extra) = self.lens2.output(dist, initPoint, initDir, n1, n2, intensity, boxedge, False, False)[0]
            distNew = self.lens1.get_distance(nextPoint, nextRefr)
            interceptNew = nextPoint + distNew*nextRefr
            if self.h - (self.boxsize*1.5) <= interceptNew[0] <= self.h + (self.boxsize*1.5) and self.k - (self.boxsize*1.5) <= interceptNew[1] <= self.k + (self.boxsize*1.5) and distNew != -1:
                (nextPoint, nextRefl, nextRefr, intensity, extra) = self.lens1.output(self.lens1.get_distance(nextPoint, nextRefr), nextPoint, nextRefr, n2, n1, intensity, boxedge, False, False)[0]
            else:
                nextPoint, nextRefl, nextRefr = intercept, initDir, initDir
        elif -1e-11 <= (intercept[0] - self.center1[0])**2 + (intercept[1] - self.center1[1])**2 - (self.r1)**2 <= 1e-11:
            (nextPoint, nextRefl, nextRefr, intensity, extra) = self.lens1.output(dist, initPoint, initDir, n1, n2, intensity, boxedge, False, False)[0]
            distNew = self.lens2.get_distance(nextPoint, nextRefr)
            interceptNew = nextPoint + distNew*nextRefr
            if self.h - (self.boxsize*1.5) <= interceptNew[0] <= self.h + (self.boxsize*1.5) and self.k - (self.boxsize*1.5) <= interceptNew[1] <= self.k + (self.boxsize*1.5) and distNew != -1:
                    (nextPoint, nextRefl, nextRefr, intensity, extra) = self.lens2.output(self.lens2.get_distance(nextPoint, nextRefr), nextPoint, nextRefr, n2, n1, intensity, boxedge, False, False)[0]
            else:
                nextPoint, nextRefl, nextRefr = intercept, initDir, initDir
        else:
            nextPoint, nextRefl, nextRefr = intercept, initDir, initDir

        if boxedge:
            nextPoint, isBoxEdge = super().box_edge(nextPoint, nextRefr)
        else:
            isBoxEdge = 0

        focalLength = 1 / ((n2 - 1) * ((1 / self.r1) + (1 / self.r2) - (((n2 - 1) * self.height) / (n2 * self.r1 * self.r2))))
        # print("Focal length: ", focalLength)

        if isPlot:
            t = np.linspace(0, 10, 500)
            plt.plot(nextPoint[0] + t*nextRefr[0], nextPoint[1] + t*nextRefr[1], 'green')
            t3 = np.linspace(-15, 6, 500)
            plt.plot(self.point1[0] + t3*self.slope[0], self.point1[1] + t3*self.slope[1], 'orange')
            plt.plot(self.point2[0] + t3*self.slope[0], self.point2[1] + t3*self.slope[1], 'orange')
            self.show_curve(boxshow)
            plt.grid(color='lightgray', linestyle='--')
            plt.xlim(self.h-20, self.h+20)
            plt.ylim(self.k-20, self.k+20)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.xlabel('x')
            plt.ylabel('y')
            lines = [Line2D([0], [0], color=c, linewidth=3) for c in ['red', 'black', 'green', 'orange']]
            labels = ['Lens', 'Input/Intermediate Rays', 'Output', 'Center Line']
            if boxshow:
                lines.append(Line2D([0], [0], color='blue', linewidth=3))
                lines.append(Line2D([0], [0], color='brown', linewidth=3))
                labels.append('Local Bounding Box')
                labels.append('Object Limit Boundary')
            plt.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=2)
            plt.show()

        return [[nextPoint, initDir, nextRefr, intensity, isBoxEdge]]


class Linear_Lens(ParentLens): 
    def __init__(self, h, k, height, boxsize, theta=0):
        self.h = h
        self.k = k
        if height <= 0:
            raise Exception("Lens thickness must be positive")
        else:
            self.height = height # thickness of lens
        self.boxsize = boxsize
        if theta >= 0:
            self.theta = theta%math.pi
        else:
            while theta < 0:
                theta += math.pi
            self.theta = theta
        self.center1 = np.array([h + np.cos(self.theta) * (height / 2), k + np.sin(self.theta) * (height / 2)])
        self.center2 = np.array([h - np.cos(self.theta) * (height / 2), k - np.sin(self.theta) * (height / 2)])
        self.slope = np.array([-math.sin(theta), math.cos(theta)])
        self.lens1 = Linear(self.center1[0], self.center1[1], self.slope[0], self.slope[1], self.boxsize, "refraction", notLens = False)
        self.lens2 = Linear(self.center2[0], self.center2[1], self.slope[0], self.slope[1], self.boxsize, "refraction", notLens = False)
    
    def get_type(self):
        return "refraction"

    def get_center(self):
        return self.h, self.k

    def show_curve(self, boxshow):
        x = np.linspace(self.h - self.boxsize * 1.5, self.h + self.boxsize * 1.5, 1000)
        y = np.linspace(self.k - self.boxsize * 1.5, self.k + self.boxsize * 1.5, 1000)
        hypx, hypy = np.meshgrid(x, y)

        # super().show_box(self.h, self.k)

        ##Equation of curve
        plt.contour(hypx, hypy, (self.func(self.center1, [hypx, hypy])), [0], colors='red')
        plt.contour(hypx, hypy, (self.func(self.center2, [hypx, hypy])), [0], colors='red')

        if boxshow:
            super().show_box(self.h, self.k)
        
    def get_distance(self, initPoint, initDir):
        dist1 = self.lens1.get_distance(initPoint, initDir)
        dist2 = self.lens2.get_distance(initPoint, initDir)

        if dist1 != -1 and dist2 == -1:
            return dist1
        elif dist1 == -1 and dist2 != -1:
            return dist2
        else:
            return min(dist1, dist2)

    def func(self, center, input):
        return (self.slope[0] * (center[1] - input[1]) + self.slope[1] * (input[0] - center[0]))

    def output(self, dist, initPoint, initDir, n1, n2, intensity, boxedge, isPlot, boxshow):
        intercept = initPoint + dist*initDir

        #Testing which order to refract through lenses, depending on location of intercept
        if (self.theta > 0 and intercept[1] < self.k + (intercept[0]-self.h) * self.slope[1]/self.slope[0]) or (self.theta == 0 and intercept[0] < self.h):
            (nextPoint, nextRefl, nextRefr, intensity, extra) = self.lens2.output(dist, initPoint, initDir, n1, n2, intensity, boxedge, False, False)[0]
            if self.lens1.get_distance(nextPoint, nextRefr) != -1:
                (nextPoint, nextRefl, nextRefr, intensity, extra) = self.lens1.output(self.lens1.get_distance(nextPoint, nextRefr), nextPoint, nextRefr, n2, n1, intensity, boxedge, False, False)[0]
            else:
                nextPoint, nextRefl, nextRefr = intercept, initDir, initDir
        elif (self.theta > 0 and intercept[1] > self.k + (intercept[0]-self.h) * self.slope[1]/self.slope[0]) or (self.theta == 0 and intercept[0] > self.h):
            (nextPoint, nextRefl, nextRefr, intensity, extra) = self.lens1.output(dist, initPoint, initDir, n1, n2, intensity, boxedge, False, False)[0]
            if self.lens2.get_distance(nextPoint, nextRefr) != -1:
                (nextPoint, nextRefl, nextRefr, intensity, extra) = self.lens2.output(self.lens2.get_distance(nextPoint, nextRefr), nextPoint, nextRefr, n2, n1, intensity, boxedge, False, False)[0]
            else:
                nextPoint, nextRefl, nextRefr = intercept, initDir, initDir
        else:
            nextPoint, nextRefl, nextRefr = intercept, initDir, initDir

        if boxedge:
            nextPoint, isBoxEdge = super().box_edge(nextPoint, nextRefr)
        else:
            isBoxEdge = 0

        if isPlot:
            t = np.linspace(0, 10, 500)
            plt.plot(nextPoint[0] + t*nextRefr[0], nextPoint[1] + t*nextRefr[1], 'green')
            self.show_curve(boxshow)
            plt.grid(color='lightgray', linestyle='--')
            plt.xlim(self.h-20, self.h+20)
            plt.ylim(self.k-20, self.k+20)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.xlabel('x')
            plt.ylabel('y')
            lines = [Line2D([0], [0], color=c, linewidth=3) for c in ['red', 'black', 'green']]
            labels = ['Lens', 'Input/Intermediate Rays', 'Output']
            if boxshow:
                lines.append(Line2D([0], [0], color='blue', linewidth=3))
                lines.append(Line2D([0], [0], color='brown', linewidth=3))
                labels.append('Local Bounding Box')
                labels.append('Object Limit Boundary')
            plt.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=2)
            plt.show()

        return [[nextPoint, initDir, nextRefr, intensity, isBoxEdge]]
