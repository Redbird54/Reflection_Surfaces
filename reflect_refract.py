import numpy as np 
import pandas as pd
import math
import matplotlib.pyplot as plt
# pip3 install numdifftools
import numdifftools as nd

class Object:
    def __init__(self):
        pass

    def show_curve(self):
        x = np.linspace(self.h - self.boxsize * 1.5, self.h + self.boxsize * 1.5, 1000)
        y = np.linspace(self.k - self.boxsize * 1.5, self.k + self.boxsize * 1.5, 1000)
        hypx, hypy = np.meshgrid(x, y)

        self.show_box(self.h, self.k)

        ##Equation of curve
        plt.contour(hypx, hypy, (self.func([hypx, hypy])), [0], colors='red')

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

    def output(self, dist, initPoint, initDir, n1, n2, isPlot):
        return initPoint, initDir, initDir

    def output_procedure(self, dist, initPoint, initDir, n1, n2, isPlot):
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
        outRefr = (mu*initNorm) + (normNorm * np.sqrt(1 - ((mu**2) * (1 - ((np.dot(normNorm, initNorm))**2))))) - (mu * np.dot(normNorm, np.dot(normNorm, initNorm)))

        ##Find output intercept with box
        if self.notLens:
            if self.objType == "reflection":
                outPoint, outDist = self.findBoxIntercept(intercept, outRefl)
                t2 = np.linspace(0, outDist, 100)
                plt.plot(intercept[0] + t2*outRefl[0], intercept[1] + t2*outRefl[1], 'black')
            elif self.objType == "refraction":
                outPoint, outDist = self.findBoxIntercept(intercept, outRefr)
                t2 = np.linspace(0, outDist, 100)
                plt.plot(intercept[0] + t2*outRefr[0], intercept[1] + t2*outRefr[1], 'black')
            else:
                outPoint, outDist = self.findBoxIntercept(intercept, outRefl)
                t2 = np.linspace(0, outDist, 100)
                plt.plot(intercept[0] + t2*outRefl[0], intercept[1] + t2*outRefl[1], 'black')
                outPoint, outDist = self.findBoxIntercept(intercept, outRefr)
                t2 = np.linspace(0, outDist, 100)
                plt.plot(intercept[0] + t2*outRefr[0], intercept[1] + t2*outRefr[1], 'black')
        else:
            outPoint = intercept
            
        ##Show plot
        if isPlot:
            plt.grid(color='lightgray', linestyle='--')
            plt.xlim(self.h-10, self.h+10)
            plt.ylim(self.k-10, self.k+10)
            plt.gca().set_aspect('equal', adjustable='box')
            self.show_curve()
            t3 = np.linspace(0, self.boxsize, 500)
            if self.notLens:
                plt.plot(intercept[0] + t3*normDir[0], intercept[1] + t3*normDir[1], 'orange')
                if self.objType == "reflection" or self.objType == "both":
                    plt.plot(outPoint[0] + t3*outRefl[0], outPoint[1] + t3*outRefl[1], 'green')
                elif self.objType == "refraction" or self.objType == "both":
                    plt.plot(outPoint[0] + t3*outRefr[0], outPoint[1] + t3*outRefr[1], 'green')
            plt.show()

        if self.objType == "reflection":
            return outPoint, outRefl, initDir
        elif self.objType == "refraction":
            return outPoint, initDir, outRefr
        else:
            return outPoint, outRefl, outRefr

    def findBoxIntercept(self, initPoint, dir):
        edge = self.boxsize / 2
        if not(self.h + edge > initPoint[0] and self.h - edge < initPoint[0] and
            self.k + edge > initPoint[1] and self.k - edge < initPoint[1]): #Intercept not in box
            return initPoint, 0
        
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
            outDist = 0
        return initPoint + outDist*dir, outDist


class Parabola(Object):
    def __init__(self, a, h, k, boxsize, outputType, theta=0, notLens="True"):
        if outputType == "reflection" or outputType == "refraction": # or outputType == "both":
            self.objType = outputType
        else:
            raise Exception("Type must be \"reflection\", \"refraction\", or \"both\"")
        self.a = a
        self.h = h
        self.k = k
        self.boxsize = boxsize
        self.theta = theta
        self.notLens = notLens

    def get_center(self):
        return self.h, self.k

    def get_distance(self, initPoint, initDir):
        a = self.a
        h = self.h
        k = self.k
        sin = np.sin(self.theta)
        cos = np.cos(self.theta)

        ##Find where ray and curve intersect
        a_term = a * ((initDir[0]*cos + initDir[1]*sin)**2)
        b_term = ((2 * a * (initDir[0]*cos + initDir[1]*sin) * ((initPoint[0] - h)*cos + (initPoint[1] - k)*sin)) 
            + (initDir[0]*sin) - (initDir[1]*cos))
        c_term = (a * (((initPoint[0] - h)*cos + (initPoint[1] - k)*sin)**2)
            + ((initPoint[0] - h)*sin) - ((initPoint[1] - k)*cos))
        if b_term**2 - 4 * a_term * c_term < 0:
            return -1
        distance_add = (-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
        
        if distance_add <= 1e-12 and distance_min <= 1e-12:
            return -1
        elif distance_add <= 1e-12:
            distance_true = distance_min
        elif distance_min <= 1e-12:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)
            intercept1 = np.array([initPoint[0] + min(distance_add, distance_min)*initDir[0], initPoint[1] + min(distance_add, distance_min)*initDir[1]])
            if any(abs(intercept1 - np.array([h, k])) > (self.boxsize * 1.5)):
                distance_true = max(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])
        if any(abs(intercept - np.array([h, k])) > (self.boxsize * 1.5)):
            return -1

        return distance_true

    def crosses_box_boundary(self, boxcenter):
        a = self.a
        h = self.h
        k = self.k
        sin = np.sin(self.theta)
        cos = np.cos(self.theta)
        x_tests = [boxcenter[0] - (self.boxsize / 2), boxcenter[0] + (self.boxsize / 2)]
        y_tests = [boxcenter[1] - (self.boxsize / 2), boxcenter[1] + (self.boxsize / 2)]

        for x in x_tests:
            a_term = a * (sin**2)
            b_term = 2 * a * sin * ((x - h)*cos - k*sin) - cos
            c_term = a * (((x - h)*cos - k*sin)**2) + (x - h)*sin + k*cos
            if not(b_term**2 - 4 * a_term * c_term < 0):
                y_add = (-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
                y_min = (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
                if (y_add >= y_tests[0] and y_add <= y_tests[1]) or (y_min >= y_tests[0] and y_min <= y_tests[1]):
                    return True
            
        for y in y_tests:
            a_term = a * (cos**2)
            b_term = 2 * a * cos * (-h*cos + (y - k)*sin) + sin
            c_term = a * ((-h*cos + (y - k)*sin)**2) - h*sin - (y - k)*cos 
            if not(b_term**2 - 4 * a_term * c_term < 0):
                x_add = (-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
                x_min = (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)    
                if (x_add >= x_tests[0] and x_add <= x_tests[1]) or (x_min >= x_tests[0] and x_min <= x_tests[1]):
                    return True
        return False

    def func(self, x):
        sin = np.sin(self.theta)
        cos = np.cos(self.theta)
        return (self.a * (((x[0] - self.h)*cos + (x[1] - self.k)*sin)**2) + (x[0] - self.h)*sin - (x[1] - self.k)*cos)

    def output(self, dist, initPoint, initDir, n1, n2, isPlot):
        print('PARABOLA')
        return super().output_procedure(dist, initPoint, initDir, n1, n2, isPlot)


class Linear(Object):
    def __init__(self, h, k, dx, dy, boxsize, outputType, notLens="True"):
        if outputType == "reflection" or outputType == "refraction": # or outputType == "both":
            self.objType = outputType
        else:
            raise Exception("Type must be \"reflection\", \"refraction\", or \"both\"")
        self.h = h
        self.k = k
        self.surfPoint = np.array([h, k])
        self.surfDir = np.array([dx, dy])
        self.boxsize = boxsize
        self.notLens = notLens

    def get_center(self):
        return self.h, self.k

    def get_distance(self, initPoint, initDir):
        surfPoint = self.surfPoint
        surfDir = self.surfDir

        ##Find where ray and curve intersect
        if np.array_equal(surfDir / np.linalg.norm(surfDir), initDir / np.linalg.norm(initDir)):
            return -1
        distance = (surfDir[0] * (initPoint[1] - surfPoint[1]) + surfDir[1] * (surfPoint[0] - initPoint[0])) / (surfDir[1]*initDir[0] - surfDir[0]*initDir[1])

        if distance <= 1e-12:
            return -1

        intercept = np.array([initPoint[0] + distance*initDir[0], initPoint[1] + distance*initDir[1]])
        if any(abs(intercept - surfPoint) > (self.boxsize * 1.5)):
            return -1

        return distance

    def crosses_box_boundary(self, boxcenter):
        x_tests = [boxcenter[0] - (self.boxsize / 2), boxcenter[0] + (self.boxsize / 2)]
        y_tests = [boxcenter[1] - (self.boxsize / 2), boxcenter[1] + (self.boxsize / 2)]

        for x in x_tests:
            if not(self.surfDir[0] == 0):
                y = self.surfPoint[1] + (x - self.surfPoint[0])*(self.surfDir[1] / self.surfDir[0])
            else:
                y = self.surfPoint[1]    
            if (y >= y_tests[0] and y <= y_tests[1]):
                return True
        for y in y_tests:
            if not(self.surfDir[1] == 0):
                x = self.surfPoint[0] + (y - self.surfPoint[1])*(self.surfDir[0] / self.surfDir[1])
            else:
                x = self.surfPoint[0]   
            if (x >= x_tests[0] and x <= x_tests[1]):
                return True
        return False

    def func(self, x):
        return (self.surfDir[0] * (self.surfPoint[1] - x[1]) + self.surfDir[1] * (x[0] - self.surfPoint[0]))

    def output(self, dist, initPoint, initDir, n1, n2, isPlot):
        print('LINEAR')
        return super().output_procedure(dist, initPoint, initDir, n1, n2, isPlot)


class Hyperbola(Object):
    def __init__(self, a, b, h, k, boxsize, outputType, theta=0, notLens="True"):
        if outputType == "reflection" or outputType == "refraction": # or outputType == "both":
            self.objType = outputType
        else:
            raise Exception("Type must be \"reflection\", \"refraction\", or \"both\"")
        self.a = a
        self.b = b
        self.h = h
        self.k = k
        self.boxsize = boxsize
        self.theta = theta
        self.notLens = notLens

    def get_center(self):
        return self.h, self.k

    def get_distance(self, initPoint, initDir):
        a = self.a
        b = self.b
        h = self.h
        k = self.k
        sin = np.sin(self.theta)
        cos = np.cos(self.theta)

        ##Find where ray and curve intersect
        a_term = (b**2) * (((initDir[0]*cos) + (initDir[1]*sin))**2) - (a**2) * (((initDir[0]*sin) - (initDir[1]*cos))**2)
        b_term = 2 * ((b**2) * (((initPoint[0] - h)*initDir[0]*(cos**2)) 
            + ((initPoint[1] - k)*initDir[1]*(sin**2)) 
            + ((initPoint[0]*initDir[1] + initPoint[1]*initDir[0] - k*initDir[0] - h*initDir[1])*sin*cos)) 
            - (a**2) * (((initPoint[1] - k)*initDir[1]*(cos**2))
            + ((initPoint[0] - h)*initDir[0]*(sin**2)) 
            - ((initPoint[0]*initDir[1] + initPoint[1]*initDir[0] - k*initDir[0] - h*initDir[1])*sin*cos)))
        c_term = ((b**2) * (((initPoint[0] - h)*cos + (initPoint[1] - k)*sin)**2)
            - (a**2) * (((initPoint[0] - h)*sin - (initPoint[1] - k)*cos)**2)
            - (a**2)*(b**2))
        if b_term**2 - 4 * a_term * c_term < 0:
            return -1
        distance_add = (-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
    
        if distance_add <= 1e-12 and distance_min <= 1e-12:
            return -1
        elif distance_add <= 1e-12:
            distance_true = distance_min
        elif distance_min <= 1e-12:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)
            intercept1 = np.array([initPoint[0] + min(distance_add, distance_min)*initDir[0], initPoint[1] + min(distance_add, distance_min)*initDir[1]])
            if any(abs(intercept1 - np.array([h, k])) > (self.boxsize * 1.5)):
                distance_true = max(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])
        if any(abs(intercept - np.array([h, k])) > (self.boxsize * 1.5)):
            return -1

        return distance_true

    def crosses_box_boundary(self, boxcenter):
        a = self.a
        b = self.b
        h = self.h
        k = self.k
        sin = np.sin(self.theta)
        cos = np.cos(self.theta)
        x_tests = [boxcenter[0] - (self.boxsize / 2), boxcenter[0] + (self.boxsize / 2)]
        y_tests = [boxcenter[1] - (self.boxsize / 2), boxcenter[1] + (self.boxsize / 2)]

        for x in x_tests:
            a_term = ((b**2) * (sin**2)) - ((a**2) * (cos**2)) 
            b_term = 2 * (((b**2) * sin * ((x - h)*cos - k*sin)) + ((a**2) * cos * ((x - h)*sin + k*cos)))
            c_term = ((b**2) * (((x - h)*cos - k*sin)**2) - (a**2) * (((x - h)*sin + k*cos)**2) - (a**2)*(b**2))
            if not(b_term**2 - 4 * a_term * c_term < 0):
                y_add = (-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
                y_min = (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
                if (y_add >= y_tests[0] and y_add <= y_tests[1]) or (y_min >= y_tests[0] and y_min <= y_tests[1]):
                    return True
            
        for y in y_tests:
            a_term = ((b**2) * (cos**2)) + ((a**2) * (sin**2)) 
            b_term = 2 * (((b**2) * cos * (-h*cos + (y - k)*sin)) - ((a**2) * sin * (h*sin + (y - k)*cos)))
            c_term = ((b**2) * ((-h*cos + (y - k)*sin)**2) + (a**2) * ((h*sin + (y - k)*cos)**2) - (a**2)*(b**2)) 
            if not(b_term**2 - 4 * a_term * c_term < 0):
                x_add = (-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
                x_min = (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)    
                if (x_add >= x_tests[0] and x_add <= x_tests[1]) or (x_min >= x_tests[0] and x_min <= x_tests[1]):
                    return True
        return False
    
    def func(self, x):
        sin = np.sin(self.theta)
        cos = np.cos(self.theta)
        return ((self.b**2) * ((((x[0] - self.h)*cos) + ((x[1] - self.k)*sin))**2) 
            - (self.a**2) * ((((x[0] - self.h)*sin) - ((x[1] - self.k)*cos))**2) 
            - (self.a**2) * (self.b**2))

    def output(self, dist, initPoint, initDir, n1, n2, isPlot):
        print('HYPERBOLA')
        return super().output_procedure(dist, initPoint, initDir, n1, n2, isPlot)


class Ellipse(Object):
    def __init__(self, a, b, h, k, boxsize, outputType, theta=0, notLens="True"):
        if outputType == "reflection" or outputType == "refraction": # or outputType == "both":
            self.objType = outputType
        else:
            raise Exception("Type must be \"reflection\", \"refraction\", or \"both\"")
        self.a = a
        self.b = b
        self.h = h
        self.k = k
        self.boxsize = boxsize
        self.theta = theta
        self.notLens = notLens

    def get_center(self):
        return self.h, self.k

    def get_distance(self, initPoint, initDir):
        a = self.a
        b = self.b
        h = self.h
        k = self.k
        sin = np.sin(self.theta)
        cos = np.cos(self.theta)

        ##Find where ray and curve intersect
        a_term = (b**2) * (((initDir[0]*cos) + (initDir[1]*sin))**2) + (a**2) * (((initDir[0]*sin) - (initDir[1]*cos))**2)
        b_term = 2 * ((b**2) * (((initPoint[0] - h)*initDir[0]*(cos**2)) 
            + ((initPoint[1] - k)*initDir[1]*(sin**2)) 
            + ((initPoint[0]*initDir[1] + initPoint[1]*initDir[0] - k*initDir[0] - h*initDir[1])*sin*cos)) 
            + (a**2) * (((initPoint[1] - k)*initDir[1]*(cos**2))
            + ((initPoint[0] - h)*initDir[0]*(sin**2)) 
            - ((initPoint[0]*initDir[1] + initPoint[1]*initDir[0] - k*initDir[0] - h*initDir[1])*sin*cos)))
        c_term = ((b**2) * (((initPoint[0] - h)*cos + (initPoint[1] - k)*sin)**2)
            + (a**2) * (((initPoint[0] - h)*sin - (initPoint[1] - k)*cos)**2)
            - (a**2)*(b**2))      
        if b_term**2 - 4 * a_term * c_term < 0:
            return -1
        distance_add = (-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)

        if distance_add <= 1e-12 and distance_min <= 1e-12:
            return -1
        elif distance_add <= 1e-12:
            distance_true = distance_min
        elif distance_min <= 1e-12:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)
            intercept1 = np.array([initPoint[0] + min(distance_add, distance_min)*initDir[0], initPoint[1] + min(distance_add, distance_min)*initDir[1]])
            if any(abs(intercept1 - np.array([h, k])) > (self.boxsize * 1.5)):
                distance_true = max(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])
        if any(abs(intercept - np.array([h, k])) > (self.boxsize * 1.5)):
            return -1
        
        return distance_true

    def crosses_box_boundary(self, boxcenter):
        a = self.a
        b = self.b
        h = self.h
        k = self.k
        sin = np.sin(self.theta)
        cos = np.cos(self.theta)
        x_tests = [boxcenter[0] - (self.boxsize / 2), boxcenter[0] + (self.boxsize / 2)]
        y_tests = [boxcenter[1] - (self.boxsize / 2), boxcenter[1] + (self.boxsize / 2)]

        for x in x_tests:
            a_term = ((b**2) * (sin**2)) + ((a**2) * (cos**2)) 
            b_term = 2 * (((b**2) * sin * ((x - h)*cos - k*sin)) - ((a**2) * cos * ((x - h)*sin + k*cos)))
            c_term = ((b**2) * (((x - h)*cos - k*sin)**2) + (a**2) * (((x - h)*sin + k*cos)**2) - (a**2)*(b**2))
            if not(b_term**2 - 4 * a_term * c_term < 0):
                y_add = (-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
                y_min = (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
                if (y_add >= y_tests[0] and y_add <= y_tests[1]) or (y_min >= y_tests[0] and y_min <= y_tests[1]):
                    return True
            
        for y in y_tests:
            a_term = ((b**2) * (cos**2)) + ((a**2) * (sin**2)) 
            b_term = 2 * (((b**2) * cos * (-h*cos + (y - k)*sin)) - ((a**2) * sin * (h*sin + (y - k)*cos)))
            c_term = ((b**2) * ((-h*cos + (y - k)*sin)**2) + (a**2) * ((h*sin + (y - k)*cos)**2) - (a**2)*(b**2)) 
            if not(b_term**2 - 4 * a_term * c_term < 0):
                x_add = (-b_term + np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)
                x_min = (-b_term - np.sqrt(b_term**2 - 4 * a_term * c_term)) / (2 * a_term)    
                if (x_add >= x_tests[0] and x_add <= x_tests[1]) or (x_min >= x_tests[0] and x_min <= x_tests[1]):
                    return True
        return False

    def func(self, x):
        sin = np.sin(self.theta)
        cos = np.cos(self.theta)
        return ((self.b**2) * (((x[0] - self.h)*cos + (x[1] - self.k)*sin)**2) 
            + (self.a**2) * (((x[0] - self.h)*sin - (x[1] - self.k)*cos)**2) 
            - ((self.a**2) * (self.b**2)))

    def output(self, dist, initPoint, initDir, n1, n2, isPlot):
        print('ELLIPSE')
        return super().output_procedure(dist, initPoint, initDir, n1, n2, isPlot)



#Convex lens using two curved surfaces
class Lens: 
    def __init__(self, h, k, r1, r2, height, boxsize, theta=0):
        self.h = h
        self.k = k
        self.r1 = r1   # radius of first curvature
        self.r2 = r2
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
        self.center1 = np.array([h + np.cos(self.theta) * (r1 - (height / 2)), k + np.sin(self.theta) * (r1 - (height / 2))])
        self.center2 = np.array([h + np.cos(self.theta) * ((height / 2) - r2), k + np.sin(self.theta) * ((height / 2) - r2)])
        self.lens1 = Ellipse(r1, r1, self.center1[0], self.center1[1], 10*self.boxsize, "refraction", notLens = False)
        self.lens2 = Ellipse(r2, r2, self.center2[0], self.center2[1], 10*self.boxsize, "refraction", notLens = False)
        
        d = np.linalg.norm(self.center1 - self.center2)

        self.theta1 = math.acos((self.r1**2 + d**2 - self.r2**2)/(2 * self.r1 * d))
        self.theta2 = math.acos((self.r2**2 + d**2 - self.r1**2)/(2 * self.r2 * d))

        self.slope = np.array([self.r1*(np.cos(self.theta + math.pi - self.theta1) - np.cos(self.theta + math.pi + self.theta1)), self.r1*(np.sin(self.theta + math.pi - self.theta1) - np.sin(self.theta + math.pi + self.theta1))])
        self.crossPoint = np.array([self.center1[0]+self.r1*np.cos(self.theta + math.pi - self.theta1), self.center1[1]+self.r1*np.sin(self.theta + math.pi - self.theta1)])

    def get_center(self):
        return self.h, self.k

    def show_curve(self):
        self.show_box(self.h, self.k)

        t = np.linspace(self.theta + math.pi - self.theta1, self.theta + math.pi + self.theta1, 100)
        t2 = np.linspace(self.theta + self.theta2, self.theta - self.theta2, 100)
        
        ##Equation of curve
        plt.plot(self.center1[0]+self.r1*np.cos(t), self.center1[1]+self.r1*np.sin(t), color='red')
        plt.plot(self.center2[0]+self.r2*np.cos(t2), self.center2[1]+self.r2*np.sin(t2), color='red')

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
        dist1 = self.lens1.get_distance(initPoint, initDir)
        if dist1 != -1:
            intercept1 = initPoint + dist1*initDir
        dist2 = self.lens2.get_distance(initPoint, initDir)
        if dist2 != -1:
            intercept2 = initPoint + dist2*initDir

        #Testing if intercepts are valid by testing if they are above or below center of lens
        if self.theta > 0:
            if dist1 != -1 and intercept1[1] > self.k + (intercept1[0]-self.h) * self.slope[1]/self.slope[0]:
                nextDist = self.lens1.get_distance(intercept1, initDir)
                if nextDist == -1:
                    dist1 = -1
                else:
                    dist1 = dist1 + nextDist
            if dist2 != -1 and intercept2[1] < self.k + (intercept2[0]-self.h) * self.slope[1]/self.slope[0]:
                nextDist = self.lens2.get_distance(intercept2, initDir)
                if nextDist == -1:
                    dist2 = -1
                else:
                    dist2 = dist2 + nextDist
        else: #Testing for left or right side when theta = 0
            if dist1 != -1 and intercept1[0] > self.h:
                nextDist = self.lens1.get_distance(intercept1, initDir)
                if nextDist == -1:
                    dist1 = -1
                else:
                    dist1 = dist1 + nextDist
            if dist2 != -1 and intercept2[0] < self.h:
                nextDist = self.lens2.get_distance(intercept2, initDir)
                if nextDist == -1:
                    dist2 = -1
                else:
                    dist2 = dist2 + nextDist

        if dist1 != -1 and dist2 == -1:
            return dist1
        elif dist1 == -1 and dist2 != -1:
            return dist2
        else:
            return min(dist1, dist2)

    def output(self, dist, initPoint, initDir, n1, n2, isPlot):
        intercept = initPoint + dist*initDir

        if (self.theta > 0 and intercept[1] > self.k + (intercept[0]-self.h) * self.slope[1]/self.slope[0]) or (self.theta == 0 and intercept[0] > self.crossPoint[0]):
            nextPoint, nextRefl, nextRefr = self.lens2.output(dist, initPoint, initDir, n1, n2, False)
            if self.lens1.get_distance(nextPoint, nextRefr) != -1:
                nextPoint, nextRefl, nextRefr = self.lens1.output(self.lens1.get_distance(nextPoint, nextRefr), nextPoint, nextRefr, n2, n1, False)
            else:
                nextPoint, nextRefl, nextRefr = self.lens2.output(dist, initPoint, initDir, n2, n1, False)
        elif (self.theta > 0 and intercept[1] < self.k + (intercept[0]-self.h) * self.slope[1]/self.slope[0]) or (self.theta == 0 and intercept[0] < self.crossPoint[0]):
            nextPoint, nextRefl, nextRefr = self.lens1.output(dist, initPoint, initDir, n1, n2, False)
            if self.lens2.get_distance(nextPoint, nextRefr) != -1:
                nextPoint, nextRefl, nextRefr = self.lens2.output(self.lens2.get_distance(nextPoint, nextRefr), nextPoint, nextRefr, n2, n1, False)
            else:
                nextPoint, nextRefl, nextRefr = self.lens1.output(dist, initPoint, initDir, n2, n1, False)
        else:
            nextPoint, nextRefl, nextRefr = initPoint, initDir, initDir

        focalLength = 1 / ((n2 - 1) * ((1 / self.r1) + (1 / self.r2) - (((n2 - 1) * self.height) / (n2 * self.r1 * self.r2))))
        print("Focal length: ", focalLength)

        if isPlot:
            t = np.linspace(0, 10, 500)
            plt.plot(nextPoint[0] + t*nextRefr[0], nextPoint[1] + t*nextRefr[1], 'green')
            t3 = np.linspace(-3, 3, 500)
            plt.plot(self.crossPoint[0] + t3*self.slope[0], self.crossPoint[1] + t3*self.slope[1], 'orange')
            self.show_curve()
            plt.grid(color='lightgray', linestyle='--')
            plt.xlim(self.h-20, self.h+20)
            plt.ylim(self.k-20, self.k+20)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()

        return nextPoint, initDir, nextRefr


class Linear_Lens: 
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
        self.slope = np.array([h*math.cos(theta) - (k+1)*math.sin(theta) - h, h*math.sin(theta) + (k+1)*math.cos(theta) - k])
        self.lens1 = Linear(self.center1[0], self.center1[1], self.slope[0], self.slope[1], 10*self.boxsize, "refraction", notLens = False)
        self.lens2 = Linear(self.center2[0], self.center2[1], self.slope[0], self.slope[1], 10*self.boxsize, "refraction", notLens = False)
    
    def get_center(self):
        return self.h, self.k

    def show_curve(self):
        x = np.linspace(self.h - self.boxsize * 1.5, self.h + self.boxsize * 1.5, 1000)
        y = np.linspace(self.k - self.boxsize * 1.5, self.k + self.boxsize * 1.5, 1000)
        hypx, hypy = np.meshgrid(x, y)

        self.show_box(self.h, self.k)

        ##Equation of curve
        plt.contour(hypx, hypy, (self.func(self.center1, [hypx, hypy])), [0], colors='red')
        plt.contour(hypx, hypy, (self.func(self.center2, [hypx, hypy])), [0], colors='red')

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
        dist1 = self.lens1.get_distance(initPoint, initDir)
        dist2 = self.lens2.get_distance(initPoint, initDir)

        if dist1 != -1 and dist2 == -1:
            #CHECK FOR TOP/BOTTOM OF BOX
            return dist1
        elif dist1 == -1 and dist2 != -1:
            #CHECK FOR TOP/BOTTOM OF BOX
            return dist2
        else:
            return min(dist1, dist2)

    def func(self, center, x):
        return (self.slope[0] * (center[1] - x[1]) + self.slope[1] * (x[0] - center[0]))

    def output(self, dist, initPoint, initDir, n1, n2, isPlot):
        intercept = initPoint + dist*initDir

        #Testiing which order to refract through lenses, depending on location of intercept
        if (self.theta > 0 and intercept[1] < self.k + (intercept[0]-self.h) * self.slope[1]/self.slope[0]) or (self.theta == 0 and intercept[0] < self.h):
            nextPoint, nextRefl, nextRefr = self.lens2.output(dist, initPoint, initDir, n1, n2, False)
            if self.lens1.get_distance(nextPoint, nextRefr) != -1:
                nextPoint, nextRefl, nextRefr = self.lens1.output(self.lens1.get_distance(nextPoint, nextRefr), nextPoint, nextRefr, n2, n1, False)
            else:
                nextPoint, nextRefl, nextRefr = self.lens2.output(dist, initPoint, initDir, n2, n1, False) #Redo lens as though it is leaving glass, makes reversing process easier
        elif (self.theta > 0 and intercept[1] > self.k + (intercept[0]-self.h) * self.slope[1]/self.slope[0]) or (self.theta == 0 and intercept[0] > self.h):
            nextPoint, nextRefl, nextRefr = self.lens1.output(dist, initPoint, initDir, n1, n2, False)
            if self.lens2.get_distance(nextPoint, nextRefr) != -1:
                nextPoint, nextRefl, nextRefr = self.lens2.output(self.lens2.get_distance(nextPoint, nextRefr), nextPoint, nextRefr, n2, n1, False)
            else:
                nextPoint, nextRefl, nextRefr = self.lens1.output(dist, initPoint, initDir, n2, n1, False)
        else:
            nextPoint, nextRefl, nextRefr = initPoint, initDir, initDir

        if isPlot:
            t = np.linspace(0, 10, 500)
            plt.plot(nextPoint[0] + t*nextRefr[0], nextPoint[1] + t*nextRefr[1], 'green')
            self.show_curve()
            plt.grid(color='lightgray', linestyle='--')
            plt.xlim(self.h-20, self.h+20)
            plt.ylim(self.k-20, self.k+20)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()

        return nextPoint, initDir, nextRefr