import numpy as np 
import pandas as pd
import math
import matplotlib.pyplot as plt
# pip3 install numdifftools
import numdifftools as nd

class Object:
    def __init__(self):
        pass

    def show_curve(self, isPlot):
        x = np.linspace(self.h-self.boxsize*1.5, self.h+self.boxsize*1.5, 1000)
        y = np.linspace(self.k-self.boxsize*1.5, self.k+self.boxsize*1.5, 1000)
        hypx, hypy = np.meshgrid(x, y)

        self.show_box(self.h,self.k,isPlot)

        ##Equation of curve
        plt.contour(hypx, hypy,(self.func([hypx,hypy])), [0], colors='red')

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-30, 30)
            plt.ylim(-30, 30)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()

    def show_box(self, h, k, isPlot):
        # PLOT CRYPTO RECTANGLE 
        plt.plot([h+(self.boxsize/2), h+(self.boxsize/2)], [k-(self.boxsize/2), k+(self.boxsize/2)], color = 'blue')
        plt.plot([h+(self.boxsize/2), h-(self.boxsize/2)], [k+(self.boxsize/2), k+(self.boxsize/2)], color = 'blue')
        plt.plot([h+(self.boxsize/2), h-(self.boxsize/2)], [k-(self.boxsize/2), k-(self.boxsize/2)], color = 'blue')
        plt.plot([h-(self.boxsize/2), h-(self.boxsize/2)], [k-(self.boxsize/2), k+(self.boxsize/2)], color = 'blue')

    def get_distance(self, initPoint, initDir):
        return -1

    def reflect(self, dist, initPoint, initDir, isPlot):
        return initPoint, initDir

    def reflect_procedure(self, dist, initPoint, initDir, isPlot):
        ##Variable for plotting & incident ray 
        t2 = np.linspace(0, dist, 100)
        self.show_curve(False)

        print("Ray Direction: ", initDir[0], initDir[1])

        intercept = np.array([initPoint[0] + dist*initDir[0], initPoint[1] + dist*initDir[1]])
        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')

        ##Normal vector is gradient of function
        normDir = nd.Gradient(self.func)(intercept)
        if np.dot(initDir, normDir) < 0:
            normDir = -normDir

        ##Direction of normal
        print("normDir: ", normDir)            	

        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = initDir - 2*(np.dot(initDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
            
        ##Print plot
        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-30, 30)
            plt.ylim(-30, 30)
            plt.gca().set_aspect('equal', adjustable='box')
            t = np.linspace(0, self.boxsize, 500)
            plt.plot(intercept[0] + t*normDir[0], intercept[1] + t*normDir[1], 'orange')
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')
            plt.show()

        return intercept, out


class Parabola(Object):
    def __init__(self, a, h, k, boxsize, theta=0):
        self.a = a
        self.h = h
        self.k = k
        self.boxsize = boxsize
        self.theta = theta

    def get_center(self):
        return self.h, self.k

    def get_distance(self, initPoint, initDir):
        a = self.a
        h = self.h
        k = self.k
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)

        ##Find where ray and curve intersect
        a_term = a * ((initDir[0]*cos + initDir[1]*sin)**2)
        b_term = ((2 * a * initPoint[0] * initDir[0]*(cos**2)) + (2 * a * initPoint[1] * initDir[1]*(sin**2)) 
            + (2 * a * initPoint[1] * initDir[0] * cos * sin) + (2 * a * initPoint[0] * initDir[1] * cos * sin) 
            - (2 * a * initDir[0] * cos*(h*cos + k*sin)) - (2 * a * initDir[1] * sin*(h*cos + k*sin)) + (initDir[0] * sin) - (initDir[1] * cos))
        c_term = ((a * (((initPoint[0]*cos + initPoint[1]*sin)**2) + ((h*cos + k*sin)**2)) 
            - 2*initPoint[0]*cos*(h*cos + k*sin) - 2*initPoint[1]*sin*(h*cos + k*sin)) 
            + ((initPoint[0]-h)*sin) - ((initPoint[1]-k)*cos))
        if b_term*b_term - 4 * a_term * c_term < 0:
            print('No intercept')
            return -1
        distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        
        if distance_add <= 1e-12 and distance_min <= 1e-12:
            print('No intercept')
            return -1
        elif distance_add <= 1e-12:
            distance_true = distance_min
        elif distance_min <= 1e-12:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])
        if any(abs(intercept-np.array([h,k])) > (self.boxsize*1.5)):
            return -1

        return distance_true

    def crosses_box_boundary(self, boxcenter):
        a = self.a
        h = self.h
        k = self.k
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        x_tests = [boxcenter[0]-(self.boxsize/2), boxcenter[0]+(self.boxsize/2)]
        y_tests = [boxcenter[1]-(self.boxsize/2), boxcenter[1]+(self.boxsize/2)]

        for x in x_tests:
            a_term = a*(sin**2)
            b_term = 2*a*sin*((x-h)*cos - k*sin) - cos
            c_term = a*(((x-h)*cos - k*sin)**2) + (x-h)*sin + k*cos
            if not(b_term*b_term - 4 * a_term * c_term < 0):
                y_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
                y_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
                if (y_add >= y_tests[0] and y_add <= y_tests[1]) or (y_min >= y_tests[0] and y_min <= y_tests[1]):
                    return True
            
        for y in y_tests:
            a_term = a*(cos**2)
            b_term = 2*a*cos*(-h*cos + (y-k)*sin) + sin
            c_term = a*((-h*cos + (y-k)*sin)**2) - h*sin - (y-k)*cos 
            if not(b_term*b_term - 4 * a_term * c_term < 0):
                x_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
                x_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)    
                if (x_add >= x_tests[0] and x_add <= x_tests[1]) or (x_min >= x_tests[0] and x_min <= x_tests[1]):
                    return True
        return False

    def func(self, x):
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        return (self.a*(((x[0]-self.h)*cos + (x[1]-self.k)*sin)**2) + (x[0]-self.h)*sin - (x[1]-self.k)*cos)

    def reflect(self, dist, initPoint, initDir, isPlot):
        print('PARABOLA')
        return super().reflect_procedure(dist, initPoint, initDir, isPlot)


class Linear(Object):
    def __init__(self, h, k, dx, dy, boxsize):
        self.h = h
        self.k = k
        self.surfPoint = np.array([h,k])
        self.surfDir = np.array([dx,dy])
        self.boxsize = boxsize

    def get_center(self):
        return self.h, self.k

    def get_distance(self, initPoint, initDir):
        surfPoint = self.surfPoint
        surfDir = self.surfDir

        ##Find where ray and curve intersect
        if np.array_equal(surfDir/np.linalg.norm(surfDir),initDir/np.linalg.norm(initDir)):
            print('No intercept')
            return -1
        distance = (surfDir[0] * (initPoint[1] - surfPoint[1]) + surfDir[1] * (surfPoint[0] - initPoint[0]))/(surfDir[1]*initDir[0] - surfDir[0]*initDir[1])

        if distance <= 1e-12:
            print('No intercept')
            return -1

        intercept = np.array([initPoint[0] + distance*initDir[0], initPoint[1] + distance*initDir[1]])
        if any(abs(intercept-surfPoint) > (self.boxsize*1.5)):
            return -1
        # if np.linalg.norm(intercept-surfPoint) > (self.boxsize*1.5):
        #     return -1

        return distance

    def crosses_box_boundary(self, boxcenter):
        x_tests = [boxcenter[0]-(self.boxsize/2), boxcenter[0]+(self.boxsize/2)]
        y_tests = [boxcenter[1]-(self.boxsize/2), boxcenter[1]+(self.boxsize/2)]

        for x in x_tests:
            if not(self.surfDir[0] == 0):
                y = self.surfPoint[1] + (x-self.surfPoint[0])*(self.surfDir[1]/self.surfDir[0])
            else:
                y = self.surfPoint[1]    
            if (y >= y_tests[0] and y <= y_tests[1]):
                return True
        for y in y_tests:
            if not(self.surfDir[1] == 0):
                x = self.surfPoint[0] + (y-self.surfPoint[1])*(self.surfDir[0]/self.surfDir[1])
            else:
                x = self.surfPoint[0]   
            if (x >= x_tests[0] and x <= x_tests[1]):
                return True
        return False

    def func(self, x):
        return (self.surfDir[0]*(self.surfPoint[1]-x[1]) + self.surfDir[1]*(x[0]-self.surfPoint[0]))

    def reflect(self, dist, initPoint, initDir, isPlot):
        print('LINEAR')
        return super().reflect_procedure(dist, initPoint, initDir, isPlot)


class Hyperbola(Object):
    def __init__(self, a, b, h, k, boxsize, theta=0):
        self.a = a
        self.b = b
        self.h = h
        self.k = k
        self.boxsize = boxsize
        self.theta = theta

    def get_center(self):
        return self.h, self.k

    def get_distance(self, initPoint, initDir):
        a = self.a
        b = self.b
        h = self.h
        k = self.k
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)

        ##Find where ray and curve intersect
        a_term = (((b**2)*(initDir[0]**2)*(cos**2)) - ((a**2)*(initDir[0]**2)*(sin**2)) 
            + (2*(b**2)*initDir[0]*initDir[1]*cos*sin) + (2*(a**2)*initDir[0]*initDir[1]*cos*sin) 
            + ((b**2)*(initDir[1]**2)*(sin**2)) - ((a**2)*(initDir[1]**2)*(cos**2)))
        b_term = 2 * (((b**2)*initPoint[0]*initDir[0]*(cos**2)) - ((b**2)*h*initDir[0]*(cos**2)) 
            + ((b**2)*initPoint[1]*initDir[1]*(sin**2)) - ((b**2)*k*initDir[1]*(sin**2)) 
            + ((b**2)*initPoint[0]*initDir[1]*cos*sin) + ((b**2)*initPoint[1]*initDir[0]*cos*sin) 
            - ((b**2)*k*initDir[0]*sin*cos) - ((b**2)*h*initDir[1]*sin*cos) 
            - ((a**2)*initPoint[1]*initDir[1]*(cos**2)) + ((a**2)*k*initDir[1]*(cos**2))
            - ((a**2)*initPoint[0]*initDir[0]*(sin**2))  + ((a**2)*h*initDir[0]*(sin**2)) 
            + ((a**2)*initPoint[0]*initDir[1]*cos*sin) + ((a**2)*initPoint[1]*initDir[0]*cos*sin) 
            - ((a**2)*k*initDir[0]*sin*cos) - ((a**2)*h*initDir[1]*sin*cos))
        c_term = ((b**2)*((cos**2)*((initPoint[0]**2)-2*initPoint[0]*h + (h**2)) 
            + (sin**2)*((initPoint[1]**2)-2*initPoint[1]*k + (k**2)) 
            + 2*sin*cos*(initPoint[0]*initPoint[1] - initPoint[0]*k - initPoint[1]*h + h*k)) 
            - (a**2)*((cos**2)*((initPoint[1]**2)-2*initPoint[1]*k + (k**2)) 
            + (sin**2)*((initPoint[0]**2)-2*initPoint[0]*h + (h**2)) 
            - 2*sin*cos*(initPoint[0]*initPoint[1] - initPoint[0]*k - initPoint[1]*h + h*k)) 
            - (a**2)*(b**2))
        if b_term*b_term - 4 * a_term * c_term < 0:
            print('No intercept')
            return -1
        distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
    
        if distance_add <= 1e-12 and distance_min <= 1e-12:
            print('No intercept')
            return -1
        elif distance_add <= 1e-12:
            distance_true = distance_min
        elif distance_min <= 1e-12:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])

        if any(abs(intercept-np.array([h,k])) > (self.boxsize*1.5)):
            return -1

        return distance_true

    def crosses_box_boundary(self, boxcenter):
        a = self.a
        b = self.b
        h = self.h
        k = self.k
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        x_tests = [boxcenter[0]-(self.boxsize/2), boxcenter[0]+(self.boxsize/2)]
        y_tests = [boxcenter[1]-(self.boxsize/2), boxcenter[1]+(self.boxsize/2)]

        for x in x_tests:
            a_term = ((b**2)*(sin**2)) - ((a**2)*(cos**2)) 
            b_term = 2 * (((b**2) * sin * ((x-h)*cos - k*sin)) + ((a**2) * cos * ((x-h)*sin + k*cos)))
            c_term = ((b**2) * (((x-h)*cos - k*sin)**2) - (a**2) * (((x-h)*sin + k*cos)**2) - (a**2)*(b**2))
            if not(b_term*b_term - 4 * a_term * c_term < 0):
                y_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
                y_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
                if (y_add >= y_tests[0] and y_add <= y_tests[1]) or (y_min >= y_tests[0] and y_min <= y_tests[1]):
                    return True
            
        for y in y_tests:
            a_term = ((b**2)*(cos**2)) + ((a**2)*(sin**2)) 
            b_term = 2 * (((b**2) * cos * (-h*cos + (y-k)*sin)) - ((a**2) * sin * (h*sin + (y-k)*cos)))
            c_term = ((b**2) * ((-h*cos + (y-k)*sin)**2) + (a**2) * ((h*sin + (y-k)*cos)**2) - (a**2)*(b**2)) 
            if not(b_term*b_term - 4 * a_term * c_term < 0):
                x_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
                x_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)    
                if (x_add >= x_tests[0] and x_add <= x_tests[1]) or (x_min >= x_tests[0] and x_min <= x_tests[1]):
                    return True
        return False
    
    def func(self, x):
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        return ((self.b**2) * ((((x[0] - self.h) * cos) + ((x[1] - self.k) * sin))**2) 
            - (self.a**2)*((((x[0] - self.h) * sin) - ((x[1] - self.k) * cos))**2) 
            - (self.a**2) * (self.b**2))

    def reflect(self, dist, initPoint, initDir, isPlot):
        print('HYPERBOLA')
        return super().reflect_procedure(dist, initPoint, initDir, isPlot)


class Ellipse(Object):
    def __init__(self, a, b, h, k, boxsize, theta=0):
        self.a = a
        self.b = b
        self.h = h
        self.k = k
        self.boxsize = boxsize
        self.theta = theta

    def get_center(self):
        return self.h, self.k

    def get_distance(self, initPoint, initDir):
        a = self.a
        b = self.b
        h = self.h
        k = self.k
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)

        ##Find where ray and curve intersect
        a_term = (((b**2)*(initDir[0]**2)*(cos**2)) + ((a**2)*(initDir[0]**2)*(sin**2)) 
            + (2*(b**2)*initDir[0]*initDir[1]*cos*sin) - (2*(a**2)*initDir[0]*initDir[1]*cos*sin) 
            + ((b**2)*(initDir[1]**2)*(sin**2)) + ((a**2)*(initDir[1]**2)*(cos**2)))
        b_term = 2 * ((b**2) * ((initPoint[0]*initDir[0]*(cos**2)) - (h*initDir[0]*(cos**2)) 
            + (initPoint[1]*initDir[1]*(sin**2)) - (k*initDir[1]*(sin**2)) 
            + (initPoint[0]*initDir[1]*cos*sin) + (initPoint[1]*initDir[0]*cos*sin) 
            - (k*initDir[0]*sin*cos) - (h*initDir[1]*sin*cos)) 
            + (a**2) * ((initPoint[1]*initDir[1]*(cos**2)) - (k*initDir[1]*(cos**2))
            + (initPoint[0]*initDir[0]*(sin**2)) - (h*initDir[0]*(sin**2)) 
            - (initPoint[0]*initDir[1]*cos*sin) - (initPoint[1]*initDir[0]*cos*sin) 
            + (k*initDir[0]*sin*cos) + (h*initDir[1]*sin*cos)))
        c_term = ((b**2)*((cos**2)*((initPoint[0]**2)-2*initPoint[0]*h + (h**2)) 
            + (sin**2)*((initPoint[1]**2)-2*initPoint[1]*k + (k**2)) 
            + 2*sin*cos*(initPoint[0]*initPoint[1] - initPoint[0]*k - initPoint[1]*h + h*k)) 
            + (a**2)*((cos**2)*((initPoint[1]**2)-2*initPoint[1]*k + (k**2)) 
            + (sin**2)*((initPoint[0]**2)-2*initPoint[0]*h + (h**2)) 
            - 2*sin*cos*(initPoint[0]*initPoint[1] - initPoint[0]*k - initPoint[1]*h + h*k)) 
            - (a**2)*(b**2))        
        if b_term*b_term - 4 * a_term * c_term < 0:
            print('No intercept')
            return -1
        distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)

        if distance_add <= 1e-12 and distance_min <= 1e-12:
            print('No intercept')
            return -1
        elif distance_add <= 1e-12:
            distance_true = distance_min
        elif distance_min <= 1e-12:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])
        if any(abs(intercept-np.array([h,k])) > (self.boxsize*1.5)):
            return -1
        
        return distance_true

    def crosses_box_boundary(self, boxcenter):
        a = self.a
        b = self.b
        h = self.h
        k = self.k
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        x_tests = [boxcenter[0]-(self.boxsize/2), boxcenter[0]+(self.boxsize/2)]
        y_tests = [boxcenter[1]-(self.boxsize/2), boxcenter[1]+(self.boxsize/2)]

        for x in x_tests:
            a_term = ((b**2)*(sin**2)) + ((a**2)*(cos**2)) 
            b_term = 2 * (((b**2) * sin * ((x-h)*cos - k*sin)) - ((a**2) * cos * ((x-h)*sin + k*cos)))
            c_term = ((b**2) * (((x-h)*cos - k*sin)**2) + (a**2) * (((x-h)*sin + k*cos)**2) - (a**2)*(b**2))
            if not(b_term*b_term - 4 * a_term * c_term < 0):
                y_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
                y_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
                if (y_add >= y_tests[0] and y_add <= y_tests[1]) or (y_min >= y_tests[0] and y_min <= y_tests[1]):
                    return True
            
        for y in y_tests:
            a_term = ((b**2)*(cos**2)) + ((a**2)*(sin**2)) 
            b_term = 2 * (((b**2) * cos * (-h*cos + (y-k)*sin)) - ((a**2) * sin * (h*sin + (y-k)*cos)))
            c_term = ((b**2) * ((-h*cos + (y-k)*sin)**2) + (a**2) * ((h*sin + (y-k)*cos)**2) - (a**2)*(b**2)) 
            if not(b_term*b_term - 4 * a_term * c_term < 0):
                x_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
                x_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)    
                if (x_add >= x_tests[0] and x_add <= x_tests[1]) or (x_min >= x_tests[0] and x_min <= x_tests[1]):
                    return True
        return False

    def func(self, x):
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        return ((self.b**2) * (((x[0] - self.h) * cos + (x[1] - self.k) * sin)**2) 
            + (self.a**2) * (((x[0] - self.h) * sin - (x[1] - self.k) * cos)**2) 
            - ((self.a**2) * (self.b**2)))

    def reflect(self, dist, initPoint, initDir, isPlot):
        print('ELLIPSE')
        return super().reflect_procedure(dist, initPoint, initDir, isPlot)
