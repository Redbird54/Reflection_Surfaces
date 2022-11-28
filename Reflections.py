import numpy as np 
import pandas as pd
import math
import matplotlib.pyplot as plt
# pip3 install numdifftools
import numdifftools as nd

class Object:
    def __init__(self):
        pass

    def get_distance(self, initPoint, initDir):
        return -1

    def reflect(self, dist, initPoint, initDir, ray_only, isPlot):
        return initPoint, initDir


class Parabola(Object):
    def __init__(self,a,h,k,theta=0):
        self.a = a
        self.h = h
        self.k = k
        self.theta = theta
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)

        ##Equation of curve
        print("Equation of curve,  y =", a,"(x - ", h,")^2 +", k)
        x = np.linspace(-10, 10, 1000)
        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)
        plt.contour(hypx, hypy,((a*(hypx**2)*(cos**2)) + (a*(hypy**2)*(sin**2)) + (2*a*hypx*hypy*cos*sin) + (hypx*cos * (-2*a*h*cos - 2*a*k*sin)) + (hypy*sin * (-2*a*h*cos - 2*a*k*sin)) + (hypx*sin) - (hypy*cos) + (a*h*h*cos*cos + 2*a*h*k*cos*sin + a*k*k*sin*sin - h*sin + k*cos)), [0], colors='red')

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
            - (2*a * initDir[0] * cos*(h*cos + k*sin)) - (2*a * initDir[1] * sin*(h*cos + k*sin)) + (initDir[0] * sin) - (initDir[1] * cos))
        c_term = ((a * (((initPoint[0]*cos + initPoint[1]*sin)**2) + ((h*cos + k*sin)**2)) 
            - 2*initPoint[0]*cos*(h*cos + k*sin) - 2*initPoint[1]*sin*(h*cos + k*sin)) 
            + ((initPoint[0]-h)*sin) - ((initPoint[1]-k)*cos))
        if b_term*b_term - 4 * a_term * c_term < 0:
            print('No intercept')
            return -1
        distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        
        if distance_add <= 1e-14 and distance_min <= 1e-14:
            print('No intercept (neg)')
            return -1
        elif distance_add <= 1e-14:
            distance_true = distance_min
        elif distance_min <= 1e-14:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)

        return distance_true

    def func(self, x):
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        return ((self.a*(x[0]**2)*(cos**2)) + (self.a*(x[1]**2)*(sin**2)) 
            + (2*self.a*x[0]*x[1]*cos*sin) + (x[0]*cos * -2 * (self.a*self.h*cos + self.a*self.k*sin)) 
            + (x[1]*sin * -2*(self.a*self.h*cos + self.a*self.k*sin)) + (x[0]*sin) - (x[1]*cos) 
            + (self.a*(self.h**2)*(cos**2) + 2*self.a*self.h*self.k*cos*sin + self.a*(self.k**2)*(sin**2) 
            - self.h*sin + self.k*cos))

    def reflect(self, dist, initPoint, initDir, ray_only, isPlot):

        print('PARABOLA')
        
        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variable for plotting & incident ray
        t = np.linspace(0, 10, 500)

        print("Ray Direction: ", initDir[0], initDir[1])


        intercept = np.array([initPoint[0] + dist*initDir[0], initPoint[1] + dist*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')


        ##Plot where the ray actually goes (initial point until intercept)
        t2 = np.linspace(0, dist, 100)
        plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')


        ##Normal vector is gradient of function
        normDir = nd.Gradient(self.func)(intercept)

        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only):
            plt.plot(intercept[0] + t*normDir[0], intercept[1] + t*normDir[1], 'orange')


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = initDir - 2*(np.dot(initDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out


class Linear(Object):

    def __init__(self,*args, **kwargs):
        self.surfPoint = np.array([0,kwargs.get('b')])
        self.surfDir = np.array([1,kwargs.get('a')])
        if all(isinstance(arg, int) for arg in args) and len(args) == 2 and len(kwargs) == 0:
            self.surfPoint = np.array([0,args[1]])
            self.surfDir = np.array([1,args[0]])
        elif all(isinstance(arg, np.ndarray) for arg in args) and len(args) == 2 and len(kwargs) == 0:
            self.surfPoint = args[0]
            self.surfDir = args[1]
        elif kwargs.get('a') and kwargs.get('b') and len(args) == 0 and len(kwargs) == 2:
            self.surfPoint = np.array([0,kwargs.get('b')])
            self.surfDir = np.array([1,kwargs.get('a')])
        elif isinstance(kwargs.get('point'), np.ndarray) and isinstance(kwargs.get('dir'), np.ndarray) and len(args) == 0 and len(kwargs) == 2:
            self.surfPoint = kwargs.get('point')
            self.surfDir = kwargs.get('dir')
        else:
            raise ValueError("Unsupported linear format")

        ##Equation of curve
        print("Equation of curve:", self.surfPoint,"+ t*", self.surfDir)
        x = np.linspace(-10, 10, 1000)
        plt.plot(self.surfPoint[0] + x*self.surfDir[0], self.surfPoint[1] + x*self.surfDir[1],'red')

    def get_distance(self, initPoint, initDir):
        surfPoint = self.surfPoint
        surfDir = self.surfDir

        ##Find where ray and curve intersect
        a = self.surfDir[1]/self.surfDir[0]
        if (initDir[1] - (a * initDir[0])) == 0:
            print('No intercept')
            return -1
        distance = (surfPoint[1] - initPoint[1] + (a * (initPoint[0] - surfPoint[0]))) / (initDir[1] - (a * initDir[0]))

        if distance <= 1e-14:
            print('No intercept')
            return -1

        return distance

    def reflect(self, dist, initPoint, initDir, ray_only, isPlot):

        print('LINEAR')

        surfPoint = self.surfPoint
        surfDir = self.surfDir

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variable for plotting & incident ray
        t = np.linspace(0, 10, 500)

        print("Ray Direction: ", initDir[0], initDir[1])
            
        intercept = np.array([initPoint[0] + dist*initDir[0], initPoint[1] + dist*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        t2 = np.linspace(0, dist, 100)
        plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')


        ##Normal vector is gradient of function
        normDir = np.array([surfDir[1], -surfDir[0]])
        
        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only):
            plt.plot(intercept[0] + t*normDir[0], intercept[1] + t*normDir[1], 'orange')


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = initDir - 2*(np.dot(initDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out


class Hyperbola(Object):

    def __init__(self,a,b,h,k,theta=0):
        self.a = a
        self.b = b
        self.h = h
        self.k = k
        self.theta = theta
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)

        ##Equation of curve
        print("Equation of curve,  0 =", b*b,"x^2 +", -a*a,"y^2 +", -2*b*b*h,"x +", 2*a*a*k,"y +", b*b*h*h - a*a*k*k - a*a*b*b)
        x = np.linspace(-10, 10, 1000)
        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)
        plt.contour(hypx, hypy,((b**2)*((((hypx-h) * cos) + ((hypy-k) * sin))**2) - (a**2)*((((hypy-k) * cos) - ((hypx-h) * sin))**2) - (a**2)*(b**2)), [0], colors='red')

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
    
        if distance_add <= 1e-14 and distance_min <= 1e-14:
            print('No intercept')
            return -1
        elif distance_add <= 1e-14:
            distance_true = distance_min
        elif distance_min <= 1e-14:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)

        return distance_true

    def func(self, x):
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        return ((self.b**2) * ((((x[0] - self.h) * cos) + ((x[1] - self.k) * sin))**2) 
            - (self.a**2)*((((x[1] - self.k) * cos) - ((x[0] - self.h) * sin))**2) 
            - (self.a**2) * (self.b**2))

    def reflect(self, dist, initPoint, initDir, ray_only, isPlot):

        print('HYPERBOLA')

        a = self.a
        b = self.b
        h = self.h
        k = self.k
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variable for plotting & incident ray
        t = np.linspace(0, 10, 500)

        print("Ray Direction: ", initDir[0], initDir[1])


        intercept = np.array([initPoint[0] + dist*initDir[0], initPoint[1] + dist*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        t2 = np.linspace(0, dist, 100)
        plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')


        ##Normal vector is gradient of function
        normDir = nd.Gradient(self.func)(intercept)
    
        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only):
            plt.plot(intercept[0] + t*normDir[0], intercept[1] + t*normDir[1], 'orange')
        

        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1
        out = initDir - 2*(np.dot(initDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out


class Ellipse(Object):

    def __init__(self,a,b,h,k,theta=0):
        self.a = a
        self.b = b
        self.h = h
        self.k = k
        self.theta = theta
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)

        ##Equation of curve
        print("Equation of curve,  0 =", b*b,"x^2 +", a*a,"y^2 +", -2*b*b*h,"x +", -2*a*a*k,"y +", b*b*h*h + a*a*k*k - a*a*b*b)
        x = np.linspace(-10, 10, 100)
        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)
        plt.contour(hypx, hypy,((b**2) * (((hypx - h) * cos + (hypy - k) * sin)**2) + (a**2) * (((hypx - h) * sin - (hypy - k) * cos)**2) - ((a**2) * (b**2))), [0], colors='red')
    

    def get_distance(self, initPoint, initDir):
        a = self.a
        b = self.b
        h = self.h
        k = self.k
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)

        ##Find where ray and curve intersect
        #Can self be generalized?
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

        if distance_add <= 1e-14 and distance_min <= 1e-14:
            print('No intercept')
            return -1
        elif distance_add <= 1e-14:
            distance_true = distance_min
        elif distance_min <= 1e-14:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)
        
        return distance_true

    def func(self, x):
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        return ((self.b**2) * (((x[0] - self.h) * cos + (x[1] - self.k) * sin)**2) 
            + (self.a**2) * (((x[0] - self.h) * sin - (x[1] - self.k) * cos)**2) 
            - ((self.a**2) * (self.b**2)))


    def reflect(self, dist, initPoint, initDir, ray_only, isPlot):

        print('ELLIPSE')

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variable for plotting & incident ray
        t = np.linspace(0, 10, 500)

        print("Ray Direction: ", initDir[0], initDir[1])

        intercept = np.array([initPoint[0] + dist*initDir[0], initPoint[1] + dist*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        t2 = np.linspace(0, dist, 100)
        plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')


        ##Normal vector is gradient of function
        normDir = nd.Gradient(self.func)(intercept)

        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only):
            plt.plot(intercept[0] + t*normDir[0], intercept[1] + t*normDir[1], 'orange')	


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = initDir - 2*(np.dot(initDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')
        

        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out