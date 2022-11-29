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

    def reflect(self, dist, initPoint, initDir, isPlot):
        return initPoint, initDir


class Parabola(Object):
    def __init__(self, a, h, k, theta=0):
        self.a = a
        self.h = h
        self.k = k
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

    def crosses_box_boundary(self, boxcenter):
        a = self.a
        h = self.h
        k = self.k
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        x_tests = [boxcenter[0]-2.5, boxcenter[0]+2.5]
        y_tests = [boxcenter[1]-2.5, boxcenter[1]+2.5]

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

        a = self.a
        h = self.h
        k = self.k
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        
        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')
            t = np.linspace(0, 10, 500)

        ##Variables for plotting
        t2 = np.linspace(0, dist, 100)
        x = np.linspace(-10, 10, 1000)
        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)

        ##Equation of curve
        print("Equation of curve,  y =", a,"(x - ", h,")^2 +", k)
        plt.contour(hypx, hypy,(a*(((hypx-h)*cos + (hypy-k)*sin)**2) + (hypx-h)*sin - (hypy-k)*cos), [0], colors='red')


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
        if isPlot:
            plt.plot(intercept[0] + t*normDir[0], intercept[1] + t*normDir[1], 'orange')


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = initDir - 2*(np.dot(initDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if isPlot:
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')

        # PLOT CRYPTO RECTANGLE 
        point1 = [h+2.5, k-2.5] # br
        point2 = [h+2.5, k+2.5] # tr
        point3 = [h-2.5, k-2.5] # bl 
        point4 = [h-2.5, k+2.5] #tl
       
        x_values = [point1[0], point2[0]] #gather x-values.
        y_values = [point1[1], point2[1]] #gather y-values.
        x2_values = [point2[0], point4[0]] #gather x-values.
        y2_values = [point2[1], point4[1]] #gather y-values.
        x3_values = [point1[0], point3[0]] #gather x-values.
        y3_values = [point1[1], point3[1]] #gather y-values.
        x4_values = [point3[0], point4[0]] #gather x-values.
        y4_values = [point3[1], point4[1]] #gather y-values.
        plt.plot(x_values, y_values, color = 'blue')
        plt.plot(x2_values, y2_values, color = 'blue')
        plt.plot(x3_values, y3_values, color = 'blue')
        plt.plot(x4_values, y4_values, color = 'blue')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out


class Linear(Object):

    def __init__(self, h, k, dx, dy):
        self.h = h
        self.k = k
        self.surfPoint = np.array([h,k])
        self.surfDir = np.array([dx,dy])

    def get_center(self):
        return self.h, self.k

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

    def crosses_box_boundary(self, boxcenter):
        x_tests = [boxcenter[0]-2.5, boxcenter[0]+2.5]
        y_tests = [boxcenter[1]-2.5, boxcenter[1]+2.5]

        for x in x_tests:
            if not(self.surfDir[0] == 0):
                y = self.surfPoint[1] + (x-self.surfPoint[0])*(self.surfDir[1]/self.surfDir[0])       
                if (y >= y_tests[0] and y <= y_tests[1]):
                    return True
        for y in y_tests:
            if not(self.surfDir[1] == 0):
                x = self.surfPoint[0] + (y-self.surfPoint[1])*(self.surfDir[0]/self.surfDir[1])      
                if (x >= x_tests[0] and x <= x_tests[1]):
                    return True
        return False


    def reflect(self, dist, initPoint, initDir, isPlot):

        print('LINEAR')

        h = self.h
        k = self.k

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')
            t = np.linspace(0, 10, 500)

        ##Variables for plotting
        t2 = np.linspace(0, dist, 100)
        x = np.linspace(-10, 10, 1000)

        ##Equation of curve
        print("Equation of curve:", self.surfPoint,"+ t*", self.surfDir)
        plt.plot(self.surfPoint[0] + x*self.surfDir[0], self.surfPoint[1] + x*self.surfDir[1],'red')


        print("Ray Direction: ", initDir[0], initDir[1])
            
        intercept = np.array([initPoint[0] + dist*initDir[0], initPoint[1] + dist*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')


        ##Normal vector is gradient of function
        normDir = np.array([self.surfDir[1], -self.surfDir[0]])
        if np.dot(initDir, normDir) < 0:
            normDir = -normDir
        
        ##Direction of normal
        print("normDir: ", normDir)
        if isPlot:
            plt.plot(intercept[0] + t*normDir[0], intercept[1] + t*normDir[1], 'orange')


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = initDir - 2*(np.dot(initDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if isPlot:
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')

        # PLOT CRYPTO RECTANGLE 
        point1 = [h+2.5, k-2.5] # br
        point2 = [h+2.5, k+2.5] # tr
        point3 = [h-2.5, k-2.5] # bl 
        point4 = [h-2.5, k+2.5] #tl
       
        x_values = [point1[0], point2[0]] #gather x-values.
        y_values = [point1[1], point2[1]] #gather y-values.
        x2_values = [point2[0], point4[0]] #gather x-values.
        y2_values = [point2[1], point4[1]] #gather y-values.
        x3_values = [point1[0], point3[0]] #gather x-values.
        y3_values = [point1[1], point3[1]] #gather y-values.
        x4_values = [point3[0], point4[0]] #gather x-values.
        y4_values = [point3[1], point4[1]] #gather y-values.
        plt.plot(x_values, y_values, color = 'blue')
        plt.plot(x2_values, y2_values, color = 'blue')
        plt.plot(x3_values, y3_values, color = 'blue')
        plt.plot(x4_values, y4_values, color = 'blue')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out


class Hyperbola(Object):

    def __init__(self, a, b, h, k, theta=0):
        self.a = a
        self.b = b
        self.h = h
        self.k = k
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

    def crosses_box_boundary(self, boxcenter):
        a = self.a
        b = self.b
        h = self.h
        k = self.k
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        x_tests = [boxcenter[0]-2.5, boxcenter[0]+2.5]
        y_tests = [boxcenter[1]-2.5, boxcenter[1]+2.5]

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
            t = np.linspace(0, 10, 500)

        ##Variables for plotting  
        t2 = np.linspace(0, dist, 100)
        x = np.linspace(-10, 10, 1000)
        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)

        ##Equation of curve
        print("Equation of curve,  0 =", b*b,"x^2 +", -a*a,"y^2 +", -2*b*b*h,"x +", 2*a*a*k,"y +", b*b*h*h - a*a*k*k - a*a*b*b)
        plt.contour(hypx, hypy,((b**2)*((((hypx-h) * cos) + ((hypy-k) * sin))**2) - (a**2)*(((((hypx-h) * sin) - (hypy-k) * cos))**2) - (a**2)*(b**2)), [0], colors='red')


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
        if isPlot:
            plt.plot(intercept[0] + t*normDir[0], intercept[1] + t*normDir[1], 'orange')
        

        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1
        out = initDir - 2*(np.dot(initDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if isPlot:
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')

        # PLOT CRYPTO RECTANGLE 
        point1 = [h+2.5, k-2.5] # br
        point2 = [h+2.5, k+2.5] # tr
        point3 = [h-2.5, k-2.5] # bl 
        point4 = [h-2.5, k+2.5] #tl
       
        x_values = [point1[0], point2[0]] #gather x-values.
        y_values = [point1[1], point2[1]] #gather y-values.
        x2_values = [point2[0], point4[0]] #gather x-values.
        y2_values = [point2[1], point4[1]] #gather y-values.
        x3_values = [point1[0], point3[0]] #gather x-values.
        y3_values = [point1[1], point3[1]] #gather y-values.
        x4_values = [point3[0], point4[0]] #gather x-values.
        y4_values = [point3[1], point4[1]] #gather y-values.
        plt.plot(x_values, y_values, color = 'blue')
        plt.plot(x2_values, y2_values, color = 'blue')
        plt.plot(x3_values, y3_values, color = 'blue')
        plt.plot(x4_values, y4_values, color = 'blue')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out


class Ellipse(Object):

    def __init__(self, a, b, h, k, theta=0):
        self.a = a
        self.b = b
        self.h = h
        self.k = k
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

    def crosses_box_boundary(self, boxcenter):
        a = self.a
        b = self.b
        h = self.h
        k = self.k
        sin = math.sin(self.theta)
        cos = math.cos(self.theta)
        x_tests = [boxcenter[0]-2.5, boxcenter[0]+2.5]
        y_tests = [boxcenter[1]-2.5, boxcenter[1]+2.5]

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
            t = np.linspace(0, 10, 500)

        ##Variable for plotting & incident ray 
        t2 = np.linspace(0, dist, 100)
        x = np.linspace(-10, 10, 100)
        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)

        ##Equation of curve
        print("Equation of curve,  0 =", b*b,"x^2 +", a*a,"y^2 +", -2*b*b*h,"x +", -2*a*a*k,"y +", b*b*h*h + a*a*k*k - a*a*b*b)
        plt.contour(hypx, hypy,((b**2) * (((hypx - h) * cos + (hypy - k) * sin)**2) + (a**2) * (((hypx - h) * sin - (hypy - k) * cos)**2) - ((a**2) * (b**2))), [0], colors='red')


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
        if isPlot:
            plt.plot(intercept[0] + t*normDir[0], intercept[1] + t*normDir[1], 'orange')	


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = initDir - 2*(np.dot(initDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if isPlot:
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')
        


        # PLOT CRYPTO RECTANGLE 
        point1 = [h+2.5, k-2.5] # br
        point2 = [h+2.5, k+2.5] # tr
        point3 = [h-2.5, k-2.5] # bl 
        point4 = [h-2.5, k+2.5] #tl
       
        x_values = [point1[0], point2[0]] #gather x-values.
        y_values = [point1[1], point2[1]] #gather y-values.
        x2_values = [point2[0], point4[0]] #gather x-values.
        y2_values = [point2[1], point4[1]] #gather y-values.
        x3_values = [point1[0], point3[0]] #gather x-values.
        y3_values = [point1[1], point3[1]] #gather y-values.
        x4_values = [point3[0], point4[0]] #gather x-values.
        y4_values = [point3[1], point4[1]] #gather y-values.
        plt.plot(x_values, y_values, color = 'blue')
        plt.plot(x2_values, y2_values, color = 'blue')
        plt.plot(x3_values, y3_values, color = 'blue')
        plt.plot(x4_values, y4_values, color = 'blue')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out