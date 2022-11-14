import numpy as np 
import pandas as pd
import math
import matplotlib.pyplot as plt
# pip3 install numdifftools
import numdifftools as nd


class Parabola:
    def __init__(someobj,a,b,c,initialPoint,directionVec):
        someobj.a = a
        someobj.b = b
        someobj.c = c
        someobj.initPoint = initialPoint
        someobj.initDir = directionVec

    def reflect_vert(this, ray_only, simpleView, isPlot): # parabola 1

        print('VERTICAL PARABOLA')

        a = this.a
        b = this.b
        c = this.c
        initPoint = this.initPoint
        initDir = this.initDir

        def func(x): 
            return a*(x[0]**2) + b*x[0] + c - x[1]

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variable for plotting & incident ray
        x = np.linspace(-10, 10, 1000)
        if not(ray_only):
            t = np.linspace(0, 10, 500)
            plt.plot(initPoint[0] + t*initDir[0], initPoint[1] + t*initDir[1],'black')

        print("Ray Direction: ", initDir[0], initDir[1])
        print('The coefficient of x^2', a)
        print('The coefficient of x', b)
        print('The constant term', c)


        ##Equation of curve
        print("Equation of curve,  y =", a,"x^2 +", b,"x +", c)
        plt.plot(x, a*x*x + b*x + c, 'red')
        

        ##Find where ray and curve intersect
        a_term = a * (initDir[0]**2)
        b_term = (2 * a * initPoint[0] * initDir[0]) + (b * initDir[0]) - initDir[1]
        c_term = (a * (initPoint[0]**2)) + (b * initPoint[0]) + c - initPoint[1]
        if b_term*b_term - 4 * a_term * c_term < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        
        if distance_add < 0 and distance_min < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        elif distance_add < 0:
            distance_true = distance_min
        elif distance_min < 0:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        if ray_only:
            t2 = np.linspace(0, distance_true, 100)
            plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')
        

        ##Differentiate above equation to find tangent dy/dx
        # tangDir = np.array([1, 2*a*intercept[0] + b])

        # ##Equation of tangent
        # print("tangDir: ", tangDir)
        # if not(ray_only) and not(simpleView):
        # 	plt.plot(intercept[0] + x*tangDir[0], intercept[1] + x*tangDir[1], 'blue')


        ##Slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
        # normDir = np.array([tangDir[1], -tangDir[0]])

        normDir = nd.Gradient(func)(intercept)
        # print("Gradient of a(x^2)+bx+c-y at intercept is ", grad2)

        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only) and not(simpleView):
            plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        rayDir = np.array([initDir[0], initDir[1]])
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out
    
    def reflect_horiz(this, ray_only, simpleView, isPlot):  # parabola 2

        print('HORIZONTAL PARABOLA')

        a = this.a
        b = this.b
        c = this.c
        initPoint = this.initPoint
        initDir = this.initDir

        def func(x): 
            return a*(x[1]**2) + b*x[1] + c - x[0]

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variable for plotting & incident ray
        x = np.linspace(-10, 10, 1000)
        if not(ray_only):
            t = np.linspace(0, 10, 500)
            plt.plot(initPoint[0] + t*initDir[0], initPoint[1] + t*initDir[1],'black')
        
        print("Ray Direction: ", initDir[0], initDir[1])
        print('The coefficient of y^2', a)
        print('The coefficient of y', b)
        print('The constant term', c)


        ##Equation of curve
        print("Equation of curve,  x =", a,"y^2 +", b,"y +", c)
        plt.plot(a*x*x + b*x + c, x, 'red')
        

        ##Find where ray and curve intersect
        a_term = a * (initDir[1]**2)
        b_term = (2 * a * initPoint[1] * initDir[1]) + (b * initDir[1]) - initDir[0]
        c_term = (a * (initPoint[1]**2)) + (b * initPoint[1]) + c - initPoint[0]
        if b_term*b_term - 4 * a_term * c_term < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        
        if distance_add < 0 and distance_min < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        elif distance_add < 0:
            distance_true = distance_min
        elif distance_min < 0:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        if ray_only:
            t2 = np.linspace(0, distance_true, 100)
            plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')
            

        ##Differentiate above equation to find tangent dy/dx
        # tangDir = np.array([2*a*intercept[1] + b, 1])

        # ##Direction of tangent
        # print("tangDir: ", tangDir)
        # if not(ray_only) and not(simpleView):
        # 	plt.plot(intercept[0] + x*tangDir[0], intercept[1] + x*tangDir[1], 'blue')


        ##Slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
        # normDir = np.array([tangDir[1], -tangDir[0]])

        normDir = nd.Gradient(func)(intercept)
        
        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only) and not(simpleView):
            plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        rayDir = np.array([initDir[0], initDir[1]])
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out


class Rotated_Parabola:
    def __init__(someobj,a,b,c,initialPoint,directionVec):
        someobj.a = a
        someobj.b = b
        someobj.c = c
        someobj.initPoint = initialPoint
        someobj.initDir = directionVec

    def reflect(this, ray_only, simpleView, isPlot):

        print('ROTATED PARABOLA')

        a = this.a
        b = this.b
        c = this.c
        initPoint = this.initPoint
        initDir = this.initDir
        theta = math.pi/2
        sin = math.sin(theta)
        cos = math.cos(theta)

        def func(x): 
            return (a*(x[0]**2)*(cos**2)) + (a*(x[1]**2)*(sin**2)) - (2*a*x[0]*x[1]*cos*sin) + (b*x[0]*cos) - (b*x[1]*sin) - (x[0]*sin) - (x[1]*cos) + c

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variable for plotting & incident ray
        x = np.linspace(-10, 10, 1000)
        if not(ray_only):
            t = np.linspace(0, 10, 500)
            plt.plot(initPoint[0] + t*initDir[0], initPoint[1] + t*initDir[1],'black')

        print("Ray Direction: ", initDir[0], initDir[1])
        print('The coefficient of x^2', a)
        print('The coefficient of x', b)
        print('The constant term', c)


        ##Equation of curve
        print("Equation of curve,  y =", a,"x^2 +", b,"x +", c)
        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)
        plt.contour(hypx, hypy,((a*(hypx**2)*(cos(theta)**2)) + (a*(hypy**2)*(sin(theta)**2)) - (2*a*hypx*hypy*cos(theta)*sin(theta)) + (b*hypx*cos(theta)) - (b*hypy*sin(theta)) - (hypx*sin(theta)) - (hypy*cos(theta)) + c), [0], colors='red')
        
        ##Find where ray and curve intersect
        a_term = a * (((initDir[0]**2) * (cos**2)) + ((initDir[1]**2) * (sin**2)) - (2 * initDir[0] * initDir[1] * cos * sin))
        b_term = ((2 * a * initPoint[0] * initDir[0]*(cos**2)) + (2 * a * initPoint[1] * initDir[1]*(sin**2)) 
            - (2 * a * initPoint[1] * initDir[0] * cos * sin) - (2 * a * initPoint[0] * initDir[1] * cos * sin) 
            + (b * initDir[0] * cos) - (b * initDir[1] * sin) - (initDir[0] * sin) - (initDir[1] * cos))
        c_term = ((a * (initPoint[0]**2) * (cos**2)) + (a * (initPoint[1]**2) * (sin**2)) 
            - (2 * a * initPoint[0] * initPoint[1] * cos * sin) + (b * initPoint[0] * cos) 
            - (b * initPoint[1] * sin) - (initPoint[0] * sin) - (initPoint[1] * cos) + c)
        if b_term*b_term - 4 * a_term * c_term < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        
        if distance_add < 0 and distance_min < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        elif distance_add < 0:
            distance_true = distance_min
        elif distance_min < 0:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        if ray_only:
            t2 = np.linspace(0, distance_true, 100)
            plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')


        ##Normal vector is gradient of function
        normDir = nd.Gradient(func)(intercept)

        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only) and not(simpleView):
            plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        rayDir = np.array([initDir[0], initDir[1]])
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out


class Linear:

    def __init__(someobj,a,b,initialPoint,directionVec):
        someobj.a = a
        someobj.b = b
        someobj.initPoint = initialPoint
        someobj.initDir = directionVec

    def reflect(this, ray_only, simpleView, isPlot):

        print('PLANE')

        a = this.a
        b = this.b
        initPoint = this.initPoint
        initDir = this.initDir

        def func(x): 
            return a*x[0] + b - x[1]

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variable for plotting & incident ray
        x = np.linspace(-10, 10, 1000)
        if not(ray_only):
            t = np.linspace(0, 10, 500)
            plt.plot(initPoint[0] + t*initDir[0], initPoint[1] + t*initDir[1],'black')

        print("Ray Direction: ", initDir[0], initDir[1])
        print('The coefficient of x', a)
        print('The constant term', b)


        ##Equation of curve
        print("Equation of curve,  y =", a,"x +", b)
        plt.plot(x, a*x + b, 'red')


        ##Find where ray and curve intersect
        if (initDir[1] - (a * initDir[0])) == 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        distance = ((a * initPoint[0]) + b - initPoint[1]) / (initDir[1] - (a * initDir[0]))

        if distance < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
            
        intercept = np.array([initPoint[0] + distance*initDir[0], initPoint[1] + distance*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        if ray_only:
            t2 = np.linspace(0, distance, 100)
            plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')
        

        ##Normal vector is gradient of function
        normDir = nd.Gradient(func)(intercept)
        
        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only) and not(simpleView):
            plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        rayDir = np.array([initDir[0], initDir[1]])
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out


class Vector_Linear:

    def __init__(someobj,initialPoint,directionVec,surfPoint,surfDir):
        someobj.initPoint = initialPoint
        someobj.initDir = directionVec
        someobj.surfPoint = surfPoint
        someobj.surfDir = surfDir   

    def reflect(this, ray_only, simpleView, isPlot):

        print('PLANE')

        surfPoint = this.surfPoint
        surfDir = this.surfDir
        initPoint = this.initPoint
        initDir = this.initDir

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variable for plotting & incident ray
        x = np.linspace(-10, 10, 1000)
        if not(ray_only):
            t = np.linspace(0, 10, 500)
            plt.plot(initPoint[0] + t*initDir[0], initPoint[1] + t*initDir[1],'black')

        print("Ray Direction: ", initDir[0], initDir[1])


        ##Equation of curve
        print("Equation of curve:", this.surfPoint,"+ t*", this.surfDir)
        plt.plot(this.surfPoint[0] + x*this.surfDir[0], this.surfPoint[1] + x*this.surfDir[1],'red')


        ##Find where ray and curve intersect
        a = this.surfDir[1]/this.surfDir[0]
        if (initDir[1] - (a * initDir[0])) == 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        distance = (surfPoint[1] - initPoint[1] + (a * (initPoint[0] - surfPoint[0]))) / (initDir[1] - (a * initDir[0]))

        if distance < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
            
        intercept = np.array([initPoint[0] + distance*initDir[0], initPoint[1] + distance*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        if ray_only:
            t2 = np.linspace(0, distance, 100)
            plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')


        ##Normal vector is gradient of function
        normDir = np.array([surfDir[1], -surfDir[0]])
        
        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only) and not(simpleView):
            plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        rayDir = np.array([initDir[0], initDir[1]])
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out

    
class Hyperbola:

    def __init__(someobj,a,b,h,k,initialPoint,directionVec):
        someobj.a = a
        someobj.b = b
        someobj.h = h
        someobj.k = k
        someobj.initPoint = initialPoint
        someobj.initDir = directionVec

    def reflect(this, ray_only, simpleView, isPlot):

        print('HYPERBOLA')

        a = this.a
        b = this.b
        h = this.h
        k = this.k
        initPoint = this.initPoint
        initDir = this.initDir

        def func(x): 
            return (b**2)*((x[0]-h)**2) - (a**2)*((x[1]-k)**2) -(a**2)*(b**2)

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variable for plotting & incident ray
        x = np.linspace(-10, 10, 1000)
        if not(ray_only):
            t = np.linspace(0, 10, 500)
            plt.plot(initPoint[0] + t*initDir[0], initPoint[1] + t*initDir[1],'black')

        print("Ray Direction: ", initDir[0], initDir[1])
        print('The coefficient of x^2', b*b)
        print('The coefficient of y^2', -a*a)
        print('The coefficient of x', -2*b*b*h)
        print('The coefficient of y', 2*a*a*k)
        print('The constant term', b*b*h*h - a*a*k*k - a*a*b*b)


        ##Equation of curve
        print("Equation of curve,  0 =", b*b,"x^2 +", -a*a,"y^2 +", -2*b*b*h,"x +", 2*a*a*k,"y +", b*b*h*h - a*a*k*k - a*a*b*b)	
        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)
        plt.contour(hypx, hypy,(((hypx - h)**2)/a**2 - ((hypy - k)**2)/b**2), [1], colors='red')

        ##Find where ray and curve intersect
        a_term = (b**2) * (initDir[0]**2) - (a**2) * (initDir[1]**2)
        b_term = 2 * ((b**2) * initPoint[0] * initDir[0] - (b**2) * h * initDir[0] 
            - (a**2) * initPoint[1] * initDir[1] + (a**2) * k * initDir[1])
        c_term = ((b**2) * ((initPoint[0] - h)**2) - (a**2) * ((initPoint[1] - k)**2) - (a**2) * (b**2))
        if b_term*b_term - 4 * a_term * c_term < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
    
        if distance_add < 0 and distance_min < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        elif distance_add < 0:
            distance_true = distance_min
        elif distance_min < 0:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        if ray_only:
            t2 = np.linspace(0, distance_true, 100)
            plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')


        ##Normal vector is gradient of function
        normDir = nd.Gradient(func)(intercept)
    
        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only) and not(simpleView):
            plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')
        

        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        rayDir = np.array([initDir[0], initDir[1]])
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1
        out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out


class Rotated_Hyperbola:

    def __init__(someobj,a,b,h,k,initialPoint,directionVec):
        someobj.a = a
        someobj.b = b
        someobj.h = h
        someobj.k = k
        someobj.initPoint = initialPoint
        someobj.initDir = directionVec

    def reflect(this, ray_only, simpleView, isPlot):

        print('ROTATED HYPERBOLA')

        a = this.a
        b = this.b
        h = this.h
        k = this.k
        initPoint = this.initPoint
        initDir = this.initDir
        theta = 4*math.pi/3
        sin = math.sin(theta)
        cos = math.cos(theta)

        def func(x): 
            return (b**2) * (((x[1] * sin) - (x[0] * cos) + h)**2) - (a**2) * (((x[1] * cos) + (x[0] * sin) - k)**2) - ((a**2) * (b**2))

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variable for plotting & incident ray
        x = np.linspace(-10, 10, 1000)
        if not(ray_only):
            t = np.linspace(0, 10, 500)
            plt.plot(initPoint[0] + t*initDir[0], initPoint[1] + t*initDir[1],'black')

        print("Ray Direction: ", initDir[0], initDir[1])
        print('The coefficient of x^2', b*b)
        print('The coefficient of y^2', -a*a)
        print('The coefficient of x', -2*b*b*h)
        print('The coefficient of y', 2*a*a*k)
        print('The constant term', b*b*h*h - a*a*k*k - a*a*b*b)


        ##Equation of curve
        print("Equation of curve,  0 =", b*b,"x^2 +", -a*a,"y^2 +", -2*b*b*h,"x +", 2*a*a*k,"y +", b*b*h*h - a*a*k*k - a*a*b*b)	
        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)
        plt.contour(hypx, hypy,((b**2) * (((hypy * sin) - (hypx * cos) + h)**2) - (a**2) * (((hypy * cos) + (hypx * sin) - k)**2) - ((a**2) * (b**2))), [0], colors='red')

        ##Find where ray and curve intersect
        a_term = (((b**2)*(initDir[0]**2)*(cos**2)) - ((a**2)*(initDir[0]**2)*(sin**2)) 
            - (2*(b**2)*initDir[0]*initDir[1]*cos*sin) - (2*(a**2)*initDir[0]*initDir[1]*cos*sin) 
            + ((b**2)*(initDir[1]**2)*(sin**2)) - ((a**2)*(initDir[1]**2)*(cos**2)))
        b_term = 2 * (((b**2)*initPoint[0]*initDir[0]*(cos**2)) - ((a**2)*initPoint[0]*initDir[0]*(sin**2)) 
            - ((b**2)*initPoint[1]*initDir[0]*cos*sin) - ((b**2)*initPoint[0]*initDir[1]*cos*sin) 
            - ((a**2)*initPoint[1]*initDir[0]*cos*sin) - ((a**2)*initPoint[0]*initDir[1]*cos*sin)
            + ((b**2)*initPoint[1]*initDir[1]*(sin**2)) - ((a**2)*initPoint[1]*initDir[1]*(cos**2))
            - ((b**2)*h*initDir[0]*cos) + ((a**2)*k*initDir[0]*sin) + ((b**2)*h*initDir[1]*sin) + ((a**2)*k*initDir[1]*cos))
        c_term = (((b**2)*(initPoint[0]**2)*(cos**2)) - ((a**2)*(initPoint[0]**2)*(sin**2)) 
            - (2*(b**2)*initPoint[0]*initPoint[1]*cos*sin) - (2*(a**2)*initPoint[0]*initPoint[1]*cos*sin) 
            + ((b**2)*(initPoint[1]**2)*(sin**2)) - ((a**2)*(initPoint[1]**2)*(cos**2)) 
            - (2*(b**2)*h*initPoint[0]*cos) + (2*(a**2)*k*initPoint[0]*sin) + (2*(b**2)*h*initPoint[1]*sin) 
            + (2*(a**2)*k*initPoint[1]*cos) + ((b**2)*(h**2)) - ((a**2)*(k**2)) - ((a**2)*(b**2)))
        if b_term*b_term - 4 * a_term * c_term < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
    
        if distance_add < 0 and distance_min < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        elif distance_add < 0:
            distance_true = distance_min
        elif distance_min < 0:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        if ray_only:
            t2 = np.linspace(0, distance_true, 100)
            plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')


        ##Normal vector is gradient of function
        normDir = nd.Gradient(func)(intercept)
    
        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only) and not(simpleView):
            plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')
        

        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        rayDir = np.array([initDir[0], initDir[1]])
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1
        out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out


class Ellipse:

    def __init__(someobj,a,b,h,k,initialPoint,directionVec):
        someobj.a = a
        someobj.b = b
        someobj.h = h
        someobj.k = k
        someobj.initPoint = initialPoint
        someobj.initDir = directionVec

    def reflect(this, ray_only, simpleView, isPlot):

        print('ELLIPSE')

        a = this.a
        b = this.b
        h = this.h
        k = this.k
        initPoint = this.initPoint
        initDir = this.initDir

        def func(x): 
            return (b**2) * ((x[0] - h)**2) + (a**2) * ((x[1] - k)**2) - ((a**2) * (b**2))

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variables for plotting & incident ray
        x = np.linspace(-5, 5, 100)
        if not(ray_only):
            t = np.linspace(0, 5, 50)
            plt.plot(initPoint[0] + t*initDir[0], initPoint[1] + t*initDir[1],'black')

        print("Ray Direction: ", initDir[0], initDir[1])
        print('The coefficient of x^2: ', b*b)
        print('The coefficient of y^2: ', a*a)
        print('The coefficient of x: ', -2*b*b*h)
        print('The coefficient of y: ', -2*a*a*k)
        print('The constant term: ', b*b*h*h + a*a*k*k - a*a*b*b)


        ##Equation of curve
        print("Equation of curve,  0 =", b*b,"x^2 +", a*a,"y^2 +", -2*b*b*h,"x +", -2*a*a*k,"y +", b*b*h*h + a*a*k*k - a*a*b*b)
        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)
        plt.contour(hypx, hypy,(((hypx - h)**2)/a**2 + ((hypy - k)**2)/b**2), [1], colors='red')

        # t = np.linspace(0, 2*math.pi, 100)
        # plt.plot(h + a*np.cos(t), k + b*np.sin(t), 'red')


        ##Find where ray and curve intersect
        a_term = (b**2) * (initDir[0]**2) + (a**2) * (initDir[1]**2)
        b_term = 2 * ((b**2) * initPoint[0] * initDir[0] - (b**2) * h * initDir[0] 
            + (a**2) * initPoint[1] * initDir[1] - (a**2) * k * initDir[1])
        c_term = ((b**2) * ((initPoint[0] - h)**2) + (a**2) * ((initPoint[1] - k)**2) - (a**2) * (b**2))
        if b_term*b_term - 4 * a_term * c_term < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)

        if distance_add < 0 and distance_min < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        elif distance_add < 0:
            distance_true = distance_min
        elif distance_min < 0:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        if ray_only:
            t2 = np.linspace(0, distance_true, 100)
            plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')


        ##Normal vector is gradient of function
        normDir = nd.Gradient(func)(intercept)

        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only) and not(simpleView):
            plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')	


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        rayDir = np.array([initDir[0], initDir[1]])
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')
        

        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out


class Rotated_Ellipse:

    def __init__(someobj,a,b,h,k,initialPoint,directionVec):
        someobj.a = a
        someobj.b = b
        someobj.h = h
        someobj.k = k
        someobj.initPoint = initialPoint
        someobj.initDir = directionVec

    def reflect(this, ray_only, simpleView, isPlot):

        print('ROTATED ELLIPSE')

        a = this.a
        b = this.b
        h = this.h
        k = this.k
        initPoint = this.initPoint
        initDir = this.initDir
        theta = math.pi/3
        sin = math.sin(theta)
        cos = math.cos(theta)

        def func(x): 
            return (b**2) * (((x[1] * sin) - (x[0] * cos) + h)**2) + (a**2) * (((x[1] * cos) + (x[0] * sin) - k)**2) - ((a**2) * (b**2))

        def func1(x): # Equation I found online, rotates in other direction? Which is correct? (Mine matches the parabola)
            return (b**2) * (((x[0] - h) * cos + (x[1] - k) * sin)**2) + (a**2) * (((x[0] - h) * sin - (x[1] - k) * cos)**2) - ((a**2) * (b**2))

        if isPlot:
            plt.grid(color='lightgray',linestyle='--')
            plt.xlim(-10, 10)
            plt.ylim(-10, 10)
            plt.gca().set_aspect('equal', adjustable='box')

        ##Variables for plotting & incident ray
        x = np.linspace(-5, 5, 100)
        if not(ray_only):
            t = np.linspace(0, 5, 50)
            plt.plot(initPoint[0] + t*initDir[0], initPoint[1] + t*initDir[1],'black')

        print("Ray Direction: ", initDir[0], initDir[1])
        print('The coefficient of x^2: ', b*b)
        print('The coefficient of y^2: ', a*a)
        print('The coefficient of x: ', -2*b*b*h)
        print('The coefficient of y: ', -2*a*a*k)
        print('The constant term: ', b*b*h*h + a*a*k*k - a*a*b*b)


        ##Equation of curve
        print("Equation of curve,  0 =", b*b,"x^2 +", a*a,"y^2 +", -2*b*b*h,"x +", -2*a*a*k,"y +", b*b*h*h + a*a*k*k - a*a*b*b)
        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)
        plt.contour(hypx, hypy,((b**2) * (((hypy * sin) - (hypx * cos) + h)**2) + (a**2) * (((hypy * cos) + (hypx * sin) - k)**2) - ((a**2) * (b**2))), [0], colors='red')
        
        #Ellipse rotated in 'correct' direction
        plt.contour(hypx, hypy,((b**2) * (((hypx - h) * cos + (hypy - k) * sin)**2) + (a**2) * (((hypx - h) * sin - (hypy - k) * cos)**2) - ((a**2) * (b**2))), [0], colors='purple')
        # t = np.linspace(0, 2*math.pi, 100)
        # plt.plot(h + a*np.cos(t), k + b*np.sin(t), 'red')


        ##Find where ray and curve intersect
        #Can this be generalized?
        a_term = (((b**2)*(initDir[0]**2)*(cos**2)) + ((a**2)*(initDir[0]**2)*(sin**2)) 
            - (2*(b**2)*initDir[0]*initDir[1]*cos*sin) + (2*(a**2)*initDir[0]*initDir[1]*cos*sin) 
            + ((b**2)*(initDir[1]**2)*(sin**2)) + ((a**2)*(initDir[1]**2)*(cos**2)))
        b_term = 2 * (((b**2)*initPoint[0]*initDir[0]*(cos**2)) + ((a**2)*initPoint[0]*initDir[0]*(sin**2)) 
            - ((b**2)*initPoint[1]*initDir[0]*cos*sin) - ((b**2)*initPoint[0]*initDir[1]*cos*sin) 
            + ((a**2)*initPoint[1]*initDir[0]*cos*sin) + ((a**2)*initPoint[0]*initDir[1]*cos*sin)
            + ((b**2)*initPoint[1]*initDir[1]*(sin**2)) + ((a**2)*initPoint[1]*initDir[1]*(cos**2))
            - ((b**2)*h*initDir[0]*cos) - ((a**2)*k*initDir[0]*sin) + ((b**2)*h*initDir[1]*sin) - ((a**2)*k*initDir[1]*cos))
        c_term = (((b**2)*(initPoint[0]**2)*(cos**2)) + ((a**2)*(initPoint[0]**2)*(sin**2)) 
            - (2*(b**2)*initPoint[0]*initPoint[1]*cos*sin) + (2*(a**2)*initPoint[0]*initPoint[1]*cos*sin) 
            + ((b**2)*(initPoint[1]**2)*(sin**2)) + ((a**2)*(initPoint[1]**2)*(cos**2)) 
            - (2*(b**2)*h*initPoint[0]*cos) - (2*(a**2)*k*initPoint[0]*sin) + (2*(b**2)*h*initPoint[1]*sin) 
            - (2*(a**2)*k*initPoint[1]*cos) + ((b**2)*(h**2)) + ((a**2)*(k**2)) - ((a**2)*(b**2)))
        if b_term*b_term - 4 * a_term * c_term < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)

        if distance_add < 0 and distance_min < 0:
            print('No intercept')
            if isPlot:
                plt.show()
            return initPoint, initDir
        elif distance_add < 0:
            distance_true = distance_min
        elif distance_min < 0:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)

        intercept = np.array([initPoint[0] + distance_true*initDir[0], initPoint[1] + distance_true*initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        ##Plot where the ray actually goes (initial point until intercept)
        if ray_only:
            t2 = np.linspace(0, distance_true, 100)
            plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')


        ##Normal vector is gradient of function
        normDir = nd.Gradient(func)(intercept)

        ##Direction of normal
        print("normDir: ", normDir)
        if not(ray_only) and not(simpleView):
            plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')	


        ##Finding equation of reflected line by projecting onto norm and subtracting vectors
        rayDir = np.array([initDir[0], initDir[1]])
        normNorm = normDir/np.linalg.norm(normDir)
        #Citation 1 
        out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

        ##Direction of reflected line
        print("Output Direction: ", out)
        if not(ray_only):
            plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')
        

        ##Print plot
        if isPlot:
            plt.show()

        return intercept, out