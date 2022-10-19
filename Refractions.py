import numpy as np 
import pandas as pd
import math
import matplotlib.pyplot as plt


#Previous implementation: block with flat surfaces
class Glass_objects: # glass block or otherwise
    def __init__(someobj,m,c,a1,b1,m1,h,n1,n2):	
        someobj.m = m
        someobj.c = c
        someobj.a1 = a1 # x coord of point from where incident
        someobj.b1 = b1 # y coord of point from where incident
        someobj.m1 = m1
        someobj.h = h   # height of glass
        someobj.n1 = n1 # medium of input
        someobj.n2 = n2 # medium of output

    def block_type1(this): # type 1 glass block

        x = np.linspace(-5, 5, 100)

        print("Equation of glass block starting side, y =",this.m,"x +",+ this.c)

        c1 = this.b1 - this.a1*this.m1  # calculating constant intercept of incident ray
        print("Equation of incident ray, y =",this.m1,"x +",+ c1)

        x_inc = (c1 - this.c)/(this.m - this.m1)
        y_inc = (this.m*(c1 - this.c))/(this.m - this.m1)  + this.c 

        print("Coordinates of point at which incident ray hits glass,  x coord , ",x_inc,"y coord", y_inc)

        theta1 = math.atan(this.m1)

        # using Snell's law, n1sin(theta1) = n2sin(theta2)
        # let slope of refracted ray be m2 , then m2 is -

        m2 = math.tan(np.arcsin((this.n1*math.sin(theta1))/this.n2))


        #constant intercept of refracted ray

        c2 = y_inc - m2*x_inc

        print("Equation of refracted ray, y =",m2,"x +",+ c2)

        # equation of line parallel to the glass starting line, with height difference h (so, other side of the planar glass)

        c3 = this.h*math.sqrt(this.m*this.m + 1) + this.c


        print("Equation of glass block ending side, y =",this.m,"x +",+ c3)

        # point coordinates of where the glass ends and so does the refracted ray
        # there are two possibilities, given the absolute sign of the equation on the numerator

        x_ref_exit = (c3 - c2)/(m2 - this.m)

        y_ref_exit = (this.m*(c3 - c2))/(m2 - this.m) + c3

        print("Coordinates of point at which refracted ray exits glass,  x coord , ",x_ref_exit,"y coord", y_ref_exit)

#Flat surface, no longer a block. Slightly vectorized. Theta incorrect
class Medium_surfaces: # flat surface between two mediums
    def __init__(someobj,n1,n2,initDir,initPoint,surfDir,surfPoint):	
        someobj.n1 = n1 # medium of input
        someobj.n2 = n2 # medium of output
        someobj.initDir = initDir
        someobj.initPoint = initPoint
        someobj.surfDir = surfDir
        someobj.surfPoint = surfPoint

    def block_type1(this): # type 1 glass block

        #m1 -> initDir
        #m -> surfDir
        #c -> surfPoint

        m = this.surfDir[1]/this.surfDir[0]
        m1 = this.initDir[1]/this.initDir[0]
        c = this.surfPoint[1]

        print("Equation of glass block starting side, y =", m,"x +",+ c)

        c1 = this.initPoint[1] - this.initPoint[0]*m1  # calculating constant intercept of incident ray
        print("Equation of incident ray, y =",m1,"x +",+ c1)

        x_inc = (c1 - c)/(m - m1)
        y_inc = (m*(c1 - c))/(m - m1)  + c 

        print("Coordinates of point at which incident ray hits glass,  x coord , ",x_inc,"y coord", y_inc)

        theta1 = math.atan(m1)

        # using Snell's law, n1sin(theta1) = n2sin(theta2)
        # let slope of refracted ray be m2 , then m2 is -

        m2 = math.tan(np.arcsin((this.n1*math.sin(theta1))/this.n2))


        #constant intercept of refracted ray

        c2 = y_inc - m2*x_inc

        print("Equation of refracted ray, y =",m2,"x +",+ c2)

        return np.array([x_inc,y_inc]),np.array([1,m2])


#Flat surface, completely vectorized with new Snell's Law for theta
class Medium_surfaces2: # flat surface between two mediums
    def __init__(someobj,n1,n2,initPoint,initDir,surfPoint,surfDir):	
        someobj.n1 = n1 # medium of input
        someobj.n2 = n2 # medium of output
        someobj.initPoint = initPoint
        someobj.initDir = initDir
        someobj.surfPoint = surfPoint
        someobj.surfDir = surfDir   

    def block_type1(this): # type 1 glass block

        x = np.linspace(-5, 5, 100)
        t = np.linspace(0, 10, 500)

        #m1 -> initDir
        #m -> surfDir
        #c -> surfPoint

        m = this.surfDir[1]/this.surfDir[0]
        m1 = this.initDir[1]/this.initDir[0]
        c = this.surfPoint[1]

        print("Vector equation of surface:", this.surfPoint,"+ t*", this.surfDir)
        plt.plot(this.surfPoint[0] + x*this.surfDir[0], this.surfPoint[1] + x*this.surfDir[1],'blue')
    
        print("Vector incident ray equation:", this.initPoint,"+ s*", this.initDir)
        plt.plot(this.initPoint[0] + t*this.initDir[0], this.initPoint[1] + t*this.initDir[1],'black')

        a = this.surfDir[1]/this.surfDir[0]
        # if (initDir[1] - (a * initDir[0])) == 0:
        # 	print('No intercept')
        # 	if isPlot:
        # 		plt.show()
        # 	return initPoint, initDir
        distance = (this.surfPoint[1] - this.initPoint[1] + (a * (this.initPoint[0] - this.surfPoint[0]))) / (this.initDir[1] - (a * this.initDir[0]))

        intercept = np.array([this.initPoint[0] + distance*this.initDir[0], this.initPoint[1] + distance*this.initDir[1]])

        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        normDir = np.array([this.surfDir[1], -this.surfDir[0]])
        if np.dot(this.initDir, normDir) < 0:
            normDir = -normDir
        print("Norm direction:", normDir)

        normNorm = normDir/np.linalg.norm(normDir)

        initNorm = this.initDir/np.linalg.norm(this.initDir)

        mu = this.n1/this.n2
        
        outDir = (mu*initNorm) + (normNorm*np.sqrt(1-(mu*mu*(1-((np.dot(normNorm, initNorm))**2))))) - (mu*np.dot(normNorm, np.dot(normNorm, initNorm)))


        print("Output direction:", outDir)
        plt.plot(intercept[0] + t*outDir[0], intercept[1] + t*outDir[1],'green')

        theta1 = math.atan(m1)

        # using Snell's law, n1sin(theta1) = n2sin(theta2)
        # let slope of refracted ray be m2 , then m2 is -

        m2 = math.tan(np.arcsin((this.n1*math.sin(theta1))/this.n2))


        #constant intercept of refracted ray

        c1 = this.initPoint[1] - this.initPoint[0]*m1
        x_inc = (c1 - c)/(m - m1)
        y_inc = (m*(c1 - c))/(m - m1)  + c 
        c2 = y_inc - m2*x_inc

        print("Equation of refracted ray, y =",m2,"x +",+ c2)
        plt.plot(x, m2*x + c2, 'red')

        

        return intercept,outDir

#Curved surface, vectorized
class Circular_Lens:
    def __init__(someobj,n1,n2,r,initPoint,initDir,centerPoint):
        someobj.n1 = n1 # refractive index of medium 1
        someobj.n2 = n2 # refractive index of medium 2
        someobj.r = r# radius of curvature
        someobj.initPoint = initPoint
        someobj.initDir = initDir
        someobj.center = centerPoint # center of circle for lens
		
    def block_type2(this): # type 2 convex surface

        x = np.linspace(-10, 10, 1000)
        t = np.linspace(0, 10, 500)

        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)
        plt.contour(hypx, hypy,(((hypx - this.center[0])**2) + ((hypy - this.center[1])**2)), [this.r**2], colors='blue')
    
        print("Vector incident ray equation:", this.initPoint,"+ s*", this.initDir)
        plt.plot(this.initPoint[0] + t*this.initDir[0], this.initPoint[1] + t*this.initDir[1],'black')

        #find intercept
        a_term = this.initDir[0] * this.initDir[0] + this.initDir[1] * this.initDir[1]
        b_term = 2 * (this.initPoint[0] * this.initDir[0] - this.center[0] * this.initDir[0] 
            + this.initPoint[1] * this.initDir[1] - this.center[1] * this.initDir[1])
        c_term = ((this.initPoint[0] - this.center[0]) * (this.initPoint[0] - this.center[0]) 
            + (this.initPoint[1] - this.center[1]) * (this.initPoint[1] - this.center[1]) - this.r**2)

        distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
        distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)

        if distance_add < 0 and distance_min < 0:
            print('No intercept')
            # if isPlot:
            # 	plt.show()
            return initPoint, initDir
        elif distance_add < 0:
            distance_true = distance_min
        elif distance_min < 0:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)
	

        intercept = np.array([this.initPoint[0] + distance_true*this.initDir[0], this.initPoint[1] + distance_true*this.initDir[1]])
        print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        #find tangent
        ##Differentiate above equation to find tangent dy/dx
        tangDir = np.array([intercept[1] - this.center[1], -(intercept[0] - this.center[0])])
		
        ##Direction of tangent
        print("tangDir: ", tangDir)
        # plt.plot(intercept[0] + x*tangDir[0], intercept[1] + x*tangDir[1], 'red')


        #find norm
        normDir = np.array([tangDir[1], -tangDir[0]])
        if np.dot(this.initDir, normDir) < 0:
            normDir = -normDir

        print("NormDot:", np.dot(this.initDir, normDir))

        ##Direction of normal
        print("normDir: ", normDir)
        # plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')	

        
        normNorm = normDir/np.linalg.norm(normDir)
        initNorm = this.initDir/np.linalg.norm(this.initDir)
        mu = this.n1/this.n2
        outDir = (mu*initNorm) + (normNorm*np.sqrt(1-(mu*mu*(1-((np.dot(normNorm, initNorm))**2))))) - (mu*np.dot(normNorm, np.dot(normNorm, initNorm)))

        print("Output direction:", outDir)
        plt.plot(intercept[0] + t*outDir[0], intercept[1] + t*outDir[1],'green')

        return intercept, outDir

#Convex lens using two of same curved surfaces
class Lens: 
    def __init__(someobj,n1,n2,r,h,initPoint,initDir,centerPoint):
        someobj.n1 = n1 # refractive index of medium 1
        someobj.n2 = n2 # refractive index of medium 2
        someobj.r = r   # radius of curvature
        someobj.h = h
        someobj.initPoint = initPoint
        someobj.initDir = initDir
        someobj.center = centerPoint # center of lens
		
    def block_type2(this): # type 2 convex surface
        conv1 = Circular_Lens(n1,n2,this.r,this.initPoint,this.initDir,np.array([this.center[0]-(this.h/2)+this.r,this.center[1]]))
        nextPoint, nextDir = conv1.block_type2()
        print("HERE", nextPoint, nextDir,np.array([this.center[0]+(this.h/2)-this.r,this.center[1]]))
        conv1 = Circular_Lens(n2,n1,this.r,nextPoint,nextDir,np.array([this.center[0]+(this.h/2)-this.r,this.center[1]]))
        nextPoint, nextDir = conv1.block_type2()

        plt.show()


#let the medium 1 be air, the medium 2 be glass

n1 = 1.0003
n2 = 1.52

# glass1 = Glass_objects(1,1,-1,-1,2,3,n1,n2)
# glass1.block_type1()	

# glass2 = Medium_surfaces(n1,n2, np.array([1,2]), np.array([-1,-1]),np.array([1,1]), np.array([0,1]))
# refPt, refDir = glass2.block_type1()

# glassOut = Medium_surfaces(n2,n1,refDir,refPt,np.array([1,1]),np.array([0,5.242640687119286]))
# refPt, refDir = glassOut.block_type1()

plt.grid(color='lightgray',linestyle='--')
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.gca().set_aspect('equal', adjustable='box')

# glass2 = Medium_surfaces2(n1,n2, np.array([-1,-1]),np.array([1,2]), np.array([0,1]), np.array([1,1]))
# refPt, refDir = glass2.block_type1()

# glass2 = Medium_surfaces2(n2,n1, refPt, refDir, np.array([0,5.242640687119286]),np.array([1,1]))
# refPt, refDir = glass2.block_type1()

# conv1 = Circular_Lens(n1,n2,3,np.array([-2,2.5]),np.array([1,0]),np.array([3,3]))
# conv1.block_type2()

# plt.show()

lens = Lens(n1,n2,5,1,np.array([-2,4]),np.array([1,0]),np.array([0,3]))
lens.block_type2()
