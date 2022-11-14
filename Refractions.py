import numpy as np 
import pandas as pd
import math
import matplotlib.pyplot as plt
# pip3 install numdifftools
import numdifftools as nd

#Flat surface, completely vectorized with new Snell's Law for theta
class Medium_flat_surfaces: # flat surface between two mediums
    def __init__(someobj,n1,n2,initPoint,initDir,surfPoint,surfDir, notLens="True"):	
        someobj.n1 = n1 # medium of input
        someobj.n2 = n2 # medium of output
        someobj.initPoint = initPoint
        someobj.initDir = initDir
        someobj.surfPoint = surfPoint
        someobj.surfDir = surfDir 
        someobj.notLens = notLens

    def block_type1(this): # type 1 glass block

        x = np.linspace(-5, 5, 100)
        t = np.linspace(0, 10, 500)

        print("Vector equation of surface:", this.surfPoint,"+ t*", this.surfDir)
        plt.plot(this.surfPoint[0] + x*this.surfDir[0], this.surfPoint[1] + x*this.surfDir[1],'blue')
    
        print("Vector incident ray equation:", this.initPoint,"+ s*", this.initDir)
        plt.plot(this.initPoint[0] + t*this.initDir[0], this.initPoint[1] + t*this.initDir[1],'black')

        a = this.surfDir[1]/this.surfDir[0]
        if (initDir[1] - (a * initDir[0])) == 0:
        	print('No intercept')
        # 	if isPlot:
        # 		plt.show()
        	return initPoint, initDir
        distance = (this.surfPoint[1] - this.initPoint[1] + (a * (this.initPoint[0] - this.surfPoint[0]))) / (this.initDir[1] - (a * this.initDir[0]))

        intercept = np.array([this.initPoint[0] + distance*this.initDir[0], this.initPoint[1] + distance*this.initDir[1]])

        if this.notLens:
            print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

        normDir = nd.Gradient(func)(intercept)
        if np.dot(this.initDir, normDir) < 0:
            normDir = -normDir
        if this.notLens:
            print("Norm direction:", normDir)

        normNorm = normDir/np.linalg.norm(normDir)

        initNorm = this.initDir/np.linalg.norm(this.initDir)

        mu = this.n1/this.n2
        
        #Citation 1
        outDir = (mu*initNorm) + (normNorm*np.sqrt(1-(mu*mu*(1-((np.dot(normNorm, initNorm))**2))))) - (mu*np.dot(normNorm, np.dot(normNorm, initNorm)))


        if this.notLens:
            print("Output direction:", outDir)
            plt.plot(intercept[0] + t*outDir[0], intercept[1] + t*outDir[1],'green')        

        return intercept,outDir


#Curved surface, vectorized
class Circular_Lens:
    def __init__(someobj,n1,n2,r,initPoint,initDir,centerPoint, notLens="True"):
        someobj.n1 = n1 # refractive index of medium 1
        someobj.n2 = n2 # refractive index of medium 2
        someobj.r = r# radius of curvature
        someobj.initPoint = initPoint
        someobj.initDir = initDir
        someobj.center = centerPoint # center of circle for lens
        someobj.notLens = notLens
		
    def block_type2(this): # type 2 convex surface

        x = np.linspace(-10, 10, 1000)
        t = np.linspace(0, 10, 500)

        def func(x): 
            return ((x[0] - this.center[0])**2) + ((x[1] - this.center[1])**2) - (this.r**2)

        hyp = np.linspace(-10, 10, 1000)
        hypx, hypy = np.meshgrid(x, hyp)
        plt.contour(hypx, hypy,(((hypx - this.center[0])**2) + ((hypy - this.center[1])**2)), [this.r**2], colors='blue')
    
        if this.notLens:
            print("Vector incident ray equation:", this.initPoint,"+ s*", this.initDir)

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
            return this.initPoint, this.initDir
        elif distance_add < 0:
            distance_true = distance_min
        elif distance_min < 0:
            distance_true = distance_add
        else:
            distance_true = min(distance_add, distance_min)
	

        intercept = np.array([this.initPoint[0] + distance_true*this.initDir[0], this.initPoint[1] + distance_true*this.initDir[1]])
        if this.notLens:
            print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')
        t2 = np.linspace(0, distance_true, 100)
        plt.plot(this.initPoint[0] + t2*this.initDir[0], this.initPoint[1] + t2*this.initDir[1],'black')

        ##Normal vector is gradient of function
        normDir = nd.Gradient(func)(intercept)
        if np.dot(this.initDir, normDir) < 0:
            normDir = -normDir

        ##Direction of normal
        if this.notLens:
            print("normDir: ", normDir)
            plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')	

        
        normNorm = normDir/np.linalg.norm(normDir)
        initNorm = this.initDir/np.linalg.norm(this.initDir)
        mu = this.n1/this.n2
        #Citation 1
        outDir = (mu*initNorm) + (normNorm*np.sqrt(1-(mu*mu*(1-((np.dot(normNorm, initNorm))**2))))) - (mu*np.dot(normNorm, np.dot(normNorm, initNorm)))

        if this.notLens:
            print("Output direction:", outDir)
            plt.plot(intercept[0] + t*outDir[0], intercept[1] + t*outDir[1],'green')

        return intercept, outDir


#Convex lens using two of same curved surfaces
class Lens: 
    def __init__(someobj,n1,n2,r1,r2,h,initPoint,initDir,centerPoint):
        someobj.n1 = n1 # refractive index of medium 1
        someobj.n2 = n2 # refractive index of medium 2
        someobj.r1 = r1   # radius of first curvature
        someobj.r2 = r2
        someobj.h = h
        someobj.initPoint = initPoint
        someobj.initDir = initDir
        someobj.center = centerPoint # center of lens
		
    def block_type2(this): # type 2 convex surface
        t = np.linspace(0, 10, 500)

        for test in range(len(this.initPoint)):
            print("Initial Point: ", this.initPoint[test])
            conv1 = Circular_Lens(n1,n2,this.r1,this.initPoint[test],this.initDir,np.array([this.center[0]-(this.h/2)+this.r1,this.center[1]]), False)
            nextPoint, nextDir = conv1.block_type2()
            conv2 = Circular_Lens(n2,n1,this.r2,nextPoint,nextDir,np.array([this.center[0]+(this.h/2)-this.r2,this.center[1]]), False)
            nextPoint, nextDir = conv2.block_type2()
            print("Output: ", nextPoint, "+ s*", nextDir)
            focalLength = 1/((this.n2 - 1)*((1/this.r1) + (1/this.r2)))
            print("Focal length: ", focalLength)

            plt.plot(nextPoint[0] + t*nextDir[0], nextPoint[1] + t*nextDir[1],'green')
        plt.show()



#let the medium 1 be air, the medium 2 be glass

n1 = 1.0003
n2 = 1.52

plt.grid(color='lightgray',linestyle='--')
plt.xlim(6, 10)
plt.ylim(0, 4)
plt.gca().set_aspect('equal', adjustable='box')

# glass2 = Medium_flat_surfaces(n1,n2, np.array([-1,-1]),np.array([1,2]), np.array([0,1]), np.array([1,1]))
# refPt, refDir = glass2.block_type1()

# glass2 = Medium_flat_surfaces(n2,n1, refPt, refDir, np.array([0,5.242640687119286]),np.array([1,1])) #Fix finding y-intercept
# refPt, refDir = glass2.block_type1()

# conv1 = Circular_Lens(n1,n2,3,np.array([-2,1.5]),np.array([1,0]),np.array([3,3]))
# conv1.block_type2()

# plt.show()

# lens = Lens(n1,n2,3,3,1,np.array([-2,4]),np.array([1,0]),np.array([0,3]))
# lens.block_type2()

# lens = Lens(n1,n2,8,8,2,np.array([np.array([-2,1]),np.array([-2,1.5]),np.array([-2,2]),np.array([-2,2.5]),np.array([-2,3]),np.array([-2,3.5]),np.array([-2,4]),np.array([-2,4.5]),np.array([-2,5])]),np.array([1,0]),np.array([0,3]))
lens = Lens(n1,n2,8,8,2,np.array([np.array([-2,1]),np.array([-2,1.5]),np.array([-2,2]),np.array([-2,2.5]),np.array([-2,3])]),np.array([1,0]),np.array([0,3]))
lens.block_type2()

