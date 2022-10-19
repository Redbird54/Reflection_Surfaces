import numpy as np 
import pandas as pd
import math
import matplotlib.pyplot as plt


class Curved_Surface:
	def __init__(someobj,a,b,c,initialPoint,directionVec):
		someobj.a = a
		someobj.b = b
		someobj.c = c
		someobj.initPoint = initialPoint
		someobj.initDir = directionVec

	def tang_normal_type1(this, ray_only, simpleView, isPlot): # parabola 1

		print('VERTICAL PARABOLA')

		a = this.a
		b = this.b
		c = this.c
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
		print('The coefficient of x^2', a)
		print('The coefficient of x', b)
		print('The constant term', c)


	    ##Equation of curve
		print("Equation of curve,  y =", a,"x^2 +", b,"x +", c)
		plt.plot(x, a*x*x + b*x + c, 'red')
		

		##Find where ray and curve intersect
		a_term = a * initDir[0] * initDir[0]
		b_term = (2 * a * initPoint[0] * initDir[0]) + (b * initDir[0]) - initDir[1]
		c_term = (a * initPoint[0] * initPoint[0]) + (b * initPoint[0]) + c - initPoint[1]
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
		tangDir = np.array([1, 2*a*intercept[0] + b])

	    ##Equation of tangent
		print("tangDir: ", tangDir)
		if not(ray_only) and not(simpleView):
			plt.plot(intercept[0] + x*tangDir[0], intercept[1] + x*tangDir[1], 'blue')


		##Slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		normDir = np.array([tangDir[1], -tangDir[0]])

		##Direction of normal
		print("normDir: ", normDir)
		if not(ray_only) and not(simpleView):
			plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')


		##Finding equation of reflected line by projecting onto norm and subtracting vectors
		rayDir = np.array([initDir[0], initDir[1]])
		normNorm = normDir/np.linalg.norm(normDir)
		out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

		##Direction of reflected line
		print("Output Direction: ", out)
		if not(ray_only):
			plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


		##Print plot
		if isPlot:
			plt.show()

		return intercept, out
	
	def tang_normal_type2(this, ray_only, simpleView, isPlot):  # parabola 2

		print('HORIZONTAL PARABOLA')

		a = this.a
		b = this.b
		c = this.c
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
		print('The coefficient of y^2', a)
		print('The coefficient of y', b)
		print('The constant term', c)


	    ##Equation of curve
		print("Equation of curve,  x =", a,"y^2 +", b,"y +", c)
		plt.plot(a*x*x + b*x + c, x, 'red')
		

		##Find where ray and curve intersect
		a_term = a * initDir[1] * initDir[1]
		b_term = (2 * a * initPoint[1] * initDir[1]) + (b * initDir[1]) - initDir[0]
		c_term = (a * initPoint[1] * initPoint[1]) + (b * initPoint[1]) + c - initPoint[0]
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
		tangDir = np.array([2*a*intercept[1] + b, 1])

		##Direction of tangent
		print("tangDir: ", tangDir)
		if not(ray_only) and not(simpleView):
			plt.plot(intercept[0] + x*tangDir[0], intercept[1] + x*tangDir[1], 'blue')


		##Slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		normDir = np.array([tangDir[1], -tangDir[0]])
		
		##Direction of normal
		print("normDir: ", normDir)
		if not(ray_only) and not(simpleView):
			plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')


		##Finding equation of reflected line by projecting onto norm and subtracting vectors
		rayDir = np.array([initDir[0], initDir[1]])
		normNorm = normDir/np.linalg.norm(normDir)
		out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

		##Direction of reflected line
		print("Output Direction: ", out)
		if not(ray_only):
			plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


		##Print plot
		if isPlot:
			plt.show()

		return intercept, out



class Plane_Surface:

	def __init__(someobj,a,b,initialPoint,directionVec):
		someobj.a = a
		someobj.b = b
		someobj.initPoint = initialPoint
		someobj.initDir = directionVec

	def tang_normal_type3(this, ray_only, simpleView, isPlot):

		print('PLANE')

		a = this.a
		b = this.b
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

		intercept = np.array([initPoint[0] + distance*initDir[0], initPoint[1] + distance*initDir[1]])

		print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

		##Plot where the ray actually goes (initial point until intercept)
		if ray_only:
			t2 = np.linspace(0, distance, 100)
			plt.plot(initPoint[0] + t2*initDir[0], initPoint[1] + t2*initDir[1],'black')
		

		##Differentiate above equation to find tangent dy/dx
		tangDir = np.array([1, a])

		##Direction of tangent
		print("tangDir: ", tangDir)
		if not(ray_only) and not(simpleView):
			plt.plot(intercept[0] + x*tangDir[0], intercept[1] + x*tangDir[1], 'blue')


		##Slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		normDir = np.array([tangDir[1], -tangDir[0]])
		
		##Direction of normal
		print("normDir: ", normDir)
		if not(ray_only) and not(simpleView):
			plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')


		##Finding equation of reflected line by projecting onto norm and subtracting vectors
		rayDir = np.array([initDir[0], initDir[1]])
		normNorm = normDir/np.linalg.norm(normDir)
		out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

		##Direction of reflected line
		print("Output Direction: ", out)
		if not(ray_only):
			plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')


		##Print plot
		if isPlot:
			plt.show()

		return intercept, out


	
class Hyperbola_Curve:

	def __init__(someobj,a,b,h,k,initialPoint,directionVec):
		someobj.a = a
		someobj.b = b
		someobj.h = h
		someobj.k = k
		someobj.initPoint = initialPoint
		someobj.initDir = directionVec

	def tang_normal_type4(this, ray_only, simpleView, isPlot): # parabola 1

		print('HYPERBOLA')

		a = this.a
		b = this.b
		h = this.h
		k = this.k
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
		a_term = b * b * initDir[0] * initDir[0] - a * a * initDir[1] * initDir[1]
		b_term = 2 * (b * b * initPoint[0] * initDir[0] - b * b * h * initDir[0] 
			- a * a * initPoint[1] * initDir[1] + a * a * k * initDir[1])
		c_term = (b * b * (initPoint[0] - h) * (initPoint[0] - h) 
			- a * a * (initPoint[1] - k) * (initPoint[1] - k) 
			- a * a * b * b)
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
		tangDir = np.array([a*a*intercept[1] - a*a*k, intercept[0]*b*b - h*b*b])
		
		##Direction of tangent
		print("tangDir: ", tangDir)
		if not(ray_only) and not(simpleView):
			plt.plot(intercept[0] + x*tangDir[0], intercept[1] + x*tangDir[1], 'blue')


		##Slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		normDir = np.array([tangDir[1], -tangDir[0]])
	
		##Direction of normal
		print("normDir: ", normDir)
		if not(ray_only) and not(simpleView):
			plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')
		

		##Finding equation of reflected line by projecting onto norm and subtracting vectors
		rayDir = np.array([initDir[0], initDir[1]])
		normNorm = normDir/np.linalg.norm(normDir)
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

	def tang_normal_type5(this, ray_only, simpleView, isPlot):

		print('ELLIPSE')

		a = this.a
		b = this.b
		h = this.h
		k = this.k
		initPoint = this.initPoint
		initDir = this.initDir

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


		##Find where ray and curve intersect
		a_term = b * b * initDir[0] * initDir[0] + a * a * initDir[1] * initDir[1]
		b_term = 2 * (b * b * initPoint[0] * initDir[0] - b * b * h * initDir[0] 
			+ a * a * initPoint[1] * initDir[1] - a * a * k * initDir[1])
		c_term = (b * b * (initPoint[0] - h) * (initPoint[0] - h) 
			+ a * a * (initPoint[1] - k) * (initPoint[1] - k) 
			- a * a * b * b)
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
		tangDir = np.array([a*a*(intercept[1] - k), -(b * b * (intercept[0] - h))])

		##Direction of tangent
		print("tangDir: ", tangDir)
		if not(ray_only) and not(simpleView):
			plt.plot(intercept[0] + x*tangDir[0], intercept[1] + x*tangDir[1], 'blue')


		##Slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		normDir = np.array([tangDir[1], -tangDir[0]])

		##Direction of normal
		print("normDir: ", normDir)
		if not(ray_only) and not(simpleView):
			plt.plot(intercept[0] + x*normDir[0], intercept[1] + x*normDir[1], 'orange')	


		##Finding equation of reflected line by projecting onto norm and subtracting vectors
		rayDir = np.array([initDir[0], initDir[1]])
		normNorm = normDir/np.linalg.norm(normDir)
		out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

		##Direction of reflected line
		print("Output Direction: ", out)
		if not(ray_only):
			plt.plot(intercept[0] + t*out[0], intercept[1] + t*out[1],'green')
		

		##Print plot
		if isPlot:
			plt.show()

		return intercept, out


##Settings for how to show plot(s)
raysOnly = True
simpleView = False
indivPlots = False

# curv1 = Curved_Surface(1,0,1,np.array([0,2]),np.array([1,2])) #a,b,c,initPoint,initDir
# nextPoint,nextDir = curv1.tang_normal_type1(raysOnly, simpleView, indivPlots)	

# curv2 = Curved_Surface(1,2,4,np.array([12,0]),np.array([-9,1])) 
# nextPoint,nextDir = curv2.tang_normal_type2(raysOnly, simpleView, indivPlots)

# curv3 = Plane_Surface(1,0,np.array([4,3]),np.array([1,2])) #a,b,initPoint,initDir
# nextPoint,nextDir = curv3.tang_normal_type3(raysOnly, simpleView, indivPlots)

# curv4 = Hyperbola_Curve(1,2,1,3,np.array([2,2]),np.array([1,4])) #a,b,h,k,initPoint,initDir
# nextPoint,nextDir = curv4.tang_normal_type4(raysOnly, simpleView, indivPlots)

curv5 = Ellipse(1,2,1,3, np.array([0,0]),np.array([1,4])) 
nextPoint,nextDir = curv5.tang_normal_type5(raysOnly, simpleView, indivPlots)

curv4 = Hyperbola_Curve(3,2,1,3,nextPoint,nextDir)
nextPoint,nextDir = curv4.tang_normal_type4(raysOnly, simpleView, indivPlots)


if not(indivPlots):
	t = np.linspace(0, 10, 500)
	plt.plot(nextPoint[0] + t*nextDir[0], nextPoint[1] + t*nextDir[1],'green')

	##Print plot
	plt.grid(color='lightgray',linestyle='--')
	plt.xlim(-10, 10)
	plt.ylim(-10, 10)
	plt.gca().set_aspect('equal', adjustable='box')
	plt.show()