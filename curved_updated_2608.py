import numpy as np 
import pandas as pd
import math
import matplotlib.pyplot as plt


class Curved_surfaces:
	def __init__(someobj,a,b,c,initialPoint,directionVec):
		someobj.a = a
		someobj.b = b
		someobj.c = c
		someobj.a1 = initialPoint[0]
		someobj.b1 = initialPoint[1]
		someobj.dx = directionVec[0]
		someobj.dy = directionVec[1]

	def tang_normal_type1(surface1): # parabola 1

		print('VERTICAL PARABOLA')

		##Variable for plotting & incident ray
		x = np.linspace(-10, 10, 1000)
		plt.plot(x, (surface1.dy/surface1.dx)*(x-surface1.a1) + surface1.b1, 'black')

		print('The coefficient of x^2', surface1.a)
		print('The coefficient of x', surface1.b)
		print('The constant term', surface1.c)


	    ##Equation of curve
		print("Equation of curve,  y =",surface1.a,"x^2 +",+ surface1.b,"x +",surface1.c)
		plt.plot(x, surface1.a*x*x + surface1.b*x + surface1.c, 'red')
		

		##Find where ray and curve intersect
		a_term = surface1.a * surface1.dx * surface1.dx
		b_term = (2 * surface1.a * surface1.a1 * surface1.dx) + (surface1.b * surface1.dx) - surface1.dy
		c_term = (surface1.a * surface1.a1 * surface1.a1) + (surface1.b * surface1.a1) + surface1.c - surface1.b1
		distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
		distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
		distance_true = min(abs(distance_add), abs(distance_min))

		intercept = np.array([surface1.a1 + distance_true*surface1.dx, surface1.b1 + distance_true*surface1.dy])

		print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

		surface1.x = intercept[0]
		surface1.y = intercept[1]
		

		##Differentiate above equation to find tangent dy/dx
		tangDir = np.array([1, 2*surface1.a*surface1.x + surface1.b])

		diff = 2*surface1.a*surface1.x + surface1.b
		c1 = surface1.y - diff*surface1.x

	    ##Equation of tangent
		print("Equations of tangent at given point, y =",diff,"x +",c1)
		print("tangDir: ", tangDir)
		plt.plot(surface1.x + x*tangDir[0], surface1.y + x*tangDir[1], 'blue')


		##Slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		normDir = np.array([tangDir[1], -tangDir[0]])

		diff_norm = -1/(diff)
		c2 = surface1.y - diff_norm*surface1.x

		#equation of normal
		print("Equations of normal at given point, y =",diff_norm,"x +",c2)
		print("normDir: ", normDir)
		plt.plot(surface1.x + x*normDir[0], surface1.y + x*normDir[1], 'orange')


		##Equation of incident line using two points formula (Should match our input):
		print("Equation of incident ray, y =",((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x)),"x +",(surface1.y - (((surface1.b1 - surface1.y)*surface1.x)/(surface1.a1 - surface1.x))))


		##Finding equation of reflected line by projecting onto norm and subtracting vectors
		rayDir = np.array([surface1.dx, surface1.dy])
		normNorm = normDir/np.linalg.norm(normDir)
		out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

		##Finding equation of reflected line using tan_theta
		diff_line = ((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x))

		tan_theta = (diff_line - diff_norm) / (1 + (diff_line * diff_norm))

		if tan_theta < 0:
			tan_theta = -tan_theta

		slope_reflected = (diff_norm - tan_theta) / (1 + (diff_norm * tan_theta))

		c3 = surface1.y - slope_reflected*surface1.x

		##Equation of reflected line
		print("Equation of reflected ray using tan_theta, y =",slope_reflected,"x +",c3)
		plt.plot(surface1.x + x*out[0], surface1.y + x*out[1],'green')
		

		##Using second set of equations, with mirroring technique:
		x_t = ((c1 - surface1.b1 + surface1.a1*diff_norm)/(diff_norm-diff)) 
		y_t = x_t*diff + c1

		x_2t = 2*x_t - surface1.a1
		y_2t = 2*y_t - surface1.b1

		c_ref = y_2t - (((surface1.y - y_2t)*x_2t)/(surface1.x - x_2t))

		print("Equation of reflected ray without using tan_theta , y =",((surface1.y - y_2t)/(surface1.x - x_2t)),"x +",c_ref)
		#plt.plot(x,((surface1.y - y_2t)/(surface1.x - x_2t))*x +c_ref,'black')


		##Print plot
		plt.grid(color='lightgray',linestyle='--')
		plt.xlim(-10, 10)
		plt.ylim(-10, 10)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()
	
	#tang_normal(-1,0,1,1,2,0,2)
	def tang_normal_type2(surface1):  # parabola 2

		print('HORIZONTAL PARABOLA')

		##Variable for plotting & incident ray
		x = np.linspace(-10, 100, 1000)
		plt.plot(x, (surface1.dy/surface1.dx)*(x-surface1.a1) + surface1.b1, 'black')
		print("Ray Direction: ", surface1.dx, surface1.dy)

		print('The coefficient of y^2', surface1.a)
		print('The coefficient of y', surface1.b)
		print('The constant term', surface1.c)


	    ##Equation of curve
		print("Equation of curve,  x =",surface1.a,"y^2 +",+ surface1.b,"y +",surface1.c)
		plt.plot(surface1.a*x*x + surface1.b*x + surface1.c, x, 'red')
		

		##Find where ray and curve intersect
		a_term = surface1.a * surface1.dy * surface1.dy
		b_term = (2 * surface1.a * surface1.b1 * surface1.dy) + (surface1.b * surface1.dy) - surface1.dx
		c_term = (surface1.a * surface1.b1 * surface1.b1) + (surface1.b * surface1.b1) + surface1.c - surface1.a1
		distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
		distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
		distance_true = min(abs(distance_add), abs(distance_min))

		intercept = np.array([surface1.a1 + distance_true*surface1.dx, surface1.b1 + distance_true*surface1.dy])

		print('Point of intersection 1: (', intercept[0], ', ', intercept[1], ')')

		surface1.x = intercept[0]
		surface1.y = intercept[1]
			

		##Differentiate above equation to find tangent dy/dx
		tangDir = np.array([2*surface1.a*surface1.y + surface1.b, 1])
		if tangDir[0] == 0:
			diff = np.inf

			##Equation of tangent
			print("Equations of tangent at given point, x =",surface1.x)
		else:
			diff = 1/(2*surface1.a*surface1.y + surface1.b)
			c1 = surface1.y - diff*surface1.x

			##Equation of tangent
			print("Equations of tangent at given point, y =",diff,"x +",c1)
		print("tangDir: ", tangDir)
		plt.plot(surface1.x + x*tangDir[0], surface1.y + x*tangDir[1], 'blue')


		##Slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		normDir = np.array([tangDir[1], -tangDir[0]])
		diff_norm = (-1 * tangDir[0])/tangDir[1]
		c2 = surface1.y - diff_norm*surface1.x

		##Equation of normal
		print("Equations of normal at given point, y =",diff_norm,"x +",c2)
		print("normDir: ", normDir)
		plt.plot(surface1.x + x*normDir[0], surface1.y + x*normDir[1], 'orange')


		##Equation of incident line using two points formula (Should match our input):
		print("Equation of incident ray, y =",((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x)),"x +",(surface1.y - (((surface1.b1 - surface1.y)*surface1.x)/(surface1.a1 - surface1.x))))


		##Finding equation of reflected line by projecting onto norm and subtracting vectors
		rayDir = np.array([surface1.dx, surface1.dy])
		normNorm = normDir/np.linalg.norm(normDir)
		out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

		##Finding equation of reflected line using tan_theta
		diff_line = ((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x))	

		tan_theta = (diff_line - diff_norm) / (1 + (diff_line * diff_norm))

		if tan_theta < 0:
			tan_theta = -tan_theta

		slope_reflected = (diff_norm - tan_theta) / (1 + (diff_norm * tan_theta))

		c3 = surface1.y - slope_reflected*surface1.x

		##Equation of reflected line
		print("Equation of reflected ray, y =",slope_reflected,"x +",c3)
		plt.plot(surface1.x + x*out[0], surface1.y + x*out[1],'green')


		##Using second set of equations, with mirroring technique:
		# x_t = ((c1 - surface1.b1 + surface1.a1*diff_norm)/(diff_norm-diff)) 
		# y_t = x_t*diff + c1

		# x_2t = 2*x_t - surface1.a1
		# y_2t = 2*y_t - surface1.b1

		# c_ref = y_2t - (((surface1.y - y_2t)*x_2t)/(surface1.x - x_2t))

		# print("Equation of reflected ray without using tan_theta , y =",((surface1.y - y_2t)/(surface1.x - x_2t)),"x +",c_ref)

		##Print plot
		plt.grid(color='lightgray',linestyle='--')
		plt.xlim(-10, 10)
		plt.ylim(-10, 10)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()



class Plane_surfaces:

	def __init__(someobj,a,b,initialPoint,directionVec):
		someobj.a = a
		someobj.b = b
		someobj.a1 = initialPoint[0]
		someobj.b1 = initialPoint[1]
		someobj.dx = directionVec[0]
		someobj.dy = directionVec[1]

	def tang_normal_type3(surface1):

		print('PLANE')

		##Variable for plotting & incident ray
		x = np.linspace(-10, 10, 1000)
		plt.plot(x, (surface1.dy/surface1.dx)*(x-surface1.a1) + surface1.b1, 'black')

		print('The coefficient of x',surface1.a)
		print('The constant term',surface1.b)


	    ##Equation of curve
		print("Equation of line,  y =",surface1.a,"x +",surface1.b)
		#plt.plot(x, surface1.a*x + surface1.b, 'red')


		##Find where ray and curve intersect
		distance = ((surface1.a * surface1.a1) + surface1.b - surface1.b1) / (surface1.dy - (surface1.a * surface1.dx))

		intercept = np.array([surface1.a1 + distance_true*surface1.dx, surface1.b1 + distance_true*surface1.dy])

		print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

		surface1.x = intercept[0]
		surface1.y = intercept[1]
		

		##Differentiate above equation to find tangent dy/dx
		tangDir = np.array([1, surface1.a])

		diff = surface1.a
		c1 = surface1.y - diff*surface1.x

	    ##Equation of tangent
		print("Equations of tangent at given point, y =",diff,"x +",c1)
		#plt.plot(x, diff*x + c1, 'blue')
		print("tangDir: ", tangDir)
		plt.plot(surface1.x + x*tangDir[0], surface1.y + x*tangDir[1], 'blue')


		##Slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		normDir = np.array([tangDir[1], -tangDir[0]])
		diff_norm = -1/(diff)
		c2 = surface1.y - diff_norm*surface1.x

		##Equation of normal
		print("Equations of normal at given point, y =",diff_norm,"x +",c2)
		#plt.plot(x, diff_norm*x + c2, 'orange')
		print("normDir: ", normDir)
		plt.plot(surface1.x + x*normDir[0], surface1.y + x*normDir[1], 'orange')


		##Equation of incident line using two points formula (Should match our input):
		print("Equation of incident ray, y =",((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x)),"x +",(surface1.y - (((surface1.b1 - surface1.y)*surface1.x)/(surface1.a1 - surface1.x))))


		##Finding equation of reflected line by projecting onto norm and subtracting vectors
		rayDir = np.array([surface1.dx, surface1.dy])
		normNorm = normDir/np.linalg.norm(normDir)
		out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

		##Finding equation of reflected line using tan_theta
		diff_line = ((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x))

		tan_theta = ((diff_line - diff_norm)/(1+(diff_line)*(diff_norm)))

		if tan_theta < 0:
			tan_theta = -tan_theta

		slope_reflected = (diff_norm - tan_theta) / (1 + (diff_norm * tan_theta))

		c3 = surface1.y - slope_reflected*surface1.x

		##Equation of reflected line
		print("Equation of reflected ray using tan theta, y =",slope_reflected,"x +",c3)
		#plt.plot(x,slope_reflected*x + c3,'green')
		plt.plot(surface1.x + x*out[0], surface1.y + x*out[1],'green')


		##Using second set of equations, with mirroring technique:
		x_t = ((c1 - surface1.b1 + surface1.a1*diff_norm)/(diff_norm-diff)) 
		y_t = x_t*diff + c1

		x_2t = 2*x_t - surface1.a1
		y_2t = 2*y_t - surface1.b1

		c_ref = y_2t - (((surface1.y - y_2t)*x_2t)/(surface1.x - x_2t))

		print("Equation of reflected ray without using tan_theta , y =",((surface1.y - y_2t)/(surface1.x - x_2t)),"x +",c_ref)

		##Print plot
		plt.grid(color='lightgray',linestyle='--')
		plt.xlim(-10, 10)
		plt.ylim(-10, 10)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()


	
class Hyperbola_curve:

	def __init__(someobj,a,b,h,k,initialPoint,directionVec):
		someobj.a = a
		someobj.b = b
		someobj.h = h
		someobj.k = k
		someobj.a1 = initialPoint[0]
		someobj.b1 = initialPoint[1]
		someobj.dx = directionVec[0]
		someobj.dy = directionVec[1]

	def tang_normal_type4(surface1): # parabola 1

		print('HYPERBOLA')

		##Variable for plotting & incident ray
		x = np.linspace(-10, 10, 1000)
		hyp = np.linspace(-10, 10, 1000)
		plt.plot(x, (surface1.dy/surface1.dx)*(x-surface1.a1) + surface1.b1, 'black')

		##Variables for calculations
		x2coeff = surface1.b*surface1.b
		y2coeff = -surface1.a*surface1.a
		xcoeff = -2*surface1.b*surface1.b*surface1.h
		ycoeff = 2*surface1.a*surface1.a*surface1.k
		constant = surface1.b*surface1.b*surface1.h*surface1.h - surface1.a*surface1.a*surface1.k*surface1.k - surface1.a*surface1.a*surface1.b*surface1.b

		print('The coefficient of x^2', x2coeff)
		print('The coefficient of y^2', y2coeff)
		print('The coefficient of x', xcoeff)
		print('The coefficient of y', ycoeff)
		print('The constant term', constant)


	    ##Equation of curve
		print("Equation of curve,  0 =",x2coeff,"x^2 +",y2coeff,"y^2 +", xcoeff,"x +",ycoeff,"y +",constant)
		#y = someothing...
		#plt.plot(x, y, 'red')		
		hypx, hypy = np.meshgrid(x, hyp)
		plt.contour(hypx, hypy,(((hypx - surface1.h)**2)/surface1.a**2 - ((hypy - surface1.k)**2)/surface1.b**2), [1], colors='red')

		##Find where ray and curve intersect
		a_term = surface1.b * surface1.b * surface1.dx * surface1.dx - surface1.a * surface1.a * surface1.dy * surface1.dy
		b_term = 2 * (surface1.b * surface1.b * surface1.a1 * surface1.dx - surface1.b * surface1.b * surface1.h * surface1.dx 
			- surface1.a * surface1.a * surface1.b1 * surface1.dy + surface1.a * surface1.a * surface1.k * surface1.dy)
		c_term = (surface1.b * surface1.b * (surface1.a1 - surface1.h) * (surface1.a1 - surface1.h) 
			- surface1.a * surface1.a * (surface1.b1 - surface1.k) * (surface1.b1 - surface1.k) 
			- surface1.a * surface1.a * surface1.b * surface1.b)
		distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
		distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
		distance_true = min(abs(distance_add), abs(distance_min))

		intercept = np.array([surface1.a1 + distance_true*surface1.dx, surface1.b1 + distance_true*surface1.dy])

		print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

		surface1.x = intercept[0]
		surface1.y = intercept[1]


		##Differentiate above equation to find tangent dy/dx
		tangDir = np.array([surface1.a*surface1.a*surface1.y - surface1.a*surface1.a*surface1.k, surface1.x*surface1.b*surface1.b - surface1.h*surface1.b*surface1.b])
		diff = (surface1.x*surface1.b*surface1.b - surface1.h*surface1.b*surface1.b)/(surface1.a*surface1.a*surface1.y - surface1.a*surface1.a*surface1.k)
		c1 = surface1.y - diff*surface1.x

		##Equation of tangent
		print("Equations of tangent at given point, y =",diff,"x +",c1)
		# plt.plot(x, diff*x + c1, 'blue')
		print("tangDir: ", tangDir)
		plt.plot(surface1.x + x*tangDir[0], surface1.y + x*tangDir[1], 'blue')


		##Slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		normDir = np.array([tangDir[1], -tangDir[0]])
		diff_norm = -1/(diff)
		c2 = surface1.y - diff_norm*surface1.x

		##Equation of normal
		print("Equations of normal at given point, y =",diff_norm,"x +",c2)
		# plt.plot(x, diff_norm*x + c2, 'orange')
		print("normDir: ", normDir)
		plt.plot(surface1.x + x*normDir[0], surface1.y + x*normDir[1], 'orange')


		##Equation of incident line using two points formula (Should match our input):
		print("Equation of incident ray, y =",((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x)),"x +",(surface1.y - (((surface1.b1 - surface1.y)*surface1.x)/(surface1.a1 - surface1.x))))
		

		##Finding equation of reflected line by projecting onto norm and subtracting vectors
		rayDir = np.array([surface1.dx, surface1.dy])
		normNorm = normDir/np.linalg.norm(normDir)
		out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

		##Finding equation of reflected line using tan_theta
		diff_line = (surface1.dy/surface1.dx)

		tan_theta = (diff_line - diff_norm) / (1 + (diff_line * diff_norm))

		if tan_theta < 0:
			tan_theta = -tan_theta

		slope_reflected = (diff_norm - tan_theta) / (1 + (diff_norm * tan_theta))

		c3 = surface1.y - slope_reflected*surface1.x

		##Equation of reflected line
		print("Equation of reflected ray using tan_theta, y =",slope_reflected,"x +",c3)
		# plt.plot(x,slope_reflected*x + c3,'green')
		plt.plot(surface1.x + x*out[0], surface1.y + x*out[1],'green')

		
		##Using second set of equations, with mirroring technique:
		x_t = ((c1 - surface1.b1 + surface1.a1*diff_norm)/(diff_norm-diff)) 
		y_t = x_t*diff + c1

		x_2t = 2*x_t - surface1.a1
		y_2t = 2*y_t - surface1.b1

		c_ref = y_2t - (((surface1.y - y_2t)*x_2t)/(surface1.x - x_2t))

		print("Equation of reflected ray without using tan_theta , y =",((surface1.y - y_2t)/(surface1.x - x_2t)),"x +",c_ref)
		#plt.plot(x,((surface1.y - y_2t)/(surface1.x - x_2t))*x +c_ref,'black')

		##Print plot
		plt.grid(color='lightgray',linestyle='--')
		plt.xlim(-10, 10)
		plt.ylim(-10, 10)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()



class Ellipse:

	def __init__(someobj,a,b,h,k,initialPoint,directionVec):
		someobj.a = a
		someobj.b = b
		someobj.h = h
		someobj.k = k
		someobj.a1 = initialPoint[0]
		someobj.b1 = initialPoint[1]
		someobj.dx = directionVec[0]
		someobj.dy = directionVec[1]

	def tang_normal_type5(surface1):

		print('ELLIPSE')

		##Variables for plotting & incident ray
		x = np.linspace(-5, 5, 100)
		t = np.linspace(0, 2*math.pi, 100)
		plt.plot(x, (surface1.dy/surface1.dx)*(x-surface1.a1) + surface1.b1, 'black')

		##Variables for calculations
		x2coeff = surface1.b*surface1.b
		y2coeff = surface1.a*surface1.a
		xcoeff = -2*surface1.b*surface1.b*surface1.h
		ycoeff = -2*surface1.a*surface1.a*surface1.k
		constant = surface1.b*surface1.b*surface1.h*surface1.h + surface1.a*surface1.a*surface1.k*surface1.k - surface1.a*surface1.a*surface1.b*surface1.b

		print('The coefficient of x^2: ', x2coeff)
		print('The coefficient of y^2: ', y2coeff)
		print('The coefficient of x: ', xcoeff)
		print('The coefficient of y: ', ycoeff)
		print('The constant term: ', constant)


	    ##Equation of curve
		print("Equation of curve,  0 =",x2coeff,"x^2 +", y2coeff,"y^2 +", xcoeff,"x +", ycoeff,"y +", constant)
		plt.plot( surface1.h+surface1.a*np.cos(t) , surface1.k+surface1.b*np.sin(t) )


		##Find where ray and curve intersect
		# vec = np.array([(surface1.a1 - surface1.x), (surface1.b1 - surface1.y)])
		a_term = surface1.b * surface1.b * surface1.dx * surface1.dx + surface1.a * surface1.a * surface1.dy * surface1.dy
		b_term = 2 * (surface1.b * surface1.b * surface1.a1 * surface1.dx - surface1.b * surface1.b * surface1.h * surface1.dx 
			+ surface1.a * surface1.a * surface1.b1 * surface1.dy - surface1.a * surface1.a * surface1.k * surface1.dy)
		c_term = (surface1.b * surface1.b * (surface1.a1 - surface1.h) * (surface1.a1 - surface1.h) 
			+ surface1.a * surface1.a * (surface1.b1 - surface1.k) * (surface1.b1 - surface1.k) 
			- surface1.a * surface1.a * surface1.b * surface1.b)
		distance_add = (-b_term + np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
		distance_min = (-b_term - np.sqrt(b_term*b_term - 4 * a_term * c_term)) / (2 * a_term)
		distance_true = min(abs(distance_add), abs(distance_min))

		intercept = np.array([surface1.a1 + distance_true*surface1.dx, surface1.b1 + distance_true*surface1.dy])

		print('Point of intersection: (', intercept[0], ', ', intercept[1], ')')

		surface1.x = intercept[0]
		surface1.y = intercept[1]


		##Differentiate above equation to find tangent dy/dx
		tangDir = np.array([surface1.a*surface1.a*(surface1.y - surface1.k), -(surface1.b * surface1.b * (surface1.x - surface1.h))])

		diff = -(surface1.b * surface1.b * (surface1.x - surface1.h)) / (surface1.a*surface1.a*(surface1.y - surface1.k))
		c1 = surface1.y - diff*surface1.x

		##Equation of tangent
		print("Equations of tangent at given point, y =",diff,"x +",c1)
		# plt.plot(x, diff*x + c1, 'blue')
		print("tangDir: ", tangDir)
		plt.plot(surface1.x + x*tangDir[0], surface1.y + x*tangDir[1], 'blue')


		##Slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		normDir = np.array([tangDir[1], -tangDir[0]])
		diff_norm = -1/(diff)
		c2 = surface1.y - diff_norm*surface1.x

		##Equation of normal
		print("Equations of normal at given point, y =",diff_norm,"x +",c2)
		# plt.plot(x, diff_norm*x + c2, 'orange')
		print("normDir: ", normDir)
		plt.plot(surface1.x + x*normDir[0], surface1.y + x*normDir[1], 'orange')


		##Equation of incident line using two points formula (Should match our input):
		print("Equation of incident ray, y =",((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x)),"x +", (surface1.y - surface1.x * ((surface1.b1 - surface1.y) / (surface1.a1 - surface1.x))))		


		##Finding equation of reflected line by projecting onto norm and subtracting vectors
		rayDir = np.array([surface1.dx, surface1.dy])
		normNorm = normDir/np.linalg.norm(normDir)
		out = rayDir - 2*(np.dot(rayDir, normNorm)*normNorm)

		##Finding equation of reflected line using tan_theta
		diff_line = (surface1.dy/surface1.dx)

		tan_theta = (diff_line - diff_norm) / (1 + (diff_line * diff_norm))

		if tan_theta < 0:
			tan_theta = -tan_theta

		slope_reflected = (diff_norm - tan_theta) / (1 + (diff_norm * tan_theta))

		c3 = surface1.y - slope_reflected*surface1.x

		##Equation of reflected line
		print("Equation of reflected ray using tan_theta, y =",slope_reflected,"x +",c3)
		# plt.plot(x,slope_reflected*x + c3,'green')
		plt.plot(surface1.x + x*out[0], surface1.y + x*out[1],'green')
		

		##Using second set of equations, with mirroring technique:
		x_t = ((c1 - surface1.b1 + surface1.a1*diff_norm)/(diff_norm-diff)) 
		y_t = x_t*diff + c1

		x_2t = 2*x_t - surface1.a1
		y_2t = 2*y_t - surface1.b1

		c_ref = y_2t - (((surface1.y - y_2t)*x_2t)/(surface1.x - x_2t))

		print("Equation of reflected ray without using tan_theta , y =",((surface1.y - y_2t)/(surface1.x - x_2t)),"x +",c_ref)
		##plt.plot(x,((surface1.y - y_2t)/(surface1.x - x_2t))*x +c_ref,'black')

		##Testing other method
		c4 = surface1.b1 - (diff * surface1.a1)

		x_norm = (c2 - c4) / (diff - diff_norm)
		y_norm = (diff_norm * x_norm) + c2

		x_out = x_norm - (surface1.a1 - x_norm)
		y_out = y_norm - (surface1.b1 - y_norm)

		slope = (y_out - surface1.y) / (x_out - surface1.x)
		intercept = (-slope * x_out) + y_out

		print("Equation of reflected ray without using tan_theta 2, y =",slope,"x +",intercept)

		##Print plot
		plt.grid(color='lightgray',linestyle='--')
		plt.xlim(-10, 10)
		plt.ylim(-10, 10)
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()




# curv1 = Curved_surfaces(1,0,1,np.array([0,2]),np.array([1,2])) #a,b,c,a1,b1,dx,dy
# curv1.tang_normal_type1()	


# curv2 = Curved_surfaces(1,2,4,np.array([-6,-2]),np.array([9,1])) 
# curv2.tang_normal_type2()

# plane1 = Plane_surfaces(1,0,np.array([4,3]),np.array([2,2])) #a,b,a1,b1,dx,dy
# plane1.tang_normal_type3()

curv4 = Hyperbola_curve(1,2,1,3,np.array([2,2]),np.array([1,4])) #a,b,h,k,a1,b1,dx,dy
curv4.tang_normal_type4()

# curv5 = Ellipse(1,2,1,3, np.array([0,0]),np.array([1,4])) 
# curv5.tang_normal_type5()
