import numpy as np 
import pandas as pd
import math
import matplotlib.pyplot as plt


class Curved_surfaces:
	def __init__(someobj,a,b,c,x,y,a1,b1):
		someobj.a = a
		someobj.b = b
		someobj.c = c
		someobj.x = x
		someobj.y = y
		someobj.a1 = a1
		someobj.b1 = b1

	def tang_normal_type1(surface1): # parabola 1

		print('The coefficient of x^2',+surface1.a)
		print('The coefficient of x', +surface1.b)
		print('The constant term', +surface1.c)


	    #equation of curve
		print("Equation of curve,  y =",surface1.a,"x^2 +",+ surface1.b,"x +", +surface1.c)
		

		x = np.linspace(-10, 10, 1000)
		#x = np.arange(-100,100,1)
		#x = np.linspace(0,2*math.pi,400)
		y = surface1.a*x*x+ surface1.b*x+surface1.c
		plt.plot(x, y, 'r')

		#surface1.y = surface1.a*surface1.x*surface1.x + surface1.b*surface1.x + surface1.c
		#plot(surface1.x,surface1.y)

		

	#differentiate above equation to find tangent dy/dx

		diff = 2*surface1.a*surface1.x + surface1.b
		c1 = surface1.y - diff*surface1.x

	    #equation of tangent
		print("Equations of tangent at given point, y =",diff,"x +",+ c1)
		y = diff*x + c1
		plt.plot(x, y, 'b')



		#y = diff*x + c1
	    #plt.plot(x,y)
	    # Show the plot
		#plt.show()

		#the slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		diff_norm = -1/(diff)
		c2 = surface1.y - diff_norm*surface1.x

		#equation of normal
		print("Equations of normal at given point, y =",diff_norm,"x +",+ c2)
		y = diff_norm*x + c2
		plt.plot(x, y, 'orange')



		#y = diff*x + c1
	    #plt.plot(x,y)
	    # Show the plot
		#plt.show()

		# set the incident light ray
		#y = a1x + b1

		#a1=0
		#b1=2
		### equation of incident line using two points formula :
		print("Equation of incident ray, y =",((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x)),"x +",+ (surface1.y - (((surface1.b1 - surface1.y)*surface1.x)/(surface1.a1 - surface1.x))))

		#equation of incident line
		#print("Equation of incident ray, y =",surface1.a1,"x +",+ surface1.b1)
		

		y = ((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x))*x + (surface1.y - (((surface1.b1 - surface1.y)*surface1.x)/(surface1.a1 - surface1.x)))
		#y = diff*x + c1
	    
	    # Show the plot
	    #plt.plot(x, y, 'green')
		plt.plot(x,y,'green')
	    # Show the plot
		#plt.show()

		#plt.show()

		diff_line = ((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x))

		tan_theta = ((diff_line - diff_norm)/(1+(diff_line)*(diff_norm)))

		#print(tan_theta)

		slope_reflected = ((-tan_theta + diff_norm)/(1+(diff_norm)*(tan_theta)))


		#print(slope_reflected)

		c3 = surface1.y - slope_reflected*surface1.x

		#equation of reflected line
		print("Equation of reflected ray using tan_theta, y =",slope_reflected,"x +",+ c3)

		y = slope_reflected*x + c3
		plt.plot(x,y,'yellow')
		#plt.show()
		
		# To load the display window
		#using second set of equations, with mirroring technique:

		x_t = ((c1 - surface1.b1 + surface1.a1*diff_norm)/(diff_norm-diff)) 
		y_t = x_t*diff + c1

		x_2t = 2*x_t - surface1.a1
		y_2t = 2*y_t - surface1.b1

		c_ref = y_2t - (((surface1.y - y_2t)*x_2t)/(surface1.x - x_2t))

		print("Equation of reflected ray without using tan_theta , y =",((surface1.y - y_2t)/(surface1.x - x_2t)),"x +",+c_ref)

		y = ((surface1.y - y_2t)/(surface1.x - x_2t))*x +c_ref

		plt.plot(x,y,'black')
		plt.show()
	#tang_normal(-1,0,1,1,2,0,2)

	def tang_normal_type2(surface1):  # parabola 2

		print('The coefficient of y^2',+surface1.a)
		print('The coefficient of y', +surface1.b)
		print('The constant term', +surface1.c)


	    #equation of curve
		print("Equation of curve,  x =",surface1.a,"y^2 +",+ surface1.b,"y +", +surface1.c)
		

			

	#differentiate above equation to find tangent dy/dx

		diff = 1/(2*surface1.a*surface1.y + surface1.b)
		c1 = surface1.y - diff*surface1.x

	    #equation of tangent
		print("Equations of tangent at given point, y =",diff,"x +",+ c1)
		#y = diff*x + c1
	    #plt.plot(x,y)
	    # Show the plot
		#plt.show()

		#the slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		diff_norm = -1/(diff)
		c2 = surface1.y - diff_norm*surface1.x

		#equation of normal
		print("Equations of normal at given point, y =",diff_norm,"x +",+ c2)


		# set the incident light ray
		#y = a1x + b1

		#a1=0
		#b1=2
		### equation of incident line using two points formula :
		print("Equation of incident ray, y =",((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x)),"x +",+ (surface1.y - (((surface1.b1 - surface1.y)*surface1.x)/(surface1.a1 - surface1.x))))
		#equation of incident line
		#print("Equation of incident ray, y =",surface1.a1,"x +",+ surface1.b1)


		diff_line = ((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x))

		tan_theta = ((diff_line - diff_norm)/(1+(diff_line)*(diff_norm)))

		#print(tan_theta)

		slope_reflected = ((-tan_theta + diff_norm)/(1+(diff_norm)*(tan_theta)))


		#print(slope_reflected)

		c3 = surface1.y - slope_reflected*surface1.x

		#equation of reflected line
		print("Equation of reflected ray, y =",slope_reflected,"x +",+ c3)

		####x = np.linspace(-100, 100, 100)
		#y = a*x*x + b*x + c

		# Create the plot
		#plt.plot(x,y,label='y = x**2')
		#plt.plot(x,(1/2)* y + 7,label='y = (1/2) * (x**2) + 7')
		#plt.plot(x,y + 3,label='y = x**2 + 3')
		#plt.plot(x,y - 5,label='y = x**2 - 5')
		#plt.plot(x,y - 3,label='y = x**2 - 3')


		x_t = ((c1 - surface1.b1 + surface1.a1*diff_norm)/(diff_norm-diff)) 
		y_t = x_t*diff + c1

		x_2t = 2*x_t - surface1.a1
		y_2t = 2*y_t - surface1.b1

		c_ref = y_2t - (((surface1.y - y_2t)*x_2t)/(surface1.x - x_2t))

		print("Equation of reflected ray without using tan_theta , y =",((surface1.y - y_2t)/(surface1.x - x_2t)),"x +",+c_ref)


class Plane_surfaces:

	def __init__(someobj,a,b,x,y,a1,b1): # line equation
		someobj.a = a
		someobj.b = b
		someobj.x = x
		someobj.y = y
		someobj.a1 = a1
		someobj.b1 = b1

	def tang_normal_type3(surface1):

		print('The coefficient of x',+surface1.a)
		print('The constant term', +surface1.b)


	    #equation of curve
		print("Equation of line,  y =",surface1.a,"x +",+ surface1.b)
		

	#differentiate above equation to find tangent dy/dx

		diff = surface1.a
		surface1.b = surface1.y - diff*surface1.x

	    #equation of tangent
		print("Equations of tangent at given point, y =",diff,"x +",+ surface1.b)
	
		c1 = surface1.b
		#the slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		diff_norm = -1/(diff)
		c2 = surface1.y - diff_norm*surface1.x

		#equation of normal
		print("Equations of normal at given point, y =",diff_norm,"x +",+ c2)


		# set the incident light ray
		#y = a1x + b1

		#a1=0
		#b1=2
		### equation of incident line using two points formula :
		print("Equation of incident ray, y =",((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x)),"x +",+ (surface1.y - (((surface1.b1 - surface1.y)*surface1.x)/(surface1.a1 - surface1.x))))
		#equation of incident line
		#print("Equation of incident ray, y =",surface1.a1,"x +",+ surface1.b1)


		diff_line = ((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x))

		tan_theta = ((diff_line - diff_norm)/(1+(diff_line)*(diff_norm)))

		#print(tan_theta)

		slope_reflected = ((-tan_theta + diff_norm)/(1+(diff_norm)*(tan_theta)))
 

		#print(slope_reflected)

		c3 = surface1.y - slope_reflected*surface1.x

		#equation of reflected line
		print("Equation of reflected ray using tan theta, y =",slope_reflected,"x +",+ c3)

		
		
		x_t = ((c1 - surface1.b1 + surface1.a1*diff_norm)/(diff_norm-diff)) 
		y_t = x_t*diff + c1

		x_2t = 2*x_t - surface1.a1
		y_2t = 2*y_t - surface1.b1

		c_ref = y_2t - (((surface1.y - y_2t)*x_2t)/(surface1.x - x_2t))

		print("Equation of reflected ray without using tan_theta , y =",((surface1.y - y_2t)/(surface1.x - x_2t)),"x +",+c_ref)
		  
		# Adding legend, which helps us recognize the curve according to it's color
		
		# To load the display window

	
class Hyperbola_curve:

	def __init__(someobj,a,b,x,y,h,k,a1,b1):
		someobj.a = a
		someobj.b = b
		someobj.x = x
		someobj.y = y
		someobj.h = h
		someobj.k = k
		someobj.a1 = a1
		someobj.b1 = b1

	def tang_normal_type4(surface1): # parabola 1

		print('HYPERBOLA')

		print('The coefficient of x^2',+surface1.b*surface1.b)
		print('The coefficient of y^2', + -surface1.a*surface1.a)
		print('The coefficient of x', + -2*surface1.b*surface1.b*surface1.h)
		print('The coefficient of y', + 2*surface1.a*surface1.a*surface1.k)
		print('The constant term', surface1.b*surface1.b*surface1.h*surface1.h - surface1.a*surface1.a*surface1.k*surface1.k - surface1.a*surface1.a*surface1.b*surface1.b)


	    #equation of curve
		print("Equation of curve,  0 =",surface1.b*surface1.b,"x^2 +",+ -surface1.a*surface1.a,"y^2 +", + -2*surface1.b*surface1.b*surface1.h,"x +", +2*surface1.a*surface1.a*surface1.k,"y +",+surface1.b*surface1.b*surface1.h*surface1.h - surface1.a*surface1.a*surface1.k*surface1.k - surface1.a*surface1.a*surface1.b*surface1.b)
		

		x = np.linspace(-10, 10, 1000)
		#x = np.arange(-100,100,1)
		#x = np.linspace(0,2*math.pi,400)
		#y = surface1.a*x*x+ surface1.b*x+surface1.c
		#plt.plot(x, y, 'r')

	#differentiate above equation to find tangent dy/dx

		diff = (surface1.x*surface1.b*surface1.b - surface1.h*surface1.b*surface1.b)/(surface1.a*surface1.a*surface1.y - surface1.a*surface1.a*surface1.k)
		c1 = surface1.y - diff*surface1.x

		 #equation of tangent
		print("Equations of tangent at given point, y =",diff,"x +",+ c1)
		#y = diff*x + c1
		#plt.plot(x, y, 'b')



		#y = diff*x + c1
	    #plt.plot(x,y)
	    # Show the plot
		#plt.show()

		#the slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		diff_norm = -1/(diff)
		c2 = surface1.y - diff_norm*surface1.x

		#equation of normal
		print("Equations of normal at given point, y =",diff_norm,"x +",+ c2)
		#y = diff_norm*x + c2
		#plt.plot(x, y, 'orange')



		#y = diff*x + c1
	    #plt.plot(x,y)
	    # Show the plot
		#plt.show()

		# set the incident light ray
		#y = a1x + b1

		#a1=0
		#b1=2
		### equation of incident line using two points formula :
		print("Equation of incident ray, y =",((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x)),"x +",+ (surface1.y - (((surface1.b1 - surface1.y)*surface1.x)/(surface1.a1 - surface1.x))))

		#equation of incident line
		#print("Equation of incident ray, y =",surface1.a1,"x +",+ surface1.b1)
		

		#y = ((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x))*x + (surface1.y - (((surface1.b1 - surface1.y)*surface1.x)/(surface1.a1 - surface1.x)))
		#y = diff*x + c1
	    
	    # Show the plot
	    #plt.plot(x, y, 'green')
		#plt.plot(x,y,'green')
	    # Show the plot
		#plt.show()

		#plt.show()

		diff_line = ((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x))

		tan_theta = ((diff_line - diff_norm)/(1+(diff_line)*(diff_norm)))

		#print(tan_theta)

		slope_reflected = ((-tan_theta + diff_norm)/(1+(diff_norm)*(tan_theta)))


		#print(slope_reflected)

		c3 = surface1.y - slope_reflected*surface1.x

		#equation of reflected line
		print("Equation of reflected ray using tan_theta, y =",slope_reflected,"x +",+ c3)

		#y = slope_reflected*x + c3
		#plt.plot(x,y,'yellow')
		#plt.show()
		
		# To load the display window
		#using second set of equations, with mirroring technique:

		x_t = ((c1 - surface1.b1 + surface1.a1*diff_norm)/(diff_norm-diff)) 
		y_t = x_t*diff + c1

		x_2t = 2*x_t - surface1.a1
		y_2t = 2*y_t - surface1.b1

		c_ref = y_2t - (((surface1.y - y_2t)*x_2t)/(surface1.x - x_2t))

		print("Equation of reflected ray without using tan_theta , y =",((surface1.y - y_2t)/(surface1.x - x_2t)),"x +",+c_ref)

		#y = ((surface1.y - y_2t)/(surface1.x - x_2t))*x +c_ref

		#plt.plot(x,y,'black')
		#plt.show()



class Ellipse:

	def __init__(someobj,a,b,x,y,h,k,a1,b1):
		someobj.a = a
		someobj.b = b
		someobj.x = x
		someobj.y = y
		someobj.h = h
		someobj.k = k
		someobj.a1 = a1
		someobj.b1 = b1

	def tang_normal_type5(surface1):

		print('ELLIPSE')

		x2coeff = surface1.b*surface1.b
		y2coeff = surface1.a*surface1.a
		xcoeff = -2*surface1.b*surface1.b*surface1.h
		ycoeff = -2*surface1.a*surface1.a*surface1.k
		constant = surface1.b*surface1.b*surface1.h*surface1.h + surface1.a*surface1.a*surface1.k*surface1.k - surface1.a*surface1.a*surface1.b*surface1.b

		print('The coefficient of x^2', + x2coeff)
		print('The coefficient of y^2', + y2coeff)
		print('The coefficient of x', + xcoeff)
		print('The coefficient of y', + ycoeff)
		print('The constant term', constant)


	    #equation of curve
		print("Equation of curve,  0 =",x2coeff,"x^2 +",+ y2coeff,"y^2 +", + xcoeff,"x +", + ycoeff,"y +", + constant)
		

		x = np.linspace(-10, 10, 1000)
		#x = np.arange(-100,100,1)
		#x = np.linspace(0,2*math.pi,400)
		#y = surface1.a*x*x+ surface1.b*x+surface1.c
		#plt.plot(x, y, 'r')

	#differentiate above equation to find tangent dy/dx

		diff = -(surface1.b * surface1.b * (surface1.x - surface1.h)) / (surface1.a*surface1.a*(surface1.y - surface1.k))
		c1 = surface1.y - diff*surface1.x

		 #equation of tangent
		print("Equations of tangent at given point, y =",diff,"x +",+ c1)
		#y = diff*x + c1
		#plt.plot(x, y, 'b')



		#y = diff*x + c1
	    #plt.plot(x,y)
	    # Show the plot
		#plt.show()

		#the slope of the normal is (-1/slope of tangent), m1m2 = -1 perpendicular lines
		diff_norm = -1/(diff)
		c2 = surface1.y - diff_norm*surface1.x

		#equation of normal
		print("Equations of normal at given point, y =",diff_norm,"x +",+ c2)
		#y = diff_norm*x + c2
		#plt.plot(x, y, 'orange')



		#y = diff*x + c1
	    #plt.plot(x,y)
	    # Show the plot
		#plt.show()

		# set the incident light ray
		#y = a1x + b1

		#a1=0
		#b1=2
		### equation of incident line using two points formula :
		print("Equation of incident ray, y =",((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x)),"x +",+ (surface1.y - surface1.x * ((surface1.b1 - surface1.y) / (surface1.a1 - surface1.x))))

		#equation of incident line
		#print("Equation of incident ray, y =",surface1.a1,"x +",+ surface1.b1)
		

		#y = ((surface1.b1 - surface1.y)/(surface1.a1 - surface1.x))*x + (surface1.y - (((surface1.b1 - surface1.y)*surface1.x)/(surface1.a1 - surface1.x)))
		#y = diff*x + c1
	    
	    # Show the plot
	    #plt.plot(x, y, 'green')
		#plt.plot(x,y,'green')
	    # Show the plot
		#plt.show()

		diff_line = (surface1.b1 - surface1.y)/(surface1.a1 - surface1.x)

		tan_theta = (diff_line - diff_norm) / (1 + diff_line * diff_norm)

		#print(tan_theta)

		if tan_theta < 0:
			tan_theta = -tan_theta

		slope_reflected = (diff_norm - tan_theta) / (1 + (diff_norm * tan_theta))


		#print(slope_reflected)

		c3 = surface1.y - slope_reflected*surface1.x

		#equation of reflected line
		print("Equation of reflected ray using tan_theta, y =",slope_reflected,"x +",+ c3)

		#y = slope_reflected*x + c3
		#plt.plot(x,y,'yellow')
		#plt.show()
		
		# To load the display window
		#using second set of equations, with mirroring technique:

		x_t = ((c1 - surface1.b1 + surface1.a1*diff_norm)/(diff_norm-diff)) 
		y_t = x_t*diff + c1

		x_2t = 2*x_t - surface1.a1
		y_2t = 2*y_t - surface1.b1

		c_ref = y_2t - (((surface1.y - y_2t)*x_2t)/(surface1.x - x_2t))

		print("Equation of reflected ray without using tan_theta , y =",((surface1.y - y_2t)/(surface1.x - x_2t)),"x +",+c_ref)

		#Testing other method
		c4 = surface1.b1 - (diff * surface1.a1)

		x_norm = (c2 - c4) / (diff - diff_norm)
		y_norm = (diff_norm * x_norm) + c2

		x_out = x_norm - (surface1.a1 - x_norm)
		y_out = y_norm - (surface1.b1 - y_norm)

		slope = (y_out - surface1.y) / (x_out - surface1.x)
		intercept = (-slope * x_out) + y_out

		print("Equation of reflected ray without using tan_theta 2, y =",slope,"x +", + intercept)

		#y = ((surface1.y - y_2t)/(surface1.x - x_2t))*x +c_ref

		#plt.plot(x,y,'black')
		#plt.show()




# curv1 = Curved_surfaces(1,0,1,1,2,0,2)
# curv1.tang_normal_type1()	


# curv2 = Curved_surfaces(1,2,4,9,1,3,2)
# curv2.tang_normal_type2()

# plane1 = Plane_surfaces(1,0,2,2,4,3)
# plane1.tang_normal_type3()

# curv4 = Hyperbola_curve(1,2,4,9,1,3,2,2)
# curv4.tang_normal_type4()

curv5 = Ellipse(1,2,4,9,1,3,2,2)
curv5.tang_normal_type5()




#figure out point where ray meets curve, given the direction is known (currently a random point on the curve is chosen and the slope is found)
#get hyperbola from github, then create the ellipse
