# Reflection_Surfaces
<!-- 
#figure out point where ray meets curve, given the direction is known (currently a random point on the curve is chosen and the slope is found)

Week 1: Added in the Ellipse in same fashion as previously created Hyperbole.  Now have
equation, tangent, normal line, and tan_theta between two intersecting lines. Slight fix to 
Hyperbola constant term.

Week 2: Goal: Turn previous equations "slopes" into vector functions. 
Question: If ray is a vector, how to find where ray meets ellipse? -> Quadratic formula to solve for distance variable
Discovered that given "points" on hyperbola and ellipse were not actually on curves, so rays were not hitting the shapes before
Has been fixed to start with a point and a direction for the ray, finding the intercept, and keeping everything else the same for now
Thought for later: how to find the distance variable of ray in future in a more automatic fashion?

Week 3: Goal: Complete finding ray intercept points for easier shapes; add in refractions 
Cleaned up comments and code
Added in plots of curves and rays (except hyperbola)
Must still convert tangent, normal, and output line into vector notation

Week 4: Goal: Finish converting lines to vector equations & add refractions
Hyperbola plot added using contour
Special case found that breaks derivative: curv2 = Curved_surfaces(1,2,4,9,1,3,-1) [derivative is 1/0] -> fixed
Want to convert dx and dy inputs into a single input array -> will allow for higher dimensions
Completed turning lines into vectors and using vectors to calculate reflected ray -> will need to adjust some things for higher dimensions
Functions now output intercept point and direction of output ray. Can now chain functions together

Week 5: Goals: Fix plots so they show curves as the rays they represent; one plot for multiple curves; add in basic refractions
Continued to clean up and remove old code; realized current method to find output direction is same as "mirroring technique" used before, not the tan_theta method
Rays now shown as such in plots
Plots can now show chain of rays, although can be very confusing to look at
Added control parameters to control what kind of plot is shown (only 1 vs multiple; show normal and tangent; show only where light rays actually travel)
--> 