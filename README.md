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
Began to work on refractions. Found that previous method was incorrect (theta is not simply arctan of slope)
Began implementing flat surface; replacing the glass "chunk" (now only one boundary, not two)
Created graph for illustration
Found equation for using Snell's Law with vectors but is very complicated, must figure out how it works
Begun work to implement circular lens, must still fix tangent
Tangent fix, lens appears to work

Week 6: Goals: Check that "lens" is correctly showing the physical properties (make sure physics equations work), quick check on polar coordinates.
Lens appears to be somewhat working, but not following physics calculations exactly. More investigation required.
Spent most of week working on intermediate presentation.

Week 7: Goals: Add in citations, work on polar coordinates, continue to check physics of lens.
Mainly worked on polar coordinates mathematics. Having trouble implementing them and not sure if possible in all cases.

Week 8: Goals: Fix the derivatives, continue to check physics of lens, start working on rotations, polar coordinates.
Added in gradients of curve at the intercept points. This is done using an import called numdifftools which calculated the derivative direction. This gradient represents the *normal* to the tangent and thus replaces the normDir.
Rotation of parabola, ellipse, and hyperbola completed. Are they rotating the "correct" direction? See ellipse for example

Week 9: Goals: Confirm different parabola forms are same, implement paper, continue to check physics of lens.
Added a flag in Refractions for when surface is part of lens or not to make more readable.
Different parabola forms are the same. In fact, the current form is *more* useful for rotation (for ellipse we change to form Ax^x + BXY + Cy^2 + Dx + Ey + F=0).
Added in Sneha's 3D models.
--> 