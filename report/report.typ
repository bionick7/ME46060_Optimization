#show figure.where(kind: table): set figure.caption(position: top)

#let vv(x) = $arrow(#x)$
#set math.vec(delim: "[")
#set math.mat(delim: "[")
#set heading(numbering: "1.1")

// TODO: title
#set align(center)
#text(size: 32pt)[
  Optimal placement of reaction control thrusters conforming to a spacecraft's surface
]

#set align(left)

#pagebreak()
#set page(numbering: "1/1", number-align: right, header: grid(columns: (100%-3em, auto), [Nick Felten], [5223679]))
#counter(page).update(1)
#set text(size: 11pt)

#outline(
  title: [Table of contents],
  indent: auto
)

#heading(level: 1, outlined: false)[Table of contents]
#figure(table(columns: 2, stroke: none, align: left,
  $vv(M)$, [Moments exserted by the thruster system],
  $vv(F)$, [Forces exserted by the thruster system],
  $vv(T)$, [Thrusts of each thruster (can be negative)],
  $vv(r)_i$, [Position of thruster $i$],
  $vv(d)_i$, [Normalized direction of thruster $i$],
  $vv(u)$, [x-position of the thrusters on the 2d map],
  $vv(v)$, [y-position of the thrusters on the 2d map],
  $vv(a)$, [Angle of the thrusters on the surface],
  $A$, [Thruster matrix],
  $vv(x)$, [Design variables],
  $f(vv(x))$, [Objective function],
  $g(vv(x))$, [Inequality constraint],
  $vv(phi.alt)(x)$, [2D Path representing the surface],
))

= Introduction
In order to control attitude and translation, spacecraft can make use of reaction control thrusters. Reaction control thrusters are small thrusters that can expel gas in a certain direction to impart a controlled force at a specific point of the vehicle into a specific direction. Ideally, a combination of a sufficient amount of these thrusters on specific points of the craft can provide any combination of torque and force within a certain magnitude. Some thruster arrangements may be more or less efficient at converting their thrust into certain torques or moments than others.

The purpose of this report is to find the optimal arrangement of thrusters on a certain spacecraft. The spacecraft in question is shown in @fig:spacecraft-render. This is a spacecraft design I came up casually some years ago. It is chosen as it provides a geometry which is not to complex, but where the optimal solution is not immediately obvious. The geometry also contains holes, and some sharp angles, which might make the optimization more interesting. The thruster design used is as shown in @fig:thruster-design. The reason for this design is the bi-directional thrust, which allows the thrust to be negative in the mathematical model. This significantly simplifies the calculations. The challenge of the optimization is mostly the restrictions of the thruster placement.
The thrusters need to lie on the surface of the vehicle and the thrust direction needs to be parallel to the surface.

A more rigorous formulation of the optimization problem and the model is given in @ch:problem-formulation. Some of the mathematical properties of the objective function, like bounds/continuity and noise are discussed in @ch:initial-problem-investigation. The optimization of a simplified 2D problem is discussed @ch:simp-problem while the actual optimization is discussed in @ch:true-problem. In @ch:results, the results of the optimization are commented on and their sensitivity is analyzed. The final conclusion is presented in @ch:conclusion.

#grid(columns: (30%, auto), align: bottom,
  [#figure(image("thruster.png", width: 100%), caption: [The concept of the thruster design, thrust can be applied in the positive or negative x-direction])<fig:thruster-design>],
  [#figure(image("spacecraft_render.png", width: 100%), caption: [The spacecraft in question. Only the surface where thruster placement is allowable is visible])<fig:spacecraft-render>]
)

= Problem formulation <ch:problem-formulation>
The moment ($M_i$) and force ($F_i$) a single thruster imparts on the spacecraft can be written as:
$
  vv(M)_i &= T_i (vv(r)_i times vv(d)_i) \
  vv(F)_i &= T_i vv(d)_i
$
where $vv(r)_i$ is the position of the thruster, $vv(d)_i$ is the direction and $T_i$ the thrust, that thruster exerts.
The final forces and moments are the sum of the forces and moments from the individual thrusters. 
$
  vec(vv(M), vv(F)) &= mat(
    vv(r)_1 times vv(d)_1, dots.c, vv(r)_n times vv(d)_n;
    vv(d)_1, dots.c, vv(d)_n
  ) vec(T_1, dots.v, T_n) \
  &= A vv(T)
$
If $n = 6$, $A$ is a square matrix and, if $A$ is invertible, the thrust of the thrusters can be calculated as
$
  vv(T) = A^(-1) vec(vv(M), vv(F))
$
If $A$ is not invertible, there are forces/moments that the design cannot provide.
The goal of this project is to find a combination of thrusters on a specific spacecraft that most 'efficiently' provide the thrusts for certain input thrusts/torques, i.e. where the average magnitude of $vv(T)$ is minimal.

== Simplifications
Several assumptions are made in the problem statement that are be stated below.
- The geometry of the spacecraft does not influence the thrust of the thrusters. Thruster are considered to always give constant thrust no matter where they are placed
- Thrusters can be placed anywhere on the surface. Some of the surface has been omitted from the model to indicate places where no placement is possible. The remaining surface is all considered an equally valid destination. In reality, the placement of some thrusters might restrict the valid placement of others.
- All forces and torques are weighed equally. For the purpose of this optimization, 1 unit of force has the same utility as one unit of torque. Since all points in this problem lie within $[-1, 1]^3$, force and torque have similar scales. The importance of different forces/moments can easily be weighted by scaling the rows of the matrix $A$ to account for e.g. the mass moment of inertia.

== Optimization problem
// (also present mathematical problem statement)
The initial attempt to represent the 'average magnitude of $vv(T)$' mathematically is to solve $vv(T)$ for multiple predefined values of $[vv(M), vv(F)]^T$ and consider the (smooth) maximum absolute value between all the solutions the value to minimize. However, this approach brought multiple problems:
- Only a finite subset of all possible inputs can ever be represented
- The maximum thrust is considered, where in practice the average thrust is more important as it relates to fuel consumption
- The function evaluation involved solving multiple 6x6 system of equations, which is rather slow.
- There are many cases where $A$ is not invertible.

Instead, I settled on minimizing of $abs(det(A^(-1)))$, since the determinant is representative of how much a vector is scaled. Minimizing $abs(det(A^(-1)))$ is equivalent to maximizing $abs(det(A))$, since $det(A^(-1)) = det(A)^(-1)$. To archive a smoother version of the absolute value, the final optimization function is.
$
  f(x) = -sqrt(det(A)^2 + 0.1^2) \
  A = mat(
    vv(r)_1 times vv(d)_1, dots.c, vv(r)_6 times vv(d)_6;
    vv(d)_1, dots.c, vv(d)_6
  )
$

== Modelling aspects
// (incl. discussion of simplification choices)
Since the matrix $A$ is constructed from the vectors $vv(r)_1 ... vv(r)_6$ and $vv(d)_1 ... vv(d)_6$, it is easy to conclude to model the problem as a 36-dimensional problem with equality constraints that force $vv(r)_i$ to conform to the surface of the vehicle and $vv(d)_i$ to be normalized and perpendicular to the the surface. However, there are 2 ways to reduce the dimensionality of the problem while simultaneously reduce the amount of constraints.

Firstly, since each thruster has to be confined to the surface, the position can be defined by a 2d vector, which has the components $u_i$ and $v_i$. To retrieve the 3d position form a 2d vector, an image as shown in @fig:position-map can be used. Each point on the image represent a vector on the surface of the object. This is  a technique common in computer graphics, called UV-mapping /*TODO: source*/. The red, green and blue components of the color in the image correspond to the x, y and z components in 3d space.

A similar mapping process can be done with the surface normal vector. This map is split up into the positive  and negative components of the normal vector (@fig:pos-normal-map and @fig:neg-normal-map), because of the limitations of the images and the tools to create them. Subtracting @fig:neg-normal-map from @fig:pos-normal-map gives the true normal vector. Since the direction is orthogonal to the surface normal, its rotation can now be expressed by a single variable, called $a_i$. The direction is retrieved by building an orthonormal basis around the normal vector via the Gram-Schmidt process, where $pi dot a_i$ defines a rotation in the plane orthogonal to the normal vector.

These two simplifications allow to reduce the problem from 6 variables per thruster to 3, yielding 18 variables in total. In addition the large amount of complicated constraints has been reduced to a single constraint: the 2d-coordinates of all thrusters must fall within the area for which a surface point is defined. For this, a signed distance field is defined, as seen in @fig:sdf-map. The signed distance field defines the distance of each point to a surface /*Source*/ , with points within the object being negative. Therefore the signed distance field can be used directly as a constraint for a single thruster. The final inequality constraint is the maximum sdf of all thrusters.
$
  g(vv(u), vv(v), vv(a)) = limits("Max")_i ("sdf"(u_i, v_i)) <= 0
$

$vv(u)$, $vv(v)$ and $vv(a)$ can take any value, but repeat outside of the (0, 1) domain.

#let grid-width = 90%
#grid(columns: 2, row-gutter: 1em,
  [#figure(image("../models/textures/ih/positions.png", width: grid-width), caption: [Position map. Image is remapped from (0, 1) to (-1, 1)])<fig:position-map>],
  [#figure(
    box(fill: black, inset:0pt, image("../models/textures/ih/sdf.png", width: grid-width)),
    caption: [Inverse SDF map. Image is remapped from (0, 1) to (1, -1) before use])<fig:sdf-map>],
  [#figure(image("../models/textures/ih/normals_pos.png", width: grid-width), caption: [Positive normal map])<fig:pos-normal-map>],
  [#figure(image("../models/textures/ih/normals_neg.png", width: grid-width), caption: [Negative normal map])<fig:neg-normal-map>],
)

Lastly, the symmetry of the spacecraft can be exploited to place pairs of thrusters instead of individual thrusters and maintain the symmetrical design. This has the risk of eliminating potential optimal solutions that are either non-symmetrical or has thrusters at the boundary, but reduces the number of design variables to 9. This is performed by representing the symmetry of the 3d model in the symmetry of the UV-map and setting the the state of the last 3 elements as written below.
$
  u_(4,5,6) &= 1 - u_(1,2,3) \
  v_(4,5,6) &= v_(1,2,3) \
  a_(4,5,6) &= -a_(1,2,3) \
$

= Initial problem investigation <ch:initial-problem-investigation>
There is no reason to assume that the objective function is either monotonous or convex and it will be clear in the simplified problem that it isn't. It is also clear from @fig:position-map that the domain is not convex as well.
== Bounds
The optimization function is 
$
  f(vv(x)) = -sqrt(det(A)^2 + 0.1^2)
$
Whose upper bound occurs as $det(A) = 0$. So $f(x) <= -0.1$. To find the lower bound, the bounds of $abs(det(A))$ need to be known. Hadamards inequality states that $abs(det(A)) <= product_(i=1)^n norm(vv(w)_i)$, where $vv(w)_i$ are the column vectors of $A$ @hadamard. In this case, 
$
  vv(w)_i &= vec(vv(r)_i times vv(d)_i vv(d)_i) \ 
  <=> norm(vv(w)_i) &= sqrt(norm(vv(r)_i times vv(d)_i)^2 + norm(vv(d)_i)^2) \
  &<= sqrt(norm(vv(r)_i)^2 + 1) \
$
The maximum value of $norm(vv(r_i))$ across the whole surface is around 1.087. Therefore:
$
  det(A)^2 <= (1.087^2 + 1)^3 tilde.equiv 10.383 \ 
  -10.383 <= f(vv(x)) <= -0.1
$
This gives a lower bound for the optimum.

== Smoothness & Continuity
The determinant varies smoothly with respect to the components of the matrix, so $f(vv(x))$ is smooth with respect to the positions and directions of the thrusters $vv(r)_i$ and $vv(d)_i$. The direction also varies smoothly with the angle $a_i$. Thus the question of smoothness and continuity relies on the smoothness of the position/normal maps. The position map is continuous because of the way it is generated: adjacent, valid points on the 2d map are also adjacent in 3d space. Since the surface normal is related to the gradient of the surface, the normal map represents the derivatives of the position map. Therefore the normal map is discontinuities if and only if the position map is non-smooth. This corresponds with sharp angles in the geometry of the mesh. Even though the shape is mostly smooth, there are some sharp corners in the mesh and some very tight bevels (only a few pixels wide), which are represented as very steep gradients in the normal map. Therefore there are some discontinuities or near discontinuities in the objective functions.

== Noise and precision
Storing data in images comes with two fundamental restriction in resolution. Both the resolution of the data itself, as well as the spatial resolution are limited by the image format. The position and normal maps are 2048$times$2048 pixel maps with a color depth of 16 bits per channel. The image is sampled by linearly interpolating between the closest pixels and can be considered, on a small scale, to consist of facets of size $5 dot 10^(-4)$. This can influence finite difference calculations, especially around sharp corners. The 16 bit color depth makes the normals accurate to the order $~ 1.5 dot 10^(-5)$ and positions accurate to the order $~ 3 dot 10^(-5)$. This is sufficient for the problem, and a much less dominant restriction than the image resolution.

I could not get an SDF image in 16 bit per channel color depth. Therefore it is also limited to 1024$times$1024 resolution as an SDF that decreases by less than 1 unit per pixel can lead to problems. This means that he constraint has significantly lower resolution than the optimization function. Because of the type of optimizer chosen and the way the constraint is implemented, however, the noise in the constraint is much less relevant.

// TODO: Smoothness drives optimization conditions

= Initial optimization on simplified problem <ch:simp-problem>
== Simplified problem model
Before solving the full problem, a simplified version of the problem is investigated in order to gage the complexity and nature of the optimization function. An easy way to do this is to look at the 2-dimensional equivalent of the original problem.

In the 2-dimensional case, there are 2 degrees of translational motion and a single degree of rotational motion. Therefore only 3 thrusters are necessary. Instead of defining the position of the thrusters as 2d point that maps to a 3d surface, it can be defined as a single value that maps to a 2d position via evaluating a path representing the surface. The constraint that the thrust direction needs to be orthogonal to the surface allows to get the direction directly from the position, as the normalized derivative of the the path. This reduces the number of variables in the 2D case to only 3. Similarly to the actual problem, the use of symmetry allows to reduce the problem to 2 variables by assuming that $x_3 = 1 - x_1$.

The path $vv(phi.alt)(x)$is represented by interpolating between 20 points with piecewise cubic hermite spline, as seen in @fig:path-2d. This way, some of the flat surfaces and sharp corners of the shape stay conserved while still offering smooth interpolation in some places. The path starts at (0.5, 0), reaches (0.5, 1) at $x = 0.5$ and returns to (0.5, 0) at $x = 1$. It is symmetrical along $x = 0.5$.

$
  vv(r)_(i, "2D") = vv(phi.alt)(x_i) \ 
  vv(d)_(i, "2D") = (dif vv(phi.alt)(x_i)) / (dif x_i) dot norm((dif vv(phi.alt)(x_i)) / (dif x_i))^(-1)
$

$
  A_"2D" = mat(
    mat(delim: "|", vv(r)_(1, "2D"), vv(d)_(1, "2D")), mat(delim: "|", vv(r)_(2, "2D"), vv(d)_(2, "2D")), mat(delim: "|", vv(r)_(3, "2D"), vv(d)_(3, "2D"));
    vv(d)_(1, "2D"), vv(d)_(2, "2D"), vv(d)_(3, "2D")
  )
$

== Algorithm selection
The reduction to two variables allows to plot the function across the whole domain at once, as seen in @fig:plot-f-2d. Furthermore, a cross-section along the $x_1 = 0.5$ is shown in @fig:cross-section-2d. These two figures confirm two of the difficulties that where suspected earlier: The objective function has many local optima and it contains a lot of plateaux and discontinuities.

#figure(image("contour_plot_2d.png"), caption: [Plot of $f(vv(x))$ for the simplified problem over its entire domain]) <fig:plot-f-2d>

#figure(image("func_cross_section_2d.png"), caption: [Plot of a cross-section of $f(vv(x))$ at $x_1$ = 0.5 for the simplified problem]) <fig:cross-section-2d>

This leads to the following criteria to choose an optimization algorithm:

+ The optimizer must be able to escape local optima.
+ The optimizer must not rely on the first/second derivative.

The last condition eliminates almost any method that are not direct search methods. Direct search methods include stochastic methods like simulated annealing or genetic algorithms, as well as methods like cyclic coordinate search and the Nelder-mead simplex method. /*TODO: Sources perhaps?*/ My initial thought was to use a stochastic method for their ability to escape local minima. I chose to use simulated annealing for its simplicity to implement. Indeed, for this relatively simple problem, simulated annealing performs relatively well and finds the global optimum in most runs. There are multiple global optima because of the symmetry in the objective function, but a solution is shown in @fig:path-2d. 

The thrusters are as far as possible from each other to maximize the angular momentum they can provide and are 'balancing' on the corners to be able to provide thrust in both x- and y- directions. Although this makes the solution technically very sensitive to the inputs, it doesn't make for an unreasonable design. This will be discussed more in the sensitivity analysis of the main problem.

#figure(image("solution_2d.png"), caption: [Plot of the path representing the 2d shape in red. Position of the thrusters in the solution plotted in blue.]) <fig:path-2d>

// - Analytical (expected) solution
// - Optimizer must be 0th order
// - Apply simulated annealing
// - Arrives at solution (sometimes)

// Motivation of optimization approach, choices
// Investigation of obtained optimum
// Observations, interpretation of results, conclusions

= Optimization of the actual problem <ch:true-problem>
The actual 3D problem has 9 design variables, so it is not possible to construct a plot like @fig:plot-f-2d, but similar properties can be assumed, since flat surfaces and tight corners are also present in the 3D case. In this more complicated case, simulated annealing turned out to be very unreliable and often converged to results which were obviously not optimal.

The reason for this might be the ease for simulated annealing to get trapped at low temperature within shallow plateaux which are local optima. This can be mitigated by 'reheating' the solution after it is stuck for to long, but the optimizer is still very inefficient in finding the right solution. Since there are a lot of variations of simulated annealing and a lot of parameters that can be adjusted, there might be a more efficient version of simulated annealing that can deal with this problem, but there is also a simpler way.

An efficient algorithm, robust to both un-smooth functions and plateaux, can be run on many random points on the domain in hopes that it finds a good global optimum in the process. After some testing, Powell's conjugate directions (PCD) method was chosen @powell. Since PCD is based on line-searches that use golden-section search or quadratic interpolation internally @fminbnd, it does not rely on derivatives and is pretty robust.

A major benefit of PCD in this case is that the search radius can be defined exactly with the bounds of the line-searches. The optimization performs better and converges faster on a smaller search radius, centered around 0. By limiting the search radius of the algorithm, the optimizer needs to deal with multiple optima less often and is more likely to find the closest local optimum quickly. It is also less likely to cycle between different optima without converging, as less of the function is visible at a time. This comes with the disadvantage that less of the solution space is explored at a time, so more cycles might be needed to find the local minimum. A good compromise is to set the search radius to 0.25.

Since PCD, like most direct search methods, does not account explicitly for constraints, the objective function needs to be transformed to give a penalty outside the feasible domain, such that the optimal solution cannot lie outside the domain. The follow transformation is used.
$
  tilde(f) = f + p(max(0, g))^2 #v(2em) p in RR^+
$
where $p$ determines the magnitude of the penalty. $p$ can be set quite high right away, since the optimizer needs to deal with very steep gradients (almost discontinuities) anyway. The value used is $p = 1000$.

To generate an initial value to use for the PCD, a random vector is generated with a uniform distribution in the domain $[0, 1]^9$ and re-generated until it conforms to the inequality constraint. After some initial tests, it was found that all the most promising results occur with 2 sets of thruster at $v_i >= 0.5$ (the upper half of the UV map, corresponding to the rear of the spacecraft) and one set at $v_i <= 0.5$. This makes physical sense, as the thrusters would benefit to be evenly distributed over the spacecraft to gain the maximum moment arm, with the rear of the spacecraft being wider and better for roll control. Ultimately, the domain for the initial value generation was constrained to accommodate these preferences and make the total search space of initial guesses smaller. Along with the restrictions in $v_i_0$, $u_i_0 < 0.5$ can be enforced, as the $f(x)$ is symmetric with respect to $u_i$. 
/*$
  v_(1,2)_0 in [0.5, 1]\
  v_(3)_0 in [0, 0.5]\
  vv(u)_0 in [0, 0.5]^3\
$*/
This reduces the hypervolume of the search space by a factor of $2^6 = 64$.

The interval of convergence for the PCD is set to a rather large number, 0.01. The idea is to have an algorithm that converges quickly, such that it can be run many times to find a good candidate for the optimal solution, then find the precise optimum with that candidate as a starting point in a second step.

// - Apply simulated annealing
// - Problems
// - Apply Powell's conj. directions on 100 randomly sampled points
// - Constrain domain

// (Same points as above)
// (possibly including different variations of model, problem formulation, optimization approach)

= Interpretation of results <ch:results>
After running the program multiple times with different rng states, most come up with similar solutions. This is a good sign, as it is likely that this family of solutions are include the global optimum. It is impossible to say though if any solutions is the true global minimum unless a better solution is found. The result that is shown here is the best one found with a utility function of $f(x) = -2.8347$. The values in terms of map coordinates ($u$, $v$) and $a$ are presented in @fig:solution-uva, while the solutions in terms of position ($vv(r)$) and direction ($vv(d)$) are presented in @fig:solution-rd. A 3d visualization of the solution is also visible in @fig:solution-3d.
// Solution makes sense physically

#figure(table(columns: 4,
  $i$,      $1$,      $2$,      $3$,
  $u$, $0.2355$, $0.4431$, $0.2622$,
  $v$, $0.1255$, $0.9431$, $0.1081$,
  $a$, $0.7742$, $0.0647$, $0.1260$,
), caption: [Solution in terms of map coordinates of first 3 thrusters]) <fig:solution-uva>

#figure(table(columns: 7,
    $i$,       $1$,       $2$,       $3$,       $4$,       $5$,       $6$, 
  $r_x$, $-0.3320$, $-0.0473$, $-0.3587$, $ 0.3327$, $ 0.0473$, $ 0.3584$,
  $r_y$, $ 0.8432$, $-0.9758$, $ 0.8344$, $ 0.8431$, $-0.9763$, $ 0.8342$,
  $r_z$, $ 0.0342$, $ 0.0036$, $-0.0050$, $ 0.0336$, $ 0.0033$, $-0.0060$,
  $d_x$, $ 0.6213$, $-0.6935$, $ 0.2543$, $ 0.6184$, $-0.5928$, $ 0.3233$,
  $d_y$, $ 0.6593$, $ 0.0775$, $ 0.3084$, $-0.6590$, $-0.0769$, $-0.3085$,
  $d_z$, $ 0.4234$, $-0.7032$, $-0.9165$, $-0.4281$, $ 0.7854$, $ 0.8945$,
), caption: [Solution in terms of 3d positions and directions]) <fig:solution-rd>

//#figure(image("solution_map.png"), caption: [Solution in map coordinates, overlaid over the (positive) normal map])
#figure(image("solution_3d.png"), caption: [Solution rendered in 3d]) <fig:solution-3d>

It can be seen from @fig:solution-3d that the thrusters from the solution form 3 pairs, one at each extremity of the spacecraft. Each pair forms a perpendicular cross. This solution maximizes both the distance of the thrusters from each other while making sure there is minimum alignment between any two thruster. It makes sense physically that this is one of the best solution that can be archived.

== Sensitivity analysis
// 2 kinds of sensitivity: w.r.t. 3d positions and w.r.t variables
// Since there are no active constraints, response of f to x is either 0 or discontinuous.
// 2nd point is discontinuous in u and v
// Takes advantage of corners, so theoretically very sensitive
// Not an issue in practice.

For this problem, a sensitivity analysis comes in two forms. One is the sensitivity with respect to the design parameters that have been optimized, which is what is typically done. These parameters, however have been chosen to reduce the dimensionality of the problem and do not correspond to any physical dimensions. Therefore, the sensitivity with respect to the position and direction is also calculated.

The single constraint was not an active constraint in the solution, so the only response that is relevant is the response of the objective function itself. Since the relevant point is a local optimum that does not lie on any constraint, the derivatives at that point with any value are either zero or the function is not smooth/discontinuous. Discontinuities can only come from the position/normal maps themselves. After looking at the discrete gradient of the normal map, it is found that points 2 and 4 (the pair in the front of the spacecraft) lie on a crease, so they have discontinuous derivatives and very large sensitivity. The other points lie on a smoother surface, so the derivative with respect to them is zero. 
The thruster pair is, in a way, exploiting the interpolation between pixels in the normal map to find a very narrow solution that results in a right angle for the thrusters.
This makes the solution very sensitive with respect to the placement of thruster 2/4 on the 2d map, but this sensitivity does not have any physical meaning.

/*#figure(table(columns: 4,
  $i$, $1$, $2$, $3$,
  $u$, $0$, $infinity$, $0$,
  $v$, $0$, $infinity$, $0$,
  $a$, $0$, $0$, $0$,
), caption: [Solution in terms of map coordinates of first 3 thrusters]) <fig:sensitivities-uva>*/

The sensitivity with respect to the positions and directions of the thrusters has a more direct physical meaning. It represents how much the objective function changes with manufacturing defects, design adjustments and physical deformation (due to e.g. temperature or stress).

The derivative with respect to each position $vv(r)_i$ and direction $vv(d)_i$ can be calculated analytically. 

$
  (dif f) / (dif mat(vv(r), vv(d))^T) = (dif f) / (dif det(A)) 
  (dif det(A)) / (dif A)
  (dif A) / (dif mat(vv(r), vv(d))^T)
$

$
  (dif f) / (dif det(A)) = det(A) / f approx -1
$
For the second derivative, Jacobi's formula can be used
$
  (dif det(A)) / (dif A) = "adj"(A) = det(A) A^(-1)
$

To get $(diff A) / (diff mat(vv(r), vv(d))^T)$, both matrices can be written in vector form.
$
  diff vec(A_11, A_21, dots.v, A_56, A_66) = mat(
    Omega_r_1, Omega_d_1, "", dots.c, 0;
    "", I, "", "", dots.v;
    "", "", dots.down, "", "";
    dots.v, "", "", Omega_r_6, Omega_d_6; 
    0, dots.c, "", "", I; 
  ) diff vec(vv(r)_1, vv(d)_1, dots.v, vv(r)_6, vv(d)_6)
$
where
$
  vv(r)_i times vv(d)_i = Omega_r_i vv(r)_i = Omega_d_i vv(d)_i
$

Finally, it needs to be considered that the direction vectors $vv(d)_i$ are always normalized. Therefore a change in one dimension of $vv(d)$ will have indirect changes, as the other components adjust.

Since, for any pair of vector components
$
  &d_i^2 + d_j^2 + ... = 1 \
  <=>& (diff d_i) / (diff d_j) = - d_j / d_i
$
The total derivative can be calculated as:
$
  (dif f) / (dif d_i) &= (diff f) / (diff d_i) + (diff f) / (diff d_j) (diff d_j) / (diff d_i) + (diff f) / (diff d_k) (diff d_k) / (diff d_i) \
  &= (diff f) / (diff d_i) - (diff f) / (diff d_j) d_i / d_j - (diff f) / (diff d_k) d_i / d_k\
  (dif f) / (dif vv(d)) &= 2 I - vv(d) (1/vv(d))^T (diff f) / (diff vv(d)) #h(2em) 1/vv(d) " refers here to the \n element-wise inversion"
$

The final sensitivities are presented in @tab:sensitivity-rd. Since both position and direction have a similar scaling (from -1 ... 1), it makes sense to present the sensitivity in linear derivatives, instead of logarithmic ones.
The highest sensitivity is in $r_x_3$ with an absolute value of 3.2. For reference: If the vehicle has a length of 10m, then a displacement of $r_x_3$ of 10 mm in the most sensitive direction results of a change in the utility function of 0.2 %, which I consider sufficiently stable.

Furthermore, it can be seen that some thrusters are more sensitive to design changes than other ones. Thruster pair 2/5, for example is much less sensitive than pair 3/6. The dimension at which the highest sensitivity for a thruster occurs can also give insight into the purpose of the unit. Since thrusters 2/5 are more sensitive to changes in $r_y$ and $d_z$, it can be deduced that is most important in pitch maneuvers (torque around the x-axis). With this information, the solution can be altered manually to better integrate it with the rest of the design, while maintaining its performance as much as possible.

//$
//  (dif mat(vv(a), vv(u), vv(v))^T) / (dif vv(x)) = mat(I_(9 times 9), -I_(9 times 9))
//$
/*
#figure(table(columns: 4,
  $i$,       $1$,        $2$,       $3$,
  $u$, $ 0.1407$, $-13.4541$, $-0.3078$,
  $v$, $-0.1660$, $ 11.0299$, $-0.1425$,
  $a$, $-0.1663$, $ -0.9259$, $-0.0001$
))
*/

#let log_sens = (
   (-1.2418, -0.0696, -3.1924,  1.2660,  0.0599,  3.1779),
   ( 0.6157, -1.3783,  0.7497,  0.6488, -1.4026,  0.7667),
   ( 0.8635, -0.0832, -0.6336,  0.8300, -0.0922, -0.8844),
   (-1.9138, -0.0374, -0.6434, -2.0050,  0.0287, -1.0214),
   (-1.2718, -0.3768, -0.8854,  1.4866,  0.3834,  0.7276),
   (-1.2939,  1.1169,  3.0216,  1.1591, -1.3929, -2.3653),
)
#figure(
  // Written for maximum obfuscation
  table(columns: 7, fill: (i, j) => if (i*j == 0) {white} else { gradient.linear(green, orange, space: oklab).sample(calc.abs(log_sens.at(j - 1).at(i - 1)) / 2.0 * 100%)}, $i$, ..range(1,7).map(x => $#x$), ..($r_x$, $r_y$, $r_z$, $d_x$, $d_y$, $d_z$).zip(log_sens).map(x => (x.at(0), ) + x.at(1).map(x => $#x$)).flatten(),), 

  caption: [The sensitivity with respect to different $vv(r)_i$ and $vv(d)_i$, colored by magnitude]
) <tab:sensitivity-rd>

= Conclusion & Recommendations <ch:conclusion>
The objective of this optimization problem was to find the optimum thruster placement for a specific spacecraft. Nominally the 6 necessary thrusters define a 36-dimensional design space with many quite complex constraints, as the thrusters need to conform to the surface of the spacecraft. With usage of position/normal maps however, the problem can be restated with 18 design variables, and a single (non active) constraint. With the use of symmetry, the amount of design variables can be further reduced to 9. This restatement comes with the cost of introducing noise, discontinuities and plateaux into the objective function, which limits the choice of optimization function. Ultimately, running Powell's conjugate direction method with a limited search radius from multiple randomly chosen starting points yields a good result most of the time. 

The best solution found is presented in @fig:solution-3d. The sensitivity within the restated problem is very high, as one thruster lies on a sharp edge, corresponding to a near discontinuity in the objective function. The solution is much less sensitive in the original problem, which is much more representative of the physical reality.

The success of this methodology relies heavily on the geometry and the UV-map that maps the 3d surface to a 2d plane. Geometry with less sharp corners and less flat planes will be easier to optimize. The choice of mapping also strongly affects the objective function. The reduction of design variables by employing symmetry required a specific mapping to work. For some geometries, the original 36-dimensional problem formulation might be more advantageous, since the objective function is much smoother, and can be looked into for future research. The main issue with that approach is to define the problem constraints mathematically.

#bibliography("bibliography.bib")
#pagebreak()

#set heading(numbering: "A")
#counter(heading).update(0)
= Aknowledgments
normal/position maps (@fig:position-map @fig:pos-normal-map @fig:neg-normal-map), as well as the 3D renders where created in Blender (https://www.blender.org/download/)\
signed distance field from @fig:sdf-map was generated in SDF maker by job talle (https://jobtalle.com/SDFMaker/) \
All computations made in matlab

= Source code
Source code can be found on github (https://github.com/bionick7/ME46060_Optimization)
// (preferably, include this digitally)
//== Additional graphs/data
// when applicable
