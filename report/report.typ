// TODO: title

#let vv(x) = $arrow(#x)$
#set math.vec(delim: "[")
#set math.mat(delim: "[")

#set text(size: 11pt)

#outline()

= Introduction
In order to control attitude and translation, spacecraft can make use of attitude control thrusters. Attitude control thrusters are small thrusters that can expell gas in a certain direction to impart a certian force at a specific point of the vehicle into a specific direction. Ideally, a combination of a sufficient ammount of these thrusters on specific points of the craft can provide any combination of torque and force within a certain magnitude. 

The purpose of this report is to find the optimal arrangement of thrusters on a certain spacecraft. The spacecraft in question is shown in @fig:spacecraft-render. This is a spacecraft design I came up causually some years ago, chosen as it provides a shape which is not to complex, but where the optimal solution is not immediatly obvious. A thruster as shown in @fig:thruster-design is used. The reason for this design is the bi-directional thrust, which allows the thrust to be negative in the mathematical model, dramatically simplifying the model.
The challenges of the optimiyation are mostly the restrictions of the thruster placement.
The thrusters need to lie on the surface and the thrust direction needs to be parallel to the surface.

#figure(image("thruster.png", width: 40%), caption: [The concept of the thruster design])<fig:thruster-design>
#figure(image("spacecraft_render.png", width: 75%), caption: [The spacecraft in question. Only the surface where thruster placement is allowable is visible])<fig:spacecraft-render>

// TODO: descriptive of optimization + Introduce spacecraft


= Problem formulation
The moment ($M_i$) and force ($F_i$) a single thruster imparts on the spacecraft can be written as:
$
  vv(M)_i &= T_i (vv(r)_i times vv(d)_i) \
  vv(F)_i &= T_i vv(d)_i
$
where $vv(r)_i$ is the position of the thruster, $vv(d)_i$ is the direction and $T_i$ the thrust, that thruster exerts.
The final forces and moments are the sum of the forces and moments from the individual thrusters. 
$
  vec(vv(M), vv(F)) &= mat(
    vv(r)_1 times vv(d)_1, ..., vv(r)_n times vv(d)_n;
    vv(d)_1, ..., vv(d)_n
  ) vec(T_1, dots.v, T_n) \
  &= A vv(T)
$
If $n = 6$, $A$ is a square matrix and, if $A$ is invertible, the thrust of the thrusters can be calculated as
$
  vv(T) = A^(-1) vec(vv(M), vv(F))
$
The goal of this project is to find a combination of thrusters on a specific spacecraft that most 'efficiently' provide the thrusts for certain input thrusts/torques, i.e. where the average magnitude of $vv(T)$ is minimal.


== Optimization problem
// (also present mathematical problem statement)
My initial attempt to represent the 'average magnitude of $vv(T)$' mathematically is to solve $vv(T)$ for multiple predefined values of $[vv(M), vv(F)]^T$ and consider the (smooth) maximum absolute value between all the solutions the value to minimize. However, this approach brought multiple problems:
- Only a small subset of all possible inputs was represented
- The maximum thrust was considered, where in practice the average thrust is more important as it relates to fuel consumption
- The function evaluation involved solving multiple 6x6 system of equations
- There are many cases where $A$ is not invertible.

Instead, I settled on minimizing of $abs(det(A^(-1)))$, since the determinant is representative of how much a vector is scaled. Minimizing $abs(det(A^(-1)))$ is equivalent to maximizing $abs(det(A))$ (Since $det(A^(-1)) = det(A)^(-1)$). To archieve a smoother version of the absolute value, the final optimization function is.
$
  f(x) = -sqrt(det(A)^2 + 0.1^2) \
  A = mat(
    vv(r)_1 times vv(d)_1, ..., vv(r)_6 times vv(d)_6;
    vv(d)_1, ..., vv(d)_6
  )
$

== Modelling aspects
// (incl. discussion of simplification choices)
Since the matrix $A$ is constructed from the vectors $vv(r)_1 ... vv(r)_6$ and $vv(d)_1 ... vv(d)_6$, it is easy to conclude to model the problem as a 36-dimensional problem with equality constraints that force $vv(r)_i$ to conform to the surface of the vehicle and $vv(d)_i$ to be normalized and perpendicular to the the surface. However, there are 2 ways to reduce the dimensionality of the problem while simultaniously reduce the ammount of constraints.

Firstly, since each thruster has to be confined to the surface, the position can be defined by a 2d vector, which will have the components $u_i$ and $v_i$. To retrieve the 3d position form a 2d vector, an image as shown in @fig:position-map can be used. Each point on the image represent a vector on the surface of the object. This is  a technique common in computer graphics, called UV-mapping /*TODO: source*/. The red, green and blue components of the color in the image correspond to the x, y and z components in 3d space.

A similar mapping process can be done with the surface normal vector. This map is split up into the positive  and negative components of the normal vector (@fig:pos-normal-map and @fig:neg-norml-map), because of the limitations of the images and the tools to create them. Subtracting @fig:neg-norml-map from @fig:pos-normal-map gives the true normal vector. Since the direction is orthogonal to the surface normal, its rotation can now be expressed by a single variable, called $a_i$.

These two simplifications allow to reduce the problem from 6 variables per thruster to 3, yielding 18 variables in total. In addition the large amount of complicated constraints has been reduced to a single constraint: the 2d-coordinates of all thrusters must fall within the area for which a surface point is defined. For this, a signed distance field is defined, as seen in @fig:sdf-map. The signed distance field defines the distance of each point to a surface, with points within the object beeing negative. Therefore the signed distance field can be used directly as a constraint for a single thruster. The final inequality constraint is the maximum sdf of all thrusters.
$
  g(vv(u), vv(v), vv(a)) = limits("Max")_i ("sdf"(u_i, v_i)) <= 0
$

$vv(u)$, $vv(v)$ and $vv(a)$ can take any value, but repeat outside of the (0, 1) domain. It was found that letting the values grow to large numbers with this behaviour makes it significantly more likely to find the global minimum than otherwise. I have no explenation for this behaviour.

#let grid-width = 90%
#grid(columns: 2, row-gutter: 1em,
  [#figure(image("../models/textures/ih/positions.png", width: grid-width), caption: [Position map. Image is remapped from (0, 1) to (-1, 1)])<fig:position-map>],
  [#figure(
    box(fill: black, inset:0pt, image("../models/textures/ih/sdf.png", width: grid-width)),
    caption: [Inverse SDF map. Image is remapped from (0, 1) to (1, -1) before use])<fig:sdf-map>],
  [#figure(image("../models/textures/ih/normals_pos.png", width: grid-width), caption: [Positive normal map])<fig:pos-normal-map>],
  [#figure(image("../models/textures/ih/normals_neg.png", width: grid-width), caption: [Negative normal map])<fig:neg-norml-map>],
)

Lastly, the symmetry of the spacecraft can be exploited to place pairs of thrusters instead of individual thrusters and maintain the symmetrical design. This has the risk of eliminating potential optimal solutions but reduces the number of design variables to 9. This is performed by representing the symmetry of the 3d model in the symmetry of the UV-map and setting the the state of the last 3 elements as written below.
$
  u_(4,5,6) &= 1 - u_(1,2,3) \
  v_(4,5,6) &= v_(1,2,3) \
  a_(4,5,6) &= -a_(1,2,3) \
$

= Initial problem investigation
It will be clear from the simplified problem, the function is neither monotonous, nor convex. It is also easy to see from @fig:position-map that the domain is non-convex as well.
== Boundedness
The optimization function is 
$
  f(vv(x)) = -det(A)^2
$
Whose upper bound occurs as $det(A) = 0$. So $f(x) <= 0.0$. To find the lower bound, the bounds of $abs(det(A))$ need to be known. Hadarmads inequality states that $abs(det(A)) <= product_(i=1)^n norm(vv(w)_i)$, where $vv(w)_i$ are the column vectors of $A$ @hadamard. In this case, 
$
  vv(w)_i &= vec(vv(r)_i times vv(d)_i vv(d)_i) \ 
  <=> norm(vv(w)_i) &= sqrt(norm(vv(r)_i times vv(d)_i)^2 + norm(vv(d)_i)^2) \
  &<= sqrt(norm(vv(r)_i)^2 + 1) \
$
The maximum value of $norm(vv(r_i))$ accross the whole surface was found to be around 1.087. Therefore:
$
  det(A)^2 <= (1.087^2 + 1)^6 tilde.equiv 107.799 \ 
  -107.799 <= f(vv(x)) <= 0.0
$
Giving the exact bounds of the problem.

== Smoothness
The determinant varies smoothly with respect to the components of the matrix, so $f(vv(x))$ is smooth with respect to the positions and directions of the thrusters $vv(r)_i$ and $vv(d)_i$. The direction also varies smoothly with the angle $a_i$. Thus the question of smoothness relies on the smoothness of the position/normal maps. The position map is continuous because of the way it is generated: adjacent, valid points on the 2d map are also adjacent in 3d space. The smoothness of the position map is slightly non-smooth 

// E.g. boundedness, monotonicity, convexity
// Numerical noise
// Sensitivity analysis

= Initial optimization on simplified problem (e.g. 2 variables)
// Motivation of optimization approach, choices
// Investigation of obtained optimum
// Observations, interpretation of results, conclusions
= Optimization of actual problem
// (Same points as above)
// (possibly including different variations of model, problem formulation, optimization approach)

= Interpretation of results
== Sensitivity analysis
= Conclusions and recommendations

#bibliography("bibliography.bib")

#pagebreak()
= Appendices
== Source code
// (preferably, include this digitally)
== Additional graphs/data
// when applicable
