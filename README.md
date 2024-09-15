# Ray-Tracing
This is the third assignment of the computer graphics course. The assignment comes with a template and my tasks were to implement ray tracing to generate realistic image of 3D scene.

### Implementation Workflow (My Tasks) 
#### To deal with diffuse, specular, ambient color shading for the Phong interpolation, shadow and reflection in main.cpp:

1. Calculate the color at the nearest intersection point for the given ray in TraceRay(…);
2. Calculate the corresponding color on the given the ray-surface intersection information in Shade(…).

#### To deal with various ray-intersections in hittable.cpp:
1. Perform intersection calculation between a ray and a triangle in Triangle::Hit(…);
2. Perform intersection calculation between a ray and a sphere in Sphere::Hit(…);
3. Perform intersection calculation between a ray and a quadric surface in Quadric::Hit(…).
