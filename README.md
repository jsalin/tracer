# tracer
A very simple software 3D ray tracing engine

(C) 2020 Jussi Salin

Just programmed by myself as a code sample as I had never used Golang before or developed any graphics code in ray tracing style, so two personal interests were pursued at once. Not intented for any other use or collaboration at this point.

The graphics style quickly resembles simple rasterization from early 90's 3D games, but actually rays (vectors) and intersections between those and polygons are calculated to ultimately produce each pixel. It is a reason this is so slow even on a modern PC. Relies on math.Atan2 function heavily. No GPU acceleration. The shading in original revision was simple flat shading, but now it's somewhere between gouraud and phong. The sides of polygons are sometimes clearly visible, but there is a circular highlight coming from light source, which changes in size based on distance and gives the mesh a reflective plastic or metal like look.

Only external dependency is Shiny, used for getting the calculated pixels to a window on the screen. By default, 16 concurrent threads are used to accelerate rendering for multi core CPUs. Supports one 3D mesh made of faces (polygons) of three vertexes (points in 3D space). Has one moveable RGB omni light source. More could be added to the code easily. Camera (the viewpoint) is static at (0,0,0), so the scene should be moved and rotated towards that coordinate for "moving inside a world" style effect in the future.

The included example STL files are the duck from 3D Studio examples and a simple sphere made with Blender. The file and scaling can be choosed by editting constants near beginning of the code file.

Ideas for future improvements:
- add more optimizations than the per-vertex atan2 pre-calc (which was already huge improvement alone, like 10x)
- port to C++ to test if there are any speed increases compared to Golang
- texture mapping

Directory structure:
- "go" has the code and data files
- "videos" has some videos of the tracer in action (few of earlier stages before pushing to GitHub)
- "screenshots" has some still screen shots, including interesting looking clitches that occurred during development

There is also a few higher quality sample videos outside GitHub, that were not rendered in real-time:
- https://youtu.be/D9EQBfdgmTY
- https://youtu.be/7mgVlcHT6SM