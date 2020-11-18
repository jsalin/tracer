/* * * * * * * * * * * * * * * * * * * * * * * * *
 * A very simple software 3D ray tracing engine  *
 * (C) 2020 Jussi Salin - See the README.md file *
 * * * * * * * * * * * * * * * * * * * * * * * * */

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This C++ version is about 10x faster (per thread) than the original Golang version.
// Maybe the math functions or std::vector is much faster then Golang versions?
//
// Windows threads are used as there was difficulty in getting either pthreads-win32 or
// SDL_threads working. There is an experimental #ifdef for trying out pthreads.
//
// * NOTE * Use only x64 Release profile when building, because Debug is unuseably slow!
//
// Uses only SDL2 as an external reference (a nuget package), to draw graphics on a window.
//
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <SDL.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <windows.h>

using namespace std;

/// <summary>
/// Choose data type we use for floating points (double, float, ...)
/// </summary>
typedef float tFloat;

/// <summary>
/// Vertex structure
/// </summary>
typedef struct
{
    tFloat x, y, z;     // Coordinate of vertex
    tFloat r, g, b;     // Color of vertex
    tFloat xa, ya;      // Speed up, store pre-calculated angles to camera vertex for each vertex (once before frame)
} Vertex;

/// <summary>
/// Face structure
/// </summary>
typedef struct
{
    int a, b, c;        // Index of each vertex of polygon
} Face;

/// <summary>
/// Pixel structure
/// </summary>
typedef struct {
    tFloat depth;       // Depth the color was calculated from
    tFloat r, g, b;     // Color
} Pixel;

/// <summary>
/// Mesh structure, which holds all vertex and face data
/// </summary>
typedef struct
{
	Face f[1000];
	Vertex v[1000];
	int fCount;
	int vCount;
} Mesh;

// Consts (adjustable)
#define THREADS 16							// How many threads to use in parallel rendering
#undef USE_PTHREADS							// Define on POSIX threads system, undefine on Windows
#undef SHOW_PROGRESS						// Define to show each row instantly as it gets rendered
const int width = 640;						// Render area width in pixels
const int height = 480;						// Render area height in pixels
const tFloat contrastDiffuse = 0.05f;		// A multiplier to shading, to increase or decrease contrast of the overall image
const tFloat contrastSpecular = 0.5f;		// A multiplier to shading, to increase or decrease contrast of the overall image
const tFloat ambient = 0.0;					// Ambient brightness, between 0.0 (black is black) and 1.0 ("fullbright")
const tFloat fov = 90.0;					// Field of vision (in degrees)
const char* filename = "../go/ball.stl";    // File to load (the 3d mesh)
const tFloat filescale = 100.0f;			// Scale for the file, as some meshes can be really big or really small (use 0.07 for duck2.stl and 100.0 for ball.stl)
const tFloat speed = 0.25;					// Speed multiplier for animation

// Globals
SDL_Window* window = NULL;
SDL_Surface* screenSurface = NULL;

/// <summary>
/// FindVertex returns index of a vertex in a list by coordinates
/// </summary>
int FindVertex(tFloat x, tFloat y, tFloat z, Mesh &mesh) {
	for (int i = 0; i < mesh.vCount; i++) {
		if ((mesh.v[i].x == x) && (mesh.v[i].y == y) && (mesh.v[i].z == z)) return i;
	}
	return -1;
}

/// <summary>
/// LoadSTL loads ascii STL file. Works with a single mesh with three vertex faces only.
/// </summary>
bool LoadSTL(string filename, Mesh &mesh)
{
	tFloat xmin = 0;
	tFloat ymin = 0;
	tFloat zmin = 0;
	tFloat xmax = 0;
	tFloat ymax = 0;
	tFloat zmax = 0;
	float x, y, z;

	// Initialize empty lists for vertexes and faces
	mesh.vCount = 0;
	mesh.fCount = 0;
	vector<int> vertexIndexes;

	// Load the lines from file to a string array
	ifstream file(filename);
	if (!file) {
		cout << "Could not open STL file " << filename << endl;
		return false;
	}

	vector<string> lines;
	string line;
	while (std::getline(file, line))
	{
		// Not trimming (C++ strings) and also not needed (C++ sscanf)
		lines.push_back(line);
	}

	// Parse distinct vertexes from the text, don't add duplicates. Also figure out size of the mesh.
	for (vector<string>::iterator i = lines.begin(); i != lines.end(); ++i)
	{
		int n = sscanf_s(i->c_str(), "vertex %f %f %f\n", &x, &y, &z);
		
		if (n == 3) 
		{
			if (x < xmin) xmin = x;
			if (x > xmax) xmax = x;
			if (y < ymin) ymin = y;
			if (y > ymax) ymax = y;
			if (z < zmin) zmin = z;
			if (z > zmax) zmax = z;

			bool duplicate = false;
			for (int j = 0; j < mesh.vCount; j++)
			{
				if ((mesh.v[j].x == x) && (mesh.v[j].y == y) && (mesh.v[j].z == z)) duplicate = true;
			}

			if (!duplicate) {
				// Not duplicate - add the vertex
				mesh.v[mesh.vCount].x = x;
				mesh.v[mesh.vCount].y = y;
				mesh.v[mesh.vCount].z = z;
				mesh.v[mesh.vCount].xa = 0;
				mesh.v[mesh.vCount].ya = 0;
				// Apply some color as the file has none
				mesh.v[mesh.vCount].r = 1.0;
				mesh.v[mesh.vCount].g = 1.0;
				mesh.v[mesh.vCount].b = 1.0;
				mesh.vCount++;
			}
		}
	}
	
	// Parse faces, which will be structured as three vertex polygons which vertex to each distinct vertex index
	for (vector<string>::iterator i = lines.begin(); i != lines.end(); ++i)	{
		int n = sscanf_s(i->c_str(), "vertex %f %f %f\n", &x, &y, &z);
		if (n == 3) {
			int va = FindVertex(x, y, z, mesh);
			if (va >= 0) {
				vertexIndexes.push_back(va);
			} else {
				cout << "Couldn't find vertex!" << endl;
			}
		}
	}
	
	for (unsigned int i = 0; i < vertexIndexes.size() / 3; i++) {		
		// Add face
		mesh.f[mesh.fCount].a = vertexIndexes[i * 3];
		mesh.f[mesh.fCount].b = vertexIndexes[i * 3 + 1];
		mesh.f[mesh.fCount].c = vertexIndexes[i * 3 + 2];
		mesh.fCount++;
	}

	// Report
	cout << "Read " << mesh.vCount << " unique vertexes in " << mesh.fCount << " faces." << std::endl;
	printf("Mesh dimensions: X: %.2f to %.2f, Y: %.2f to %.2f, Z: %.2f to %.2f\n", xmin, xmax, ymin, ymax, zmin, zmax);

	return true;
}

/// <summary>
/// MoveVertexes moves vertexes around in space by (x,y,z) offset
/// </summary>
void MoveVertexes(tFloat x, tFloat y, tFloat z, Mesh &mesh) {
	for (int i = 0; i < mesh.vCount; i++) {
		mesh.v[i].x += x;
		mesh.v[i].y += y;
		mesh.v[i].z += z;
	}
}

/// <summary>
/// RotateVertexes rotates vertexes by angle specified for each axis in radians
/// </summary>
void RotateVertexes(tFloat xa, tFloat ya, tFloat za, Mesh &mesh) {

	// Rotate each vertex with help of a rotation matrix
	for (int i = 0; i < mesh.vCount; i++) {
		tFloat x = mesh.v[i].x;
		tFloat y = mesh.v[i].y;
		tFloat z = mesh.v[i].z;

		// x angle
		tFloat tx = x;
		tFloat ty = cos(xa) * y - sin(xa) * z;
		tFloat tz = sin(xa) * y + cos(xa) * z;

		// y angle
		tFloat tx2 = cos(ya) * tx - sin(ya) * ty;
		tFloat ty2 = sin(ya) * tx + cos(ya) * ty;
		tx = tx2;
		ty = ty2;

		// z angle
		tx2 = cos(za) * tx + sin(za) * tz;
		tFloat tz2 = -sin(za) * tx + cos(za) * tz;
		tx = tx2;
		tz = tz2;

		mesh.v[i].x = tx;
		mesh.v[i].y = ty;
		mesh.v[i].z = tz;
	}
}

/// <summary>
/// ScaleVertexes scales size of a mesh by a factor for each axis
/// </summary>
void ScaleVertexes(tFloat x, tFloat y, tFloat z, Mesh &mesh) {
	for (int i = 0; i < mesh.vCount; i++) {
		mesh.v[i].x *= x;
		mesh.v[i].y *= y;
		mesh.v[i].z *= z;
	}
}

/// <summary>
/// CopyVertexes makes a copy of a mesh
/// </summary>
void CopyMesh(Mesh &source, Mesh &target) {
	memcpy(&target, &source, sizeof(source));
}

/// <summary>
/// InsidePolygon checks if coordinates (a,b) are inside or outside a 2D polygon defined by a list of coordinates
/// </summary>
bool InsidePolygon(tFloat a, tFloat b, tFloat corners[3][2]) {
	bool inside = false;
	int size = 3;

	// Go through every line of the polygon
	int j = size - 1;
	for (int i = 0; i < size; i++) {
		tFloat ai = corners[i][0];
		tFloat bi = corners[i][1];
		tFloat aj = corners[j][0];
		tFloat bj = corners[j][1];

		// Calculate if intersecting and flip result bool till we reach a final result from all lines
		if (((bi > b) != (bj > b)) && (a < (aj - ai) * (b - bi) / (bj - bi) + ai)) inside = !inside;
		
		j = i;
	}

	return inside;
}

/// <summary>
/// Limit a value between 0.0 and 1.0
/// </summary>
tFloat Limit(tFloat value) {
	if (value > 1.0) return 1.0;
	if (value < 0.0) return 0.0;
	return value;
}

/// <summary>
/// Container for data passed for render thread (SDL thread)
/// </summary>
typedef struct
{
	int y1;
	int y2;
	Mesh mesh;
	Vertex light;
} RenderBag;

/// <summary>
/// Render rows between y1 and y2 to an image. Parameters passed as void* because this is a thread.
/// </summary>
#ifdef USE_PTHREADS
void *Render(void* data) {
#else
DWORD WINAPI Render(LPVOID data) {
#endif
	// Unpack parameters
	RenderBag* bag = (RenderBag*)data;
	int y1 = bag->y1;
	int y2 = bag->y2;
	Mesh *mesh = &(bag->mesh);
	Vertex light = bag->light;
	
	// Array of pixels where to render to (32bpp RGBA)
	Uint32* pixels = (Uint32*)screenSurface->pixels;

	for (tFloat y = (tFloat)y1; y < (tFloat)y2; y++) {
		for (tFloat x = 0; x < width; x++) {

			// Background is a gradient
			//pixels[((int)y * screenSurface->w + (int)x)] = (Uint32)(y / height * 64);
			Pixel foremostPixel;
			foremostPixel.r = 0;
			foremostPixel.g = 0;
			foremostPixel.b = y / height / 4;
			foremostPixel.depth = 10000;
			
			// Angle of the camera ray. They open up outwards from a single 3d-vertex (0,0,0) at the center of the image, which creates illusion of perspective.
			// Also convert FOV in degrees to radians
			tFloat cax = (x - width / 2) * fov * 0.01745329f / height;
			tFloat cay = (y - height / 2) * fov * 0.01745329f / height;			

			// Loop through every face in the scene, or mesh
			for (int f = 0; f < mesh->fCount; f++) {

				// If the camera ray is within a 2D polygon formed by the angles of the 3D polygon towards the camera (0,0,0), we found a pixel to fill with some color (not background)
				tFloat polygonPoints[3][2] = {
					{mesh->v[mesh->f[f].a].xa, mesh->v[mesh->f[f].a].ya},
					{mesh->v[mesh->f[f].b].xa, mesh->v[mesh->f[f].b].ya},
					{mesh->v[mesh->f[f].c].xa, mesh->v[mesh->f[f].c].ya}
				};
				
				if (InsidePolygon(cax, cay, polygonPoints)) {

					// Calculate at with depth the polygon we just hit is at (for depth sorting, as we might hit more polygons at different depths with the same camera ray)
					tFloat avgdepth = (mesh->v[mesh->f[f].a].z + mesh->v[mesh->f[f].b].z + mesh->v[mesh->f[f].c].z) / 3;

					if (avgdepth < foremostPixel.depth) {

						// Calculate angle between the omni light source and the polygon
						tFloat laxa = -atan2(mesh->v[mesh->f[f].a].x - light.x, mesh->v[mesh->f[f].a].z - light.z);
						tFloat laya = -atan2(mesh->v[mesh->f[f].a].y - light.y, mesh->v[mesh->f[f].a].z - light.z);
						tFloat laxb = -atan2(mesh->v[mesh->f[f].b].x - light.x, mesh->v[mesh->f[f].b].z - light.z);
						tFloat layb = -atan2(mesh->v[mesh->f[f].b].y - light.y, mesh->v[mesh->f[f].b].z - light.z);
						tFloat laxc = -atan2(mesh->v[mesh->f[f].c].x - light.x, mesh->v[mesh->f[f].c].z - light.z);
						tFloat layc = -atan2(mesh->v[mesh->f[f].c].y - light.y, mesh->v[mesh->f[f].c].z - light.z);

						// Average one brightness factor from the two angle pairs per vertex
						tFloat laa = sqrt(pow(laxa, 2.0f) + pow(laya, 2.0f));
						tFloat lab = sqrt(pow(laxb, 2.0f) + pow(layb, 2.0f));
						tFloat lac = sqrt(pow(laxc, 2.0f) + pow(layc, 2.0f));

						// Calculate colour for each vertex for gouraud shading, taking into account angle to light source and diffuse contrast
						tFloat ra = mesh->v[mesh->f[f].a].r / laa * light.r * contrastDiffuse;
						tFloat ga = mesh->v[mesh->f[f].a].g / laa * light.g * contrastDiffuse;
						tFloat ba = mesh->v[mesh->f[f].a].b / laa * light.b * contrastDiffuse;

						tFloat rb = mesh->v[mesh->f[f].b].r / lab * light.r * contrastDiffuse;
						tFloat gb = mesh->v[mesh->f[f].b].g / lab * light.g * contrastDiffuse;
						tFloat bb = mesh->v[mesh->f[f].b].b / lab * light.b * contrastDiffuse;

						tFloat rc = mesh->v[mesh->f[f].c].r / lac * light.r * contrastDiffuse;
						tFloat gc = mesh->v[mesh->f[f].c].g / lac * light.g * contrastDiffuse;
						tFloat bc = mesh->v[mesh->f[f].c].b / lac * light.b * contrastDiffuse;

						// Distance along surface from vertexes to where the camera ray hit the surface
						tFloat da = sqrt(pow(laxa - cax, 2.0f) + pow(laya - cay, 2.0f)) / contrastSpecular;
						tFloat db = sqrt(pow(laxb - cax, 2.0f) + pow(layb - cay, 2.0f)) / contrastSpecular;
						tFloat dc = sqrt(pow(laxc - cax, 2.0f) + pow(layc - cay, 2.0f)) / contrastSpecular;

						// Let color be determined by each corner of polygon and how far the camera ray is from those corners
						// Also apply ambient light level
						tFloat r = ((ra / da) + (rb / db) + (rc / dc)) + ambient;
						tFloat g = ((ga / da) + (gb / db) + (gc / dc)) + ambient;
						tFloat b = ((ba / da) + (bb / db) + (bc / dc)) + ambient;

						// Keep color below 1.0 to avoid clitches when going overbright
						foremostPixel = Pixel{ avgdepth, Limit(r), Limit(g), Limit(b) };
					}
				}
			}

			// We have now determined the foremost face that hit the camera ray, and therefore the color of the pixel, if not background
			// Convert color from three 0.0 to 1.0 values into 32-bit RGBA
			pixels[((int)y * screenSurface->w + (int)x)] = ((Uint32)(foremostPixel.r * 255) << 16) | ((Uint32)(foremostPixel.g * 255) << 8) | (Uint32)(foremostPixel.b * 255);
		}
#ifdef SHOWPROGRESS
		SDL_UpdateWindowSurface(window); // Update window for each row to show thread process independently
#endif
	}
	
	return 0;
}

// Engine thread, which keeps animating and rendering the scene, updating the frames to window
void Engine(Mesh &mesh) {
	tFloat frame = 0;

#ifdef USE_PTHREADS
	pthread_t threads[THREADS];
#else
	HANDLE threads[THREADS];
#endif

	RenderBag bags[THREADS];

	Uint32 engineStart = SDL_GetTicks();
	for (;;) {
		printf("Frame %.0f", frame / speed);

		// Transform vertexes to different position, angle etc. for each frame to archieve animation
		Mesh transformed;
		CopyMesh(mesh, transformed);
		RotateVertexes(frame / 20, frame / 25, frame / 30, transformed);
		MoveVertexes(0, 0, 200 + sin(frame / 10) * 50, transformed);

		// A pre-calculation per each frame that speeds up hugely, as two simple atan2's of each vertex are requested so often in rendering
		for (int i = 0; i < transformed.vCount; i++) {
			transformed.v[i].xa = atan2(transformed.v[i].x, transformed.v[i].z);
			transformed.v[i].ya = atan2(transformed.v[i].y, transformed.v[i].z);
		}

		// Light sources location and color is also animated
		Vertex light;
		light.x = sin(frame / 10.0f) * 160.0f;
		light.y = sin(frame / 7.0f) * 160.0f;
		light.z = -500.0f;
		light.r = (sin(frame / 3.0f) + (tFloat)M_PI) / ((tFloat)M_PI * 2.0f);
		light.g = (sin(frame / 10.0f) + (tFloat)M_PI) / ((tFloat)M_PI * 2.0f);
		light.b = (sin(frame / 5.0f) + (tFloat)M_PI) / ((tFloat)M_PI * 2.0f);
		light.xa = 0;
		light.ya = 0;

		// Render in threads, dividing the image vertically to partitions. Big speed up with multi core processors.
		Uint32 frameStart = SDL_GetTicks(); // Begin taking time
		for (int t = 0; t < THREADS; t++) {
			int y1 = height / THREADS * t;
			int y2 = height / THREADS * (t + 1);

			bags[t].y1 = y1;
			bags[t].y2 = y2;
			bags[t].mesh = transformed;
			bags[t].light = light;

			string threadName = "Render" + to_string(t);

#ifdef USE_PTHREADS
			int ret = pthread_create(&threads[t], NULL, Render, &bags[t]);
#else
			DWORD myThreadID;
			threads[t] = CreateThread(0, 0, Render, &bags[t], 0, &myThreadID);
#endif
		}
		for (int t = 0; t < THREADS; t++) {
#ifdef USE_PTHREADS
			int ret = pthread_join(threads[t], NULL);
#else
			WaitForSingleObject(threads[t], 1000);
			CloseHandle(threads[t]);
#endif
		}
		Uint32 frameDuration = SDL_GetTicks() - frameStart; // End taking time
		tFloat avgFPS = (frame / speed) / (SDL_GetTicks() - engineStart) * 1000;
		printf(" in %d ms. FPS: %.2f\n", frameDuration, avgFPS);

		// Copy the rendered image on the window
		SDL_UpdateWindowSurface(window);

		frame += speed;

		// Handle SDL events
		SDL_Event input;
		while (SDL_PollEvent(&input) > 0)
		{
			// Exit at quit request or a key press
			if ((input.type == SDL_QUIT) || (input.type == SDL_KEYDOWN)) return;
		}
	}
}

/// <summary>
/// Main function
/// </summary>
int main(int argc, char* args[])
{
	Mesh mesh;

	// Load 3D mesh data
	// ball.stl: X: -1.115109 to 1.115109, Y: -1.115109 to 1.115109, Z: -1.115109 to 1.115109, 144 vertexes in 284 faces.
	// duck2.stl: X: -1144.435059 to 672.233276, Y: -2355.597900 to 1428.072998, Z: -1920.546021 to 339.157715, 268 vertexes in 516 faces.
	if (!LoadSTL(filename, mesh)) {
		return 0;
	}
	ScaleVertexes(filescale, filescale, filescale, mesh);

	// Initialize SDL to show graphics on window
    if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
    }
    else
    {
		string title = "Tracer (" + to_string(width) + "x" + to_string(height) + ", " + filename + ")";
        window = SDL_CreateWindow(title.c_str(), SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height, SDL_WINDOW_SHOWN);
        if (window == NULL)
        {
            printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
        }
        else
        {
            screenSurface = SDL_GetWindowSurface(window);

			// SDL initialization is a success, start main engine thread
			Engine(mesh);
        }
    }

    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
