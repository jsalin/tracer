package main

/* * * * * * * * * * * * * * * * * * * * * * * * *
 * A very simple software 3D ray tracing engine  *
 * (C) 2020 Jussi Salin - See the README.md file *
 * * * * * * * * * * * * * * * * * * * * * * * * */

import (
	"bufio"
	"fmt"
	"image"
	"image/color"
	"log"
	"math"
	"os"
	"strings"
	"sync"
	"time"

	"golang.org/x/exp/shiny/driver"
	"golang.org/x/exp/shiny/screen"
	"golang.org/x/mobile/event/key"
	"golang.org/x/mobile/event/lifecycle"
)

// Vertex structure
type Vertex struct {
	x, y, z float64 // Coordinate of vertex
	r, g, b float64 // Color of vertex
	xa, ya  float64 // Speed up, store pre-calculated angles to camera vertex for each vertex (once before frame)
}

// Face structure
type Face struct {
	a, b, c int // Index of each vertex of polygon
}

// Pixel structure
type Pixel struct {
	depth   float64 // Depth the color was calculated from
	r, g, b float64 // Color
}

// Consts (adjustable)
const threads = 16           // How many threads to use in parallel rendering
const width = 320            // Render area width in pixels
const height = 240           // Render area height in pixels
const contrast = 0.01        // A multiplier to shading, to increase or decrease contrast of the overall image
const ambient = 0.0          // Ambient brightness, between 0.0 (black is black) and 1.0 ("fullbright")
const fov = 90.0             // Field of vision (in degrees)
const filename = "duck2.stl" // File to load (the 3d mesh)
const filescale = 0.07       // Scale for the file, as some meshes can be really big or really small (use 0.07 for duck2.stl and 100.0 for ball.stl)

// LoadSTL loads ascii STL file. Works with a single mesh with three vertex faces only.
func LoadSTL(filename string) ([]Vertex, []Face) {
	var xmin float64 = 0
	var ymin float64 = 0
	var zmin float64 = 0
	var xmax float64 = 0
	var ymax float64 = 0
	var zmax float64 = 0
	var x, y, z float64

	// Initialize empty lists for vertexes and faces
	vertexes := make([]Vertex, 0)
	faces := make([]Face, 0)
	vertexIndexes := make([]int, 0)

	// Load the lines from file to a string array
	file, err := os.Open(filename)
	if err != nil {
		fmt.Println("Error opening STL file.")
		return nil, nil
	}
	var lines []string
	r := bufio.NewReader(file)
	for {
		line, err := r.ReadBytes('\n')
		if err != nil {
			break
		}

		lines = append(lines, strings.Trim(string(line), " ")) // Trim indentication some STL writers add
	}

	// Parse distinct vertexes from the text, don't add duplicates. Also figure out size of the mesh.
	for _, line := range lines {
		n, err := fmt.Sscanf(line, "vertex %f %f %f\n", &x, &y, &z)
		if err == nil && n == 3 {
			if x < xmin {
				xmin = x
			}
			if x > xmax {
				xmax = x
			}
			if y < ymin {
				ymin = y
			}
			if y > ymax {
				ymax = y
			}
			if z < zmin {
				zmin = z
			}
			if z > zmax {
				zmax = z
			}

			duplicate := false
			for _, oldVertex := range vertexes {
				if oldVertex.x == x && oldVertex.y == y && oldVertex.z == z {
					duplicate = true
				}
			}

			if duplicate == false {
				// Apply some color as the file has none
				// r := rand.Float64()
				r := 1.0
				g := 1.0
				b := 1.0
				vertexes = append(vertexes, Vertex{x, y, z, r, g, b, 0, 0})
			}
		}
	}

	// Parse faces, which will be structured as three vertex polygons which vertex to each distinct vertex index
	for _, line := range lines {
		n, err := fmt.Sscanf(strings.Trim(line, " "), "vertex %f %f %f\n", &x, &y, &z)
		if err == nil && n == 3 {
			va := FindVertex(x, y, z, &vertexes)
			if va >= 0 {
				vertexIndexes = append(vertexIndexes, va)
			} else {
				fmt.Println("Couldn't find vertex!")
			}
		}
	}
	for i := 0; i < len(vertexIndexes)/3; i++ {
		faces = append(faces, Face{vertexIndexes[i*3], vertexIndexes[i*3+1], vertexIndexes[i*3+2]})
	}

	// Report
	fmt.Println("Read", len(vertexes), "unique vertexes in", len(faces), "faces.")
	fmt.Printf("Mesh dimensions: X: %f to %f, Y: %f to %f, Z: %f to %f\n", xmin, xmax, ymin, ymax, zmin, zmax)

	return vertexes, faces
}

// FindVertex returns index of a vertex in a list by coordinates
func FindVertex(x float64, y float64, z float64, vertexes *[]Vertex) int {
	for i := 0; i < len(*vertexes); i++ {
		if ((*vertexes)[i].x == x) && ((*vertexes)[i].y == y) && ((*vertexes)[i].z == z) {
			return i
		}
	}
	return -1
}

// MoveVertexes moves vertexes around in space by (x,y,z) offset
func MoveVertexes(x float64, y float64, z float64, vertexes *[]Vertex) {
	for i := 0; i < len(*vertexes); i++ {
		(*vertexes)[i].x += x
		(*vertexes)[i].y += y
		(*vertexes)[i].z += z
	}
}

// RotateVertexes rotates vertexes by angle specified for each axis in radians
func RotateVertexes(xa float64, ya float64, za float64, vertexes *[]Vertex) {

	// Rotate each vertex with help of a rotation matrix
	for i := 0; i < len(*vertexes); i++ {
		var x = (*vertexes)[i].x
		var y = (*vertexes)[i].y
		var z = (*vertexes)[i].z

		// x angle
		var tx = x
		var ty = math.Cos(xa)*y - math.Sin(xa)*z
		var tz = math.Sin(xa)*y + math.Cos(xa)*z

		// y angle
		var tx2 = math.Cos(ya)*tx - math.Sin(ya)*ty
		var ty2 = math.Sin(ya)*tx + math.Cos(ya)*ty
		tx = tx2
		ty = ty2

		// z angle
		tx2 = math.Cos(za)*tx + math.Sin(za)*tz
		var tz2 = -math.Sin(za)*tx + math.Cos(za)*tz
		tx = tx2
		tz = tz2

		(*vertexes)[i].x = tx
		(*vertexes)[i].y = ty
		(*vertexes)[i].z = tz
	}
}

// ScaleVertexes scales size of a mesh by a factor for each axis
func ScaleVertexes(x float64, y float64, z float64, vertexes *[]Vertex) {
	for i := 0; i < len(*vertexes); i++ {
		(*vertexes)[i].x *= x
		(*vertexes)[i].y *= y
		(*vertexes)[i].z *= z
	}
}

// CopyVertexes makes a copy of a vertex array (to manipulate data while keeping the original in memory)
func CopyVertexes(source *[]Vertex) []Vertex {
	target := make([]Vertex, 0)
	for i := 0; i < len(*source); i++ {
		target = append(target, (*source)[i])
	}
	return target
}

// InsidePolygon checks if coordinates (a,b) are inside or outside a 2D polygon defined by a list of coordinates
func InsidePolygon(a float64, b float64, corners [][2]float64) bool {
	inside := false

	// Go through every line of the polygon
	j := len(corners) - 1
	for i := 0; i < len(corners); i++ {
		ai := corners[i][0]
		bi := corners[i][1]
		aj := corners[j][0]
		bj := corners[j][1]

		// Calculate if intersecting and flip result bool till we reach a final result from all lines
		if ((bi > b) != (bj > b)) && (a < (aj-ai)*(b-bi)/(bj-bi)+ai) {
			inside = !inside
		}

		j = i
	}

	return inside
}

// Limit a value between 0.0 and 1.0
func Limit(value float64) float64 {
	if value > 1.0 {
		return 1.0
	}
	if value < 0.0 {
		return 0.0
	}
	return value
}

// Render rows between y1 and y2 to an image
func Render(y1 int, y2 int, img *image.RGBA, vertexes *[]Vertex, faces *[]Face, light Vertex, wg *sync.WaitGroup) {
	defer wg.Done()

	for y := float64(y1); y < float64(y2); y++ {
		for x := float64(0); x < float64(width); x++ {

			// Background is a gradient
			img.Set(int(x), int(y), color.RGBA{0, 0, uint8(y / height * 64), 255})

			// Angle of the camera ray. They open up outwards from a single 3d-vertex (0,0,0) at the center of the image, which creates illusion of perspective.
			// Also convert FOV in degrees to radians
			cax := (x - width/2) * fov * 0.01745329 / height
			cay := (y - height/2) * fov * 0.01745329 / height

			foremostPixel := Pixel{0, 0, 0, 0}
			pixelFound := false

			// Loop through every face in the scene, or mesh
			for _, f := range *faces {

				// If the camera ray is within a 2D polygon formed by the angles of the 3D polygon towards the camera (0,0,0), we found a pixel to fill with some color (not background)
				polygonPoints := [][2]float64{{(*vertexes)[f.a].xa, (*vertexes)[f.a].ya}, {(*vertexes)[f.b].xa, (*vertexes)[f.b].ya}, {(*vertexes)[f.c].xa, (*vertexes)[f.c].ya}}
				if InsidePolygon(cax, cay, polygonPoints) == true {

					// Calculate at with depth the polygon we just hit is at (for depth sorting, as we might hit more polygons at different depths with the same camera ray)
					avgdepth := ((*vertexes)[f.a].z + (*vertexes)[f.b].z + (*vertexes)[f.c].z) / 3

					if (pixelFound == false) || (avgdepth < foremostPixel.depth) {

						// Calculate angle between the omni light source and the polygon
						laxa := math.Atan2((*vertexes)[f.a].x-light.x, (*vertexes)[f.a].z-light.z)
						laya := math.Atan2((*vertexes)[f.a].y-light.y, (*vertexes)[f.a].z-light.z)
						laxb := math.Atan2((*vertexes)[f.b].x-light.x, (*vertexes)[f.b].z-light.z)
						layb := math.Atan2((*vertexes)[f.b].y-light.y, (*vertexes)[f.b].z-light.z)
						laxc := math.Atan2((*vertexes)[f.c].x-light.x, (*vertexes)[f.c].z-light.z)
						layc := math.Atan2((*vertexes)[f.c].y-light.y, (*vertexes)[f.c].z-light.z)

						// Average one brightness factor from the two angle pairs / vertex
						// Need absolute, otherwise negative half gets cut into black color
						laa := math.Abs(laxa+laya) / 2
						lab := math.Abs(laxb+layb) / 2
						lac := math.Abs(laxc+layc) / 2

						// Calculate average color of the polygon, taking into account angle and contrast value
						ra := (*vertexes)[f.a].r / laa * light.r * contrast
						ga := (*vertexes)[f.a].g / laa * light.g * contrast
						ba := (*vertexes)[f.a].b / laa * light.b * contrast

						rb := (*vertexes)[f.b].r / lab * light.r * contrast
						gb := (*vertexes)[f.b].g / lab * light.g * contrast
						bb := (*vertexes)[f.b].b / lab * light.b * contrast

						rc := (*vertexes)[f.c].r / lac * light.r * contrast
						gc := (*vertexes)[f.c].g / lac * light.g * contrast
						bc := (*vertexes)[f.c].b / lac * light.b * contrast

						// Distance along surface from vertexes to where the camera ray hit the surface
						da := math.Sqrt(math.Pow(laxa-cax, 2)+math.Pow(laya-cay, 2)) / 3
						db := math.Sqrt(math.Pow(laxb-cax, 2)+math.Pow(layb-cay, 2)) / 3
						dc := math.Sqrt(math.Pow(laxc-cax, 2)+math.Pow(layc-cay, 2)) / 3

						// Let color be determined by each corner of polygon and how far the camera ray is from those corners
						// Also apply contrast adjustment and ambient light level
						r := ((ra/da)+(rb/db)+(rc/dc))/3 + ambient
						g := ((ga/da)+(gb/db)+(gc/dc))/3 + ambient
						b := ((ba/da)+(bb/db)+(bc/dc))/3 + ambient

						// Keep color below 1.0 to avoid clitches when going overbright
						foremostPixel = Pixel{avgdepth, Limit(r), Limit(g), Limit(b)}
						pixelFound = true
					}
				}

			}

			// We have now determined the foremost face that hit the camera ray, and therefore the color of the pixel, if not background
			if pixelFound == true {
				img.Set(int(x), int(y), color.RGBA{uint8(foremostPixel.r * 255), uint8(foremostPixel.g * 255), uint8(foremostPixel.b * 255), 255})
			}
		}

	}
}

// Engine thread, which keeps animating and rendering the scene, updating the frames to window
func Engine(s screen.Screen, w screen.Window, vertexes *[]Vertex, faces *[]Face) {
	var frame float64 = 0

	for {
		fmt.Print("Frame ", frame)

		// Transform vertexes to different position, angle etc. for each frame to archieve animation
		transformed := CopyVertexes(vertexes)
		RotateVertexes(frame/20, frame/25, frame/30, &transformed)
		MoveVertexes(0, 0, 200+math.Sin(frame/5)*20, &transformed)

		// Create a shiny buffer for the new rendered image
		size := image.Point{width, height}
		b, err := s.NewBuffer(size)
		if err != nil {
			log.Fatal(err)
		}
		defer b.Release()

		// A pre-calculation per each frame that speeds up hugely, as two simple atan2's of each vertex are requested so often in rendering
		for i := 0; i < len(transformed); i++ {
			transformed[i].xa = math.Atan2(transformed[i].x, transformed[i].z)
			transformed[i].ya = math.Atan2(transformed[i].y, transformed[i].z)
		}

		// Light sources location and color is also animated
		lightR := (math.Sin(frame/3) + math.Pi) / (math.Pi * 2)
		lightG := (math.Sin(frame/10) + math.Pi) / (math.Pi * 2)
		lightB := (math.Sin(frame/5) + math.Pi) / (math.Pi * 2)
		//lightR := 1.0
		//lightG := 1.0
		//lightB := 1.0

		light := Vertex{math.Sin(frame/10) * 160, math.Sin(frame/7) * 120, -500, lightR, lightG, lightB, 0, 0}

		// Render in threads, dividing the image vertically to partitions. Big speed up with multi core processors.
		start := time.Now() // Begin taking time
		wg := new(sync.WaitGroup)
		for t := 0; t < threads; t++ {
			wg.Add(1)
			y1 := height / threads * t
			y2 := height / threads * (t + 1)
			go Render(y1, y2, b.RGBA(), &transformed, faces, light, wg)
		}
		wg.Wait()
		duration := time.Now().Sub(start) // End taking time
		fmt.Println(" in", duration)

		// Copy the rendered image in shiny buffer to a texture and then draw it on the window
		t0, err := s.NewTexture(size)
		if err != nil {
			log.Fatal(err)
		}
		defer t0.Release()
		t0.Upload(image.Point{}, b, b.Bounds())
		w.Copy(image.Point{0, 0}, t0, t0.Bounds(), screen.Src, nil)

		// These shiny buffers have to be released manually after each frame or we leak memory at great rate
		t0.Release()
		b.Release()

		frame++
	}
}

// Main function
func main() {
	// Load 3D mesh data
	// ball.stl: X: -1.115109 to 1.115109, Y: -1.115109 to 1.115109, Z: -1.115109 to 1.115109, 144 vertexes in 284 faces.
	// duck2.stl: X: -1144.435059 to 672.233276, Y: -2355.597900 to 1428.072998, Z: -1920.546021 to 339.157715, 268 vertexes in 516 faces.
	vertexes, faces := LoadSTL(filename)
	ScaleVertexes(filescale, filescale, filescale, &vertexes)

	// Create shiny window, start drawing to it and listen to events
	driver.Main(func(s screen.Screen) {
		w, err := s.NewWindow(&screen.NewWindowOptions{
			Title:  "Tracer",
			Width:  width,
			Height: height,
		})
		if err != nil {
			fmt.Println(err)
			return
		}
		defer w.Release()

		// Start background thread, which keeps producing animated frames to the window
		go Engine(s, w, &vertexes, &faces)

		// Main event loop
		for {
			switch e := w.NextEvent().(type) {

			case lifecycle.Event:
				// Exit when window is closed
				if e.To == lifecycle.StageDead {
					return
				}

			case key.Event:
				// Exit if escape key is pressed
				if e.Code == key.CodeEscape {
					return
				}

			}
		}
	})
}
