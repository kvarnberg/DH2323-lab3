/*
 * SDL2 skeleton for lab assignments 1–3 of the KTH course DH2323,
 * Computer Graphics and Interaction (and also SU DA3001, Datalogi II).
 *
 * See README.md for details.
 */

#include <iostream>
#include </usr/local/include/glm/glm.hpp>
#include <vector>
#include "SDL2Auxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::ivec2;
using glm::mat3;
using glm::vec2;
using glm::vec3;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
int t;
vector<Triangle> triangles;
SDL2Aux *sdlAux;

const float focalLength = 500;
vec3 cameraPos(0, 0, -3.001);

mat3 R = mat3(vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1)); // rotation matrix for camera
float yaw = 0.0f;
float pi = 3.14159265359;

vec3 currentColor;

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
struct Pixel
{
  int x;
  int y;
  float zinv;
  vec3 pos3d;
};

struct Vertex
{
  vec3 position;
};

vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 14.1f * vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f * vec3(1, 1, 1);

vec3 currentNormal;
vec3 currentReflectance;

// --------------------------------------------------------
// FUNCTION DECLARATIONS

void Draw();
void Rotate();
void Update();
// void Interpolate(ivec2 a, ivec2 b, vector<ivec2> &result);
void VertexShader(const vec3 &v, Pixel &p);
void DrawLineSDL(Pixel a, Pixel b, vec3 color);
void DrawPolygonEdges(const vector<Vertex> &vertices);
void ComputePolygonRows(const vector<Pixel> &vertexPixels,
                        vector<Pixel> &leftPixels,
                        vector<Pixel> &rightPixels);
void DrawPolygonRows(const vector<Pixel> &leftPixels,
                     const vector<Pixel> &rightPixels);
void DrawPolygon(const vector<Vertex> &vertices);
void Interpolate(Pixel a, Pixel b, vector<Pixel> &result);
void PixelShader(const Pixel &p);

// --------------------------------------------------------
// FUNCTION DEFINITIONS

int main(int argc, char *argv[])
{
  LoadTestModel(triangles);
  sdlAux = new SDL2Aux(SCREEN_WIDTH, SCREEN_HEIGHT);
  t = SDL_GetTicks(); // Set start value for timer.

  while (!sdlAux->quitEvent())
  {
    Update();
    Draw();
  }

  sdlAux->saveBMP("screenshot.bmp");
  return 0;
}

void Rotate()
{
  R = mat3(glm::cos(yaw), 0, glm::sin(yaw), 0, 1, 0, -glm::sin(yaw), 0, glm::cos(yaw));
}

void Update()
{
  // Compute frame time:
  int t2 = SDL_GetTicks();
  float dt = float(t2 - t); // delta time
  t = t2;
  cout << "Render time: " << dt << " ms." << endl;

  const Uint8 *state = SDL_GetKeyboardState(NULL);
  float anglechange = ((2.0f * pi) / 1000.0f) * dt;

  // from instructions: vectors representing the axes directions that we update R for
  vec3 right(R[0][0], R[0][1], R[0][2]);   // x axis
  vec3 down(R[1][0], R[1][1], R[1][2]);    // y axis
  vec3 forward(R[2][0], R[2][1], R[2][2]); // z axis

  // copied the code from lab2 for these interactions

  if (state[SDL_SCANCODE_UP])
  {
    cameraPos += anglechange * forward;
  }
  if (state[SDL_SCANCODE_DOWN])
  {
    cameraPos -= anglechange * forward;
  }
  if (state[SDL_SCANCODE_RIGHT])
  {
    cameraPos += anglechange * right;
  }
  if (state[SDL_SCANCODE_LEFT])
  {
    cameraPos -= anglechange * right;
  }
  if (state[SDL_SCANCODE_Z])
  {
    yaw += anglechange;
    Rotate();
  }
  if (state[SDL_SCANCODE_X])
  {
    yaw -= anglechange;
    Rotate();
  }

  if (state[SDL_SCANCODE_W])
  {
    lightPos += anglechange * forward;
  }

  if (state[SDL_SCANCODE_S])
  {
    lightPos -= anglechange * forward;
  }

  if (state[SDL_SCANCODE_A])
  {
    lightPos += anglechange * right;
  }

  if (state[SDL_SCANCODE_D])
  {
    lightPos -= anglechange * right;
  }

  if (state[SDL_SCANCODE_Q])
  {
    lightPos += anglechange * down;
  }

  if (state[SDL_SCANCODE_E])
  {
    lightPos -= anglechange * down;
  }
}

void Draw()
{
  sdlAux->clearPixels();

  for (int y = 0; y < SCREEN_HEIGHT; ++y)
    for (int x = 0; x < SCREEN_WIDTH; ++x)
      depthBuffer[y][x] = 0;

  for (int i = 0; i < triangles.size(); ++i)
  {
    vector<Vertex> vertices(3);
    vertices[0].position = triangles[i].v0;
    vertices[1].position = triangles[i].v1;
    vertices[2].position = triangles[i].v2;
    currentColor = triangles[i].color;
    currentReflectance = triangles[i].color;
    currentNormal = triangles[i].normal;
    DrawPolygon(vertices);
  }
  sdlAux->render();
}

void VertexShader(const Vertex &v, Pixel &p)
{
  // total transformation, v is point in world. equation 8, P' = (P - C) R

  vec3 transform = (v.position - cameraPos) * R;
  p.zinv = 1.0f / transform.z;

  // equations 3 and 4, x = f X/Z+W/2, y = f Y/Z+H/2

  p.x = focalLength * (transform.x / transform.z) + (SCREEN_WIDTH / 2);
  p.y = focalLength * (transform.y / transform.z) + (SCREEN_HEIGHT / 2);
  p.pos3d = v.position;
};

void DrawLineSDL(Pixel a, Pixel b, vec3 color)
{
  // we first need to know how many pixels the line should consist of. a and b is start and end of segment.
  // Depending on whether the line is mostly horizontal or vertical we will use one pixel per column or one pixel per row
  Pixel delta;
  delta.x = glm::abs(b.x - a.x);
  delta.y = glm::abs(b.y - a.y);
  int pixels = glm::max(delta.x, delta.y) + 1;

  // get the pixel positions of the line by calling the Interpolation function, which then puts pixel values
  // onto vector Pixel line
  vector<Pixel> line(pixels);
  Interpolate(a, b, line);
  for (int i = 0; i < line.size(); ++i)
  {
    // bounds checking, this is to prohibit it from crashing when moving camera
    if (line[i].x >= 0 && line[i].x < SCREEN_WIDTH && line[i].y >= 0 && line[i].y < SCREEN_HEIGHT)
    {
      PixelShader(line[i]);
    }
  }
};

void DrawPolygonEdges(const vector<Vertex> &vertices)
{
  int V = vertices.size();
  // Transform each vertex from 3D world position to 2D image position:
  vector<Pixel> projectedVertices(V);
  for (int i = 0; i < V; ++i)
  {
    VertexShader(vertices[i], projectedVertices[i]);
  }
  // Loop over all vertices and draw the edge from it to the next vertex:
  for (int i = 0; i < V; ++i)
  {
    int j = (i + 1) % V; // The next vertex
    vec3 color(1, 1, 1);
    DrawLineSDL(projectedVertices[i], projectedVertices[j],
                color);
  }
}

void ComputePolygonRows(const vector<Pixel> &vertexPixels,
                        vector<Pixel> &leftPixels,
                        vector<Pixel> &rightPixels)
{
  // vertexPixels looks like this: vertexPixels[0] = ivec2(10, 5), vertexPixels[1] = ivec2( 5,10), vertexPixels[2] = ivec2(15,15)
  // 1. Find max and min y-value of the polygon
  // and compute the number of rows it occupies. max-min +1
  int min = std::min({vertexPixels[0].y, vertexPixels[1].y, vertexPixels[2].y});
  int max = std::max({vertexPixels[0].y, vertexPixels[1].y, vertexPixels[2].y});

  int ROWS = (max - min) + 1;

  // 2. Resize leftPixels and rightPixels
  // so that they have an element for each row.
  leftPixels.resize(ROWS);
  rightPixels.resize(ROWS);

  // 3. Initialize the x-coordinates in leftPixels to some really large value and the x-coordinatesin
  // rightPixels to some really small value. I give y-coordinates the min y-value of the polygon.
  for (int i = 0; i < ROWS; ++i)
  {
    leftPixels[i].x = +numeric_limits<int>::max();
    leftPixels[i].y = i + min;
    rightPixels[i].x = -numeric_limits<int>::max();
    rightPixels[i].y = i + min;
  }

  // 4. Loop through all edges of the polygon and use
  // linear interpolation to find the x-coordinate for
  // each row it occupies. Update the corresponding
  // values in rightPixels and leftPixels.

  for (int i = 0; i < vertexPixels.size(); ++i)
  {
    // this is the next vertice, same as in DrawPolygonEdges
    int i_next = (i + 1) % vertexPixels.size();
    // the line/vector between the two vertices (as ivec2 vector)
    ivec2 d_pixels;
    d_pixels.x = glm::abs(vertexPixels[i_next].x - vertexPixels[i].x);
    d_pixels.y = glm::abs(vertexPixels[i_next].y - vertexPixels[i].y);
    // int pixels = glm::max(d_pixels.x,d_pixels) + 1;
    int pixels = glm::max(d_pixels.x, d_pixels.y) + 1;
    // make an edge(or vector) and have it as the size of pixels
    vector<Pixel> edge(pixels);
    // making edges, by interpolating between two vertices and returning resulting interpolated edge
    Interpolate(vertexPixels[i], vertexPixels[i_next], edge);

    // From instructions: ... for each y-coordinate of this line we check
    // the corresponding elements in the left and right arrays.
    // If the current x-value in the left array at that place is
    // larger than the x-value of the line we replace it.
    // If the current x-value in the right array at that place is
    // smaller than the x-value of the line we replace it.
    // After we have looped through all edges we then have the
    // smallest x-value for each row in the left array and the largest value for each row in the right array.
    for (int j = 0; j < edge.size(); ++j)
    {
      for (int k = 0; k < ROWS; ++k)
      {
        if (leftPixels[k].y == edge[j].y)
        {
          if (leftPixels[k].x > edge[j].x)
          {
            leftPixels[k].x = edge[j].x;
            leftPixels[k].pos3d = edge[j].pos3d;
            leftPixels[k].zinv = edge[j].zinv;
          }
          if (rightPixels[k].x <= edge[j].x)
          {
            rightPixels[k].x = edge[j].x;
            rightPixels[k].pos3d = edge[j].pos3d;
            rightPixels[k].zinv = edge[j].zinv;
          }
        }
      }
    };
  }
}

void DrawPolygonRows(const vector<Pixel> &leftPixels,
                     const vector<Pixel> &rightPixels)
{
  // This function should call PutPixelSDL for each pixel between the start and end for each row.
  for (int i = 0; i < leftPixels.size(); ++i)
  {
    // make use of DrawLineSDL to send in vectors that can then be looped through and put on screen
    DrawLineSDL(leftPixels[i], rightPixels[i], currentColor);
  }
};

void DrawPolygon(const vector<Vertex> &vertices)
{
  int V = vertices.size();
  vector<Pixel> vertexPixels(V);
  for (int i = 0; i < V; ++i)
    VertexShader(vertices[i], vertexPixels[i]);
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  DrawPolygonRows(leftPixels, rightPixels);
}

void Interpolate(Pixel a, Pixel b, vector<Pixel> &result)
{
  // Previously we just interpolated the position between the vertices when we a triangle. Now we also want to
  // interpolate the inverse depth. Thus, you need to implement an Interpolation function that interpolates our Pixel
  // struct linearly instead of glm::ivec2:

  // inverse depth 1/z
  int N = result.size();

  // steps, same as in previous interpolation
  float step_x = (b.x - a.x) / float(max(N - 1, 1));
  float step_y = (b.y - a.y) / float(max(N - 1, 1));
  float step_zinv = (b.zinv - a.zinv) / float(max(N - 1, 1));
  // vec3 step_illumination = (b.illumination - a.illumination) / float(max(N - 1, 1));

  // these, then divide current.posd3 with zinv further down.
  a.pos3d = a.pos3d * a.zinv;
  b.pos3d = b.pos3d * b.zinv;

  vec3 step_pos3d = (b.pos3d - a.pos3d) / float(max(N - 1, 1));
  // same as before but with pixel
  Pixel current(a);
  for (int i = 0; i < N; ++i)
  {
    // interpolation as done in lab1, (1-t)a + t*b, or startValue + fraction * (endValue - startValue), t ∈ [0, 1]
    current.x = a.x + (float)i * step_x;
    current.y = a.y + (float)i * step_y;
    current.zinv = a.zinv + (float)i * step_zinv;
    current.pos3d = (a.pos3d + (float)i * step_pos3d) / current.zinv;
    result[i] = current;
  };
};

void PixelShader(const Pixel &p)
{
  int x = p.x;
  int y = p.y;
  // check for depth
  if (p.zinv > depthBuffer[y][x])
  {
    depthBuffer[y][x] = p.zinv;

    // same as in lab2
    vec3 n = currentNormal;
    vec3 r = lightPos - p.pos3d;
    vec3 rnorm = glm::normalize(r);
    vec3 D;

    float r_length = glm::length(r);

    D = vec3(lightPower * glm::max(glm::dot(rnorm, n), 0.0f)) / (4.0f * glm::pow(r_length, 2.0f) * pi);
    vec3 illumination = currentReflectance * (D + indirectLightPowerPerArea);
    sdlAux->putPixel(x, y, currentColor * illumination);
  }
};
