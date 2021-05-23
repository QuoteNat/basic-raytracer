#include <SDL.h>
#include <iostream>
#include <armadillo>
#include <vector>
#include "functions.h"

//Screen dimension constants
const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 640;
const double VIEWPORT_WIDTH = 1;
const double VIEWPORT_HEIGHT = 1;
const double VIEWPORT_DISTANCE = 1;
const COLOR BACKGROUND = {{0, 0, 0}};

// prototypes
void initSpheres(std::vector<SPHERE> &);
void initLights(std::vector<LIGHT> &);

int main () {
    // variables
    SDL_Event event;
    bool quit = false;
    // origin of the camera
    arma::Row<double> O = arma::rowvec({0, 0, 0});

    // spheres
    std::vector<SPHERE> spheres;
    initSpheres(spheres);

    // lights
    std::vector<LIGHT> lights;
    initLights(lights);
    // create window and renderer
    SDL_Window *window;
    SDL_Renderer *renderer;
    SDL_CreateWindowAndRenderer(SCREEN_WIDTH, SCREEN_HEIGHT, 0, &window, &renderer);


    // render the image
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);
    for (int x = -SCREEN_WIDTH/2; x <= SCREEN_WIDTH/2; x++) {
        for (int y = -SCREEN_HEIGHT/2; y <= SCREEN_HEIGHT/2; y++) {
            //std::cout << "Rendering (" << x << ", " << y << ")\n";
            // get distance
            arma::Row<double> D = CanvasToViewport(x, y, SCREEN_WIDTH, SCREEN_HEIGHT, VIEWPORT_WIDTH, VIEWPORT_HEIGHT, VIEWPORT_DISTANCE);
            // ray trace
            COLOR color = TraceRay(O, D, 1, INFINITY, spheres, lights, BACKGROUND);
            // draw the point
            SDL_SetRenderDrawColor(renderer, color.rgb[0], color.rgb[1], color.rgb[2], 255);
            SDL_RenderDrawPoint(renderer, (SCREEN_WIDTH/2) + x, (SCREEN_HEIGHT / 2) - y);
            //SDL_RenderPresent(renderer);
        }
    }
    std::cout << "Done\n";

    // keep window open until it is quit
    while (!quit) {
        SDL_WaitEvent(&event);

        switch(event.type) {
            case SDL_QUIT:
                quit = true;
                break;
        }
        SDL_RenderPresent(renderer);

    }
    //cleanup
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
}

/**
 * @brief Initializes the spheres in spheres.
 *
 * @param spheres A vector of spheres.
 */
void initSpheres(std::vector<SPHERE> &spheres) {
    SPHERE temp, temp2, temp3;
    // red
    temp.center = {0, -1, 3};
    temp.radius = 1;
    temp.color = {255, 0, 0};
    temp.name = "Red";
    spheres.push_back(temp);
    // blue
    temp2.center = {2, 0, 4};
    temp2.radius = 1;
    temp2.color = {0, 0, 255};
    temp2.name = "Blue";
    spheres.push_back(temp2);

    // green
    temp3.center = {-2, 0, 4};
    temp3.radius = 1;
    temp3.color = {0, 255, 0};
    temp3.name = "Green";
    spheres.push_back(temp3);

    // yellow
    SPHERE temp4;
    temp4.center = {0, -5001, 0};
    temp4.radius = 5000;
    temp4.color = {255, 255, 0};
    spheres.push_back(temp4);

    std::cout << spheres.size();
}

void initLights(std::vector<LIGHT> &lights) {
    // ambient light
    LIGHT ambientLight = {AMBIENT, .2};
    lights.push_back(ambientLight);

    // point light
    LIGHT pointLight;
    pointLight.type = POINT;
    pointLight.intensity = 0.6;
    pointLight.position = {2, 1, 0};
    lights.push_back(pointLight);

    // directional light
    LIGHT directionalLight;
    directionalLight.type = DIRECTIONAL;
    directionalLight.intensity = 0.2;
    directionalLight.direction = {1, 4, 4};
    lights.push_back(directionalLight);
}

