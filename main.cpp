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
const COLOR BACKGROUND = {{255, 255, 255}};

// prototypes
void initSpheres(std::vector<SPHERE> &);

int main () {
    // variables
    SDL_Event event;
    bool quit = false;
    // origin of the camera
    arma::Row<double> O = arma::rowvec({0, 0, 0});

    //spheres
    std::vector<SPHERE> spheres;
    initSpheres(spheres);
    // create window and renderer
    SDL_Window *window;
    SDL_Renderer *renderer;
    SDL_CreateWindowAndRenderer(SCREEN_WIDTH, SCREEN_HEIGHT, 0, &window, &renderer);


    // render the image
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);
    for (int x = -SCREEN_WIDTH/2; x <= SCREEN_WIDTH/2; x++) {
        for (int y = -SCREEN_HEIGHT/2; y <= SCREEN_HEIGHT/2; y++) {
            // get distance
            arma::Row<double> D = CanvasToViewport(x, y, SCREEN_WIDTH, SCREEN_HEIGHT, VIEWPORT_WIDTH, VIEWPORT_HEIGHT, VIEWPORT_DISTANCE);
            // ray trace
            COLOR color = TraceRay(O, D, 1, INFINITY, spheres, BACKGROUND);
            // draw the point
            SDL_SetRenderDrawColor(renderer, color.rgb[0], color.rgb[1], color.rgb[2], 255);
            SDL_RenderDrawPoint(renderer, (SCREEN_WIDTH/2) + x, (SCREEN_HEIGHT / 2) - y);
            //SDL_RenderPresent(renderer);
        }
    }

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
    // circle 1
    temp.center = {0, -1, 3};
    temp.radius = 1;
    temp.color.rgb[0] = 255;
    temp.name = "Red";
    spheres.push_back(temp);
    // circle 2
    temp2.center = {2, 0, 4};
    temp2.radius = 1;
    temp2.color.rgb[2] = 255;
    temp2.name = "Blue";
    spheres.push_back(temp2);

    // circle 3
    temp3.center = {-2, 0, 4};
    temp3.radius = 1;
    temp3.color.rgb[1] = 255;
    temp3.name = "Green";
    spheres.push_back(temp3);

    std::cout << spheres.size();
}
