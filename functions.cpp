#include <armadillo>
#include <cmath>
#include <iterator>
#include <iostream>
#include "functions.h"

/**
 * @brief Determines which pixel on the viewport corresponds with the canvas pixel (x, y).
 *
 * @param x x position on the canvas
 * @param y y postition on the canvas
 * @param Cw screen width
 * @param Ch screen height
 * @param Vw viewport width
 * @param Vh viewport height
 * @param d distance from canvas to viewport
 * @return arma::Row<double>
 */
arma::Row<double> CanvasToViewport(double x, double y, double Cw, double Ch,  double Vw, double Vh, double d)  {
    return {x*Vw/Cw, y*Vh/Ch, d};
}

/**
 * @brief Returns the intersections between a ray and a sphere, or [inf, inf] if there are none
 *
 * @param O Origin of the camera
 * @param D Direction vector for the ray
 * @param sphere Sphere object
 * @param intersects Array of 2 ints for the intersections to be written to. INFINITY is written to both indexes if there are no intersections.
 */
void IntersectRaySphere(arma::Row<double> O, arma::Row<double> D, SPHERE sphere, double intersects[2]) {
    double r = sphere.radius;
    arma::Row<double> CO = O - sphere.center;
    double a = arma::dot(D, D);
    double b = 2 * arma::dot(CO, D);
    double c = arma::dot(CO, CO) - r*r;
    double discriminant = b*b - 4*a*c;
    if (discriminant < 0) {
        intersects[0] = INFINITY;
        intersects[1] = INFINITY;
    } else {
        intersects[0] = (-b + sqrt(discriminant)) / (2*a);
        intersects[1] = (-b - sqrt(discriminant)) / (2*a);
    }
}

double magnitude (arma::Row<double> vec) {
    return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

double ComputeLighting(std::vector<LIGHT> &lights, arma::Row<double> P, arma::Row<double> N) {
    double intensity = 0.0;
    for (int i=0; i < lights.size(); i++) {
        // if light type is ambient
        if (lights[i].type == AMBIENT) {
            // add intensity to i
            intensity += lights[i].intensity;
        } else {
            arma::Row<double> L;
            // set L to the position of the point light, or the direction of the DIRECTIONAL light.
            if (lights[i].type == POINT) {
                L = lights[i].position - P;
            } else {
                L = lights[i].direction;
            }

            // get he dot product of N dot L.
            double n_dot_l = arma::dot(N, L);
            if (n_dot_l > 0) {
                // intensity + the normalized value of n_dot_l
                intensity += lights[i].intensity * n_dot_l / (magnitude(L) * magnitude(N));
            }
        }
    }
    return intensity;
}


COLOR adjustColor(COLOR color, double intensity) {
    // change the intensity of the color value by intensity
    for (int i=0; i < 3; i++) {
        color.rgb[i] = int(color.rgb[i] * intensity);
    }
    return color;
}

/**
 * @brief Performs a raytrace and returns the color of the collision if there are any.
 *
 * @param O Origin vector.
 * @param D Direction vector.
 * @param t_min Minimum distance of the ray.
 * @param t_max Maximum distance of the ray.
 * @param spheres An array containing the SPHERE objects in the scene.
 * @param background The background color.
 * @return COLOR The color of the object the ray intersects, or background if there are no intersections.
 */
COLOR TraceRay(arma::Row<double> O, arma::Row<double> D, double t_min, double t_max, std::vector<SPHERE>& spheres, std::vector<LIGHT>& lights,COLOR background) {
    // closest intersectioin
    double closest_t = INFINITY;
    SPHERE* closest_sphere = NULL;
    SPHERE temp;
    // color for when no intersects are found/background color
    //closest_sphere.color = background;
    // for each sphere in the scene
    for(int i=0; i < spheres.size(); i++) {
        SPHERE sphere = spheres[i];
        double intersects[2];
        // check for intersects between sphere and the ray
        IntersectRaySphere(O, D, sphere, intersects);
        // if an intersect is within the acceptable range and is less than closest_t, set closest_t to the intersect and the closest_sphere to sphere.
        if (intersects[0] > t_min && intersects[0] < t_max && intersects[0] < closest_t) {
            closest_t = intersects[0];
            temp = sphere;
            closest_sphere = &temp;
        }
        if (intersects[1] > t_min && intersects[1] < t_max && intersects[1] < closest_t) {
            closest_t = intersects[1];
            temp = sphere;
            closest_sphere = &temp;
        }
    }

    if (closest_sphere == NULL) return background;

    arma::Row<double> P = O + arma::as_scalar(closest_t) * D; // Compute intersection [broken]
    arma::Row<double> N = P - closest_sphere->center; // Compute sphere normal at intersection
    N = N / magnitude(N); // normalize
    return adjustColor(closest_sphere->color, ComputeLighting(lights, P, N));
}


