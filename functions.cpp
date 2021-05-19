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
COLOR TraceRay(arma::Row<double> O, arma::Row<double> D, double t_min, double t_max, SPHERE spheres[3], COLOR background) {
    double closest_t = INFINITY;
    SPHERE closest_sphere;
    closest_sphere.color = background;
    for(int i=0; i < 3; i++) {
        SPHERE sphere = spheres[i];
        double intersects[2];
        IntersectRaySphere(O, D, sphere, intersects);
        if (intersects[0] > t_min && intersects[0] < t_max && intersects[0] < closest_t) {
            closest_t = intersects[0];
            closest_sphere = sphere;
        }
        if (intersects[1] > t_min && intersects[1] < t_max && intersects[1] < closest_t) {
            closest_t = intersects[1];
            closest_sphere = sphere;
        }
    }
    return closest_sphere.color;
}
