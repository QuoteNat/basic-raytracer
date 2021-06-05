#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

enum LightType{AMBIENT, POINT, DIRECTIONAL};

/**
 * @brief Defines an RGB value
 *
 */
struct COLOR {
    arma::Row<double> rgb;
};

/**
 * @brief Defines a sphere object in a 3d space.
 *
 */
struct SPHERE {
    std::string name;
    // the center of the sphere in 3d space
    arma::Row<double> center = arma::rowvec({0, 0, 0});
    // the radius of the sphere
    double radius = 0;
    // the color of the sphere
    struct COLOR color;
    // how "shiny" the sphere is
    int specular = -1;
    // how reflective the sphere is from 0 to 1.
    double reflective = 0;
};

struct LIGHT {
    LightType type;
    double intensity = 0.0;
    arma::Row<double> position = arma::rowvec({0, 0, 0});
    arma::Row<double> direction = arma::rowvec({0, 0, 0});
};

arma::Row<double> CanvasToViewport(double x, double y, double Cw, double Ch,  double Vw, double Vh, double d);
COLOR TraceRay(arma::Row<double> O, arma::Row<double> D, double t_min, double t_max, std::vector<SPHERE>& spheres, std::vector<LIGHT>& lights, COLOR background, int recursionDepth);

#endif // FUNCTIONS_H_INCLUDED
