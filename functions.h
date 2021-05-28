#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

enum LightType{AMBIENT, POINT, DIRECTIONAL};

/**
 * @brief Defines an RGB value
 *
 */
struct COLOR {
    int rgb[3];
};

/**
 * @brief Defines a sphere object in a 3d space.
 *
 */
struct SPHERE {
    std::string name;
    arma::Row<double> center = arma::rowvec({0, 0, 0});
    double radius = 0;
    struct COLOR color = {{0, 0, 0}};
    int specular = -1;
};

struct LIGHT {
    LightType type;
    double intensity = 0.0;
    arma::Row<double> position = arma::rowvec({0, 0, 0});
    arma::Row<double> direction = arma::rowvec({0, 0, 0});
};

arma::Row<double> CanvasToViewport(double x, double y, double Cw, double Ch,  double Vw, double Vh, double d);
COLOR TraceRay(arma::Row<double> O, arma::Row<double> D, double t_min, double t_max, std::vector<SPHERE>& spheres, std::vector<LIGHT>& lights,COLOR background);

#endif // FUNCTIONS_H_INCLUDED
