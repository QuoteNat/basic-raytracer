#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

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
};

arma::Row<double> CanvasToViewport(double x, double y, double Cw, double Ch,  double Vw, double Vh, double d);
COLOR TraceRay(arma::Row<double> O, arma::Row<double> D, double t_min, double t_max, SPHERE spheres[3], COLOR background);

#endif // FUNCTIONS_H_INCLUDED
