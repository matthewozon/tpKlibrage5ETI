#ifndef C_KLIBRATION_H
#define C_KLIBRATION_H

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>

class C_klibration
{
public:
    C_klibration(std::string fileNameSceneScreen);

    //void assessParameters(void); //the method that gathers all steps of the estimation
    //void saveParameters(std::string fileNameParmeters);

private:
    //parameters to be estimated
    //double Tx, Ty, Tz;
    //double omega, phi, gamma;
    //double f, sx, sy, alpha;

    //working methods
    //void calculateLVector(void); //a method to calucultate the L vector
    ///here you must add all other method to compute the calibration.

    //contains all points of the scene in object and screen frame
    //std::vector< C_vector<double> > m_Scene;
    //std::string m_fileNameScene;
    //void loadObjectAndScreenPoints();
    //void coordinatesFromLine(std::string line, double& x, double& y, double& z, double& u, double& v);
};

#endif // C_KLIBRATION_H
