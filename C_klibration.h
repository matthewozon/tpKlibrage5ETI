#ifndef C_KLIBRATION_H
#define C_KLIBRATION_H

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>

#include <C_vector.h>
#include <C_matrix.h>

#define CX 256
#define CY 256

class C_klibration
{
public:
    C_klibration(std::string fileNameSceneScreen);
    virtual ~C_klibration();

    void assessParameters(void); //the method that gathers all steps of the estimation
    bool saveParameters(std::string fileNameParmeters);


private:
    //parameters to be estimated
    double Tx, Ty, Tz;
    double omega, phi, gamma;
    double f, sx, sy, alpha;
    double r11, r12, r13;
    double r21, r22, r23;
    double r31, r32, r33;
    double f2, beta/*f1/f2*/;

    //working methods
    C_vector<double>* calculateLVector(void); //a method to calucultate the L vector
    ///here you must add all other method to compute the calibration.
    C_matrix<double>* createMatrixA(void);
    C_vector<double>* createVectorU(void);

    C_vector<double> * calculateTzAndf2(void);
    C_matrix<double>* createMatrixB(void);
    C_vector<double>* createVectorR(void);

    //contains all points of the scene in object and screen frame
    std::vector< C_vector<double> > m_Scene;
    std::string m_fileNameScene;
    bool loadObjectAndScreenPoints();
    void coordinatesFromLine(std::string line, double& x, double& y, double& z, double& u, double& v);
};

#endif // C_KLIBRATION_H
