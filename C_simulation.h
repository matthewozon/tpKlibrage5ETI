#ifndef C_SIMULATION_H
#define C_SIMULATION_H

#include <C_matrix.h>
#include <C_vector.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>

class C_simulation
{
public:
    C_simulation(std::string fileNameSceneObject); /*loads the points of the scene in objetScene from a file*/

    //model the camera Mint & Mext
    void setRotationAndTranslation(double omega, double phi, double gamma, double dx, double dy, double dz);
    void setIntrinsecParameters(double f, double sx, double sy, double cx=0.0, double cy=0.0);
    //project points on screen
    void projectOntoScreen(void);
    //save simulated points
    bool saveVector(std::string fileNameSceneScreen);

    std::vector< C_vector<double> > m_objetScreen;

private:
    //matrices for the model
    C_matrix<double>* m_Mext;
    C_matrix<double>* m_Mint;

    //contains all points of the scene in object frame
    std::vector< C_vector<double> > m_objetScene;

    std::string m_fileNameScene; //name of file that contains all scene points
    bool loadScenePointsFromFile(void); /*this method streams the file "m_fileNameScene",
                                      reads the number of points (#row), stores the points
                                      in a vector of C_vectors (coordinates are seperated by space)*/
    void coordinatesFromLine(std::string line, double& x, double& y, double& z);// a tool that gives back 3 numbers (x,y,z) from a line a characters
    C_matrix<double>* setFocal(double f); //initialize the P matrix (made out of the focal distance)
    C_matrix<double>* setPixelSizeAndCenter(double sx, double sy, double cx, double cy); //the S matrix: scaling factors in x and y dimensions and the offset to the center

};

#endif // C_SIMULATION_H
