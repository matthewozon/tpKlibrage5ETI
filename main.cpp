#include <iostream>
#include <stdio.h>
#include <sstream>

#include <C_matrix.h>
#include <C_vector.h>
#include <C_simulation.h>
#include <C_klibration.h>

#define PI 3.14159265359L


int main()
{
    ///simulation of the acquisition
    std::string fileNameSceneObject = "points";
    C_simulation theScene(fileNameSceneObject);
    theScene.setRotationAndTranslation(PI/4,PI/8.0,PI/8,100.0,100.0, 1000.0);
    theScene.setIntrinsecParameters(5.0, 0.01, 0.01, CX, CY);
    theScene.projectOntoScreen();
    std::string fileNameSceneScreen = "pointsOntoScreen";
    theScene.saveVector(fileNameSceneScreen);

    ///calibration
    C_klibration theKLIB(fileNameSceneScreen);
    theKLIB.assessParameters();
    theKLIB.saveParameters("parameters");
    return 0;
}







/*C_matrix<double> M(3,4);
C_vector<double> tmp(4), tmpCalc(3);
tmp = 1;
tmp.show();
M=1;
M.show();
tmpCalc = M * tmp;
tmpCalc.show();

return 6;*/
