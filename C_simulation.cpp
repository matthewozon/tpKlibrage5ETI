#include <C_simulation.h>

C_simulation::C_simulation(std::string fileNameSceneObject)
{
    m_fileNameScene = fileNameSceneObject;
    loadScenePointsFromFile();

    //init matrices
    m_Mext = new C_matrix<double>(4,4);
    m_Mint = new C_matrix<double>(3,4);
}

C_simulation::~C_simulation()
{
    //init matrices
    delete m_Mext;
    delete m_Mint;
}


bool C_simulation::loadScenePointsFromFile(void)
{
    std::fstream fhdr;
    fhdr.open (m_fileNameScene.data(), std::ifstream::in); // | ifstream::binary
    if (fhdr.is_open())
    {
        std::string line;
        double x=0.0, y=0.0, z=0.0;
        C_vector<double> tmp(4);
        while(getline (fhdr,line)!=0)
        {
            coordinatesFromLine(line, x, y, z);
            tmp.set(0,x);
            tmp.set(1,y);
            tmp.set(2,z);
            tmp.set(3,1.0);
            m_objetScene.push_back(tmp);
        }
        fhdr.close();
    }
    else
    {
        return false;
    }
    return true;
}

void C_simulation::coordinatesFromLine(std::string line, double& x, double& y, double& z)
{
    ///we assume the the line contains exactly three numbers separated by space
    std::size_t found = line.find_first_of(" ");
    x = atof(line.substr(0,found).data());
    std::size_t found2 = line.find_first_of(" ", found+1);
    y = atof(line.substr(found+1,found2).data());
    z = atof(line.substr(found2+1,std::string::npos).data());
    return;
}


void C_simulation::setRotationAndTranslation(double omega, double phi, double gamma, double dx, double dy, double dz)
{
    //

    m_Mext->set(0,0,cos(omega)*cos(gamma)); m_Mext->set(0,1,-cos(omega)*sin(gamma)); m_Mext->set(0,2,sin(omega)); m_Mext->set(0,3,dx);
    m_Mext->set(1,0,sin(phi)*sin(omega)*cos(gamma)+cos(phi)*sin(gamma)); m_Mext->set(1,1,-sin(phi)*sin(omega)*sin(gamma)+cos(phi)*cos(gamma)); m_Mext->set(1,2,-sin(phi)*cos(omega)); m_Mext->set(1,3,dy);
    m_Mext->set(2,0,-cos(phi)*sin(omega)*cos(gamma) + sin(phi)*sin(gamma)); m_Mext->set(2,1,cos(phi)*sin(omega)*sin(gamma) + sin(phi)*cos(gamma)); m_Mext->set(2,2,cos(phi)*cos(omega)); m_Mext->set(2,3,dz);
    m_Mext->set(3,0,0.0); m_Mext->set(3,1,0.0); m_Mext->set(3,2,0.0); m_Mext->set(3,3,1.0);
    return;
}


C_matrix<double>* C_simulation::setFocal(double f)
{
    C_matrix<double>* m_Focal = new C_matrix<double>(4,4);
    *m_Focal = 0.0;
    m_Focal->set(0,0,f);
    m_Focal->set(1,1,f);
    m_Focal->set(2,2,f);
    m_Focal->set(3,2,1.0);
    return m_Focal;
}

C_matrix<double>* C_simulation::setPixelSizeAndCenter(double sx, double sy, double cx, double cy)
{
    C_matrix<double>* m_Scaling = new C_matrix<double>(3,4);
    *m_Scaling = 0.0;
    m_Scaling->set(0,0,1.0/sx);
    m_Scaling->set(1,1,1.0/sy);
    m_Scaling->set(0,3,cx);
    m_Scaling->set(1,3,cy);
    m_Scaling->set(2,3,1.0);
    return m_Scaling;
}

void C_simulation::setIntrinsecParameters(double f, double sx, double sy, double cx, double cy)
{
    C_matrix<double>* m_Scaling = setPixelSizeAndCenter(sx, sy, cx, cy);
    C_matrix<double>* m_Focal = setFocal(f);
    //m_Scaling->show();
    //m_Focal->show();
    *m_Mint = *m_Scaling * *m_Focal;
    delete m_Focal;
    delete m_Scaling;
    return;
}



void C_simulation::projectOntoScreen(void)
{
    C_vector<double> tmp(2), tmpCalc(3);
    C_matrix<double> M(3,4);
    m_Mext->show();
    //std::getchar();
    //m_Mint->show();
    M = *m_Mint * *m_Mext;
    //M.show();

    for(unsigned int i=0 ; i<m_objetScene.size() ; i++)
    {
        M * m_objetScene.at(i);
        tmpCalc = M * m_objetScene.at(i);
        tmp.set(0,tmpCalc[0]/tmpCalc[2]);
        tmp.set(1,tmpCalc[1]/tmpCalc[2]);
        m_objetScreen.push_back(tmp);
    }
    return;
}


bool C_simulation::saveVector(std::string fileNameSceneScreen)
{
    std::ofstream myfile;
    myfile.open (fileNameSceneScreen.data());
    for(unsigned int i=0 ; i<m_objetScreen.size() ; i++)
    {
        myfile << i << " " << m_objetScene.at(i)[0] << " " << m_objetScene.at(i)[1] << " " << m_objetScene.at(i)[2]\
                   << " " << m_objetScreen.at(i)[0] << " " << m_objetScreen.at(i)[1] << "\n";
    }
    myfile.close();
    return true;
}
