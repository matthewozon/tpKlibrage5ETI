#include <C_klibration.h>

C_klibration::C_klibration(std::string fileNameSceneScreen)
{
    m_fileNameScene = fileNameSceneScreen;
    loadObjectAndScreenPoints();
}

C_klibration::~C_klibration()
{
}

bool C_klibration::loadObjectAndScreenPoints()
{
    std::fstream fhdr;
    fhdr.open (m_fileNameScene.data(), std::ifstream::in); // | ifstream::binary
    if (fhdr.is_open())
    {
        std::string line;
        double x=0.0, y=0.0, z=0.0, u=0.0, v=0.0;
        C_vector<double> tmp(5);
        while(getline (fhdr,line)!=0)
        {
            coordinatesFromLine(line, x, y, z, u, v);
            tmp.set(0,x);
            tmp.set(1,y);
            tmp.set(2,z);
            tmp.set(3,u);
            tmp.set(4,v);
            m_Scene.push_back(tmp);
        }
        fhdr.close();
    }
    else
    {
        return false;
    }
    return true;
}

void C_klibration::coordinatesFromLine(std::string line, double& x, double& y, double& z, double& u, double& v)
{
    std::size_t found = line.find_first_of(" ");
    std::size_t found1 = line.find_first_of(" ", found+1); //skip first term that happen to be the index of the point
    x = atof(line.substr(found+1,found1).data());
    std::size_t found2 = line.find_first_of(" ", found1+1);
    y = atof(line.substr(found1+1,found2).data());
    std::size_t found3 = line.find_first_of(" ", found2+1);
    z = atof(line.substr(found2+1,found3).data());
    std::size_t found4 = line.find_first_of(" ", found3+1);
    u = atof(line.substr(found3+1,found4).data());
    v = atof(line.substr(found4+1,std::string::npos).data());
    return;
}

void C_klibration::assessParameters(void)
{
    unsigned short I = 2;
    //first try to estimate L
    C_vector<double> * L = calculateLVector();

    //estimate parameter that can be calculated
    Ty = 1.0/sqrt((*L)[4]*(*L)[4] + (*L)[5]*(*L)[5] + (*L)[6]*(*L)[6]); //must determine the sign of Ty

    //assume Ty positive at first
    beta = Ty*sqrt((*L)[0]*(*L)[0] + (*L)[1]*(*L)[1] + (*L)[2]*(*L)[2]);
    Tx = Ty * (*L)[3]/beta;

    //rotation matrix
    r11 = (*L)[0]*Ty/beta; r12 = (*L)[1]*Ty/beta; r13 = (*L)[2]*Ty/beta;
    r21 = (*L)[4]*Ty; r22 = (*L)[5]*Ty; r23 = (*L)[6]*Ty;
    r31 = r12*r23 - r13*r22;
    r32 = r13*r21 - r11*r23;
    r33 = r11*r22 - r12*r21;


    //focal and focal ratio
    C_vector<double> * TzF2 = calculateTzAndf2();
    Tz = (*TzF2)[0];
    f2 = (*TzF2)[1];
    delete TzF2;

    ///OK til here

    ////////////////////////////////////////////////////////////////////////
    //calculate the projection of the first point onto the screen and store it in Vy1
    double Vy1 = CY + f2*(r21*(m_Scene.at(I))[0] + r22*(m_Scene.at(I))[1] + r23*(m_Scene.at(I))[2] + Ty)/(r31*(m_Scene.at(I))[0] + r32*(m_Scene.at(I))[1] + r33*(m_Scene.at(I))[2] + Tz);
    double Vx1 = CX + f2*(r11*(m_Scene.at(I))[0] + r12*(m_Scene.at(I))[1] + r13*(m_Scene.at(I))[2] + Tx)/(r31*(m_Scene.at(I))[0] + r32*(m_Scene.at(I))[1] + r33*(m_Scene.at(I))[2] + Tz);
    double D1 = ABS(Vy1-(m_Scene.at(I))[4]) + ABS(Vx1-(m_Scene.at(I))[3]);



    ///////////////////////////////////////////////////////////////////////
    //assume Ty negative
    beta = Ty*sqrt((*L)[0]*(*L)[0] + (*L)[1]*(*L)[1] + (*L)[2]*(*L)[2]);
    Ty=-Ty;
    Tx = Ty * (*L)[3]/beta;

    //rotation matrix
    r11 = (*L)[0]*Ty/beta; r12 = (*L)[1]*Ty/beta; r13 = (*L)[2]*Ty/beta;
    r21 = (*L)[4]*Ty; r22 = (*L)[5]*Ty; r23 = (*L)[6]*Ty;
    r31 = r12*r23 - r13*r22;
    r32 = r13*r21 - r11*r23;
    r33 = r11*r22 - r12*r21;

    //focal and focal ratio
    TzF2 = calculateTzAndf2();
    Tz = (*TzF2)[0];
    f2 = (*TzF2)[1];
    delete TzF2;

    double Vy2 = CY + f2*(r21*(m_Scene.at(I))[0] + r22*(m_Scene.at(I))[1] + r23*(m_Scene.at(I))[2] + Ty)/(r31*(m_Scene.at(I))[0] + r32*(m_Scene.at(I))[1] + r33*(m_Scene.at(I))[2] + Tz);
    double Vx2 = CX + f2*(r11*(m_Scene.at(I))[0] + r12*(m_Scene.at(I))[1] + r13*(m_Scene.at(I))[2] + Tx)/(r31*(m_Scene.at(I))[0] + r32*(m_Scene.at(I))[1] + r33*(m_Scene.at(I))[2] + Tz);
    double D2 = ABS(Vy2-(m_Scene.at(I))[4]) + ABS(Vx2-(m_Scene.at(I))[3]);

    if(D1<D2)
    {
        //back to previous configuration
        Ty = ABS(Ty);
        Tx = Ty * (*L)[3]/beta;

        //rotation matrix
        r11 = (*L)[0]*Ty/beta; r12 = (*L)[1]*Ty/beta; r13 = (*L)[2]*Ty/beta;
        r21 = (*L)[4]*Ty; r22 = (*L)[5]*Ty; r23 = (*L)[6]*Ty;
        r31 = r12*r23 - r13*r22;
        r32 = r13*r21 - r11*r23;
        r33 = r11*r22 - r12*r21;

        //focal and focal ratio
        TzF2 = calculateTzAndf2();
        Tz = (*TzF2)[0];
        f2 = (*TzF2)[1];
        delete TzF2;
    }

    delete L;
    return;
}

bool C_klibration::saveParameters(std::string fileNameParmeters)
{
    std::ofstream myfile;
    myfile.open (fileNameParmeters.data(), std::ifstream::out);
    if(!myfile.is_open()) return false;

    myfile << "Ext param : translation & rotation \n";
    myfile << "\tTranslation: \n";
    myfile << "(Tx, Ty, Tz) = (" << Tx << ", " << Ty << ", " << Tz  << ") mm\n";
    myfile << "\tRotation: \n";
    myfile << "r11 = " << r11 << "\t r12 = " << r12 << "\t r13 = " << r13 <<"\n";
    myfile << "r21 = " << r21 << "\t r22 = " << r22 << "\t r23 = " << r23 <<"\n";
    myfile << "r31 = " << r31 << "\t r32 = " << r32 << "\t r33 = " << r33 <<"\n";
    myfile << "Rotation angles: \n";
    myfile << "phi = " << -atan(r23/r33) << "\n";
    myfile << "gamma = " << -atan(r12/r11) << "\n";
    myfile << "omega = " << -atan(r13/(r33*cos(phi)-r23*sin(phi))) << "\n";
    myfile << "\nInt param : f1/f2 & f2\n";
    myfile << "beta " << beta << " mm\n";
    myfile << "f2 " << f2 << " mm\n";


    myfile.close();
    return true;
}


C_vector<double> * C_klibration::calculateLVector(void)
{
    C_matrix<double>* AtA = new C_matrix<double>(7,7);
    C_matrix<double>* A = createMatrixA();
    C_vector<double>* U = createVectorU();
    *AtA =  (A->Transpose()) * *A;
    C_vector<double>* L = AtA->LineAlgEq_LU((A->Transpose()) * *U);
    delete AtA;
    delete U;
    delete A;
    return L;
}

C_matrix<double>* C_klibration::createMatrixA(void)
{
    C_matrix<double>* A = new C_matrix<double>(m_Scene.size(),7);
    for(unsigned int i=0 ; i<m_Scene.size() ; i++)
    {
        double x = (m_Scene.at(i))[0];
        double y = (m_Scene.at(i))[1];
        double z = (m_Scene.at(i))[2];
        double u = (m_Scene.at(i))[3]-CX;
        double v = (m_Scene.at(i))[4]-CY;
        A->set(i,0, v * x );
        A->set(i,1, v * y );
        A->set(i,2, v * z );

        A->set(i,3, v );

        A->set(i,4, -u * x );
        A->set(i,5, -u * y );
        A->set(i,6, -u * z );
    }
    return A;
}

C_vector<double>* C_klibration::createVectorU(void)
{
    C_vector<double>* U = new C_vector<double>(m_Scene.size());
    for(unsigned int i=0 ; i<m_Scene.size() ; i++)
    {
        U->set(i,m_Scene.at(i)[3]-CX);
    }
    return U;
}

C_vector<double> * C_klibration::calculateTzAndf2(void)
{
    C_matrix<double>* BtB = new C_matrix<double>(2,2);
    C_matrix<double>* B = createMatrixB();
    C_vector<double>* R = createVectorR();
    *BtB =  (B->Transpose()) * *B;
    C_vector<double>* TzF2 = BtB->LineAlgEq_LU((B->Transpose()) * *R);
    delete BtB;
    delete R;
    delete B;
    return TzF2;
}

C_matrix<double>* C_klibration::createMatrixB(void)
{
    C_matrix<double>* B = new C_matrix<double>(m_Scene.size(),2);
    for(unsigned int i=0 ; i<m_Scene.size() ; i++)
    {
        B->set(i,0, m_Scene.at(i)[4]-CY );
        B->set(i,1, -(r21*m_Scene.at(i)[0] + r22*m_Scene.at(i)[1] + r23*m_Scene.at(i)[2] + Ty) );
    }
    return B;
}

C_vector<double>* C_klibration::createVectorR(void)
{
    C_vector<double>* R = new C_vector<double>(m_Scene.size());
    for(unsigned int i=0 ; i<m_Scene.size() ; i++)
    {
        R->set(i,-(m_Scene.at(i)[4]-CY)*(r31*m_Scene.at(i)[0] + r32*m_Scene.at(i)[1] + r33*m_Scene.at(i)[2]));
    }
    return R;
}
