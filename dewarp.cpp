#include <iostream>
#include <string>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
using namespace std;

class Pose
{
public:
    double x;
    double y;
    double theta;
    Pose(double x, double y, double theta) : x(x), y(y), theta(theta) {
    }

    Eigen::Vector2d Pose2World(Eigen::Vector2d p) {
        Eigen::Vector2d rv;
        rv(0) = p(0) * cos(theta) - p(1) * sin(theta) + x;
        rv(1) = p(0) * sin(theta) + p(1) * cos(theta) + y;

        return rv;
    }

    Eigen::Transform<double,2,Eigen::Affine> Pose2WorldTransform() {
        Eigen::Rotation2D<double> rotate(theta);
        Eigen::Translation<double,2> translate(x,y);
        Eigen::Transform<double,2,Eigen::Affine> t;
        t.setIdentity();
        t.rotate(theta);
        t.translate(Eigen::Vector2d(x,y));
        //(rotate);
        return t;
    }
};

std::string to_string(Pose & pose) {
    stringstream ss;
    ss << pose.x << ", " << pose.y << ", " << pose.theta;
    return ss.str();
}

int main(int, char**)
{
    cout << "start" << endl;
    
    Pose pose(0,0,M_PI);
    Eigen::Vector2d p;
    p << 0,1;
    cout << "pose :" << to_string(pose) << endl;
    cout << "p :" << p << endl;
    cout << "world_p :" << pose.Pose2World(p) << endl;
    auto t = pose.Pose2WorldTransform();
    cout << "transform_p :" << t*p << endl;
    //cout << "tranform :" << pose.Pose2WorldTransform() << endl;
    cout << "done" << endl;
    return 0;
}