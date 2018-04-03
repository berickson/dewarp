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
        Eigen::Transform<double,2,Eigen::Affine> t;
        t.setIdentity();
        t.rotate(theta);
        t.translate(Eigen::Vector2d(x,y));
        return t;
    }
    Eigen::Transform<double,2,Eigen::Affine> World2PoseTransform() {
        return Pose2WorldTransform().inverse();
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
    
    Pose pose(0,0,EIGEN_PI);
    Eigen::Vector2d p,w,wp;
    p << 0,1;
    cout << "pose :" << to_string(pose) << endl;
    cout << "p :" << p << endl;
    cout << "world_p :" << pose.Pose2World(p) << endl;
    auto t = pose.Pose2WorldTransform();
    w = t*p;
    cout << "transform_p :" << w << endl;
    wp = pose.World2PoseTransform() * w;
    cout << "untransform_p :" << wp << endl;
    //cout << "tranform :" << pose.Pose2WorldTransform() << endl;
    cout << "done" << endl;
    return 0;
}