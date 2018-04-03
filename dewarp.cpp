#include <iostream>
#include <string>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iomanip>
#include <vector>

using namespace std;

template <class T=double>
class Pose
{
public:
    typedef Eigen::Matrix<T,2,1> Vector2T;
    T x;
    T y;
    T theta;

    Pose() : x(0), y(0), theta(0) {

    }
    Pose(Vector2T p, T theta) : x(p(0)), y(p(1)), theta(theta) {
        x=p(0);
        y=p(1);
        theta=theta;

    }

    Pose move(Vector2T p, T dtheta) {
        auto w = Pose2World(p);
        x = w(0);
        y = w(1);
        theta += dtheta;
    }

    Vector2T Pose2World(Eigen::Vector2d p) const {
        Vector2T rv;
        rv(0) = p(0) * cos(theta) - p(1) * sin(theta) + x;
        rv(1) = p(0) * sin(theta) + p(1) * cos(theta) + y;

        return rv;
    }

    Eigen::Transform<T,2,Eigen::Affine> Pose2WorldTransform() const {
        Eigen::Transform<T,2,Eigen::Affine> t;
        t.setIdentity();
        t.rotate(theta);
        t.translate(Vector2T(x,y));
        return t;
    }

    Eigen::Transform<T,2,Eigen::Affine> World2PoseTransform() const {
        return Pose2WorldTransform().inverse();
    }
};

std::string to_string(Pose<> & pose) {
    stringstream ss;
    ss << std::fixed << std::setprecision(2) << pose.x << ", " << pose.y << ", " << pose.theta;
    return ss.str();
}

template <class T = double>
bool is_between(T a, T b, T c) {
    return b <= a && a <=c || c <= a && a <= b;
}

template <class T = double>
class Line {
    typedef Eigen::Matrix<T,2,1> Vector2T;
    typedef Eigen::Matrix<T,3,1> Vector3T;
 public:
    T a;
    T b;
    T c;

   static Line<T> from_points(Vector2T p1, Vector2T p2) {
       auto h1 = Vector3T(p1(0),p1(1),1);      
       auto h2 = Vector3T(p2(0),p2(1),1);
       auto l = h1.cross(h2);
       Line line;
       line.a = l(0);
       line.b = l(1);
       line.c = l(2);
       return line;  
    }

    Vector2T intersection(Line l2) {
        auto h = Vector3T(a,b,c).cross(Vector3T(l2.a, l2.b, l2.c));
        if(h(2)==0) {
            return Vector2T(NAN, NAN);
        }
        return Vector2T(h(0)/h(2), h(1)/h(2));
    }
};

template <class T=double>
class LineSegment {
    typedef Eigen::Matrix<T,2,1> Vector2T;
    public:
    Vector2T p1;
    Vector2T p2;
    LineSegment(Vector2T p1, Vector2T p2) : p1(p1), p2(p2) {      
    }

    Vector2T intersection(LineSegment s2) {
        auto l2 = Line<T>::from_points(s2.p1, s2.p2);
        return intersection(l2);
    }

    Vector2T intersection(Line<T> l2) {
        auto l1 = Line<T>::from_points(p1,p2);
        Vector2T p = l1.intersection(l2);
        if(is_between(p(0),p1(0),p2(0))) {
            return p;
        } else {
            return Vector2T(NAN, NAN);
        }

    }
};

template<class T=double>
T sign_of(T v) {
    if(v >= 0.){
         return 1.;
    }
    return -1.;
}

template<class T=double>
T fake_scan(T theta, vector<LineSegment<T>> & world) {
    typedef Eigen::Matrix<T,2,1> Vector2T;
    T dx = cos(theta);
    T dy = sin(theta);
    auto l = Line<T>::from_points({0,0},{dx,dy});
    T best_d = NAN;
    for( auto segment : world) {
        Vector2T p = segment.intersection(l);
        if(sign_of(p(0))==sign_of(dx) && sign_of(p(1)==sign_of(dy))) {
            double d = p.norm();
            if(isnan(best_d)) {
                best_d = d;
            } else {
                best_d = std::min<T>(d, best_d);
            }
        }
    }
    
    return best_d;
}

void test_fake_scan() {
    vector<LineSegment<double>> world;
    world.push_back(LineSegment<double>({0,1}, {1,1}));
    world.push_back(LineSegment<double>({-10,2}, {10,2}));
    for( int i = 0; i < 360; i++) {
        double theta = i * EIGEN_PI / 180.;
        cout << "degrees: " << i <<  " d:" << fake_scan(theta, world) << endl;
    }


}

void test_intersection() {
    auto l1 = LineSegment<>({1,1},{2,2});
    auto l2 = LineSegment<>({0,1},{1,0});
    auto intersection = l1.intersection(l2);
    cout << "intersection at:" << intersection << endl;
}

int main(int, char**)
{
    cout << "start" << endl;
    
    Pose<> pose({0, 0}, EIGEN_PI);
    Eigen::Vector2d p(0,1), w, wp;

    cout << "pose :" << to_string(pose) << endl;
    cout << "p :" << p << endl;
    cout << "world_p :" << pose.Pose2World(p) << endl;
    auto t = pose.Pose2WorldTransform();
    w = t*p;
    cout << "transform_p :" << w << endl;
    wp = pose.World2PoseTransform() * w;
    cout << "untransform_p :" << wp << endl;
    cout << "done" << endl;

    Pose<> p2;
    for(int i = 0; i < 11; i++) {
        cout << "i:" << i << " pose: " << to_string(p2) << endl;
        p2.move({1, 0}, 2. * EIGEN_PI / 10.);
    }

    auto line = Line<>::from_points({2,1},{3,2});
    cout << "line: " << line.a << ", " << line.b << ", " << line.c << endl;

    test_intersection();
    test_fake_scan();
    return 0;
}