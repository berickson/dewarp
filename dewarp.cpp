#include <iostream>
#include <string>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iomanip>
#include <vector>
#include <functional>
#include <chrono>

#define degrees2radians(theta) ((theta) * EIGEN_PI / 180.)
#define radians2degrees(theta) ((theta) * 180. / EIGEN_PI)

using namespace std;
using namespace std::chrono;


class Stopwatch {
public:
    time_point<system_clock> start_time;
    duration<double> elapsed_time = duration<double>::zero();
    bool started = false;
    void start() {
        start_time = system_clock::now();
        started = true;
    }

    void stop() {
        if(started) {
            elapsed_time += system_clock::now() - start_time;
            started = false;
        }
    }

    double get_elapsed_seconds() {
        if(started) {
            return (elapsed_time + system_clock::now() - start_time).count();
        }
        return elapsed_time.count();
    }
};


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


template <class T = double>
struct ScanLine {
    T theta = NAN;
    T d = NAN;
};


std::string to_string(Pose<> & pose) {
    stringstream ss;
    ss << std::fixed << std::setprecision(2) << pose.x << ", " << pose.y << ", " << pose.theta;
    return ss.str();
}

template <class T>
void print_scan(vector<T> scan) {
    double theta = 0;
    double dtheta = 2 * EIGEN_PI / scan.size();
    cout << "degrees, d, px, py" << endl;
    for(auto d : scan) {
        double px = d*cos(theta);
        double py = d*sin(theta);
        cout << radians2degrees(theta) << ", " << d << ", " << px << ", " << py << endl;
        theta += dtheta;
    }
}

template <class T = double>
bool is_between(T a, T b, T c, T e=1E-5) {
    return (b-e <= a && a <=c+e) || (c-e <= a && a <= b+e);
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
    std::string to_string() {
        stringstream ss;
        ss << "(" << p1(0) << "," << p1(1)<< ")-(" << p2(0) << "," << p2(1) << ")";
        return ss.str(); 
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
T fake_laser_reading(const Pose<T> & pose, T theta, const vector<LineSegment<T>> & world) {
    typedef Eigen::Matrix<T,2,1> Vector2T;
    T dx = cos(theta + pose.theta);
    T dy = sin(theta + pose.theta);
    auto l = Line<T>::from_points({pose.x,pose.y},{pose.x+dx,pose.y+dy});
    T best_d = NAN;
    for( auto segment : world) {
        Vector2T p = segment.intersection(l);
        if(isnan(p(0))) {
            continue;
        }
        p(0) -= pose.x;
        p(1) -= pose.y;
        if((sign_of(p(0))==sign_of(dx)) && (sign_of(p(1))==sign_of(dy))) {
            T d = p.norm();
            if(isnan(best_d)) {
                best_d = d;
            } else {
                best_d = std::min<T>(d, best_d);
            }
        }
    }

    return best_d;
}

template <class T=double>
vector<T> scan_with_twist(vector<LineSegment<T>> & world, int scan_count, T twist_x, T twist_y, T twist_theta, Pose<T> initial_pose = Pose<T>()) {
    vector<T> output;
    Pose<T> pose = initial_pose;
    for(int i=0; i < scan_count; i++) {
        T scanner_theta = i * EIGEN_PI / 180;
        T d = fake_laser_reading<>(pose, scanner_theta, world);
        output.push_back(d);
        pose.move({twist_x/scan_count, twist_y/scan_count}, twist_theta/scan_count);
    }
    return output;
}


template <class T=double>
vector<ScanLine<T>> scan_with_twist2(vector<LineSegment<T>> & world, int scan_count, Pose<T> initial_pose = Pose<T>(), T twist_x = 0, T twist_y = 0, T twist_theta = 0) {
    vector<ScanLine<T>> output;
    Pose<T> pose = initial_pose;
    for(int i=0; i < scan_count; i++) {
        T scanner_theta = i * EIGEN_PI / 180;
        T d = fake_laser_reading<>(pose, scanner_theta, world);
        ScanLine<T> scan_line;
        scan_line.theta = scanner_theta;
        scan_line.d = d;
        output.push_back(scan_line);
        pose.move({twist_x/scan_count, twist_y/scan_count}, twist_theta/scan_count);
    }
    return output;
}


void print_world(vector<LineSegment<double>> & world) {
    for(auto segment : world) {
        std::cout << segment.to_string() << endl;
    }
}


Stopwatch untwist_timer;

template <class T=double>
vector<T> untwist_scan(vector<T> &twisted_readings, T twist_x, T twist_y, T twist_theta, Pose<T> initial_pose = Pose<T>()) {
    untwist_timer.start();
    typedef Eigen::Matrix<T,2,1> Vector2T;
    int count = twisted_readings.size();
    //cout << "count: " << count << endl;
    Pose<T> pose = initial_pose;
    vector<LineSegment<T>> world;
    T inc_theta = 2.*EIGEN_PI / count;
    Vector2T p1;
    p1 << NAN,NAN;
    Vector2T p2;
    p2 << NAN,NAN;
    for(int i = 0; i < twisted_readings.size()+1; i++) {
    //for(int i = 267; i < 273; i++) {
        T scan_theta = (T) i / count * 2. * EIGEN_PI;
        T d1 = twisted_readings[i%count];
        p1 = p2;
        p2 = pose.Pose2World({cos(scan_theta)*d1, sin(scan_theta)*d1});
        if(! isnan(p1(0)) && !isnan(p2(0))) {
          world.push_back(LineSegment<T>(p1, p2));
        }
        pose.move({twist_x/count, twist_y/count}, twist_theta/count);
    }

    //cout << endl << endl <<  "untwist world" << endl;
    //print_world(world);

    vector<T> output;

    pose.x=0;
    pose.y=0;
    pose.theta=0;
    for(int i = 0; i < count; i++) {
        T scan_theta = (T) i / count * 2 * EIGEN_PI;
        output.push_back(fake_laser_reading<T>(pose, scan_theta, world));
    }
    untwist_timer.stop();
    return output;
}

void test_fake_scan() {
    vector<double> scan;
    vector<LineSegment<double>> world;
    world.push_back(LineSegment<double>({-10,2}, {10,2}));
    world.push_back(LineSegment<double>({-10,-2}, {10,-2}));
    Pose<double> pose;
    double d = fake_laser_reading<double>(pose, 260*EIGEN_PI/180., world);
    for( int i = 0; i < 360; i++) {
        double theta = i * EIGEN_PI / 180.;
        double d = fake_laser_reading<double>(pose, theta, world);
        scan.push_back(d);
    }
    print_scan(scan);
}

vector<LineSegment<double>> get_world() {
    vector<LineSegment<double>> world;
    world.push_back(LineSegment<double>({0,1}, {1,1}));
    world.push_back(LineSegment<double>({-10,2}, {10,2}));
    world.push_back(LineSegment<double>({-10,-3}, {10,-3}));
    world.push_back(LineSegment<double>({10,2}, {10,-3}));
    return world;
}


void test_scan_with_twist() {
    auto world = get_world();
    double twist_theta = 5.*EIGEN_PI/180.;
    double twist_x = .2;
    double twist_y = .1;
    double reading_count = 360;
    auto twisted = scan_with_twist(world, reading_count, twist_x, twist_y, twist_theta);
    auto untwisted = untwist_scan(twisted, twist_x, twist_y, twist_theta);
    cout << endl << endl << "twisted scan" << endl;
    print_scan(twisted);
    cout << endl << endl << "untwisted scan" << endl;
    print_scan(untwisted);
}

void test_intersection() {
    double dx = -2.89707;
    double dy = -3;
    Eigen::Vector2d v{dx,dy};
    v = v/v.norm();
    auto s = LineSegment<>({-2.89755,-3},{-2.89707,-3});
    auto l = Line<double>::from_points({0,0},v);
    auto p = s.intersection(l);
    cout << "intersection at:" << p<< endl;
    if((sign_of(p(0))==sign_of(dx)) && (sign_of(p(1))==sign_of(dy))) {
    
    cout << "signs ok";
    }
}

template<class T = double>
T scan_difference(vector<T> scan1, vector<T> scan2) {
    double total_difference = 0;
    size_t valid_count = 0;
    for(unsigned i = 0; i < scan1.size(); i++) {
        double d = scan1[i]-scan2[i];
        if(!isnan(d)) {
            total_difference += fabs(d);
            valid_count++;
        }
    }
    if(valid_count == 0) {
        return NAN;
    }
    return total_difference / valid_count;
}

struct TwiddleResult{
    vector<double> p;
    double error;
};

double abs_sum(vector<double> & v) {
    double rv = 0;
    for(auto p : v) {
        rv += fabs(p);
    }
}

// based loosely on https://martin-thoma.com/twiddle/
TwiddleResult twiddle(vector<double> guess, std::function<double(const vector<double>)> f, double threshold = 0.001) {
    
    // initialize parameters to guess
    vector<double> p = guess;

    // potential changes
    auto dp = std::vector<double>(guess.size(), 0.2);

    // change growth rate
    const double growth_rate = 1.5;

    double best_error = f(p);

    while(abs_sum(dp) > threshold) {
        for(int i = 0; i< p.size(); i++) {
            p[i] += dp[i];
            double error = f(p);

            if (error < best_error) {
                best_error = error;
                dp[i] *= 1.1;
            } else {
                // There was no improvement
                p[i] -= 2*dp[i];  // Go into the other direction
                error = f(p);

                if (error < best_error) {
                  // There was an improvement
                    best_error = error;
                    dp[i] *= growth_rate;
                } else  {
                    // There was no improvement
                    p[i] += dp[i];
                    // As there was no improvement, the step size in either
                    // direction, the step size might simply be too big.
                    dp[i] /= growth_rate;
                }
            }
        }
    }
    TwiddleResult rv;
    rv.p = p;
    rv.error = best_error;
    return rv;
}

template <class T = double>
Pose<T> match_scans(vector<double> scan1, vector<double> scan2) {
    auto error_function = [&](vector<double> params){
        Pose<T> pose;
        pose.x = params[0];
        pose.y = params[1];
        pose.theta = params[2];
        auto scan2b = untwist_scan(scan2,0.,0.,0.,pose);
        double d = scan_difference(scan1, scan2b);
        cout << "difference: " << d << " pose: " << to_string(pose) << endl;

        return d;
    };

    TwiddleResult r = twiddle({0,0,0}, error_function);
    Pose<T> match;
    match.x = r.p[0];
    match.y = r.p[1];
    match.theta = r.p[2];
    return match;
}




template <class T = double>
vector<ScanLine<T>> move_scan(vector<ScanLine<T>> scan, Pose<T> pose) {
    typedef Eigen::Matrix<T,2,1> Vector2T;
    Vector2T p, p_new;
    T radians_resolution = 2. * EIGEN_PI / scan.size();
    auto transform = pose.Pose2WorldTransform();

    vector<ScanLine<T>> moved_scan;
    
    for(auto & scan_line : scan) {
        double theta = scan_line.theta;
        double d = scan_line.d;
        ScanLine<T> moved_scan_line;
        if(!isnan(d)) {
            p << d * cos(theta), d * sin(theta);
            p_new = transform * p;
            moved_scan_line.theta = atan2(p_new(1), p_new(0));
            moved_scan_line.d = sqrt(p_new(1)*p_new(1)+p_new(0)*p_new(0));
        }
        moved_scan.push_back(moved_scan_line);
    }
    return moved_scan;
}

template <class T = double>
Pose<T> match_scans2(vector<double> scan1, vector<double> scan2) {
    auto error_function = [&](vector<double> params){
        Pose<T> pose;
        pose.x = params[0];
        pose.y = params[1];
        pose.theta = params[2];

        auto scan2b = move_scan(scan2, pose);
        double d = scan_difference2(scan1, scan2b);
        cout << "difference: " << d << " pose: " << to_string(pose) << endl;

        return d;
    };

    TwiddleResult r = twiddle({0,0,0}, error_function);
    Pose<T> match;
    match.x = r.p[0];
    match.y = r.p[1];
    match.theta = r.p[2];
    return match;
}



void test_scan_difference() {
    size_t n_points = 360;
    auto world = get_world();
    Pose<double> pose1;
    pose1.x = 0;
    Pose<double> pose2;
    pose2.x = .25;
    pose2.y = -.75;
    pose2.theta = .05;
    auto scan1 = scan_with_twist<double>(world, n_points, 0, 0, 0, pose1);
    auto scan2 = scan_with_twist<double>(world, n_points, 0, 0, 0, pose2);

    cout << scan_difference(scan1, scan2);
    auto matched_pose = match_scans(scan1, scan2);

    cout << "match_scans: " << to_string(matched_pose) << endl;
}

void test_twiddle() {
    auto lambda = [](std::vector<double> v)->double{return fabs(v[0]-3)+fabs(v[1]-5) + fabs(v[2]);};
    auto rv = twiddle({0,0,0}, lambda);
    cout << rv.p[0] << ", " << rv.p[1] << ", " << rv.p[2] << " error: " << rv.error << endl;
        
}

void test_move_scan() {
    auto world = get_world();
    Pose<double> pose1;
    pose1.x = 0;
    pose1.y = 0;
    pose1.theta = 0;
    Pose<double> pose2;
    pose2.x = .0;
    pose2.y = 0;
    pose2.theta = degrees2radians(1);
    auto scan1 = scan_with_twist2<double>(world, 360, pose1);
    auto scan2 = move_scan(scan1, pose2);
    cout << "original scan, moved scan" << endl;
    for(unsigned i = 0; i < scan2.size(); i++) {
        cout << radians2degrees(scan1[i].theta) << " " << scan1[i].d 
             << " -> "  
             << radians2degrees(scan2[i].theta) << " " << scan2[i].d << endl;
    }
}

int main(int, char**)
{
    test_move_scan();
    //test_twiddle();
    //return 0;
    //test_scan_difference();
    cout << "time untwisting: " << untwist_timer.get_elapsed_seconds() << endl;
    return 0;
    cout << "newly compiled 4" << endl;
    test_scan_with_twist();
    return 0;
    test_intersection();

    test_fake_scan();
    return 0;
    test_fake_scan();
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

    return 0;
}