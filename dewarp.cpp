#include <iostream>
#include <string>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iomanip>
#include <vector>
#include <functional>
#include <chrono>
#include <random>

#define degrees2radians(theta) ((theta) * EIGEN_PI / 180.)
#define radians2degrees(theta) ((theta) * 180. / EIGEN_PI)

using namespace std;
using namespace std::chrono;

template <class T>
struct Point2d {
    T x;
    T y;
};


class Stopwatch {
public:
    size_t start_count = 0;
    time_point<system_clock> start_time;
    duration<double> elapsed_time = duration<double>::zero();
    bool started = false;
    void start() {
        ++ start_count;
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
private:
    typedef Eigen::Matrix<T,2,1> Vector2T;
    T x;
    T y;
    T theta;
    T cos_theta;
    T sin_theta;

public:
    inline T get_x() const { return x; }
    inline T get_y()  const { return y; }
    inline T get_theta()  const { return theta; }

    void update_trig() {
        cos_theta = cos(theta);
        sin_theta = sin(theta);
    }

    Pose(T x = 0.0, T y=0.0, T theta = 0.0) : x(x), y(y), theta(theta) {
        update_trig();
    }

    Pose(Vector2T p, T theta) : x(p(0)), y(p(1)), theta(theta) {
        update_trig();
    }

    void move(Point2d<T> p, T dtheta) {
        Point2d<T> w;
        Pose2World(p, w);
        x = w.x;
        y = w.y;
        theta += dtheta;

        update_trig();
    }

    inline void  Pose2World(const Point2d<T> & p, Point2d<T> & rv) const {
        rv.x = p.x * cos_theta - p.y * sin_theta + x;
        rv.y = p.x * sin_theta + p.y * cos_theta + y;
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
    ScanLine(T theta = NAN, T d = NAN) : theta(theta), d(d) {
    }

    T theta = NAN;
    T d = NAN;
};



template <class T>
std::string to_string(Pose<T> & pose) {
    stringstream ss;
    ss << std::fixed << std::setprecision(3) << pose.get_x() 
       << ", " << pose.get_y() 
       << ", " << std::setprecision(1) << radians2degrees(pose.get_theta()) << "°";
    return ss.str();
}

template <class T>
void print_scan(vector<ScanLine<T>> scan) {
    cout << "degrees, d, px, py" << endl;
    for(auto s : scan) {
        double theta = s.theta;
        double px = s.d*cos(theta);
        double py = s.d*sin(theta);
        cout << radians2degrees(s.theta) << ", " << s.d << ", " << px << ", " << py << endl;
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

   static Line<T> from_points(const Point2d<T> & p1, const Point2d<T> & p2) {
       auto h1 = Vector3T(p1.x,p1.y,1);      
       auto h2 = Vector3T(p2.x,p2.y,1);
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
    Point2d<T> p1;
    Point2d<T> p2;
    LineSegment(Point2d<T> p1, Point2d<T> p2) : p1(p1), p2(p2) {      
    }

    Point2d<T> intersection(LineSegment & s2) {
        auto l2 = Line<T>::from_points(s2.p1, s2.p2);
        return intersection(l2);
    }

    Vector2T intersection(Line<T> l2) {
        auto l1 = Line<T>::from_points(p1,p2);
        Vector2T p = l1.intersection(l2);
        if(is_between<T>(p(0),p1.x,p2.x)) {
            return p;
        } else {
            return Vector2T(NAN, NAN);
        }
    }
    std::string to_string() {
        stringstream ss;
        ss << "(" << p1.x << "," << p1.y<< ")-(" << p2.x << "," << p2.y << ")";
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
    T dx = cos(theta + pose.get_theta());
    T dy = sin(theta + pose.get_theta());
    auto l = Line<T>::from_points({pose.get_x(),pose.get_y()},{pose.get_x()+dx,pose.get_y()+dy});
    T best_d = NAN;
    for( auto segment : world) {
        Vector2T p = segment.intersection(l);
        if(isnan(p(0))) {
            continue;
        }
        p(0) -= pose.get_x();
        p(1) -= pose.get_y();
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
vector<ScanLine<T>> scan_with_twist(vector<LineSegment<T>> & world, int scan_count, T twist_x, T twist_y, T twist_theta, Pose<T> initial_pose = Pose<T>()) {
    vector<ScanLine<T>> output;
    Pose<T> pose = initial_pose;
    for(int i=0; i < scan_count; i++) {
        T scanner_theta = 2*EIGEN_PI * i / scan_count;;
        T d = fake_laser_reading<>(pose, scanner_theta, world);
        output.push_back({scanner_theta, d});
        pose.move({twist_x/scan_count, twist_y/scan_count}, twist_theta/scan_count);
    }
    return output;
}


template <class T=double>
vector<ScanLine<T>> scan_with_twist2(vector<LineSegment<T>> & world, int scan_count, Pose<T> initial_pose = Pose<T>(), T twist_x = 0, T twist_y = 0, T twist_theta = 0) {
    vector<ScanLine<T>> output;
    Pose<T> pose = initial_pose;
    for(int i=0; i < scan_count; i++) {
        T scanner_theta = (2. * EIGEN_PI) * i / scan_count;
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
Stopwatch move_scan_timer;
Stopwatch scan_difference_timer;
Stopwatch match_scans_timer;

template <class T=double>
vector<ScanLine<T>> untwist_scan(vector<ScanLine<T>> &twisted_readings, T twist_x, T twist_y, T twist_theta, Pose<T> initial_pose = Pose<T>()) {
    untwist_timer.start();
    //typedef Eigen::Matrix<T,2,1> Vector2T;
    int count = twisted_readings.size();
    //cout << "count: " << count << endl;
    Pose<T> pose = initial_pose;
    vector<LineSegment<T>> world;
    world.reserve(count);
    Point2d<T> p1 = {NAN, NAN};
    Point2d<T> p2 = {NAN, NAN};
    for(size_t i = 0; i < twisted_readings.size()+1; i++) {
        T scan_theta = (T) i / count * 2. * EIGEN_PI;
        T d1 = twisted_readings[i%count].d;
        p1=p2;
        pose.Pose2World({cos(scan_theta)*d1, sin(scan_theta)*d1}, p2);
        if(! isnan(p1.x) && !isnan(p2.x)) {
          world.push_back(LineSegment<T>(p1, p2));
        }
        pose.move({twist_x/count, twist_y/count}, twist_theta/count);
    }


    vector<ScanLine<T>> output;
    output.reserve(count);

    Pose<T> pose2(0,0,0);
    for(int i = 0; i < count; i++) {
        T scan_theta = (T) i / count * 2 * EIGEN_PI;
        output.push_back(ScanLine<T>(scan_theta, fake_laser_reading<T>(pose2, scan_theta, world)));
    }
    untwist_timer.stop();
    return output;
}

void test_fake_scan() {
    vector<ScanLine<double>> scan;
    vector<LineSegment<double>> world;
    world.push_back(LineSegment<double>({-10,2}, {10,2}));
    world.push_back(LineSegment<double>({-10,-2}, {10,-2}));
    Pose<double> pose;
    for( int i = 0; i < 360; i++) {
        double theta = degrees2radians(i);
        double d = fake_laser_reading<double>(pose, theta, world);
        ScanLine<double> s = {theta,d};
        scan.push_back(s);
    }
    print_scan(scan);
}

template <class T = double>
vector<LineSegment<T>> get_world() {
    vector<LineSegment<T>> world;
    world.push_back(LineSegment<T>({0,2}, {1,2}));
    world.push_back(LineSegment<T>({-10,2}, {10,2}));
    world.push_back(LineSegment<T>({-10,-3}, {10,-3}));
    world.push_back(LineSegment<T>({10,2}, {10,-3}));
    return world;
}


void test_scan_with_twist() {
    auto world = get_world();
    double twist_theta = degrees2radians(5);
    double twist_x = .2;
    double twist_y = .1;
    double reading_count = 360;
    Pose<double> pose1;
    auto twisted = scan_with_twist2(world, reading_count, pose1, twist_x, twist_y, twist_theta);
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
    auto l = Line<double>::from_points({0,0},{v(0), v(1)});
    auto p = s.intersection(l);
    cout << "intersection at:" << p<< endl;
    if((sign_of(p(0))==sign_of(dx)) && (sign_of(p(1))==sign_of(dy))) {
    
    cout << "signs ok";
    }
}

template<class T = double>
T scan_difference(const vector<ScanLine<T>> & scan1, const vector<ScanLine<T>> & scan2) {
    double total_difference = 0;
    size_t valid_count = 0;
    for(unsigned i = 0; i < scan1.size(); i++) {
        double d = scan1[i].d-scan2[i].d;
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

template <class T>
struct TwiddleResult{
    vector<T> p;
    T error;
};

template <class T>
double abs_sum(vector<T> & v) {
    double rv = 0;
    for(auto p : v) {
        rv += fabs(p);
    }
    return rv;
}

// based loosely on https://martin-thoma.com/twiddle/
template <class T>
TwiddleResult<T> twiddle(vector<T> guess, std::function<T(const vector<T>)> f, T threshold = 0.003) {
    
    // initialize parameters to guess
    vector<T> p = guess;
    T best_error = f(p);

    // potential changes
    auto dp = std::vector<T>(guess.size(), 0.01);
    const T growth_rate = 1.5;


    while(abs_sum(dp) > threshold) {
        for(size_t i = 0; i< p.size(); i++) {
            p[i] += dp[i];
            T error = f(p);

            if (error < best_error) {
                best_error = error;
                dp[i] *= growth_rate;
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
    //cout << "twiddle done" << endl;
    TwiddleResult<T> rv;
    rv.p = p;
    rv.error = best_error;
    return rv;
}

template <class T = double>
Pose<T> match_scans(vector<ScanLine<double>> scan1, vector<ScanLine<double>> scan2) {
    auto error_function = [&scan1, &scan2](vector<double> params){
        Pose<T> pose(params[0], params[1], params[2]);
        auto scan2b = untwist_scan(scan2,0.,0.,0.,pose);
        //cout << "scan2b:" << endl;
        //print_scan(scan2b);
        double d = scan_difference(scan1, scan2b);
        // cout << "difference: " << d << " pose: " << to_string(pose) << endl;

        return d;
    };

    TwiddleResult<T> r = twiddle<T>({0,0,0}, error_function);
    Pose<T> match(r.p[0], r.p[1], r.p[2]);
    return match;
}



template <class T = double>
void move_scan(
    const vector<Point2d<T>> &scan_xy,  
    Pose<T> pose, 
    vector<Point2d<T>> &moved_scan) {
        
    move_scan_timer.start();
    Point2d<T> p, p_new;

    moved_scan.resize(0);

    for(auto & xy : scan_xy) {
        if(!isnan(xy.x)) {
            p.x = xy.x;
            p.y = xy.y;
            pose.Pose2World(p, p_new);
            //moved_scan_line.theta = atan2(p_new.y, p_new.x);

            //if(moved_scan_line.theta < 0) {
            //    moved_scan_line.theta += 2 * EIGEN_PI;
            //}
            //moved_scan_line.d = sqrt(p_new.y*p_new.y+p_new.x*p_new.x);
        } else {
            p_new.x = NAN;
            p_new.y = NAN;
        }
        moved_scan.emplace_back(p_new);
    }
    move_scan_timer.stop();
}

template <class T>
double prorate(T x, T x1, T x2, T y1, T y2) {
    return (x-x1) / (x2-x1) * (y2-y1) + y1;
}

template <class T>
inline bool is_ccw(const Point2d<T> & p1, const Point2d<T> & p2) {
    //return atan2(p1.y, p1.x) > atan2(p2.y, p2.x);
    return (p1.x*p2.y-p1.y*p2.x) >= 0;
}

// computes scan difference of scans without requiring matching equally spaced scan angles
// tbd whether there is a requirement for increasing scan angles
template<class T>
T scan_difference2(const vector<Point2d<T>> & scan1, const vector<Point2d<T>> & scan2) {
    typedef Eigen::Matrix<T,2,1> Vector2T;
    
    scan_difference_timer.start();
    // walk around scan1, finding correspondences in scan2
    T total_difference = 0;
    int points_compared = 0;

    // index of rays in scan2 to compare
    size_t i2a=0;
    size_t i2b=1;

    const size_t scan2_size = scan2.size();
    for(auto & p1 : scan1) {
        if(isnan(p1.x)){
            continue;
        }
        // find matching angle
        bool found = false;
        for(size_t i = 0; !found && i < scan2_size; ++i) {
            auto & l2_a = scan2[i2a];
            auto & l2_b = scan2[i2b];
            if(!isnan(l2_a.x) && !isnan(l2_b.x) && is_ccw(l2_a, p1) && is_ccw(p1, l2_b)) {
                Line<T> l1 = Line<T>::from_points(l2_a, l2_b);
                Line<T> l2 = Line<T>::from_points({0,0},p1);
                Vector2T v = l1.intersection(l2);
                Point2d<T> p;
                p.x = v(0);
                p.y = v(1);

                T dx = p.x-p1.x;
                T dy = p.y-p1.y;

                total_difference += sqrt(dx*dx+dy*dy);
                ++points_compared;
                found = true;
            } else {
                // increment line2 rays and wrap around if necessary
                if(++i2a == scan2_size) {
                    i2a = 0;
                }
                if(++i2b == scan2_size) {
                    i2b = 0;
                }
            }
        }
    }
    scan_difference_timer.stop();
    return total_difference / points_compared;
}

template <class T> vector<Point2d<T>> get_scan_xy(const vector<ScanLine<T>> & scan) {
    vector<Point2d<T>> scan_xy;
    scan_xy.reserve(scan.size());
    for(auto & scan_line : scan) {
        scan_xy.push_back({scan_line.d * cos(scan_line.theta),scan_line.d * sin(scan_line.theta) });
    }
    return scan_xy;
}

template <class T = double>
Pose<T> match_scans2(vector<ScanLine<T>> & scan1, vector<ScanLine<T>> & scan2) {
    match_scans_timer.start();
    auto scan1_xy = get_scan_xy(scan1);
    auto scan2_xy = get_scan_xy(scan2);
    vector<Point2d<T>> scan2b;
    scan2b.reserve(scan1.size());
    auto error_function = [&scan1_xy, &scan2_xy, &scan2b](vector<T> params){
        Pose<T> pose(params[0], params[1], params[2]);

        //vector<ScanLine<T>> scan2b;
        move_scan(scan2_xy, pose, scan2b);
        T d = scan_difference2(scan1_xy, scan2b);
        //cout << "difference: " << d << " pose: " << to_string(pose) << endl;

        return d;
    };

    TwiddleResult<T> r = twiddle<T>({0,0,0}, error_function);
    Pose<T> match(r.p[0], r.p[1], r.p[2]);
    match_scans_timer.stop();
    return match;
}


template<class T = double>
void test_match_scans2() {
    size_t n_points = 360;
    auto world = get_world<T>();

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);
    std::normal_distribution<T> x_value(0.0,0.5);
    std::normal_distribution<T> y_value(0.0,0.5);
    std::normal_distribution<T> theta_value(0.0,degrees2radians(3));

    Pose<T> pose1;
    Pose<T> pose2(x_value(gen), y_value(gen), (theta_value(gen)));
    auto scan1 = scan_with_twist<T>(world, n_points, 0, 0, 0, pose1);
    auto scan2 = scan_with_twist<T>(world, n_points, 0, 0, 0, pose2);

    auto matched_pose = match_scans2(scan1, scan2);

    cout << "actual pose: " << to_string(pose2)
         << " -> matched pose: " << to_string(matched_pose) << endl;
}

void test_prorate() {
    double y = prorate(15,10,20,10,40);
    cout << "prorate result: " << y << endl;
}

void test_match_scans() {
    size_t n_points = 360;
    auto world = get_world();
    Pose<double> pose1;
    Pose<double> pose2(.27, -.45, degrees2radians(45));
    auto scan1 = scan_with_twist<double>(world, n_points, 0, 0, 0, pose1);
    auto scan2 = scan_with_twist<double>(world, n_points, 0, 0, 0, pose2);

    auto matched_pose = match_scans(scan1, scan2);

    cout << "actual pose: " << to_string(pose2) << endl;
    cout << "match pose : " << to_string(matched_pose) << endl;
}

void test_twiddle() {
    auto lambda = [](std::vector<double> v)->double{return fabs(v[0]-3)+fabs(v[1]-5) + fabs(v[2]);};
    auto rv = twiddle<double>({0,0,0}, lambda, 1E-10);
    cout << rv.p[0] << ", " << rv.p[1] << ", " << rv.p[2] << " error: " << rv.error << endl;
        
}

template <class T>
void test_move_scan(bool trace = false) {
    auto world = get_world<T>();
    Pose<T> pose1;
    Pose<T> pose2(0, 5, degrees2radians(3));
    auto scan1 = scan_with_twist2<T>(world, 360, pose1);
    auto scan1_xy = get_scan_xy(scan1);
    vector<Point2d<T>> scan2;
    move_scan(scan1_xy, pose2, scan2);
    if(trace) {
        cout << "original scan, moved scan" << endl;
        for(unsigned i = 0; i < scan2.size(); i++) {
            cout << " -> "  <<scan1_xy[i].x <<  ", " << scan1_xy[i].y

                << " -> "  <<scan2[i].x <<  ", " << scan2[i].y  << endl;
                //<< radians2degrees(scan2[i].theta) << "°, " << scan2[i].d << endl;
        }
    }
}

int main(int, char**)
{
    //test_prorate();

    //test_move_scan<float>(true);
    //for(int i = 0; i < 100000; ++i) test_move_scan<float>();

    for(int i = 0; i < 100; ++i) test_match_scans2<float>();

    //test_twiddle();
    //return 0;
    //test_match_scans();
    //test_twiddle();
    //test_scan_with_twist();
    //test_intersection();
    // test_fake_scan();

    cout << "time untwisting: " << untwist_timer.get_elapsed_seconds() << endl;
    cout << "time moving: " << move_scan_timer.get_elapsed_seconds() << endl;
    cout << "time diffing: " << scan_difference_timer.get_elapsed_seconds() << " count: " << scan_difference_timer.start_count<< endl;
    cout << "total time matching: " << match_scans_timer.get_elapsed_seconds() << endl;
    return 0;

}