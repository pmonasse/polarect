#include "polarect.h"
#include <limits>
#include <iostream>
#include <cstdlib>
#include <cassert>

typedef libNumerics::matrix<double> Mat;
typedef libNumerics::vector<double> Vec;

/// Return left epipole of matrix F. The returned vector is of norm 1.
static Vec leftEpipole(const Mat& F) {
    Vec e(3), e2(3);
    double norm=-1;
    for(int i=0; i<3; i++)
        for(int j=i+1; j<3; j++) {
            e2 = cross(F.col(i),F.col(j));
            double norm2=e2.qnorm();
            if(norm < norm2) {
                norm = norm2;
                e = e2;
            }
        }
    return (e/sqrt(norm));
}

/// Normalize point:
/// - if finite, set e(2) to +1/-1 (sign unchanged)
/// - if infinite (too far from O), set e(2) to 0 and set norm to 1.
/// The point is considered infinite if its distance to a 1000x1000 square is
/// large enough to provoke an error at most 1 pixel on the opposite side if
/// sent to infinity.
static void normalize_epipole(Vec& e) {
    const double TH=1.0/(1000*1000);
    if(e(2)*e(2) < TH*TH*(e(0)*e(0)+e(1)*e(1))) {
        if(e(0)<0)
            e = -e;
        e(2) = 0;
        e /= sqrt( e.qnorm() );
    } else
        e /= fabs(e(2));
}

/// Return epipoles of F oriented such that
/// eL x m_i(L) ~+ F m_i(R) and eR x m_i(R) ~+ Ft m_i(L)
/// when m_i(L) and m_i(R) have positive 3rd coordinate (typically 1).
bool orientedEpipoles(const std::vector<Match>& m, const Mat& F,
                      Vec& eL, Vec& eR) {
    eL = leftEpipole(F);
    eR = leftEpipole(F.t());
    if(eL(2)<0) eL = -eL;
    if(eR(2)<0) eR = -eR;
    int nL=0, nR=0;
    for(std::vector<Match>::const_iterator it=m.begin(); it!=m.end(); ++it) {
        Vec xL(it->x1, it->y1, 1.0), xR(it->x2, it->y2, 1.0);
        Vec exL=cross(eL,xL), exR=cross(eR,xR);
        Vec FxR=F*xR, FtxL=F.t()*xL;
        nL += (dot(exL, FxR)  >= 0)? +1: -1;
        nR += (dot(exR, FtxL) >= 0)? +1: -1;
    }
    if(nL<0) eL=-eL;
    if(nR<0) eR=-eR;
    normalize_epipole(eL);
    normalize_epipole(eR);
    if((size_t)abs(nL)!=m.size() || (size_t)abs(nR)!=m.size())
        std::cout << "orientEpipoles: "
                  << m.size()-(size_t)abs(nL) << ',' << m.size()-(size_t)abs(nR)
                  << " (/" << m.size() << ") matches badly oriented"<<std::endl;
    return (nL!=0 && nR!=0);
}

/// Polar transform of vector (x,y). Output is (rho,theta).
/// Special case: center at infinity, direction of line (c(0),c(1)).
/// Then theta is the (signed) distance from (0,0) to the line through (x,y) and
/// rho is the abscissa along the line wrt projection of (0,0) on the line.
std::pair<double,double> Polarizer::polar(double x, double y) const {
    double rho,theta;
    if(inf()) {
        rho   =  c(0)*x + c(1)*y;
        theta = -c(1)*x + c(0)*y;
    } else {
        Vec v(x-c(0)/c(2), y-c(1)/c(2));
        rho = sqrt( v.qnorm() );
        theta = 0;
        if(rho>0) {
            v /= rho;
            theta = atan2(v(1),v(0));
        }
    }
    return std::make_pair(rho,theta);
}

/// Convert rho to abscissa of rectified image.
double Polarizer::sample_rho(double rho) const {
    rho -= r;
    if(rotate)
       rho = width()-1-rho;
    return rho;
}

/// Convert theta to ordinate of rectified image.
double Polarizer::sample_theta(double theta) const {
    if(!inf()) { // Handle theta modulo 2*pi
        if(theta<t)
            theta += 2*M_PI;
        else if(theta>T)
            theta -= 2*M_PI;
    }
    theta -= t;
    if(! inf()) // Finite center: delta_theta = 1/R
        theta *= R;
    if(rotate)
        theta = height()-1-theta;
    return theta;
}

/// Width of polar image with rho in [r,R].
int Polarizer::width() const {
    int w = 0;
    for(double rho=r; rho<R; rho++)
        ++w;
    return w;
}

/// Height of polar image with theta in [t,T].
int Polarizer::height() const {
    double deltaT = inf()? 1: 1/R;
    int h = 0;
    for(double theta=t; theta<T; theta+=deltaT)
        ++h;
    return h;
}

/// Transfer epipolar line indexed by \a theta with \a F to view \a P.
/// Index \a theta can be an angle (finite center) or signed distance to O
/// (center at infinity). The result is cos/sin of line direction (center of
/// \a P finite) or projection of O on line (center of \a P at infinity).
/// The current object is in the left view of \a F and \a P in the right view.
std::pair<double,double> Polarizer::transfer_theta(double theta,
                                                   const Polarizer& P,
                                                   const Mat* F) const {
    Vec p(3);
    if(inf()) { // center of this at infinity
        p = Vec(-theta*c(1), theta*c(0), 1.0); // Projection of O on line
        if(! F)
            return std::make_pair(p(0), p(1));
    } else { // finite center
        double dx=cos(theta), dy=sin(theta);
        if(! F)
            return std::make_pair(dx,dy);
        p = Vec(c(0)+R*dx*c(2), c(1)+R*dy*c(2), c(2));
    }
    p = (*F)*p;
    if(P.inf()) {
        theta = p(2)/(p(0)*P.c(1)-p(1)*P.c(0));
        return std::make_pair(-theta*P.c(1),theta*P.c(0));
    }
    p(2) = 0;
    if(c(2)*P.c(2)<0) p = -p;
    double rho = sqrt( p.qnorm() );
    return std::make_pair(p(1)/rho, -p(0)/rho);
}

/// Clamp value of \a x in range [m,M]
static double clamp(double x, double m, double M) {
    if(x<m) return m;
    if(x>M) return M;
    return x;
}

/// Distance from center to filled rectangle [0,w]x[0,h].
double Polarizer::dist_rect(int w, int h) const {
    if(inf())
        return 0;
    Vec v(clamp(c(0)/c(2),0,w) - c(0)/c(2), clamp(c(1)/c(2),0,h) - c(1)/c(2));
    return sqrt( v.qnorm() );
}

/// The two extreme corners of image viewed from each region.
/// The regions wrt an image in rectangle [0,w]x[0,h] are:
/// (0,0) 0 (1,0) w (2,0)
///     0 --------- 0
/// (0,1) | (1,1) | (2,1)
///     h --------- h
/// (0,2) 0 (1,2) w (2,2)
/// Corners are numbered from 0 (top-left) to 3, turning clockwise.
/// The exception is region (1,1), inside the image, where the extreme
/// corners do not exist. \sa region.
const std::pair<int,int> reg2corners[3*3] = {
    std::make_pair(1,3),std::make_pair(1,0),std::make_pair(2,0),
    std::make_pair(0,3),std::make_pair(0,0),std::make_pair(2,1),
    std::make_pair(0,2),std::make_pair(3,2),std::make_pair(3,1)};

/// Return quadrant of center wrt rectangle [0,w]x[0,h].
std::pair<int,int> Polarizer::region(int w, int h) const {
    std::pair<int,int> p = std::make_pair(1,1);
    if(inf()) {
        p.first = 0;
        p.second = (c(1)>=0)? 0: 2;
        return p;
    }
    Vec d = c(2)*c;
    if(d(0)<=0)      p.first = 0;
    if(d(0)>=w*d(2)) p.first = 2;
    if(d(1)<=0)      p.second = 0;
    if(d(1)>=h*d(2)) p.second = 2;
    return p;
}

/// Return min/max radius in \a r and \a R and min/max angle in \a t and \a T.
Polarizer::Polarizer(const Vec& center, int w, int h)
: c(center) {
    std::pair<int,int> reg = region(w, h);
    rotate = (reg.first==2) || (reg.first==1 && reg.second==2);
    std::pair<double,double> pts[4] =
        {std::make_pair(0.,0.), std::make_pair(w*1.,0.),
         std::make_pair(w*1.,h*1.), std::make_pair(0.,h*1.)};
    r = R = dist_rect(w, h);
    for(int i=0; i<4; i++) {
        std::pair<double,double> rhoTheta = polar(pts[i].first, pts[i].second);
        if(r > rhoTheta.first) r = rhoTheta.first;
        if(R < rhoTheta.first) R = rhoTheta.first;
        if(i == reg2corners[3*reg.second+reg.first].first)
            t = rhoTheta.second;
        if(i == reg2corners[3*reg.second+reg.first].second)
            T = rhoTheta.second;
    }
    if(t>=T)
        T += 2*M_PI;
    std::cout << "Region: (" <<reg.first << ',' <<reg.second << ')' <<std::endl;
}

/// Intersection of interval [t1,T1] with [t2,T2] or [T2,t2] modulo 2pi.
/// The intervals on the circle are interpreted as width less than pi. 
static void inter_mod_2pi(double& t1, double& T1, double t2, double T2,
                          bool lessPi) {
    std::cout << "[" << t1 << ',' << T1 << "] inter "
              << "[" << t2 << ',' << T2 << "] = " << std::flush;
    assert(-M_PI<=t1 && t1<=M_PI && t1<=T1);
    assert(-M_PI<=t2 && t2<=M_PI);
    if(t2>T2)
        std::swap(t2,T2);
    if(T2-t2<10*std::numeric_limits<double>::epsilon()) { // Full circle
        t2=t1; T2=T1;
    } else if(lessPi) {
        if(T2>t2+M_PI) {
            t2 += 2*M_PI;
            std::swap(t2,T2);
        }
    } else {
        if(T2<t2+M_PI) {
            T2 -= 2*M_PI;
            std::swap(t2,T2);
        }
    }
    if(T2<t1) {
        t2 += 2*M_PI;
        T2 += 2*M_PI;
    }
    if(t2>T1) {
        t2 -= 2*M_PI;
        T2 -= 2*M_PI;
    }
    if(t1<t2) t1 = t2;
    if(T1>T2) T1 = T2;
    std::cout << "[" << t1 << ',' << T1 << "]" << std::endl;
}

/// Restrict angle interval by epipolar lines in \a polR transferred from \a F.
/// If the center of \c this is at infinity, angle -> signed distance to O.
void Polarizer::restrict_angles(const Polarizer& polR, const Mat& F) {
    std::pair<double,double> p;
    p = polR.transfer_theta(polR.t, *this, &F);
    double t2 = inf()? -c(1)*p.first+c(0)*p.second: atan2(p.second,p.first);
    p = polR.transfer_theta(polR.T, *this, &F);
    double T2 = inf()? -c(1)*p.first+c(0)*p.second: atan2(p.second,p.first);
    if(inf()) {
        if(t2>T2)
            std::swap(t2,T2);
        if(t<t2) t = t2;
        if(T>T2) T = T2;
    } else
        inter_mod_2pi(t, T, t2, T2, (polR.T <= polR.t+M_PI));
}

/// Apply mirror to each line.
static void mirrorx(std::pair<double,double>* map, int w, int h) {
    for(int i=0; i<h; i++)
        for(std::pair<double,double> *p=map+i*w, *q=map+(i+1)*w-1; p<q; p++,q--)
            std::swap(*p,*q);
}

/// Apply mirror to each column.
static void mirrory(std::pair<double,double>* map, int w, int h) {
    for(std::pair<double,double> *p=map, *q=map+(h-1)*w; p<q; p+=w,q-=w)
        for(int i=0; i<w; i++)
            std::swap(p[i],q[i]);
}

/// Perform x-mirror of pullback \a map if required to preserve orientation.
static bool restore_orientation(std::pair<double,double>* map, int w, int h) {
    std::pair<double,double> pp1 = map[(  h/3)*w+(  w/3)];
    std::pair<double,double> pp2 = map[(  h/3)*w+(2*w/3)];
    std::pair<double,double> pp3 = map[(2*h/3)*w+(  w/3)];
    std::pair<double,double> v12 = std::make_pair(pp2.first-pp1.first,
                                                  pp2.second-pp1.second);
    std::pair<double,double> v13 = std::make_pair(pp3.first-pp1.first,
                                                  pp3.second-pp1.second);
    bool rotate = (v12.first*v13.second<v12.second*v13.first);
    if(rotate)
        mirrorx(map, w, h);
    return rotate;
}

/// Generate pullback map to transform image to polar (x=radius,y=theta).
/// When pol!=this and F!=0, the angles are regularly sampled in the other image
/// of the stereo pair and transferred through the fundamental matrix F.
std::pair<double,double>* Polarizer::pullback(int w, int h,
    const Polarizer& pol, const Mat* F) {
    std::pair<double,double> *pb=new std::pair<double,double>[w*h], *p=pb;
    double deltaT = pol.inf()? 1: 1/pol.R;
    std::pair<double,double> p0=std::make_pair(c(0)/c(2),c(1)/c(2));
    for(int y=0; y<h; y++) {
        double theta=pol.t+y*deltaT;
        std::pair<double,double> cossin = pol.transfer_theta(theta, *this, F);
        double dx=cossin.first, dy=cossin.second;
        if(inf()) {
            p0 = cossin;
            dx=c(0); dy=c(1);
        }
        for(int x=0; x<w; x++) {
            double rho=r+x;
            *p++ = std::make_pair(p0.first+rho*dx, p0.second+rho*dy);
        }
    }
    if(pol.rotate)
        mirrory(pb, w, h);
    rotate = restore_orientation(pb, w, h);
    return pb;
}

/// Constructor
Polarectifyer::Polarectifyer(const Mat& F0, const Vec& eL, const Vec& eR,
                             int wL, int hL, int wR, int hR)
: polL(eL,wL,hL), polR(eR,wR,hR), F(F0), pbMapL(0), pbMapR(0) {
    std::cout << "eL =" << eL << std::endl;
    std::cout << "eR =" << eR << std::endl;
    polL.restrict_angles(polR, F);
    polR.restrict_angles(polL, F.t());

    wL=widthL(); wR=widthR(); hL=height(); // Same height of polar images
    std::cout << "ImageL: " << wL << 'x' << hL << std::endl;
    std::cout << "ImageR: " << wR << 'x' << hL << std::endl;
    pbMapL = polL.pullback_map(wL, hL);
    pbMapR = polR.pullback_map(wR, hL, polL, F.t());
}

/// Destructor
Polarectifyer::~Polarectifyer() {
    delete [] pbMapL;
    delete [] pbMapR;
}

/// Return pullback map of left or right rectified view.
const std::pair<double,double>* Polarectifyer::pullback_map(bool left) const {
    return (left? pbMapL: pbMapR);
}

/// Transform point from original image to mapped point in rectified image.
void Polarectifyer::pushforward(std::pair<double,double>& xy, bool left) const {
    const Polarizer& pol = (left? polL: polR);
    double theta=0;
    if(! left) { // Handle non-uniform theta sampling in right image
        Vec p(xy.first, xy.second, 1.);
        p = F*p;
        if(polL.inf())
            theta = p(2)/(p(0)*polL.c(1)-p(1)*polL.c(0));
        else {
            if(polL.c(2)<0) p = -p;
            theta = atan2(-p(0),p(1));
        }
    }
    xy = pol.polar(xy.first, xy.second);
    if(! left)
        xy.second = theta;
    xy.first = pol.sample_rho(xy.first);
    xy.second = polL.sample_theta(xy.second);
}
