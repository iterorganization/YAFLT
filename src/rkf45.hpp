#ifndef RKF45_H
#define RKF45_H

#include <flt.hpp>

class FLT; /*Forward declaration*/
class RKF45
{
private:
    double m_abserr_save;
    int m_flag_save;
    double m_h;
    int m_init;
    int m_kflag;
    int m_kop;
    int m_nfe;
    double m_relerr_save;
    double m_remin;
    double m_r8epsilon;
public:
    RKF45();
    ~RKF45(){};
    double r8_abs(double x);
    double r8_epsilon();
    double r8_max(double x, double y);
    double r8_min(double x, double y);
    void r8_fehl(FLT *obj, double y[2], double t, double h, double yp[2],
                 double f1[2], double f2[2], double f3[2], double f4[2],
                 double f5[2], double s[2]);
    int r8_rkf45(FLT *obj, double y[2], double yp[2], double *t, double tout,
                 double *relerr, double abserr, int flag);
    double r8_sign(double x);
};

#endif /*RKF45_H*/