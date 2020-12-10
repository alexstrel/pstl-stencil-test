#pragma once

#ifndef M_PI
#define M_PI (3.1415926535897932384626)
#endif

constexpr double _pi  = 3.1415926535897932384626;
constexpr double _2pi = (2*_pi);

template<typename T>
void create_field(std::vector<T> &out,
                  const std::array<int, 3> nl,
                  const std::array<T, 3> dnl,
                  const T kappa,
                  const T length,
                  const T time) {
  const T exponent  = exp(-3.0*kappa*_pi*_pi*time/(length*length));
  T z = dnl[2];
  for (int k = 0; k < nl[2]; k++){
    T y = dnl[1];
    for (int j = 0; j < nl[1]; j++){
      T x = dnl[0];
      for (int i = 0; i < nl[0]; i++){
        int s = k*nl[0]*nl[1] + j*nl[0] + i;
        out[s] = exponent * sin(_pi*x/length) * sin(_pi*y/length)* sin(_pi*z/length);
        x += dnl[0];
      }
      y += dnl[1];
    }
    z += dnl[2];
  }
  return;
}
/*
typedef struct Double2 {
  double x;
  double y;

  Double2() : x(0.0), y(0.0) {};
  Double2(double x_, double y_) : x(x_), y(y_) {};
  Double2(const Double2 &reg) : x(reg.x), y(reg.y) {};

  Double2 operator=(Double2 &val)  {return Double2(val);}

} double2;

double2 operator+(const double2 &x, const double2 &y) {
  double2 res((x.x+y.x), (x.y+y.y));
  return res;
}
*/
template<typename T>
T accuracy(const std::vector<T> &in1, std::vector<T> &in2) {

#if 1
  double err   = 0.0;
  double norm  = 0.0;

  for (int i = 0; i < in1.size(); i++){
    err  += (in1[i] - in2[i])*(in1[i] - in2[i]);
    norm += in1[i]*in1[i];
  }
  return (T)sqrt(err/norm);
#else
  auto policy = std::execution::par_unseq;
  double2 r   = std::transform_reduce(policy,
                        in1.begin(),
                        in1.end(),
                        in2.begin(),
                        double2{0.0, 0.0},
                        [=](auto &accum, auto &local_prod ) {return accum+local_prod;},
                        [=](auto &i1, auto &i2) { return double2((i1-i2)*(i1-i2), i1*i1);}
);
  return (T)sqrt(r.x/r.y);
#endif
}
