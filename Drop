#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "log-conform.h"
#include "curvature.h"
#include "tension.h"

#define rho_r 831.7
#define rho_a 1.
#define Re 0.043
#define Re_a 66.7
#define Fr 102.
#define Wi 1.
#define We 0.019
#define LEVEL 8

scalar lambdav[], mupv[];

u.n[right] = neumann(0);
p[right]   = dirichlet(0);
u.t[left] = dirichlet(0);
tau_qq[left] = dirichlet(0);
f[left] = 0.;

int main() {
  size (3);
  init_grid (1 << LEVEL);

  rho1 = rho_r;
  rho2 = rho_a;
  mu2 = 1./Re_a;
  mu1 = 1./Re;
  mup = mupv;
  lambda = lambdav;
  f.sigma = 1./We;

  dt = 0.001;
  run();
}

event init (t = 0) {
  scalar s = tau_p.x.x;
  s[left] = dirichlet(0.);
  fraction (f, - sq(x - 0.6) - sq(y) + sq(0.5));
  foreach()
    u.x[] = - f[];
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 1./Fr*f[];
}

event properties (i++) {
  foreach() {
    mupv[] = clamp(f[],0,1) * 1./Re;
    lambdav[] = clamp(f[],0,1) * Wi;
  }
}

event adapt (i++) {
  adapt_wavelet ({f,u.x,u.y}, (double[]){1e-3,5e-4,5e-4},
maxlevel = LEVEL, minlevel = LEVEL - 2);
}

event images (t += 1./100.) {
  output_ppm (f, linear = true);
}

event printdata(i += 1; t <= 1.9) {
  FILE *fpx = fopen("data1.txt","a");
  scalar posx[] = 0;
  position (f, posx, {0,1});
  fprintf(fpx, "%g\t%g\n",t, 2.*statsf(posx).max);
}
