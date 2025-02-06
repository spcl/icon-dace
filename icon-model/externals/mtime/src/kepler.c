// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
/*
 * Description:
 *
 * Computes the solar position.
 *
 * Method:
 *
 * This routine computes three orbital parameters depending on the time
 * of the day as well as of the year (both in radians). Basic equations
 * used are the Kepler equation for the eccentric anomaly, and Lacaille's
 * formula.
 *
 * Input argument:
 *
 * pvetim - time of the year from the vernal equinox in radians
 *
 * Output arguments:
 *
 * pdisse - radius length in AU
 *          (ratio of the solar constant to its annual mean)
 * pdec   - declination of Sun
 * pra    - right ascension of Sun
 *
 * references:
 *
 * Monin, A. S.: An Introduction to the Theory of Climate
 *               D. Reidel Publishing Company, Dordrecht, 1986 (pp 10-12).
 * Meeus, J.: Astronomische Algorithmen, 2ed
 *            Johann Ambrosius Barth, Leipzig, 1994 (pp 199-222).
 *
 * S. J. Lorenz, Uni Bremen, July 1996, original version
 * S. J. Lorenz, Uni Bremen, July 1998, changed
 * U. Schlese, DKRZ, September 1998, changed
 * L. Kornblueh, MPI, December 1998, f90 rewrite
 * L. Kornblueh, MPI, February 2003, precision changes and proper commenting
 *
 */

#include <math.h>

#define _pi 3.14159265358979323846
#define _dtor _pi / 180.0

void
orbit_kepler(double pvetim, double *pdisse, double *pdec, double *pra)
{
  double cecc = 0.016715;   // Eccentricity of Earth's Orbit
  double cobld = 23.44100;  // Obliquity of Earth [Deg]
  double clonp = 282.7000;  // Longitude of Perihelion (from vernal equinox)

  // changed from 1.0e-6, due to stability problems in the solution explained in Meeus.

  double ceps = 1.0e-9;

  double zoblr, zlonpr, zsqecc, zeve, ztlonpr, zm, zeold, zenew, zeps, zcose;
  double zdisse, znu, zlambda, zsinde, zdecli;

  double za, zb, zs, zs0, zsqrt, zz;

  int iter;

  // set failed return values

  *pdisse = 0.0;
  *pdec = -999.0;
  *pra = -999.0;

  // conversion of inclination (obliquity) from degrees to rad

  zoblr = cobld * _dtor;

  // conversion of longitude of Perihelion from vernal equinox fromm degrees to rad

  zlonpr = clonp * _dtor;

  // intermediate variables

  zsqecc = sqrt((1.0 + cecc) / (1.0 - cecc));

  // calculation of eccentric anomaly of vernal equinox (Lacaille's formula)

  zeve = 2 * atan(tan(0.5 * zlonpr) / zsqecc);

  // calculation of true anomaly of vernal equinox (Kepler)

  ztlonpr = zeve - cecc * sin(zeve);

  // use Newtons method for determing the eccentric anomaly for the actual true anomaly

  // true anomaly

  zm = pvetim - ztlonpr;

  zeold = zm / (1.0 - cecc);
  zenew = zm;

  // for the iteration a first guess of the eccentric anomly is required
  // for some special cases the original assumption does not converge.
  // Following Meeus (pp. 213-214) this is covered by the following
  // calculation:

  if (cecc > 0.975)
    {

      if (fabs(zm) < 0.52359)  // M < ~30 deg in radians
        {
          za = (1.0 - cecc) / (4.0 * cecc + 0.5);
          zb = zm / (8.0 * cecc + 1.0);
          zsqrt = copysign(sqrt(zb * zb + za * za * za), zb);
          zz = (fabs(zb + zsqrt));
          zz = copysign(sqrt(zz * zz * zz), zb + zsqrt);
          zs0 = zz - 0.5 * za;
          zs = zs0 - (0.078 * zs0 * zs0 * zs0 * zs0 * zs0) / (1.0 + cecc);
          zenew = zm + cecc * (3.0 * zs - 4.0 * zs * zs * zs);
        }
    }

  // do the Newton iterations

  iter = 0;

  for (;;)
    {
      zeps = zeold - zenew;

      if (iter >= 25)
        {
          return;  // eccentric anomaly not found!
        }

      if (fabs(zeps) < ceps) break;

      iter++;

      zeold = zenew;

      zcose = cos(zenew);
      zenew = (zm + cecc * (sin(zenew) - zenew * zcose)) / (1.0 - cecc * zcose);
    }

  // calculation of distance Earth-Sun in AU

  zdisse = (1.0 / (1.0 - cecc * cos(zenew)));
  zdisse = zdisse * zdisse;

  // calculation of the true anomaly

  znu = 2.0 * atan(zsqecc * tan(zenew * 0.5));
  zlambda = znu + zlonpr;
  zsinde = sin(zoblr) * sin(zlambda);

  // calculation of the declination

  zdecli = asin(zsinde);

  // finalize return values

  *pdisse = zdisse;
  *pdec = zdecli;
  *pra = atan2(cos(zoblr) * sin(zlambda), cos(zlambda));

  return;
}
