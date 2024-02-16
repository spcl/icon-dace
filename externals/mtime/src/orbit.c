// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#include <math.h>

#include "orbit.h"

#define _pi      3.14159265358979323846
#define _dtor    _pi/180.0
#define _stor    _dtor/3600.0
#define _sectot  1.0/(86400.0*36525.0)

/*
static
int deltattable[187] = {
  1240, 1150, 1060, 980, 910, 850, 790, 740, 700, 650,
  620,  580,  550, 530, 500, 480, 460, 440, 420, 400,  
  370,  350,  330, 310, 280, 260, 240, 220, 200, 180,  
  160,  140,  130, 120, 110, 100,  90,  90,  90,  90,  
  90,   90,   90,  90, 100, 100, 100, 100, 100, 110, 
  110,  110,  110, 110, 110, 110, 110, 120, 120, 120, 
  120,  120,  130, 130, 130, 130, 140, 140, 140, 150, 
  150,  150,  150, 160, 160, 160, 160, 160, 170, 170, 
  170,  170,  170, 170, 170, 170, 160, 160, 150, 140, 
  137,  131,  127, 125, 125, 125, 125, 125, 125, 123, 
  120,  114,  106,  96,  86,  75,  66,  60,  57,  56,  
  57,   59,   62,  65,  68,  71,  73,  75,  77,  78, 
  79,   75,   64,  54,  29,  16, -10, -27, -36, -47, 
  -54,  -52,  -55, -56, -58, -59, -62, -64, -61, -47, 
  -27,    0,   26,  54,  77, 105, 134, 160, 182, 202, 
  212,  224,  235, 239, 243, 240, 239, 239, 237, 240, 
  243,  253,  262, 273, 282, 291, 300, 307, 314, 322, 
  331,  340,  350, 365, 383, 402, 422, 445, 465, 485, 
  505,  522,  538, 549, 558, 569, 580 
};
*/

static
double oblpoly[11] = {
  84381.448, -4680.93, -1.55, 1999.25, -51.38, -249.67, -39.050, 7.12, 27.87, 5.79, 2.45
};

static
double pom[4] = {
  125.0445222,  -1934.1362608,  0.00207833,  2.220e-6
};

static
double pmmoon[4] = {
  357.5277233,  35999.0503400, -0.00016030, -3.330e-6
};

static
double pmsun[4] = {
  134.9629814, 477198.8673981,  0.00869720,  1.778e-5
};

static
double pf[4] = {
  93.2719103, 483202.0175381, -0.00368250,  3.056e-6
};

static
double pd[4] = {
  297.8503631, 445267.1114800, -0.00191420,  5.278e-6
};
  
static
double abconst = 20.49552*_stor;

static
void sunlonandecc (double t, double *lon, double *ecc)
{
  double l, m, c, e, e2;

  l = (280.46646+t*(36000.76983+t*0.0003032))*_dtor;
  m = (357.52910+t*(35999.05028-t*0.0001561))*_dtor;

  e = 0.016708617-t*(0.000042040+t*0.0000001236);
  e2 = e*e;

  // equation of the center, in terms of e and m

  c = e*(2-0.25*e2)*sin(m)+1.25*e2*sin(2*m)+1.0833333333*e*e2*sin(3*m);

  *lon = l+c;
  *ecc = e;

  return;
}

static
void aberration (double t, double obl, double *ra, double *decl)
{
  double cosobl, sinobl;
  double lon, e, pir;
  double cosra, sinra;
  double cosdecl, sindecl;
  double coslon, sinlon;

  cosobl = cos(obl);
  sinobl = sin(obl);

  cosra = cos(*ra);
  sinra = sin(*ra);

  cosdecl = cos(*decl);
  sindecl = sin(*decl);

  sunlonandecc (t, &lon, &e);
      
  coslon = cos(lon);
  sinlon = sin(lon);

  // FK5 - include the e-terms
  pir = (102.93735+t*(1.71954+t*0.00046))*_dtor;
  coslon = coslon-e*cos(pir);
  sinlon = sinlon-e*sin(pir);

  *ra   = *ra-abconst*(cosra*coslon*cosobl+sinra*sinlon)/cosdecl;
  *decl = *decl-abconst*(coslon*(sinobl*cosdecl-sinra*sindecl*cosobl)+cosra*sindecl*sinlon);

  return;
}

static
void ecltoequ(double l, double b, double obl, double *ra, double *decl)
{
  double sinobl, cosobl;
  double sinl, cosl;
  double sinb, cosb;

  sinobl = sin(obl);
  cosobl = cos(obl);

  sinl = sin(l);
  cosl = cos(l);

  sinb = sin(b);
  cosb = cos(b);

  *ra = atan2(cosb*sinl*cosobl-sinb*sinobl,cosb*cosl);
  if (*ra < 0) 
    {
      *ra += 2*_pi;
    }
  *decl = asin(sinb*cosobl+cosb*sinobl*sinl);

  return;
}

static
void sphtorect(double *s, double *r)
{
  r[0] = s[2]*cos(s[0])*cos(s[1]);
  r[1] = s[2]*sin(s[0])*cos(s[1]);
  r[2] = s[2]*sin(s[1]);
  
  return;
}

static
void recttosph(double *r, double *s)
{
  s[0] = atan2(r[1],r[0]);
  if (s[0] < 0.0) 
    {
      s[0] = s[0]+2*_pi;
    }
  s[2] = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
  s[1] = asin(r[2]/s[2]);
  
  return;
}

static
void heliotogeo(double *shelio, double *searth, double *sgeo)
{
  double r[3], rearth[3];
  int i;
  
  sphtorect(shelio, r);
  sphtorect(searth, rearth);
  for (i = 0; i < 3; i++)
    {
      r[i] = r[i]-rearth[i];
    }
  recttosph(r, sgeo);

  return;
}

static
double siderealtime(double t)
{
  double sidtime, theta;

  theta = t*(360.98564736629*36525+t*(0.000387933-t/38710000));
  sidtime = fmod(((280.46061837+theta)*_dtor), 2*_pi);

  return sidtime;
};

static
double evalpoly (double *p, int n, double x)
{
  double poly;
  int i;

  poly = p[n-1];

  for (i = n-2; i >= 0; i--)
    {
      poly = poly*x+p[i];
    }

  return poly;
}

static
void nutationconst(double t, double *nutlon, double *nutobl)
{
  double om, msun, mmoon, f, d;
  double n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11;

  om    = fmod((evalpoly (pom,    4, t)*_dtor),2*_pi);
  msun  = fmod((evalpoly (pmmoon, 4, t)*_dtor),2*_pi);
  mmoon = fmod((evalpoly (pmsun,  4, t)*_dtor),2*_pi);
  f     = fmod((evalpoly (pf,     4, t)*_dtor),2*_pi);
  d     = fmod((evalpoly (pd,     4, t)*_dtor),2*_pi);

  n1 = (-0.01742 * t - 17.1996) * sin(om)            
    + (-1.3187 - 0.00016 * t) * sin(2*f-2*d+2*om) 
    + (-0.2274 - 0.00002 * t) * sin(2*f+2*om)     
    + (0.2062 + 0.00002 * t) * sin(2*om)          
    + (0.1426 - 0.00034 * t) * sin(msun)          
    + (0.0712 + 0.00001 * t) * sin(mmoon);
  n2 = + (-0.0517 + 0.00012 * t) * sin(msun+2*f-2*d+2*om) 
    + (-0.0386 - 0.00004 * t) * sin(2*f+om)            
    -0.0301 * sin(mmoon+2*f+2*om)                         
    + (0.0217 - 0.00005 * t) * sin(-msun+2*f-2*d+2*om) 
    -0.0158 * sin(mmoon-2*d)                              
    + (0.0129 + 0.00001 * t) * sin(2*f-2*d+om);
  n3 = +0.0123 * sin(-mmoon+2*f+2*om)                
    +0.0063 * sin(2*d)                            
    + (0.0063 + 0.00001 * t) * sin(mmoon+om)   
    -0.0059 * sin(-mmoon+2*f+2*d+2*om)            
    + (-0.0058 - 0.00001 * t) * sin(-mmoon+om) 
    -0.0051 * sin(mmoon+2*f+om);
  n4 = +0.0048 * sin(2*mmoon-2*d)      
    +0.0046 * sin(-2*mmoon+2*f+om)  
    -0.0038 * sin(2*f+2*d+2*om)     
    -0.0031 * sin(2*mmoon+2*f+2*om) 
    +0.0029 * sin(2*mmoon)          
    +0.0029 * sin(mmoon+2*f-2*d+2*om);
  n5 = +0.0026 * sin(2*f)                        
    -0.0022 * sin(2*f-2*d)                    
    +0.0021 * sin(-mmoon+2*f+om)              
    + (0.0017 - 0.00001 * t) * sin(2*msun) 
    +0.0016 * sin(-mmoon+2*d+om)              
    + (-0.0016 + 0.00001 * t) * sin(2*msun+2*f-2*d+2*om);
  n6 = -0.0015 * sin(msun+om)           
    -0.0013 * sin(mmoon-2*d+om)      
    -0.0012 * sin(-msun+om)          
    +0.0011 * sin(2*mmoon-2*f)       
    -0.0010 * sin(-mmoon+2*f+2*d+om) 
    -0.0008 * sin(mmoon+2*f+2*d+2*om);
  n7 = +0.0007 * sin(msun+2*f+2*om)  
    -0.0007 * sin(mmoon+msun-2*d) 
    -0.0007 * sin(-msun+2*f+2*om) 
    -0.0007 * sin(2*f+2*d+om)     
    +0.0006 * sin(mmoon+2*d)      
    +0.0006 * sin(2*mmoon+2*f-2*d+2*om);
  n8 = +0.0006 * sin(mmoon+2*f-2*d+om) 
    -0.0006 * sin(-2*mmoon+2*d+om)  
    -0.0006 * sin(2*d+om)           
    +0.0005 * sin(mmoon-msun)       
    -0.0005 * sin(-msun+2*f-2*d+om) 
    -0.0005 * sin(-2*d+om);
  n9 = -0.0005 * sin(2*mmoon+2*f+om)      
    -0.0003 * sin(-2*mmoon+2*f+2*om)   
    +0.0004 * sin(2*mmoon-2*d+om)      
    +0.0004 * sin(msun+2*f-2*d+om)     
    -0.0003 * sin(mmoon-msun+2*f+2*om) 
    -0.0003 * sin(-mmoon-msun+2*f+2*d+2*om);
  n10 =-0.0003 * sin(3*mmoon+2*f+2*om)   
    -0.0003 * sin(-msun+2*f+2*d+2*om) 
    -0.0003 * sin(mmoon-msun-d)       
    -0.0004 * sin(mmoon-d)            
    -0.0004 * sin(msun-2*d)           
    +0.0004 * sin(mmoon-2*f);
  n11 =-0.0004 * sin(d)          
    -0.0003 * sin(mmoon+msun) 
    +0.0003 * sin(mmoon+2*f);

  *nutlon = n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11;

  n1 =   (0.00089 * t + 9.2025) * cos(om)           
    + (0.5736 - 0.00031 * t) * cos(2*f-2*d+2*om) 
    + (0.0977 - 0.00005 * t) * cos(2*f+2*om)     
    + (-0.0895 + 0.00005 * t) * cos(2*om)        
    + (0.0054 - 0.00001 * t) * cos(msun)         
    -0.0007 * cos(mmoon);
  n2 = + (0.0224 - 0.00006 * t) * cos(msun+2*f-2*d+2*om)   
    +0.0200 * cos(2*f+om)                               
    + (0.0129 - 0.00001 * t) * cos(mmoon+2*f+2*om)      
    + (-0.0095 + 0.00003 * t) * cos(-msun+2*f-2*d+2*om) 
    -0.0070 * cos(2*f-2*d+om)                           
    -0.0053 * cos(-mmoon+2*f+2*om);
  n3 = -0.0033 * cos(mmoon+om)            
    +0.0026 * cos(-mmoon+2*f+2*d+2*om) 
    +0.0032 * cos(-mmoon+om)           
    +0.0027 * cos(mmoon+2*f+om)        
    -0.0024 * cos(-2*mmoon+2*f+om)     
    +0.0016 * cos(2*f+2*d+2*om);
  n4 = +0.0013 * cos(2*mmoon+2*f+2*om)    
    -0.0012 * cos(mmoon+2*f-2*d+2*om)  
    -0.0010 * cos(-mmoon+2*f+om)       
    -0.0008 * cos(-mmoon+2*d+om)       
    +0.0007 * cos(2*msun+2*f-2*d+2*om) 
    +0.0009 * cos(msun+om);
  n5 = +0.0007 * cos(mmoon-2*d+om)       
    +0.0006 * cos(-msun+om)           
    +0.0005 * cos(-mmoon+2*f+2*d+om)  
    +0.0003 * cos(mmoon+2*f+2*d+2*om) 
    -0.0003 * cos(msun+2*f+2*om)      
    +0.0003 * cos(-msun+2*f+2*om);
  n6 = +0.0003 * cos(2*f+2*d+om)           
    -0.0003 * cos(2*mmoon+2*f-2*d+2*om) 
    -0.0003 * cos(mmoon+2*f-2*d+om)     
    +0.0003 * cos(-2*mmoon+2*d+om)      
    +0.0003 * cos(2*d+om)               
    +0.0003 * cos(-msun+2*f-2*d+om);
  n7 = +0.0003 * cos(-2*d+om) 
    +0.0003 * cos(2*mmoon+2*f+om);
  
  *nutobl = n1+n2+n3+n4+n5+n6+n7;
  
  *nutlon *= _stor;
  *nutobl *= _stor;
  
  return;
}

static
double obliquity(double t)
{
  double obl;
  double u;

  u = 0.01*t;

  obl = evalpoly(oblpoly, 11, u)*_stor;

  return obl;
}

/*
static
double approxdeltat(double t)
{
  double deltat;
  double y;
  int ind;

  y = 2000+t*100;
    if (y > 2000.0) 
      {
	deltat = 102.3+t*(123.5+t*32.5);
      }
    else
      {
	if (y < 1620.0) 
	  {
	    if (y < 948.0) 
	      {
		deltat =  2715.6+t*(573.36+t*46.5);
	      }
	    else
	      {
		deltat = 50.6+t*(67.5+t*22.5);
	      }
	  }
	else
	  {
	    // interpolate from the above table
	    ind = (int)((y-1620)/2)+1;
	    if (ind > 186) 
	      {
		ind = 186;
	      }
	    y = y/2-ind-810;
	    deltat = 0.1*(deltattable[ind-1]+(deltattable[ind]-deltattable[ind-1])*y);
	  }
    }

    return deltat;
}
*/

static
double jdtot(double jd)
{
  // convert Julian Day to centuries since J2000.0.

  double jdt; // the T value corresponding to the Julian Day

  jdt = (jd-2451545.0)/36525.0;

  return jdt;
}

void orbit_vsop87 (double jde, double *ra, double *dec, double *dis, double *gha)
{
  double t, obl, nutlon, nutobl, l, b, r;

  t = jdtot(jde);
  obl = obliquity(t);
  
  nutationconst(t, &nutlon, &nutobl);

  earth_position (t, &l, &b, &r);

  double shelio[3] = { 0.0, 0.0, 0.0 };
  double searth[3] = { l , b, r };
  double sgeo[3];

  heliotogeo(shelio, searth, sgeo);

  ecltoequ(sgeo[0]+nutlon, sgeo[1], obl+nutobl, ra, dec);
  aberration (t, obl, ra, dec);

  *dis  = sgeo[2];
  *gha  = siderealtime(t);

  return;
}
