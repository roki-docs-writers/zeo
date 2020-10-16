/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_vec3d - 3D vector.
 */

#include <zeo/zeo_vec3d.h>

/* ********************************************************** */
/* CLASS: zVec3D
 * 3D vector class
 * ********************************************************** */

/* OBJECT:
 * zvec3Dzero, zvec3Dx, zvec3Dy, zvec3Dz
 * - 3D zero vector and unit vectors along (x,y,z) axis.
 */
const zVec3D zvec3Dzero = { { 0, 0, 0 } };
const zVec3D zvec3Dx    = { { 1, 0, 0 } };
const zVec3D zvec3Dy    = { { 0, 1, 0 } };
const zVec3D zvec3Dz    = { { 0, 0, 1 } };

/* zVec3DCreate
 * - create a 3D vector.
 */
zVec3D *zVec3DCreate(zVec3D *v, double x, double y, double z)
{
  zVec3DSetElem( v, zX, x );
  zVec3DSetElem( v, zY, y );
  zVec3DSetElem( v, zZ, z );
  return v;
}

/* zVec3DCreatePolar
 * - create a 3D vector by the set of value for a polar
 *   expression.
 */
zVec3D *zVec3DCreatePolar(zVec3D *v, double r, double theta, double phi)
{
  double rs;

  rs = r * sin( theta );
  return zVec3DCreate( v, rs*cos(phi), rs*sin(phi), r*cos(theta) );
}

/* zVec3DMatch
 * - check if two 3D vectors match each other.
 */
bool zVec3DMatch(zVec3D *v1, zVec3D *v2)
{
  return zVec3DElem(v1,zX) == zVec3DElem(v2,zX) &&
         zVec3DElem(v1,zY) == zVec3DElem(v2,zY) &&
         zVec3DElem(v1,zZ) == zVec3DElem(v2,zZ);
}

/* zVec3DEqual
 * - check if two 3D vectors are equal.
 */
bool zVec3DEqual(zVec3D *v1, zVec3D *v2)
{
  zVec3D err;

  return zVec3DIsTiny( zVec3DSub( v1, v2, &err ) );
}

/* zVec3DIsTol
 * - check if a 3D vector is tiny enough.
 */
bool zVec3DIsTol(zVec3D *v, double tol)
{
  return zIsTol( zVec3DElem(v,zX), tol ) &&
         zIsTol( zVec3DElem(v,zY), tol ) &&
         zIsTol( zVec3DElem(v,zZ), tol );
}

/* zVec3DIsNan
 * - check if a 3D vector includes NaN or Inf components.
 */
bool zVec3DIsNan(zVec3D *v)
{
  return zIsNan( zVec3DElem(v,zX) ) || zIsInf( zVec3DElem(v,zX) ) ||
         zIsNan( zVec3DElem(v,zY) ) || zIsInf( zVec3DElem(v,zY) ) ||
         zIsNan( zVec3DElem(v,zZ) ) || zIsInf( zVec3DElem(v,zZ) );
}

/* ********************************************************** */
/* arithmetics
 * ********************************************************** */

/* zVec3DAdd
 * - add two 3D vectors.
 */
zVec3D *zVec3DAdd(zVec3D *v1, zVec3D *v2, zVec3D *v)
{
  return zVec3DCreate( v,
    zVec3DElem(v1,zX) + zVec3DElem(v2,zX),
    zVec3DElem(v1,zY) + zVec3DElem(v2,zY),
    zVec3DElem(v1,zZ) + zVec3DElem(v2,zZ) );
}

/* zVec3DSub
 * - subtract a 3D vector from another.
 */
zVec3D *zVec3DSub(zVec3D *v1, zVec3D *v2, zVec3D *v)
{
  return zVec3DCreate( v,
    zVec3DElem(v1,zX) - zVec3DElem(v2,zX),
    zVec3DElem(v1,zY) - zVec3DElem(v2,zY),
    zVec3DElem(v1,zZ) - zVec3DElem(v2,zZ) );
}

/* zVec3DRev
 * - reverse a 3D vector.
 */
zVec3D *zVec3DRev(zVec3D *v, zVec3D *rv)
{
  return zVec3DCreate( rv,
    -zVec3DElem(v,zX), -zVec3DElem(v,zY), -zVec3DElem(v,zZ) );
}

/* zVec3DMul
 * - multiply a 3D vector by a scalar.
 */
zVec3D *zVec3DMul(zVec3D *v, double k, zVec3D *mv)
{
  return zVec3DCreate( mv,
    zVec3DElem(v,zX) * k, zVec3DElem(v,zY) * k, zVec3DElem(v,zZ) * k );
}

/* zVec3DDiv
 * - divide a 3D vector by a scalar.
 */
zVec3D *zVec3DDiv(zVec3D *v, double k, zVec3D *dv)
{
  if( k == 0 ){
    ZRUNWARN( ZEO_ERR_ZERODIV );
    return NULL;
  }
  return zVec3DMul( v, 1.0/k, dv );
}

/* zVec3DAmp
 * - amplify a 3D vector by another.
 */
zVec3D *zVec3DAmp(zVec3D *v, zVec3D *a, zVec3D *av)
{
  return zVec3DCreate( av,
    zVec3DElem(v,zX) * zVec3DElem(a,zX),
    zVec3DElem(v,zY) * zVec3DElem(a,zY),
    zVec3DElem(v,zZ) * zVec3DElem(a,zZ) );
}

/* zVec3DCat
 * - concatenate 3D vectors.
 */
zVec3D *zVec3DCat(zVec3D *v1, double k, zVec3D *v2, zVec3D *v)
{
  return zVec3DCreate( v,
    zVec3DElem(v1,zX) + zVec3DElem(v2,zX) * k,
    zVec3DElem(v1,zY) + zVec3DElem(v2,zY) * k,
    zVec3DElem(v1,zZ) + zVec3DElem(v2,zZ) * k );
}

/* zVec3DSqrNorm
 * - squared norm of a 3D vector.
 */
double zVec3DSqrNorm(zVec3D *v)
{
  return zVec3DInnerProd( v, v );
}

/* zVec3DWSqrNorm
 * - squared weighted norm of a 3D vector.
 */
double zVec3DWSqrNorm(zVec3D *v, zVec3D *w)
{
  return zSqr(zVec3DElem(v,zX)) * zVec3DElem(w,zX)
       + zSqr(zVec3DElem(v,zY)) * zVec3DElem(w,zY)
       + zSqr(zVec3DElem(v,zZ)) * zVec3DElem(w,zZ);
}

/* zVec3DSqrDist
 * - distance between two positions.
 */
double zVec3DSqrDist(zVec3D *v1, zVec3D *v2)
{
  zVec3D v;

  return zVec3DSqrNorm( zVec3DSub( v1, v2, &v ) );
}

/* zVec3DNormalizeNC
 * - normalize a 3D vector without checking vector size.
 */
double zVec3DNormalizeNC(zVec3D *v, zVec3D *nv)
{
  double l;

  l = zVec3DNorm(v);
  zVec3DDiv( v, l, nv );
  return l;
}

/* zVec3DNormalize
 * - normalize a 3D vector.
 */
double zVec3DNormalize(zVec3D *v, zVec3D *nv)
{
  if( zVec3DIsTiny( v ) ){
    ZRUNWARN( ZEO_ERR_ZERONORM );
    return -1;
  }
  return zVec3DNormalizeNC( v, nv );
}

/* zVec3DInnerProd
 * - inner product of two 3D vectors.
 */
double zVec3DInnerProd(zVec3D *v1, zVec3D *v2)
{
  return zVec3DElem(v1,zX) * zVec3DElem(v2,zX)
       + zVec3DElem(v1,zY) * zVec3DElem(v2,zY)
       + zVec3DElem(v1,zZ) * zVec3DElem(v2,zZ);
}

/* zVec3DOuterProd
 * - outer product of two 3D vectors.
 */
zVec3D *zVec3DOuterProd(zVec3D *v1, zVec3D *v2, zVec3D *v)
{
  return zVec3DCreate( v,
    zVec3DElem(v1,zY) * zVec3DElem(v2,zZ)
      - zVec3DElem(v1,zZ) * zVec3DElem(v2,zY),
    zVec3DElem(v1,zZ) * zVec3DElem(v2,zX)
      - zVec3DElem(v1,zX) * zVec3DElem(v2,zZ),
    zVec3DElem(v1,zX) * zVec3DElem(v2,zY)
      - zVec3DElem(v1,zY) * zVec3DElem(v2,zX) );
}

/* zVec3DOuterProdNorm
 * - norm of outer product of two 3D vectors.
 */
double zVec3DOuterProdNorm(zVec3D *v1, zVec3D *v2)
{
  zVec3D v;

  return zVec3DNorm( zVec3DOuterProd( v1, v2, &v ) );
}

/* zVec3DGrassmannProd
 * - scalar triple product of 3 3D vectors.
 */
double zVec3DGrassmannProd(zVec3D *v1, zVec3D *v2, zVec3D *v3)
{
  zVec3D v;

  return zVec3DInnerProd( v1, zVec3DOuterProd( v2, v3, &v ) );
}

/* zVec3DTripleProd
 * - vector triple product of 3 3D vectors.
 */
zVec3D *zVec3DTripleProd(zVec3D *v1, zVec3D *v2, zVec3D *v3, zVec3D *v)
{
  zVec3DMul( v2, zVec3DInnerProd(v1,v3), v );
  return zVec3DCatDRC( v,-zVec3DInnerProd(v1,v2), v3 );
}

/* ********************************************************** */
/* geometry
 * ********************************************************** */

/* zVec3DInterDiv
 * - interior division of two the positions.
 */
zVec3D *zVec3DInterDiv(zVec3D *v1, zVec3D *v2, double ratio, zVec3D *v)
{
  zVec3DSub( v2, v1, v );
  return zVec3DCat( v1, ratio, v, v );
}

/* zVec3DMid
 * - middle point of the two positions.
 */
zVec3D *zVec3DMid(zVec3D *v1, zVec3D *v2, zVec3D *v)
{
  return zVec3DInterDiv( v1, v2, 0.5, v );
}

/* zVec3DAngle
 * - angle between the two vectors.
 */
double zVec3DAngle(zVec3D *v1, zVec3D *v2, zVec3D *n)
{
  double c, s;
  zVec3D d;

  zVec3DOuterProd( v1, v2, &d );
  c = zVec3DInnerProd( v1, v2 );
  s = zVec3DNorm( &d );
  if( n && zVec3DInnerProd( &d, n ) < 0 ) s = -s;
  return atan2( s, c );
}

/* zVec3DProj
 * - projection of a vector.
 */
zVec3D *zVec3DProj(zVec3D *v, zVec3D *n, zVec3D *pv)
{
  zVec3D d;

  if( !zVec3DNormalize( n, &d ) ) return NULL;
  return zVec3DMul( &d, zVec3DInnerProd(&d,v), pv );
}

/* zVec3DOrthogonalize
 * - orthogonalize a vector.
 */
zVec3D *zVec3DOrthogonalize(zVec3D *v, zVec3D *n, zVec3D *ov)
{
  zVec3D tmp;

  if( !zVec3DProj( v, n, &tmp ) ) return NULL;
  return zVec3DSub( v, &tmp, ov );
}

/* zVec3DOrthoSpace
 * - create the orthogonal space of a vector.
 */
bool zVec3DOrthoSpace(zVec3D *v, zVec3D *sv1, zVec3D *sv2)
{
  zVec3DCopy( v, sv1 );
  if( zVec3DElem(v,zY) != 0 || zVec3DElem(v,zZ) != 0 )
    zVec3DElem( sv1, zX ) += 1.0;
  else
    zVec3DElem( sv1, zY ) = 1.0;
  if( !zVec3DOrthogonalize( sv1, v, sv1 ) ) return false;
  zVec3DOuterProd( v, sv1, sv2 );
  return true;
}

/* zVec3DRot
 * - rotate a 3D vector along an arbitrary axis.
 */
zVec3D *zVec3DRot(zVec3D *v, zVec3D *aa, zVec3D *rv)
{
  zVec3D v0, v1, v2, n;
  double angle;

  if( zVec3DIsTiny( aa ) ){
    zVec3DCopy( v, rv );
    return rv;
  }
  angle = zVec3DNormalize( aa, &n );
  zVec3DProj( v, &n, &v0 );
  zVec3DSub( v, &v0, &v1 );
  zVec3DOuterProd( &n, v, &v2 );
  zVec3DCat( &v0, cos(angle), &v1, rv );
  return zVec3DCatDRC( rv, sin(angle), &v2 );
}

/* ********************************************************** */
/* differential kinematics
 * ********************************************************** */

/* zVec3DDif
 * - numerical differentiation of 3D vector.
 */
zVec3D *zVec3DDif(zVec3D *v, zVec3D *vnew, double dt, zVec3D *vel)
{
  zVec3DSub( vnew, v, vel );
  return zVec3DDivDRC( vel, dt );
}

/* zVec3DZYXVel2AngVel
 * - convert z-y-x Eulerian angle differential to angular velocity.
 */
zVec3D *zVec3DZYXVel2AngVel(zVec3D *zyxvel, zVec3D *zyx, zVec3D *angvel)
{
  double sa, ca, sb, cb;

  zSinCos( zVec3DElem(zyx,0), &sa, &ca );
  zSinCos( zVec3DElem(zyx,1), &sb, &cb );
  return zVec3DZYXVel2AngVelSC( zyxvel, sa, ca, sb, cb, angvel );
}

/* zVec3DZYXVel2AngVelSC
 * - convert z-y-x Eulerian angle differential to angular velocity,
 *   directly accepting sine/cosine sets.
 */
zVec3D *zVec3DZYXVel2AngVelSC(zVec3D *zyxvel, double sa, double ca, double sb, double cb, zVec3D *angvel)
{
  return zVec3DCreate( angvel,
   -sa*zVec3DElem(zyxvel,zY) + ca*cb*zVec3DElem(zyxvel,zZ),
    ca*zVec3DElem(zyxvel,zY) + sa*cb*zVec3DElem(zyxvel,zZ),
       zVec3DElem(zyxvel,zX) -    sb*zVec3DElem(zyxvel,zZ) );
}

/* zVec3DAngVel2ZYXVel
 * - convert angular velocity to z-y-x Eulerian angle differential.
 */
zVec3D *zVec3DAngVel2ZYXVel(zVec3D *angvel, zVec3D *zyx, zVec3D *zyxvel)
{
  double sa, ca, sb, cb;

  zSinCos( zVec3DElem(zyx,0), &sa, &ca );
  zSinCos( zVec3DElem(zyx,1), &sb, &cb );
  return zVec3DAngVel2ZYXVelSC( angvel, sa, ca, sb, cb, zyxvel );
}

/* zVec3DAngVel2ZYXVelSC
 * - convert angular velocity to z-y-x Eulerian angle differential.
 *   directly accepting sine/cosine sets.
 */
zVec3D *zVec3DAngVel2ZYXVelSC(zVec3D *angvel, double sa, double ca, double sb, double cb, zVec3D *zyxvel)
{
  /* at the singular point, 'zyxvel' remains the same value. */
  if( zIsTiny(cb) ) return zyxvel;
  zVec3DSetElem( zyxvel, zZ,
    ca/cb*zVec3DElem(angvel,zX)+sa/cb*zVec3DElem(angvel,zY) );
  zVec3DSetElem( zyxvel, zX, zVec3DElem(zyxvel,zZ)*sb+zVec3DElem(angvel,zZ) );
  zVec3DSetElem( zyxvel, zY, -sa*zVec3DElem(angvel,zX)+ca*zVec3DElem(angvel,zY) );
  return zyxvel;
}

/* zVec3DZYZVel2AngVel
 * - convert z-y-z Eulerian angle differential to angular velocity.
 */
zVec3D *zVec3DZYZVel2AngVel(zVec3D *zyzvel, zVec3D *zyz, zVec3D *angvel)
{
  double sa, ca, sb, cb;

  zSinCos( zVec3DElem(zyz,0), &sa, &ca );
  zSinCos( zVec3DElem(zyz,1), &sb, &cb );
  return zVec3DZYZVel2AngVelSC( zyzvel, sa, ca, sb, cb, angvel );
}

/* zVec3DZYZVel2AngVelSC
 * - convert z-y-z Eulerian angle differential to angular velocity.
 *   directly accepting sine/cosine sets.
 */
zVec3D *zVec3DZYZVel2AngVelSC(zVec3D *zyzvel, double sa, double ca, double sb, double cb, zVec3D *angvel)
{
  return zVec3DCreate( angvel,
   -sa*zVec3DElem(zyzvel,zY) + ca*sb*zVec3DElem(zyzvel,zZ),
    ca*zVec3DElem(zyzvel,zY) + sa*sb*zVec3DElem(zyzvel,zZ),
       zVec3DElem(zyzvel,zX) +    cb*zVec3DElem(zyzvel,zZ) );
}

/* zVec3DAngVel2ZYZVel
 * - convert angular velocity to z-y-z Eulerian angle differential.
 */
zVec3D *zVec3DAngVel2ZYZVel(zVec3D *angvel, zVec3D *zyz, zVec3D *zyzvel)
{
  double sa, ca, sb, cb;

  zSinCos( zVec3DElem(zyz,0), &sa, &ca );
  zSinCos( zVec3DElem(zyz,1), &sb, &cb );
  return zVec3DAngVel2ZYZVelSC( angvel, sa, ca, sb, cb, zyzvel );
}

/* zVec3DAngVel2ZYZVelSC
 * - convert angular velocity to z-y-z Eulerian angle differential.
 *   directly accepting sine/cosine sets.
 */
zVec3D *zVec3DAngVel2ZYZVelSC(zVec3D *angvel, double sa, double ca, double sb, double cb, zVec3D *zyzvel)
{
  /* at the singular point, 'zyzvel' remains the same value. */
  if( zIsTiny(sb) ) return zyzvel;
  zVec3DSetElem( zyzvel, zZ,
    ca/sb*zVec3DElem(angvel,zX)+sa/sb*zVec3DElem(angvel,zY) );
  zVec3DSetElem( zyzvel, zX, -zVec3DElem(zyzvel,zZ)*cb+zVec3DElem(angvel,zZ) );
  zVec3DSetElem( zyzvel, zY, -sa*zVec3DElem(angvel,zX)+ca*zVec3DElem(angvel,zY) );
  return zyzvel;
}

/* ********************************************************** */
/* I/O
 * ********************************************************** */

/* zVec3DFRead
 * - input of 3D vector from file.
 */
zVec3D *zVec3DFRead(FILE *fp, zVec3D *v)
{
  zVec3DSetElem( v, zX, zFDouble( fp ) );
  zVec3DSetElem( v, zY, zFDouble( fp ) );
  zVec3DSetElem( v, zZ, zFDouble( fp ) );
  return v;
}

/* zVec3DFWrite
 * - output of 3D vector to file.
 */
void zVec3DFWrite(FILE *fp, zVec3D *v)
{
  if( !v )
    fprintf( fp, "(null 3D vector)\n" );
  else
    fprintf( fp, "{ %.10g, %.10g, %.10g }\n",
      zVec3DElem(v,zX), zVec3DElem(v,zY), zVec3DElem(v,zZ) );
}

/* zVec3DDataFWrite
 * - output of 3D vector data to file.
 */
void zVec3DDataFWrite(FILE *fp, zVec3D *v)
{
  if( !v ) return;
  fprintf( fp, "%.10g %.10g %.10g\n",
    zVec3DElem(v,zX), zVec3DElem(v,zY), zVec3DElem(v,zZ) );
}

/* METHOD:
 * zVec3DFWriteXML - xml output.
 * ... yet testing.
 */
void zVec3DFWriteXML(FILE *fp, zVec3D *v)
{
  fprintf( fp, "\"%.10g %.10g %.10g\"",
    zVec3DElem(v,zX), zVec3DElem(v,zY), zVec3DElem(v,zZ) );
}
