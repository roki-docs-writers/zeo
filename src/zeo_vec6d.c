/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_vec6d - 6D spatial vector.
 */

#include <zeo/zeo_vec6d.h>

/* ********************************************************** */
/* CLASS: zVec6D
 * 6D quasi vector class - union of linear/angular vector
 * ********************************************************** */

/* OBJECT:
 * zvec6Dzero, zvec6Dlinx, zvec6Dliny, zvec6Dlinz
 * zvec6Dangx, zvec6Dangy, zvec6Dangz
 * - 6D zero vector and unit vectors along (x,y,z) axis.
 */
const zVec6D zvec6Dzero = { { 0, 0, 0, 0, 0, 0 } };
const zVec6D zvec6Dlinx = { { 1, 0, 0, 0, 0, 0 } };
const zVec6D zvec6Dliny = { { 0, 1, 0, 0, 0, 0 } };
const zVec6D zvec6Dlinz = { { 0, 0, 1, 0, 0, 0 } };
const zVec6D zvec6Dangx = { { 0, 0, 0, 1, 0, 0 } };
const zVec6D zvec6Dangy = { { 0, 0, 0, 0, 1, 0 } };
const zVec6D zvec6Dangz = { { 0, 0, 0, 0, 0, 1 } };

/* zVec6DCreate
 * - create 6D vector.
 */
zVec6D *zVec6DCreate(zVec6D *v, double x, double y, double z, double xa, double ya, double za)
{
  zVec3DCreate( zVec6DLin(v), x , y , z  );
  zVec3DCreate( zVec6DAng(v), xa, ya, za );
  return v;
}

/* zVec3DToVec6D
 * - convert from a set of two 3D vectors to 6D vector.
 */
zVec6D *zVec3DToVec6D(zVec6D *v, zVec3D *vlin, zVec3D *vang)
{
  zVec6DSetLin( v, vlin );
  zVec6DSetAng( v, vang );
  return v;
}

/* zVec6DMatch
 * - check if the two 6D vectors match each other.
 */
bool zVec6DMatch(zVec6D *v1, zVec6D *v2)
{
  return ( zVec3DMatch( zVec6DLin(v1), zVec6DLin(v2) ) &&
           zVec3DMatch( zVec6DAng(v1), zVec6DAng(v2) ) );
}

/* zVec6DEqual
 * - check if the two 6D vectors are equal.
 */
bool zVec6DEqual(zVec6D *v1, zVec6D *v2)
{
  return ( zVec3DEqual( zVec6DLin(v1), zVec6DLin(v2) ) &&
           zVec3DEqual( zVec6DAng(v1), zVec6DAng(v2) ) );
}

/* zVec6DIsTol
 * - check if the 6D vector is tiny enough.
 */
bool zVec6DIsTol(zVec6D *v, double tol)
{
  return zVec3DIsTol( zVec6DLin(v), tol ) &&
         zVec3DIsTol( zVec6DAng(v), tol );
}

/* ********************************************************** */
/* arithmetics
 * ********************************************************** */

/* zVec6DAdd
 * - add two 6D vectors.
 */
zVec6D *zVec6DAdd(zVec6D *v1, zVec6D *v2, zVec6D *v)
{
  v->e[zX] = v1->e[zX] + v2->e[zX];
  v->e[zY] = v1->e[zY] + v2->e[zY];
  v->e[zZ] = v1->e[zZ] + v2->e[zZ];
  v->e[zXA] = v1->e[zXA] + v2->e[zXA];
  v->e[zYA] = v1->e[zYA] + v2->e[zYA];
  v->e[zZA] = v1->e[zZA] + v2->e[zZA];
  /* zVec3DAdd( zVec6DLin(v1), zVec6DLin(v2), zVec6DLin(v) ); */
  /* zVec3DAdd( zVec6DAng(v1), zVec6DAng(v2), zVec6DAng(v) ); */
  return v;
}

/* zVec6DSub
 * - subtract a 6D vector from another.
 */
zVec6D *zVec6DSub(zVec6D *v1, zVec6D *v2, zVec6D *v)
{
  v->e[zX] = v1->e[zX] - v2->e[zX];
  v->e[zY] = v1->e[zY] - v2->e[zY];
  v->e[zZ] = v1->e[zZ] - v2->e[zZ];
  v->e[zXA] = v1->e[zXA] - v2->e[zXA];
  v->e[zYA] = v1->e[zYA] - v2->e[zYA];
  v->e[zZA] = v1->e[zZA] - v2->e[zZA];

  /* zVec3DSub( zVec6DLin(v1), zVec6DLin(v2), zVec6DLin(v) ); */
  /* zVec3DSub( zVec6DAng(v1), zVec6DAng(v2), zVec6DAng(v) ); */
  return v;
}

/* zVec6DRev
 * - reverse 6D vector.
 */
zVec6D *zVec6DRev(zVec6D *v, zVec6D *rv)
{
  rv->e[0] = -v->e[0];
  rv->e[1] = -v->e[1];
  rv->e[2] = -v->e[2];
  rv->e[3] = -v->e[3];
  rv->e[4] = -v->e[4];
  rv->e[5] = -v->e[5];

  /* zVec3DRev( zVec6DLin(v), zVec6DLin(rv) ); */
  /* zVec3DRev( zVec6DAng(v), zVec6DAng(rv) ); */
  return rv;
}

/* zVec6DMul
 * - multiply 6D vector by a scalar value.
 */
zVec6D *zVec6DMul(zVec6D *v, double k, zVec6D *mv)
{
  mv->e[zX] = k * v->e[zX];
  mv->e[zY] = k * v->e[zY];
  mv->e[zZ] = k * v->e[zZ];
  mv->e[zXA] = k * v->e[zXA];
  mv->e[zYA] = k * v->e[zYA];
  mv->e[zZA] = k * v->e[zZA];

  /* zVec3DMul( zVec6DLin(v), k, zVec6DLin(mv) ); */
  /* zVec3DMul( zVec6DAng(v), k, zVec6DAng(mv) ); */
  return mv;
}

/* zVec6DDiv
 * - divide 6D vector by a scalar value.
 */
zVec6D *zVec6DDiv(zVec6D *v, double k, zVec6D *dv)
{
  if( k == 0 ){
    ZRUNWARN( ZEO_ERR_ZERODIV );
    return NULL;
  }
  k = 1.0 / k;
  dv->e[zX] = k * v->e[zX];
  dv->e[zY] = k * v->e[zY];
  dv->e[zZ] = k * v->e[zZ];
  dv->e[zXA] = k * v->e[zXA];
  dv->e[zYA] = k * v->e[zYA];
  dv->e[zZA] = k * v->e[zZA];

  /* zVec3DDiv( zVec6DLin(v), k, zVec6DLin(dv) ); */
  /* zVec3DDiv( zVec6DAng(v), k, zVec6DAng(dv) ); */
  return dv;
}

/* zVec6DCat
 * - concatenate 6D vector.
 */
zVec6D *zVec6DCat(zVec6D *v1, double k, zVec6D *v2, zVec6D *v)
{
  v->e[zX] = v1->e[zX] + k * v2->e[zX];
  v->e[zY] = v1->e[zY] + k * v2->e[zY];
  v->e[zZ] = v1->e[zZ] + k * v2->e[zZ];
  v->e[zXA] = v1->e[zXA] + k * v2->e[zXA];
  v->e[zYA] = v1->e[zYA] + k * v2->e[zYA];
  v->e[zZA] = v1->e[zZA] + k * v2->e[zZA];
  return v;
}

/* zVec6DInnerProd
 * - inner product of two 6D vectors.
 */
double zVec6DInnerProd(zVec6D *v1, zVec6D *v2)
{
  return v1->e[zX]*v2->e[zX]   + v1->e[zY]*v2->e[zY]   + v1->e[zZ]*v2->e[zZ] +
         v1->e[zXA]*v2->e[zXA] + v1->e[zYA]*v2->e[zYA] + v1->e[zZA]*v2->e[zZA];

  /* return zVec3DInnerProd( zVec6DLin(v1), zVec6DLin(v2) ) */
  /*        + zVec3DInnerProd( zVec6DAng(v1), zVec6DAng(v2) ); */
}

/* zVec6DLinShift
 * - shift velocity type of 6D vector.
 */
zVec6D *zVec6DLinShift(zVec6D *src, zVec3D *pos, zVec6D *dest)
{
  zVec3D v;

  zVec6DCopy( src, dest );
  v.e[0] = src->e[zYA] * pos->e[zZ] - src->e[zZA] * pos->e[zY];
  v.e[1] = src->e[zZA] * pos->e[zX] - src->e[zXA] * pos->e[zZ];
  v.e[2] = src->e[zXA] * pos->e[zY] - src->e[zYA] * pos->e[zX];
  dest->e[zX] += v.e[0];
  dest->e[zY] += v.e[1];
  dest->e[zZ] += v.e[2];
  return dest;

  /* zVec6DCopy( src, dest ); */
  /* return zVec6DLinShiftDRC( dest, pos ); */
}

/* zVec6DLinShiftDRC
 * - shift velocity type 6D vector directly.
 */
zVec6D *zVec6DLinShiftDRC(zVec6D *vec, zVec3D *pos)
{
  zVec3D v;

  v.e[0] = vec->e[zYA] * pos->e[zZ] - vec->e[zZA] * pos->e[zY];
  v.e[1] = vec->e[zZA] * pos->e[zX] - vec->e[zXA] * pos->e[zZ];
  v.e[2] = vec->e[zXA] * pos->e[zY] - vec->e[zYA] * pos->e[zX];
  vec->e[zX] += v.e[0];
  vec->e[zY] += v.e[1];
  vec->e[zZ] += v.e[2];

  /* zVec3DOuterProd( zVec6DAng(vec), pos, &v ); */
  /* zVec3DAddDRC( zVec6DLin(vec), &v ); */
  return vec;
}

/* zVec6DAngShift
 * - shift force type of 6D vector.
 */
zVec6D *zVec6DAngShift(zVec6D *src, zVec3D *pos, zVec6D *dest)
{
  zVec3D v;

  zVec6DCopy( src, dest );
  v.e[0] = pos->e[zY] * src->e[zZ] - pos->e[zZ] * src->e[zY];
  v.e[1] = pos->e[zZ] * src->e[zX] - pos->e[zX] * src->e[zZ];
  v.e[2] = pos->e[zX] * src->e[zY] - pos->e[zY] * src->e[zX];
  dest->e[zXA] += v.e[0];
  dest->e[zYA] += v.e[1];
  dest->e[zZA] += v.e[2];
  return dest;

  /* zVec6DCopy( src, dest ); */
  /* return zVec6DAngShiftDRC( dest, pos ); */
}

/* zVec6DAngShiftDRC
 * - shift force type of 6D vector directly.
 */
zVec6D *zVec6DAngShiftDRC(zVec6D *vec, zVec3D *pos)
{
  zVec3D v;

  v.e[0] = pos->e[zY] * vec->e[zZ] - pos->e[zZ] * vec->e[zY];
  v.e[1] = pos->e[zZ] * vec->e[zX] - pos->e[zX] * vec->e[zZ];
  v.e[2] = pos->e[zX] * vec->e[zY] - pos->e[zY] * vec->e[zX];
  vec->e[zXA] += v.e[0];
  vec->e[zYA] += v.e[1];
  vec->e[zZA] += v.e[2];

  /* zVec3D n; */
  /* zVec3DOuterProd( pos, zVec6DLin(vec), &n ); */
  /* zVec3DAddDRC( zVec6DAng(vec), &n ); */
  return vec;
}

/* ********************************************************** */
/* differential kinematics
 * ********************************************************** */

/* zVec6DDif
 * - numerical differentiation of 6D vector.
 */
zVec6D *zVec6DDif(zVec6D *v, zVec6D *vnew, double dt, zVec6D *vel)
{
  zVec6DSub( vnew, v, vel );
  return zVec6DDivDRC( vel, dt );
}

/* ********************************************************** */
/* I/O
 * ********************************************************** */

/* zVec6DFRead
 * - input of 6D vector from file.
 */
zVec6D *zVec6DFRead(FILE *fp, zVec6D *v)
{
  register int i;

  for( i=zX; i<=zZA; i++ )
    v->e[i] = zFDouble( fp );
  return v;
}

/* zVec6DDataFWrite
 * - output of 6D vector data to file.
 */
zVec6D *zVec6DDataFWrite(FILE *fp, zVec6D *v)
{
  if( !v ) return NULL;
  fprintf( fp, " %.10g %.10g %.10g %.10g %.10g %.10g",
    v->e[zX], v->e[zY], v->e[zZ], v->e[zXA], v->e[zYA], v->e[zZA] );
  return v;
}

/* zVec6DDataNLFWrite
 * - output of 6D vector data with the new line to file.
 */
zVec6D *zVec6DDataNLFWrite(FILE *fp, zVec6D *v)
{
  if( !zVec6DDataFWrite( fp, v ) ) return NULL;
  fprintf( fp, "\n" );
  return v;
}

/* zVec6DFWrite
 * - output of 6D vector to file.
 */
zVec6D *zVec6DFWrite(FILE *fp, zVec6D *v)
{
  if( !v ){
    fprintf( fp, "(null 6D vector)\n" );
    return NULL;
  }
  zVec3DFWrite( fp, zVec6DLin(v) );
  zVec3DFWrite( fp, zVec6DAng(v) );
  return v;
}

/* METHOD:
 * zVec6DFWriteXML - xml output.
 * ... yet testing.
 */
void zVec6DFWriteXML(FILE *fp, zVec6D *v)
{
  fprintf( fp, "\"" );
  zVec6DDataFWrite( fp, v );
  fprintf( fp, "\"" );
}
