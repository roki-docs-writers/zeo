/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_vec6d - 6D spatial vector.
 */

#ifndef __ZEO_VEC6D_H__
#define __ZEO_VEC6D_H__

#include <zeo/zeo_vec3d.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zVec6D
 * spatial 6D vector class = a pair of linear/angular vectors
 * ********************************************************** */

typedef union {
  double e[6];
  zVec3D v[2];
} zVec6D;

#define zVec6DArray(u)       ( (u)->e )
#define zVec6DElem(u,i)      ( zVec6DArray(u)[i] )
#define zVec6DSetElem(u,i,e) ( zVec6DElem(u,i) = (e) )
#define zVec6DLin(u)         ( &(u)->v[0] )
#define zVec6DAng(u)         ( &(u)->v[1] )
#define zVec6DSetLin(u,l)    zVec3DCopy( l, zVec6DLin(u) )
#define zVec6DSetAng(u,r)    zVec3DCopy( r, zVec6DAng(u) )

/* OBJECT:
 * zvec6Dzero, zvec6Dlinx, zvec6Dliny, zvec6Dlinz
 * zvec6Dangx, zvec6Dangy, zvec6Dangz
 * - 6D zero vector and unit vectors along (x,y,z) axis.
 */
extern const zVec6D zvec6Dzero;
extern const zVec6D zvec6Dlinx;
extern const zVec6D zvec6Dliny;
extern const zVec6D zvec6Dlinz;
extern const zVec6D zvec6Dangx;
extern const zVec6D zvec6Dangy;
extern const zVec6D zvec6Dangz;
#define ZVEC6DZERO ( (zVec6D *)&zvec6Dzero )
#define ZVEC6DLINX ( (zVec6D *)&zvec6Dlinx )
#define ZVEC6DLINY ( (zVec6D *)&zvec6Dliny )
#define ZVEC6DLINZ ( (zVec6D *)&zvec6Dlinz )
#define ZVEC6DANGX ( (zVec6D *)&zvec6Dangx )
#define ZVEC6DANGY ( (zVec6D *)&zvec6Dangy )
#define ZVEC6DANGZ ( (zVec6D *)&zvec6Dangz )

/* METHOD:
 * zVec6DCreate, zVec3DToVec6D,
 * zVec6DCopy, zVec6DClear
 * - creation, copy and cleanup of 6D vector.
 *
 * 'zVec6DCreate()' creates a 6D vector 'v' which consists
 * of 'x', 'y', 'z', 'xa', 'ya', and 'za'.
 *
 * 'zVec3DToVec6D()' creates a 6D vector 'v' from two 3D
 * vectors, 'vlin' and 'vang'.
 *
 * 'zVec6DCopy()' copies 'src' to 'dest'.
 *
 * 'zVec6DClear()' clears 'v', or sets all factors for 0.
 * [RETURN VALUE]
 * 'zVec6DCreate()', zVec3DToVec6D and
 * 'zVec6DClear()' return a pointer to 'v'.
 *
 * 'zVec6DCopy()' returns no value.
 * [NOTES]
 * It is also possible to write simply *dest = *src instead of
 * 'zVec6DCopy()'. Actually, 'zVec6DCopy()' is defined
 * as macro(see "zeo_vec_vec6d.h").
 *
 * 'zVec6DClear()' is a macro for
 * zVec6DCreate( v, 0, 0, 0, 0, 0, 0 ).
 */
__EXPORT zVec6D *zVec6DCreate(zVec6D *v, double x, double y, double z, double xa, double ya, double za);
__EXPORT zVec6D *zVec3DToVec6D(zVec6D *v, zVec3D *v1, zVec3D *v2);
#define zVec6DCopy(src,dest) ( *(dest) = *(src) )
#define zVec6DClear(v)       zVec6DCreate( v, 0, 0, 0, 0, 0, 0 )

/* METHOD:
 * zVec6DMatch, zVec6DEqual
 * - check if the two 6D vectors are equal.
 *
 * 'zVec6DMatch()' and 'zVec6DEqual()' check if the two
 * 6D vectors, 'v1' and 'v2', are equal. They return
 * a boolean value as a result.
 *
 * 'zVec6DMatch()' strictly compares the two vectors,
 * while 'zVec6DEqual()' checks if the error between
 * 'v1' and 'v2' are sufficiently small.
 * [RETURN VALUE]
 * 'zVec6DMatch()' and 'zVec6DEqual()' return the true
 * value if 'v1' and 'v2' are equal, or false otherwise.
 * [NOTES]
 * It is also possible to write just *v1 == *v2, instead
 * of calling 'zVec6DMatch( v1, v2 )'.
 */
__EXPORT bool zVec6DMatch(zVec6D *v1, zVec6D *v2);
__EXPORT bool zVec6DEqual(zVec6D *v1, zVec6D *v2);

/* METHOD:
 * zVec6DIsTol, zVec6DIsTiny
 * - check if 6D vector is tiny.
 *
 * 'zVec6DIsTol()' checks if the absolute values of every
 * components of 6D vector 'v' are smaller than 'tol'.
 *
 * 'zVec6DIsTiny()' applies zTOL (defined in "zeo_misc.h")
 * to the tolerance of 'zVec6DIsTol()'.
 * [RETURN VALUE]
 * 'zVec6DIsTol()' and 'zVec6DIsTiny()' return the
 * true value when the absolute values of every components
 * of 'v' are smaller than 'tol' and zTOL, respectively,
 * or the false value, otherwise.
 * [NOTES]
 * 'tol' must be positive.
 * [SEE ALSO]
 * zIsTol, zIsTiny
 */
__EXPORT bool zVec6DIsTol(zVec6D *v, double tol);
#define zVec6DIsTiny(v) zVec6DIsTol( v, zTOL )

/* ********************************************************** */
/* arithmetics
 * ********************************************************** */

/* METHOD:
 * zVec6DAdd, zVec6DSub, zVec6DRev, zVec6DMul,
 * zVec6DDiv, zVec6DCat,
 * zVec6DAddDRC, zVec6DSubDRC, zVec6DRevDRC,
 * zVec6DMulDRC, zVec6DDivDRC, zVec6DCatDRC
 * - four rules of the arithmetics for 6D vector.
 *
 * 'zVec6DAdd()' adds the two 6D vectors, 'v1' and 'v2'.
 * The result is put into 'v'.
 *
 * 'zVec6DSub()' subtracts the 6D vector 'v2' from
 * the other 'v1'. The result is put into 'v'.
 *
 * 'zVec6DRev()' reverses the 6D vector 'v'. The result is
 * put into 'rv'.
 *
 * 'zVec6DMul()' multiplies the 6D vector 'v' by value
 * 'k'. The result is put into 'mv'.
 *
 * 'zVec6DDiv()' divides the 6D vector 'v' by 'k'.
 * The result is put into 'dv'.
 *
 * 'zVec6DCat()' concatenates the 6D vector 'v2' to 'v1',
 * multiplied by a scalar value 'k'. The result is put into 'v'.
 *
 * 'zVec6DAddDRC()' directly adds 'v2' to 'v1'.
 *
 * 'zVec6DSubDRC()' directly subtracts 'v2' from 'v1'.
 *
 * 'zVec6DRevDRC()' directly reverses 'v'.
 *
 * 'zVec6DMulDRC()' directly multiplies 'v' by 'k'.
 *
 * 'zVec6DDivDRC()' directly divides 'v' by 'k'.
 *
 * 'zVec6DCat()' directly concatenates 'v2' multiplied 'v2'
 * by 'k' to 'v1'.
 * [RETURN VALUE]
 * Each function returns a pointer to the resultant vector.
 *
 * And, 'zVec6DDiv()' and 'zVec6DDivDRC()' return the
 * NULL pointer if the given scalar value is 0.
 */
__EXPORT zVec6D *zVec6DAdd(zVec6D *v1, zVec6D *v2, zVec6D *v);
__EXPORT zVec6D *zVec6DSub(zVec6D *v1, zVec6D *v2, zVec6D *v);
__EXPORT zVec6D *zVec6DRev(zVec6D *v, zVec6D *rv);
__EXPORT zVec6D *zVec6DMul(zVec6D *v, double k, zVec6D *mv);
__EXPORT zVec6D *zVec6DDiv(zVec6D *v, double k, zVec6D *dv);
__EXPORT zVec6D *zVec6DCat(zVec6D *v1, double k, zVec6D *v2, zVec6D *v);

#define zVec6DAddDRC(v1,v2)   zVec6DAdd(v1,v2,v1)
#define zVec6DSubDRC(v1,v2)   zVec6DSub(v1,v2,v1)
#define zVec6DRevDRC(v)       zVec6DRev(v,v)
#define zVec6DMulDRC(v,k)     zVec6DMul(v,k,v)
#define zVec6DDivDRC(v,k)     zVec6DDiv(v,k,v)
#define zVec6DCatDRC(v1,k,v2) zVec6DCat(v1,k,v2,v1)

/*! \brief inner product of two 6D vectors.
 *
 * zVec6DInnerProd() returns the inner product of two 6D vectors
 * \a v1 and \a v2.
 * \retval \a v1 ^T \a v2
 */
__EXPORT double zVec6DInnerProd(zVec6D *v1, zVec6D *v2);

/* METHOD:
 * zVec6DLinShift, zVec6DLinShiftDRC,
 * zVec6DAngShift, zVec6DAngShiftDRC
 * - shift 6D vector in 3D space.
 *
 * 'zVec6DLinShift()' shifts the velocity type of 6D vector
 * 'src' at the original point to the equivalent 6D vector
 * 'dest' at the point 'pos', namely:
 *  v1 = v0 + w0 x 'pos'
 *  w1 = w0
 * where 'src'=[v0 w0]^T and 'dest'=[v1 w1]^T.
 *
 * 'zVec6DLinShiftDRC()' directly shifts the velocity
 * type of 6D vector 'vec' at the original point to the
 * equivalent 6D vector at the point 'pos'.
 *
 * 'zVec6DAngShift()' shifts the force type of 6D vector
 * 'src' which works at 'pos' to the equivalent 6D vector
 * 'dest' which works at the original point, namely:
 *  f1 = f0
 *  n1 = n0 + 'pos' x f0
 * where 'src'=[f0 n0]^T and 'dest'=[f1 n1]^T.
 *
 * 'zVec6DAngShiftDRC()' directly shifts the force type of
 * 6D vector 'vec' which works at 'pos' to the equivalent 6D
 * vector which works at the original point.
 * [RETURN VALUE]
 * 'zVec6DLinShift()' and 'zVec6DAngShift()' return
 * pointers to 'dest'.
 * 'zVec6DLinShiftDRC()' and 'zVec6DAngShift()'
 * return pointers to 'vec'.
 */
__EXPORT zVec6D *zVec6DLinShift(zVec6D *src, zVec3D *pos, zVec6D *dest);
__EXPORT zVec6D *zVec6DLinShiftDRC(zVec6D *vec, zVec3D *pos);
__EXPORT zVec6D *zVec6DAngShift(zVec6D *src, zVec3D *pos, zVec6D *dest);
__EXPORT zVec6D *zVec6DAngShiftDRC(zVec6D *vec, zVec3D *pos);

/* ********************************************************** */
/* differential kinematics
 * ********************************************************** */

__EXPORT zVec6D *zVec6DDif(zVec6D *v, zVec6D *vnew, double dt, zVec6D *vel);

/* ********************************************************** */
/* I/O
 * ********************************************************** */

/* METHOD:
 * zVec6DFRead, zVec6DRead, zVec6DFWrite, zVec6DWrite,
 * zVec6DDataFWrite, zVec6DDataWrite
 * - input/output of 6D vector.
 *
 * 'zVec6DFRead()' reads 6 values from the current position
 * of the file 'fp', and creates a 6D vector 'v' from them.
 * 'zVec6DRead()' simply reads 6 values from the standard
 * input.
 *
 * 'zVec6DFWrite()' writes the 6D vector 'v' to the current
 * position of the file 'fp' in the following style.
 *  ( x, y, z, xa, ya, za )
 * When the NULL pointer is given, it writes the following string.
 *  (null 6D vector)
 * 'zVec6DWrite()' simply writes 'v' to the standard out.
 *
 * 'zVec6DDataFWrite()' writes the 6D vector 'v' to the current
 * position of the file 'fp' in the following style.
 *  x y z xa ya za
 * When the NULL pointer is given, it writes nothing.
 * 'zVec6DDataWrite()' simply writes 'v' to the standard out
 * in the same style with 'zVec6DDataFWrite()'.
 * [RETURN VALUE]
 * 'zVec6DFRead()' and 'zVec6DRead()' return a pointer to 'v'.
 *
 * 'zVec6DFWrite()', 'zVec6DWrite()', 'zVec6DDataFWrite()'
 * and 'zVec6DDataWrite()' return no value.
 */
__EXPORT zVec6D *zVec6DFRead(FILE *fp, zVec6D *v);
#define zVec6DRead(v) zVec6DFRead( stdin, (v) )
__EXPORT void zVec6DFWrite(FILE *fp, zVec6D *v);
#define zVec6DWrite(v) zVec6DFWrite( stdout, (v) )
__EXPORT void zVec6DDataFWrite(FILE *fp, zVec6D *v);
#define zVec6DDataWrite(v) zVec6DDataFWrite( stdout, (v) )

/* METHOD:
 * zVec6DFWriteXML - xml output.
 * ... yet testing.
 */
__EXPORT void zVec6DFWriteXML(FILE *fp, zVec6D *v);

__END_DECLS

#endif /* __ZEO_VEC6D_H__ */
