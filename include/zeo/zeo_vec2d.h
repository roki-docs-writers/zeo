/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_vec2d - 2D vector.
 */

#ifndef __ZEO_VEC2D_H__
#define __ZEO_VEC2D_H__

#include <zeo/zeo_misc.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zVec2D
 * 2D vector class
 * ********************************************************** */

typedef double zVec2D[2];

/* OBJECT:
 * zvec2Dzero, zvec2Dx, zvec2Dy
 * - 2D zero vector and unit vectors along (x,y) axis.
 */
extern const zVec2D zvec2Dzero;
extern const zVec2D zvec2Dx;
extern const zVec2D zvec2Dy;
#define ZVEC2DZERO ( (zVec2D)zvec2Dzero )
#define ZVEC2DX    ( (zVec2D)zvec2Dx )
#define ZVEC2DY    ( (zVec2D)zvec2Dy )

/* METHOD:
 * zVec2DCreate, zVec2DCopy, zVec2DClear
 * - creation, copy and cleanup of a 2D vector.
 *
 * 'zVec2DCreate()' creates a 2D vector 'v' which consists
 * of 'x' and 'y'.
 *
 * 'zVec2DCopy()' copies 'src' to 'dest'.
 *
 * 'zVec2DClear()' clears 'v', or sets all factors for 0.
 * [RETURN VALUE]
 * 'zVec2DCreate()' and 'zVec2DClear()' return a pointer
 * to 'v'.
 *
 * 'zVec2DCopy()' returns no value.
 * [NOTES]
 * 'zVec2DClear()' is a macro for zVec2DCreate( v, 0, 0 ).
 */
__EXPORT double *zVec2DCreate(zVec2D v, double a1, double a2);
#define zVec2DCopy(src,dest) do{ (dest)[0] = (src)[0]; (dest)[1] = (src)[1]; } while(0)
#define zVec2DClear(v)       zVec2DCreate( (v), 0, 0 )

/*! \brief create a 2D vector by the set of value for a polar expression.
 *
 * zVec2DCreatePolar() creates a 2D vector \a v from a set
 * of values for a polar coordinate ( \a r, \a theta ),
 * where \a r is for the radius from the original point, and
 * \a theta is for the longitudinal angle.
 * \retval \a v
 */
__EXPORT double *zVec2DCreatePolar(zVec2D v, double r, double theta);

/* METHOD:
 * zVec2DEqual - check if the two 2D vectors are equal.
 *
 * 'zVec2DEqual()' checks if the two 2D vectors, 'v1' and 'v2',
 * are equal, and returns a boolean value as a result.
 * [RETURN VALUE]
 * 'zVec2DEqual()' returns the true value if 'v1' and 'v2'
 * are equal, or false otherwise.
 */
__EXPORT bool zVec2DEqual(zVec2D v1, zVec2D v2);

/* METHOD:
 * zVec2DIsTol, zVec2DIsTiny
 * - check if 2D vector is tiny.
 *
 * 'zVec2DIsTol()' checks if the absolute values of every
 * components of 2D vector 'v' are smaller than 'tol'.
 *
 * 'zVec2DIsTiny()' applies zTOL (defined in "zm_misc.h") to
 * the tolerance of 'zVec2DIsTol()'.
 * [RETURN VALUE]
 * 'zVec2DIsTol()' and 'zVec2DIsTiny()' return the
 * true value when the absolute values of every components
 * of 'v' are smaller than 'tol' and zTOL, respectively,
 * or the false value, otherwise.
 * [NOTES]
 * 'tol' must be positive.
 * [SEE ALSO]
 * zIsTol, zIsTiny
 */
__EXPORT bool zVec2DIsTol(zVec2D v, double tol);
#define zVec2DIsTiny(v) zVec2DIsTol( v, zTOL )

/* ********************************************************** */
/* arithmetics
 * ********************************************************** */

/* METHOD:
 * zVec2DAdd, zVec2DSub, zVec2DRev, zVec2DMul,
 * zVec2DDiv, zVec2DCat,
 * zVec2DAddDRC, zVec2DSubDRC, zVec2DRevDRC,
 * zVec2DMulDRC, zVec2DDivDRC, zVec2DCatDRC
 * - the four rules of the arithmetics for 2D vector.
 *
 * 'zVec2DAdd()' adds the two 2D vectors, 'v1' and 'v2'.
 * The result is put into 'v'.
 *
 * 'zVec2DSub()' subtracts the 2D vector 'v2' from
 * the other 'v1'. The result is put into 'v'.
 *
 * 'zVec2DRev()' reverses the 2D vector 'v'. The result
 * is put into 'rv'.
 *
 * 'zVec2DMul()' multiplies the 2D vector 'v' by value
 * 'k'. The result is put into 'mv'.
 *
 * 'zVec2DDiv()' divides the 2D vector 'v' by 'k'.
 * The result is put into 'dv'.
 *
 * 'zVec2DCat()' concatenates the 2D vector 'v2' to 'v1',
 * multiplied by a scalar value 'k'. The result is put into 'v'.
 *
 * 'zVec2DAddDRC()' directly adds 'v2' to 'v1'.
 *
 * 'zVec2DSubDRC()' directly subtracts 'v2' from 'v1'.
 *
 * 'zVec2DRevDRC()' directly reverses 'v'.
 *
 * 'zVec2DMulDRC()' directly multiplies 'v' by 'k'.
 *
 * 'zVec2DDivDRC()' directly divides 'v' by 'k'.
 *
 * 'zVec2DCat()' directly concatenates 'v2' multiplied 'v2'
 * by 'k' to 'v1'.
 * [RETURN VALUE]
 * Each function returns a pointer to the resultant vector.
 *
 * And, 'zVec2DDiv()' and 'zVec2DDivDRC()' return the
 * NULL pointer if the given scalar value is 0.
 */
__EXPORT double *zVec2DAdd(zVec2D v1, zVec2D v2, zVec2D v);
__EXPORT double *zVec2DSub(zVec2D v1, zVec2D v2, zVec2D v);
__EXPORT double *zVec2DRev(zVec2D v, zVec2D rv);
__EXPORT double *zVec2DMul(zVec2D v, double k, zVec2D mv);
__EXPORT double *zVec2DDiv(zVec2D v, double k, zVec2D dv);
__EXPORT double *zVec2DCat(zVec2D v1, double k, zVec2D v2, zVec2D v);

#define zVec2DAddDRC(v1,v2)   zVec2DAdd(v1,v2,v1)
#define zVec2DSubDRC(v1,v2)   zVec2DSub(v1,v2,v1)
#define zVec2DRevDRC(v)       zVec2DRev(v,v)
#define zVec2DMulDRC(v,k)     zVec2DMul(v,k,v)
#define zVec2DDivDRC(v,k)     zVec2DDiv(v,k,v)
#define zVec2DCatDRC(v1,k,v2) zVec2DCat(v1,k,v2,v1)

/* METHOD:
 * zVec2DNorm, zVec2DSqrNorm, zVec2DSqrDist, zVec2DDist
 * - norm of the vector.
 *
 * 'zVec2DNorm()' calculates a norm of the 2D vector 'v'.
 *
 * 'zVec2DSqrNorm()' calculates a squared norm of 'v'.
 *
 * 'zVec2DDist()' calculates a distance between the two points
 * indicated by 'v1' and 'v2', which are position vectors for
 * each point.
 *
 * 'zVec2DSqrDist()' calculates a squared distance between
 * 'v1' and 'v2'.
 * [RETURN VALUE]
 * 'zVec2DNorm()' returns a norm of 'v'.
 *
 * 'zVec2DSqrNorm()' returns a squared norm of 'v'.
 *
 * 'zVec2DDist()' returns a distance between 'v1' and 'v2'.
 *
 * 'zVec2DSqrDist()' returns a squared distance between 'v1'
 * and 'v2'.
 * [NOTES]
 * 'zVec2DNorm()' is a macro for sqrt(zVec2DSqrNorm(v)) and
 * 'zVec2DDist()' is a macro for sqrt(zVec2DSqrDist(v))
 * (see "zeo_vec2d.h").
 */
__EXPORT double zVec2DSqrNorm(zVec2D v);
#define zVec2DNorm(v) sqrt(zVec2DSqrNorm((v)))

__EXPORT double zVec2DSqrDist(zVec2D v1, zVec2D v2);
#define zVec2DDist(v1,v2) sqrt(zVec2DSqrDist((v1),(v2)))

/* METHOD:
 * zVec2DNormalize, zVec2DNormalizeDRC
 * - normalization of a 2D vector.
 *
 * 'zVec2DNormalize()' normalizes the 2D vector 'v'.
 * The result is put into 'nv'.
 *
 * 'zVec2DNormalizeDRC()' normalizes the vector 'v'
 * directly.
 *
 * As a result of nomalization, the norm of a vector
 * will be 1.
 * [RETURN VALUE]
 * Both functions return a pointer to the result vector.
 * If failing normalization, i.e. the norm of 'v' is 0,
 * the value returned is the NULL pointer.
 */
__EXPORT double *zVec2DNormalize(zVec2D v, zVec2D nv);
#define zVec2DNormalizeDRC(v) zVec2DNormalize(v,v)

/* METHOD:
 * zVec2DInnerProd, zVec2DOuterProd
 * - inner/outer products.
 * [SYNOPSIS]
 * double zVec2DInnerProd(zVec2D v1, zVec2D v2);
 * double zVec2DOuterProd(zVec2D v1, zVec2D v2);
 * [DESCRIPTION]
 * 'zVec2DInnerProd()' calculates the inner product
 * of the two 2D vectors, 'v1' and 'v2'.
 *
 * 'zVec2DOuterProd()' calculates the outer product
 * of the two 2D vectors 'v1' and 'v2'. Since 'v1'
 * and 'v2' are 2-dimensional vectors, the result
 * is a scalar value.
 * [RETURN VALUE]
 * 'zVec2DInnerProd()' and 'zVec2DOuterProd()' return
 * scalar values as results.
 */
__EXPORT double zVec2DInnerProd(zVec2D v1, zVec2D v2);
__EXPORT double zVec2DOuterProd(zVec2D v1, zVec2D v2);

/* ********************************************************** */
/* geometry
 * ********************************************************** */

/* METHOD:
 * zVec2DInterDiv, zVec2DMid - interior division.
 *
 * 'zVec2DInterDiv()' calculates the interior division vector of
 * 'v1' and 'v2' with a division ratio 'ratio'. The result is
 * put into 'v'.
 *
 * i.e. 'v' = (1-'ratio')*'v1' + 'ratio'*'v2' .
 *
 * 'zVec2DMid()' calculates the middle point vector
 * of 'v1' and 'v2'. The result is put into 'v'.
 * i.e. 'v' = ( 'v1' + 'v2' ) / 2 .
 * [RETURN VALUE]
 * Both functions return a pointer to 'v'.
 */
__EXPORT double *zVec2DInterDiv(zVec2D v1, zVec2D v2, double ratio, zVec2D v);
__EXPORT double *zVec2DMid(zVec2D v1, zVec2D v2, zVec2D v);

/* METHOD:
 * zVec2DAngle - angle between the two vectors.
 *
 * 'zVec2DAngle()' calculates the angle between two
 * vectors 'v1' and 'v2'.
 * [RETURN VALUE]
 * The value returned is the angle between 'v1' and 'v2'.
 */
__EXPORT double zVec2DAngle(zVec2D v1, zVec2D v2);

/* METHOD:
 * zVec2DProject, zVec2DRot
 * - projection and rotation of a 2D vector.
 *
 * 'zVec2DProject()' projects vector 'v' onto the line
 * directed by 'n', and put the result into 'pv'; 'pv'
 * is parallel to 'n' and the subtraction vector from
 * 'pv' to 'v' is orthogonal to 'n'.
 *
 * 'zVec2DRot()' rotates 'v' with 'angle'.
 * The result is set into 'rv'.
 * [RETURN VALUE]
 * 'zVec2DProject()' returns a pointer to 'pv', or the
 * null pointer if 'n' is the zero vector.
 *
 * 'zVec2DRot()' returns a pointer to 'rv'.
 */
__EXPORT double *zVec2DProject(zVec2D v, zVec2D n, zVec2D pv);
__EXPORT double *zVec2DRot(zVec2D v, double angle, zVec2D rv);

/* ********************************************************** */
/* I/O
 * ********************************************************** */

/* METHOD:
 * zVec2DFRead, zVec2DRead, zVec2DFWrite, zVec2DWrite,
 * zVec2DDataFWrite, zVec2DDataWrite
 * - input/output of 2D vector.
 *
 * 'zVec2DFRead()' reads two values from the current
 * position of the file 'fp', and creates a 2D vector
 * 'v' from them. 'zVec2DRead()' simply reads two values
 * from the standard input.
 *
 * 'zVec2DFWrite()' writes the 2D vector 'v' to the current
 * position of the file 'fp' in the following style.
 *  ( x, y )
 * When the NULL pointer is given, it writes the following string.
 *  (null 2D vector)
 * 'zVec2DWrite()' simply writes 'v' to the standard out.
 *
 * 'zVec2DDataFWrite()' writes the 2D vector 'v' to the current
 * position of the file 'fp' in the following style.
 *  x y
 * When the NULL pointer is given, it writes nothing.
 * 'zVec2DDataWrite()' simply writes 'v' to the standard out
 * in the same style with 'zVec2DDataFWrite()'.
 * [RETURN VALUE]
 * 'zVec2DFRead()' and 'zVec2DRead()' return a pointer to 'v'.
 *
 * 'zVec2DFWrite()', 'zVec2DWrite()', 'zVec2DDataFWrite()'
 * and 'zVec2DDataWrite()' return no value.
 */
__EXPORT double *zVec2DFRead(FILE *fp, zVec2D v);
#define zVec2DRead(v) zVec2DFRead( stdin, (v) )
__EXPORT void zVec2DFWrite(FILE *fp, zVec2D v);
#define zVec2DWrite(v) zVec2DFWrite( stdout, (v) )
__EXPORT void zVec2DDataFWrite(FILE *fp, zVec2D v);
#define zVec2DDataWrite(v) zVec2DDataFWrite( stdout, (v) )

__END_DECLS

#endif /* __ZEO_VEC2D_H__ */
