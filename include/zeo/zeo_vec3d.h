/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_vec3d - 3D vector.
 */

#ifndef __ZEO_VEC3D_H__
#define __ZEO_VEC3D_H__

#include <zeo/zeo_misc.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zVec3D
 * 3D vector class
 * ********************************************************** */

typedef struct{
  double e[3];
} zVec3D;

#define zVec3DArray(v)       (v)->e
#define zVec3DElem(v,i)      ( zVec3DArray(v)[i] )
#define zVec3DSetElem(v,i,x) ( zVec3DElem(v,i) = (x) )

/* OBJECT:
 * zvec3Dzero, zvec3Dx, zvec3Dy, zvec3Dz
 * - 3D zero vector and unit vectors along (x,y,z) axis.
 */
extern const zVec3D zvec3Dzero;
extern const zVec3D zvec3Dx;
extern const zVec3D zvec3Dy;
extern const zVec3D zvec3Dz;
#define ZVEC3DZERO ( (zVec3D *)&zvec3Dzero )
#define ZVEC3DX    ( (zVec3D *)&zvec3Dx )
#define ZVEC3DY    ( (zVec3D *)&zvec3Dy )
#define ZVEC3DZ    ( (zVec3D *)&zvec3Dz )

/*! \brief create, copy and cleanup a 3D vector.
 *
 * zVec3DCreate() creates a 3D vector \a v which consists
 * of \a x, \a y and \a z.
 *
 * zVec3DCopy() copies \a src to \a dest.
 *
 * zVec3DClear() sets all components of \a v for 0.
 * \return
 * zVec3DCreate() and zVec3DClear() return a pointer \a v.
 *
 * zVec3DCopy() returns no value.
 * \notes
 * It is also possible to write simply *dest = *src instead of
 * zVec3DCopy(). Actually, zVec3DCopy() is defined as macro
 * (see "zeo_vec_vec3d.h").
 *
 * zVec3DClear() is a macro for zVec3DCreate( v, 0, 0, 0 ).
 */
__EXPORT zVec3D *zVec3DCreate(zVec3D *v, double x, double y, double z);
#define zVec3DCopy(src,dest) ( *(dest) = *(src) )
#define zVec3DClear(v)       zVec3DCreate( (v), 0, 0, 0 )

/*! \brief creation of a 3D vector by the set of value for a polar expression.
 *
 * zVec3DCreatePolar() creates a vector \a v from a set of values
 * for a polar coordinate ( \a r, \a theta, \a phi ), where
 * \a r is for the radius from the original point, \a theta is for
 * the longitudinal angle, and \a phi is for the latitudinal angle.
 * \retval \a v
 */
__EXPORT zVec3D *zVec3DCreatePolar(zVec3D *v, double r, double theta, double phi);

/* METHOD:
 * zVec3DMatch, zVec3DEqual
 * - check if the two 3D vectors are equal.
 *
 * 'zVec3DMatch()' and 'zVec3DEqual()' check if the two
 * 3D vectors, 'v1' and 'v2', are equal. They return
 * a boolean value as a result.
 *
 * 'zVec3DMatch()' strictly compares the two vectors,
 * while 'zVec3DEqual()' checks if the error between
 * 'v1' and 'v2' are sufficiently small.
 * [RETURN VALUE]
 * 'zVec3DMatch()' and 'zVec3DEqual()' return the true
 * value if 'v1' and 'v2' are equal, or false otherwise.
 * [NOTES]
 * It is also possible to write just *v1 == *v2, instead
 * of calling 'zVec3DMatch( v1, v2 )'.
 */
__EXPORT bool zVec3DMatch(zVec3D *v1, zVec3D *v2);
__EXPORT bool zVec3DEqual(zVec3D *v1, zVec3D *v2);

/* METHOD:
 * zVec3DIsTol, zVec3DIsTiny
 * - check if 3D vector is tiny.
 *
 * 'zVec3DIsTol()' checks if the absolute values of every
 * components of 3D vector 'v' are smaller than 'tol'.
 *
 * 'zVec3DIsTiny()' applies zTOL (defined in "zm_misc.h") to
 * the tolerance of 'zVec3DIsTol()'.
 * [RETURN VALUE]
 * 'zVec3DIsTol()' and 'zVec3DIsTiny()' return the
 * true value when the absolute values of every components
 * of 'v' are smaller than 'tol' and zTOL, respectively, or
 * the false value, otherwise.
 * [NOTES]
 * 'tol' must be positive.
 * [SEE ALSO]
 * zIsTol, zIsTiny
 */
__EXPORT bool zVec3DIsTol(zVec3D *v, double tol);
#define zVec3DIsTiny(v) zVec3DIsTol( v, zTOL )

/*! \brief check if a 3D vector includes NaN or Inf components. */
__EXPORT bool zVec3DIsNan(zVec3D *v);

/* ********************************************************** */
/* arithmetics
 * ********************************************************** */

/*! \brief the four rules of the arithmetics for 3D vector.
 *
 * 'zVec3DAdd()' adds the two 3D vectors, 'v1' and 'v2'.
 * The result is put into 'v'.
 *
 * 'zVec3DSub()' subtracts the 3D vector 'v2' from
 * the other 'v1'. The result is put into 'v'.
 *
 * 'zVec3DRev()' reverses the 3D vector 'v'. The result
 * is put into 'rv'.
 *
 * 'zVec3DMul()' multiplies the 3D vector 'v' by value
 * 'k'. The result is put into 'mv'.
 *
 * 'zVec3DDiv()' divides the 3D vector 'v' by 'k'.
 * The result is put into 'dv'.
 *
 * 'zVec3DCat()' concatenates the 3D vector 'v2' to 'v1',
 * multiplied by a scalar value 'k'. The result is put into 'v'.
 *
 * 'zVec3DAddDRC()' directly adds 'v2' to 'v1'.
 *
 * 'zVec3DSubDRC()' directly subtracts 'v2' from 'v1'.
 *
 * 'zVec3DRevDRC()' directly reverses 'v'.
 *
 * 'zVec3DMulDRC()' directly multiplies 'v' by 'k'.
 *
 * 'zVec3DDivDRC()' directly divides 'v' by 'k'.
 *
 * 'zVec3DCat()' directly concatenates 'v2' multiplied 'v2'
 * by 'k' to 'v1'.
 * [RETURN VALUE]
 * Each function returns a pointer to the resultant vector.
 *
 * And, 'zVec3DDiv()' and 'zVec3DDivDRC()' return the
 * NULL pointer if the given scalar value is 0.
 */
__EXPORT zVec3D *zVec3DAdd(zVec3D *v1, zVec3D *v2, zVec3D *v);
__EXPORT zVec3D *zVec3DSub(zVec3D *v1, zVec3D *v2, zVec3D *v);
__EXPORT zVec3D *zVec3DRev(zVec3D *v, zVec3D *rv);
__EXPORT zVec3D *zVec3DMul(zVec3D *v, double k, zVec3D *mv);
__EXPORT zVec3D *zVec3DDiv(zVec3D *v, double k, zVec3D *dv);
__EXPORT zVec3D *zVec3DAmp(zVec3D *v, zVec3D *a, zVec3D *av);
__EXPORT zVec3D *zVec3DCat(zVec3D *v1, double k, zVec3D *v2, zVec3D *v);

#define zVec3DAddDRC(v1,v2)   zVec3DAdd(v1,v2,v1)
#define zVec3DSubDRC(v1,v2)   zVec3DSub(v1,v2,v1)
#define zVec3DRevDRC(v)       zVec3DRev(v,v)
#define zVec3DMulDRC(v,k)     zVec3DMul(v,k,v)
#define zVec3DDivDRC(v,k)     zVec3DDiv(v,k,v)
#define zVec3DAmpDRC(v,a)     zVec3DAmp(v,a,v)
#define zVec3DCatDRC(v1,k,v2) zVec3DCat(v1,k,v2,v1)

/*! \brief norm of the vector.
 *
 * 'zVec3DNorm()' calculates a norm of the 3D vector 'v'.
 *
 * 'zVec3DSqrNorm()' calculates a squared norm of 'v'.
 *
 * 'zVec3DDist()' calculates a distance between the two points
 * indicated by 'v1' and 'v2', which are position vectors for
 * each point.
 *
 * 'zVec3DSqrDist()' calculates a squared distance between
 * 'v1' and 'v2'.
 * [RETURN VALUE]
 * 'zVec3DNorm()' returns a norm of 'v'.
 *
 * 'zVec3DSqrNorm()' returns a squared norm of 'v'.
 *
 * 'zVec3DDist()' returns a distance between 'v1' and 'v2'.
 *
 * 'zVec3DSqrDist()' returns a squared distance between 'v1'
 * and 'v2'.
 * [NOTES]
 * 'zVec3DNorm()' is a macro for sqrt(zVec3DSqrNorm(v)) and
 * 'zVec3DDist()' is a macro for sqrt(zVec3DSqrDist(v))
 * (see "zeo_vec_vec3d.h").
 */
__EXPORT double zVec3DSqrNorm(zVec3D *v);
#define zVec3DNorm(v) sqrt(zVec3DSqrNorm((v)))

__EXPORT double zVec3DWSqrNorm(zVec3D *v, zVec3D *w);
#define zVec3DWNorm(v,w) sqrt( zVec3DWSqrNorm(v,w) )

__EXPORT double zVec3DSqrDist(zVec3D *v1, zVec3D *v2);
#define zVec3DDist(v1,v2) sqrt(zVec3DSqrDist((v1),(v2)))

/* METHOD:
 * zVec3DNormalize, zVec3DNormalizeDRC
 * - normalization of a 3D vector.
 *
 * 'zVec3DNormalize()' normalizes the 3D vector 'v'.
 * The result is put into 'nv'.
 *
 * 'zVec3DNormalizeDRC()' normalizes the vector 'v' directly.
 *
 * As a result of nomalization, the norm of a vector will be 1.
 * [RETURN VALUE]
 * Both functions return a pointer to the result vector.
 * If failing normalization, i.e. the norm of 'v' is 0, the
 * value returned is the NULL pointer.
 */
__EXPORT double zVec3DNormalizeNC(zVec3D *v, zVec3D *nv);
__EXPORT double zVec3DNormalize(zVec3D *v, zVec3D *nv);
#define zVec3DNormalizeNCDRC(v) zVec3DNormalizeNC(v,v)
#define zVec3DNormalizeDRC(v)   zVec3DNormalize(v,v)

/* METHOD:
 * zVec3DInnerProd, zVec3DOuterProd, zVec3DOuterProdNorm
 * zVec3DGrassmannProd, zVec3DTripleProd
 * - inner/outer products.
 *
 * 'zVec3DInnerProd()' calculates the inner product of the
 * two 3D vectors, 'v1' and 'v2'.
 *
 * 'zVec3DOuterProd()' calculates the outer product of the
 * two 3D vectors 'v1' and 'v2'. The outer product is a 3D vector,
 * and put it into 'v', i.e. 'v' = 'v1' x 'v2'.
 *
 * 'zVec3DOuterProdNorm()' calculates only the norm of the
 * outer product of 'v1' and 'v2'.
 *
 * 'zVec3DGrassmannProd()' calculates the scalar triple product
 * of 'v1', 'v2' and 'v3', i.e. 'v1'.('v2' x 'v3'), which is
 * also expressed as [ 'v1' 'v2' 'v3' ].
 *
 * 'zVec3DTripleProd()' calculates the vector triple product
 * of 'v1', 'v2' and 'v3', and put it into 'v', i.e.
 * 'v' = 'v1' x ('v2' x 'v3').
 * [RETURN VALUE]
 * 'zVec3DInnerProd()', 'zVec3DOuterProdNorm()' and
 * 'zVec3DGrassmannProd()' return scalar values as results.
 *
 * 'zVec3DOuterProd()' and 'zVec3DTripleProd()' return
 * a pointer to 'v'.
 * [NOTES]
 * For 'zVec3DOuterProd()' and 'zVec3DTripleProd()',
 * it is allowed to let 'v' point to the same address with
 * 'v1' or 'v2'.
 */
__EXPORT double zVec3DInnerProd(zVec3D *v1, zVec3D *v2);
__EXPORT zVec3D *zVec3DOuterProd(zVec3D *v1, zVec3D *v2, zVec3D *v);
__EXPORT double zVec3DOuterProdNorm(zVec3D *v1, zVec3D *v2);
__EXPORT double zVec3DGrassmannProd(zVec3D *v1, zVec3D *v2, zVec3D *v3);
__EXPORT zVec3D *zVec3DTripleProd(zVec3D *v1, zVec3D *v2, zVec3D *v3, zVec3D *v);

/* ********************************************************** */
/* geometry
 * ********************************************************** */

/*! \brief interior division.
 *
 * zVec3DInterDiv() calculates the interior division vector
 * of \a v1 and \a v2 with a division ratio \a ratio.
 * The result is put into \a v.
 *
 * i.e. \a v = (1-\a ratio)* \a v1 + \a ratio * \a v2.
 *
 * zVec3DMid() calculates the middle point vector of \a v1
 * and \a v2. The result is put into \a v.
 * i.e. \a v = ( \a v1 + \a v2 ) / 2 .
 * \return
 * Both functions return a pointer to \a v.
 */
__EXPORT zVec3D *zVec3DInterDiv(zVec3D *v1, zVec3D *v2, double ratio, zVec3D *v);
__EXPORT zVec3D *zVec3DMid(zVec3D *v1, zVec3D *v2, zVec3D *v);

/*! \brief angle between the two vectors.
 *
 * zVec3DAngle() calculates the angle between the two vectors
 * \a v1 and \a v2. When \a n is not the null vector, signed
 * angle (i.e. the angle from \a v1 to \a v2 about the axis
 * along \a n) is computed.
 * \return
 * The value returned is the angle between \a v1 and \a v2.
 */
__EXPORT double zVec3DAngle(zVec3D *v1, zVec3D *v2, zVec3D *n);

/* METHOD:
 * zVec3DProj, zVec3DOrthogonalize, zVec3DOrthoSpace,
 * zVec3DRot
 * - projection, orthogonalization and rotation of a 3D vector.
 *
 * 'zVec3DProj()' projects vector 'v' onto the line directed by
 * 'n', and put the result into 'pv'; 'pv' is parallel to 'n' and
 * the subtraction vector from 'pv' to 'v' is orthogonal to 'n'.
 *
 * 'zVec3DOrthogonalize()' orthogonalizes 'v' against 'n', and put
 * it into 'ov'; 'ov' is orthogonal to 'n'.
 *
 * 'zVec3DOrthoSpace()' creates the orthogonal space to 'v',
 * and put them into 'sv1' and 'sv2'; 'v', 'sv1' and 'sv2' are
 * orthogonal with each other, and are normalized.
 * Note that the orthogonal is not unique in nature. This
 * function only creates "one of" them.
 *
 * 'zVec3DRot()' rotates 'v' by angle-axis vector 'aa',
 * whose direction is that of the rotation axis and norm is
 * the rotation angle.
 * The result is set into 'rv'.
 * [RETURN VALUE]
 * 'zVec3DProj()' and 'zVec3DOrthogonalize()' return a
 * pointer to 'pv' and 'ov', respectively, or the null pointer
 * if they fails because of 'n' being the zero vector.
 *
 * 'zVec3DOrthoSpace()' returns the true value if it succeeds to
 * create the orthogonal space, or the false value, otherwise.
 *
 * 'zVec3DRot()' returns a pointer to 'rv'.
 */
__EXPORT zVec3D *zVec3DProj(zVec3D *v, zVec3D *n, zVec3D *pv);
__EXPORT zVec3D *zVec3DOrthogonalize(zVec3D *v, zVec3D *n, zVec3D *ov);
__EXPORT bool zVec3DOrthoSpace(zVec3D *v, zVec3D *sv1, zVec3D *sv2);
__EXPORT zVec3D *zVec3DRot(zVec3D *v, zVec3D *aa, zVec3D *rv);

/* ********************************************************** */
/* differential kinematics
 * ********************************************************** */

/*! \brief compute average velocity of a 3D vector.
 *
 * zVec3DDif() computes the avarage velocity when a 3D vector
 * \a v changes to \a vnew in the duration \a dt.
 * \a vel is a velocity vector computed.
 * \retval \a vel
 */
__EXPORT zVec3D *zVec3DDif(zVec3D *v, zVec3D *vnew, double dt, zVec3D *vel);

/* METHOD:
 * zVec3DZYXVel2AngVel, zVec3DZYXVel2AngVelSC,
 * zVec3DAngVel2ZYXVel, zVec3DAngVel2ZYXVelSC,
 * zVec3DZYZVel2AngVel, zVec3DZYZVel2AngVelSC,
 * zVec3DAngVel2ZYZVel, zVec3DAngVel2ZYZVelSC
 * - convert from/to Eulerian angle differential to/from angular velocity.
 *
 * 'zVec3DZYXVel2AngVel()' converts a set of differential
 * values of z-y-x Eulerian angle 'zyxvel' at the attitude
 * represented by 'zyx' to an equivalent angular velocity
 * vector. The result is put into 'angvel'.
 *
 * 'zVec3DZYXVel2AngVelSC()' directly accepts sets of
 * sine/cosine values for the z-y-x Eulerian angle.
 * The set of 'sa'/'ca' is for the first angle, while
 * that of 'sb'/'cb' for the second.
 *
 * 'zVec3DAngVel2ZYXVel()' and 'zVec3DAngVel2ZYXVelSC()'
 * are their inverse conversions.
 *
 * 'zVec3DZYZVel2AngVel()' converts a set of differential
 * values of z-y-z Eulerian angle 'zyxvel' at the attitude
 * represented by 'zyz' to an equivalent angular velocity
 * vector. The result is put into 'angvel'.
 *
 * 'zVec3DZYZVel2AngVelSC()' directly accepts sets of
 * sine/cosine values for the z-y-z Eulerian angle.
 * The set of 'sa'/'ca' is for the first angle, while
 * that of 'sb'/'cb' for the second.
 *
 * 'zVec3DAngVel2ZYZVel()' and 'zVec3DAngVel2ZYZVelSC()'
 * are their inverse conversions.
 *
 * Note that conversion from angular velocity to the
 * derivatives of Eulerian angles has singular points
 * due to the mathematical representation.
 * In the case of z-y-x Eulerian angle, points where
 * cosine of the second value is zero are singular.
 * In the case of z-y-z Eulerian angle, points where
 * sine of the second value is zero are singular.
 * At such singular points, 'zVec3DAngVel2ZYXVel()',
 * 'zVec3DAngVel2ZYXVelSC()', 'zVec3DAngVel2ZYZVel()'
 * and 'zVec3DAngVel2ZYZVelSC()' do nothing.
 * [RETURN VALUE]
 * 'zVec3DZYXVel2AngVel()', 'zVec3DZYXVel2AngVelSC()',
 * 'zVec3DZYZVel2AngVel()' and 'zVec3DZYZVel2AngVelSC()'
 * return a pointer 'angvel'.
 *
 * 'zVec3DAngVel2ZYXVel()' and 'zVec3DAngVel2ZYXVelSC()',
 * return a pointer 'zyxvel', while
 * 'zVec3DAngVel2ZYZVel()' and 'zVec3DAngVel2ZYZVelSC()'
 * return a pointer 'zyzvel'.
 * [SEE ALSO]
 * zMat3DZYX, zMat3DToZYX, zMat3DZYZ, zMat3DToZYZ
 */
__EXPORT zVec3D *zVec3DZYXVel2AngVel(zVec3D *zyxvel, zVec3D *zyx, zVec3D *angvel);
__EXPORT zVec3D *zVec3DZYXVel2AngVelSC(zVec3D *zyxvel, double sa, double ca, double sb, double cb, zVec3D *angvel);
__EXPORT zVec3D *zVec3DAngVel2ZYXVel(zVec3D *angvel, zVec3D *zyx, zVec3D *zyxvel);
__EXPORT zVec3D *zVec3DAngVel2ZYXVelSC(zVec3D *angvel, double sa, double ca, double sb, double cb, zVec3D *zyxvel);
__EXPORT zVec3D *zVec3DZYZVel2AngVel(zVec3D *zyzvel, zVec3D *zyz, zVec3D *angvel);
__EXPORT zVec3D *zVec3DZYZVel2AngVelSC(zVec3D *zyzvel, double sa, double ca, double sb, double cb, zVec3D *angvel);
__EXPORT zVec3D *zVec3DAngVel2ZYZVel(zVec3D *angvel, zVec3D *zyz, zVec3D *zyzvel);
__EXPORT zVec3D *zVec3DAngVel2ZYZVelSC(zVec3D *angvel, double sa, double ca, double sb, double cb, zVec3D *zyzvel);

/* ********************************************************** */
/* I/O
 * ********************************************************** */

/* METHOD:
 * zVec3DFRead, zVec3DRead, zVec3DFWrite, zVec3DWrite,
 * zVec3DDataFWrite, zVec3DDataWrite
 * - input/output of 3D vector.
 *
 * 'zVec3DFRead()' reads three values from the current
 * position of the file 'fp', and creates a 3D vector
 * 'v' from them. 'zVec3DRead()' simply reads three
 * values from the standard input.
 *
 * 'zVec3DFWrite()' writes the 3D vector 'v' to the
 * current position of the file 'fp' in the following style.
 *  ( x, y, z )
 * When the NULL pointer is given, it writes the following string.
 *  (null 3D vector)
 * 'zVec3DWrite()' simply writes 'v' to the standard out.
 *
 * 'zVec3DDataFWrite()' writes the 3D vector 'v' to the current
 * position of the file 'fp' in the following style.
 *  x y z
 * When the NULL pointer is given, it writes nothing.
 * 'zVec3DDataWrite()' simply writes 'v' to the standard out
 * in the same style with 'zVec3DDataFWrite()'.
 * [RETURN VALUE]
 * 'zVec3DFRead()' and 'zVec3DRead()' return a pointer to 'v'.
 *
 * 'zVec3DFWrite()', 'zVec3DWrite()', 'zVec3DDataFWrite()'
 * and 'zVec3DDataWrite()' return no value.
 */
__EXPORT zVec3D *zVec3DFRead(FILE *fp, zVec3D *v);
#define zVec3DRead(v) zVec3DFRead( stdin, (v) )
__EXPORT void zVec3DFWrite(FILE *fp, zVec3D *v);
#define zVec3DWrite(v) zVec3DFWrite( stdout, (v) )
__EXPORT void zVec3DDataFWrite(FILE *fp, zVec3D *v);
#define zVec3DDataWrite(v) zVec3DDataFWrite( stdout, (v) )

/* METHOD:
 * zVec3DFWriteXML - xml output.
 * ... yet testing.
 */
__EXPORT void zVec3DFWriteXML(FILE *fp, zVec3D *v);

__END_DECLS

#include <zeo/zeo_vec3d_list.h>  /* 3D vector list */
#include <zeo/zeo_vec3d_tree.h>  /* 3D vector tree */

#endif /* __ZEO_VEC3D_H__ */
