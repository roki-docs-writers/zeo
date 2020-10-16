/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_ep - Euler parameter (unit quaternion) class.
 */

#ifndef __ZEO_EP_H__
#define __ZEO_EP_H__

#include <zeo/zeo_mat3d.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zEP
 * Euler parameter class
 * ********************************************************** */

typedef union{
  struct{
    double w;
    zVec3D v;
  } ex;
  double e[4];
} zEP;

/* METHOD:
 * zEPCreateSC, zEPCreate, zEPIdent, zEPCopy
 * - create and copy of Euler parameter.
 *
 * 'zEPCreate()' creates Euler parameter 'ep'.
 * 'ep' equals to [cos(0.5*'theta'), sin(0.5*'theta')x'axis'].
 *
 * 'zEPIdent()' creates Euler parameter equivalent to the
 * identity transformation, namely, 'ep' equals to
 * [ 1, 0, 0, 0 ].
 *
 * 'zEPCopy()' copies Euler parameter 'src' to the other
 * 'dest'.
 * [RETURN VALUE]
 * 'zEPCreate()' and 'zEPIdent()' return a pointer 'ep'.
 * 'zEPCopy()' returns no value.
 */
__EXPORT zEP *zEPCreate(zEP *ep, double theta, zVec3D *axis);
#define zEPIdent(e) zEPCreate( e, 0, ZVEC3DZERO )
#define zEPCopy(s,d) ( *(d) = *(s) )

__EXPORT bool zEPIsIdent(zEP *ep);

/* METHOD:
 * zEP2AA, zAA2EP, zMat3DEP, zMat3DToEP
 * - alternate angle-axis vector, attitude matrix and Euler parameter.
 *
 * 'zEP2AA()' converts Euler parameter 'ep' to the
 * equivalent angle-axis vector 'aa'.
 *
 * 'zAA2EP()' converts 'aa' to the equivalent 'ep'.
 *
 * 'zMat3DEP()' creates an attitude matrix 'm' from
 * Euler parameter 'ep'.
 *
 * 'zMat3DToEP()' converts 'm' to the equivalent 'ep'.
 * [RETURN VALUE]
 * 'zAA2EP()' returns a pointer 'ep'.
 * 'zEP2AA()' returns a pointer 'aa'.
 * 'zMat3DEP()' returns a pointer 'm'.
 * 'zMat3DToEP()' returns a pointer 'ep'.
 */
__EXPORT zVec3D *zEP2AA(zEP *ep, zVec3D *aa);
__EXPORT zEP *zAA2EP(zVec3D *aa, zEP *ep);
__EXPORT zMat3D *zMat3DEP(zMat3D *m, zEP *ep);
__EXPORT zEP *zMat3DToEP(zMat3D *m, zEP *ep);

/* METHOD
 * zEPRotVec - rotate vector by Euler parameter.
 *
 * 'zEPRotVec()' rotates a vector 'v' by Euler parameter
 * 'ep'. The result is put into 'rv'.
 * [RETURN VALUE]
 * 'zEPRotVec()' returns a pointer 'rv'.
 * [SEE ALSO]
 * zMat3DEP
 */
__EXPORT zVec3D *zEPRotVec(zEP *ep, zVec3D *v, zVec3D *rv);

/* METHOD
 * zEPVel2AngVel, zAngVel2EPVel
 * - convert from/to rotation velocity to/from Euler parameter derivative.
 *
 * 'zEPVel2AngVel()' converts the derivative of Euler
 * parameter 'epvel' at the attitude represented by
 * 'ep' to the equivalent rotation velocity 'angvel'.
 *
 * 'zAngVel2EPVel()' converts 'angvel' at the attitude
 * 'ep' to the equivalent 'epvel', where 'epvel' is
 * perpendicular to 'ep'.
 * [RETURN VALUE]
 * 'zEPVel2AngVel()' returns a pointer 'angvel'.
 * 'zAngVel2EPVel()' returns a pointer 'epvel'.
 */
__EXPORT zVec3D *zEPVel2AngVel(zEP *epvel, zEP *ep, zVec3D *angvel);
__EXPORT zEP *zAngVel2EPVel(zVec3D *angvel, zEP *ep, zEP *epvel);

/* METHOD:
 * zEPRev, zEPMul, zEPCat, zEPRevDRC, zEPMulDRC, zEPCatDRC
 * zEPInnerProd, zEPNorm, zEPNormalize
 * - Euler parameter arithmetics.
 *
 * 'zEPRev()' reverses Euler parameter 'ep1'. The result
 * is put into 'ep'. 'zEPRevDRC()' directly reverses
 * 'ep'.
 *
 * 'zEPMul()' multiplies 'ep1' by a scalar value 'k'.
 * The result is put into 'ep'. 'zEPMulDRC()' directly
 * multiplies 'ep' by 'k'.
 *
 * 'zEPCat()' concatenates Euler parameter 'ep2'
 * multiplied by 'k' to 'ep1'. The results is put into
 * 'ep'. 'zEPCatDRC()' directly concatenates 'ep2'
 * multiplied multiplied by 'k' to 'ep1'.
 *
 * 'zEPInnerProd()' calculates the inner products of
 * 'ep1' and 'ep2'. 'zEPNorm()' calculates the norm
 * of 'ep', equivalent to zEPInnerProd( 'ep', 'ep' ).
 *
 * 'zEPNormalize()' directly normalizes 'ep'. Since
 * Euler parameter is a class of unit quaternions,
 * it is automatically called in 'zEPCreate()' and
 * some other functions.
 * [RETURN VALUE]
 * 'zEPRev()', 'zEPMul()' and 'zEPCat()' return a
 * pointer to the resultant Euler parameter 'ep'.
 *
 * 'zEPRevDRC()', 'zEPMulDRC()', 'zEPCatDRC()' and
 * 'zEPNormalize()' return a pointer to Euler
 * parameter directly modified.
 *
 * 'zEPInnerProd()' and 'zEPNorm()' return a value
 * calculated.
 */
__EXPORT zEP *zEPSub(zEP *ep1, zEP *ep2, zEP *ep);
__EXPORT zEP *zEPRev(zEP *ep1, zEP *ep);
__EXPORT zEP *zEPMul(zEP *ep1, double k, zEP *ep);
__EXPORT zEP *zEPCat(zEP *ep1, double k, zEP *ep2, zEP *ep);
#define zEPSubDRC(e1,e)         zEPSub( e1, e, e1 )
#define zEPRevDRC(e)            zEPRev( e, e )
#define zEPMulDRC(e,k)          zEPMul( e, k, e )
#define zEPCatDRC(ep1,k,ep2)    zEPCat( ep1, k, ep2, ep1 )

__EXPORT zEP *zEPDif(zEP *ep1, zEP *ep2, double dt, zEP *ep_vel);

__EXPORT double zEPInnerProd(zEP *ep1, zEP *ep2);
__EXPORT double zEPNorm(zEP *ep);
__EXPORT zEP *zEPNormalize(zEP *ep);

/*! \brief cascade a Euler parameter to another.
 */
__EXPORT zEP *zEPCascade(zEP *e1, zEP *e2, zEP *e);

/* METHOD:
 * zEPInterDiv - interior division of Euler parameter.
 *
 * 'zEPInterDiv()' calculates the interior division of
 * two Euler parameters 'ep1' and 'ep2' in accordance
 * with the spherical interpolation.
 * 't' is the dividing ratio, which is regularly chosen
 * within the range from 0 to 1. When 't' is out of
 * the range, 'zEPInterDiv()' computes the outer
 * division for extrapolation.
 * The result is put into 'ep'.
 *
 * The computation is based on SLERP(spherical linear
 * interpolation proposed by K. Shoemake, 1985).
 * [RETURN VALUE]
 * 'zEPInterDiv()' returns a pointer 'ep'.
 */
__EXPORT zEP *zEPInterDiv(zEP *ep1, zEP *ep2, double t, zEP *ep);

/*! \brief interior division of two attitude matrices for SLERP.
 *
 * zMat3DInterDiv() calculates the interior division of
 * two attitude matrices \a m1 and \a m2 in accordance
 * with the spherical interpolation.
 * \a t is the dividing ratio, which is regularly chosen
 * within the range from 0 to 1. When \a t is out of
 * the range, zMat3DInterDiv() computes the outer
 * division for extrapolation.
 * The result is put into \a m.
 * The computation enables SLERP(spherical linear interpolation).
 *
 * \retval a pointer \a m.
 */
__EXPORT zMat3D *zMat3DInterDiv(zMat3D *m1, zMat3D *m2, double t, zMat3D *m);

/* METHOD:
 * zEPFWrite, zEPWrite
 * - output of Euler parameter.
 *
 * 'zEPFWrite()' writes Euler parameter to the current
 * position of the file 'fp' in the following form.
 *  'e0' { 'e1', 'e2', 'e3' }
 * 'zEPWrite()' writes 'ep' to the standard output.
 * [RETURN VALUE]
 * Neither 'zEPFWrite()' nor 'zEPWrite()' return
 * any values.
 */
__EXPORT void zEPFWrite(FILE *fp, zEP *ep);
#define zEPWrite(e) zEPFWrite( stdout, e )

__END_DECLS

#endif /* __ZEO_EP_H__ */
