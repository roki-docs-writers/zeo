/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_mat3d - 3x3 matrix.
 */

#ifndef __ZEO_MAT3D_H__
#define __ZEO_MAT3D_H__

#include <zeo/zeo_vec6d.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zMat3D
 * 3D matrix class
 * ********************************************************** */

typedef union{
  double e[3][3]; /*!< 3x3 matrix */
  zVec3D v[3];    /*!< 3 column vectors */
  double c[9];    /*!< 9 components */
} zMat3D;

#define zMat3DArray(m)         (m)->c
#define zMat3DElem(m,r,c)      (m)->e[(c)][(r)]
#define zMat3DVec(m,i)         ( &(m)->v[(i)] )
#define zMat3DElem9(m,i)       ( zMat3DArray(m)[i] )
#define zMat3DSetElem(m,r,c,x) ( zMat3DElem(m,r,c) = (x) )
#define zMat3DSetVec(m,i,v)    zVec3DCopy(v,zMat3DVec(m,i))
#define zMat3DSetElem9(m,i,c)  ( zMat3DElem9(m,i) = (c) )

/* OBJECT:
 * zmat3Dzero, zmat3Dident
 * - 3D zero matrix and identity matrix.
 */
extern const zMat3D zmat3Dzero;
extern const zMat3D zmat3Dident;
#define ZMAT3DZERO  ( (zMat3D *)&zmat3Dzero )
#define ZMAT3DIDENT ( (zMat3D *)&zmat3Dident )

/* METHOD:
 * zMat3DCreate, zMat3DCopy, zMat3DClear, zMat3DIdent
 * - creation, copy, cleanup and identification of 3D matrix.
 *
 * 'zMat3DCreate()' creates a 3D matrix like the following.
 *  | 'a11' 'a12' 'a13' |
 *  | 'a21' 'a22' 'a23' |
 *  | 'a31' 'a32' 'a33' |
 *
 * 'zMat3DCopy()' copies 3D matrix 'src' to 'dest'.
 *
 * 'zMat3DClear()' sets all factors of a 3D matrix 'm' as 0.
 *
 * 'zMat3DIdent()' sets a 3D matrix 'm' for the identity matrix.
 * [RETURN VALUE]
 * 'zMat3DCreate()', 'zMat3DClear()' and 'zMat3DIdent()'
 * return a pointer 'm'.
 *
 * 'zMat3DCopy()' returns no value.
 * [NOTES]
 * It is also possible to write simply *dest = *src instead of
 * 'zMat3DCopy()'. Actually, 'zMat3DCopy()' is defined
 * as macro(see "zeo_vec_mat3d.h").
 *
 * 'zMat3DClear()' is a macro for
 * zMat3DCreate( v, 0, 0, 0, 0, 0, 0, 0, 0, 0 ).
 *
 * 'zMat3DIdent()' is a macro for
 * zMat3DCreate( v, 1, 0, 0, 0, 1, 0, 0, 0, 1 ).
 */
__EXPORT zMat3D *zMat3DCreate(zMat3D *m,
  double a11, double a12, double a13,
  double a21, double a22, double a23,
  double a31, double a32, double a33);
#define zMat3DCopy(src,dest) ( *(dest) = *(src) )
#define zMat3DClear(m)       zMat3DCreate( (m), 0, 0, 0, 0, 0, 0, 0, 0, 0 )
#define zMat3DIdent(m)       zMat3DCreate( (m), 1, 0, 0, 0, 1, 0, 0, 0, 1 )

/* METHOD:
 * zMat3DT
 * - transposition of a 3D matrix.
 *
 * 'zMat3DT()' transpose a 3D matrix 'm' and put it into 'tm'.
 * [RETURN VALUE]
 * 'zMat3DT()' returns a pointer 'tm'.
 * [NOTES]
 * It is not permitted to let 'tm' point to the same address
 * with 'm'. When 'tm' is equal to 'm', anything might happen.
 */
__EXPORT zMat3D *zMat3DT(zMat3D *m, zMat3D *tm);

/* METHOD:
 * zMat3DRow, zMat3DCol
 * - abstraction of row/column vector from 3D matrix.
 * [SYNOPSIS]
 * zVec3D *zMat3DRow(zMat3D *m, int i, zVec3D *v);
 * zVec3D *zMat3DCol(zMat3D *m, int i, zVec3D *v);
 * [DESCRIPTION]
 * 'zMat3DRow()' abstract the 'i'th row from the 3D matrix
 * 'm' and put it into 'v'.
 *
 * 'zMat3DCol()' abstract the 'i'th column from the 3D matrix
 * 'm' and put it into 'v'.
 * [RETURN VALUE]
 * Both 'zMat3DRow()' and 'zMat3DCol()' return a pointer
 * to 'v'.
 */
__EXPORT zVec3D *zMat3DRow(zMat3D *m, int i, zVec3D *v);
__EXPORT zVec3D *zMat3DCol(zMat3D *m, int i, zVec3D *v);

/* ********************************************************** */
/* arithmetics
 * ********************************************************** */

/* METHOD:
 * zMat3DAdd, zMat3DSub, zMat3DRev, zMat3DMul, zMat3DDiv, zMat3DCat,
 * zMat3DAddDRC, zMat3DSubDRC, zMat3DRevDRC,
 * zMat3DMulDRC, zMat3DDivDRC, zMat3DCatDRC
 * - four rules of the arithmetics for 3D matrix.
 *
 * 'zMat3DAdd()' adds the two 3D matrices, 'm1' and 'm2'.
 * The result is put into 'm'.
 *
 * 'zMat3DSub()' subtracts 'm2' from 'm1'.
 * The result is put into 'm'.
 *
 * 'zMat3DRev()' reverses 'm'. The result is put into 'rm'.
 *
 * 'zMat3DMul()' multiplies 'm' by a scalar value 'k'.
 * The result is put into 'mm'.
 *
 * 'zMat3DDiv()' divides 'm' by 'k'.
 * The result is put into 'dm'.
 *
 * 'zMat3DCat()' concatenates multiplied 'm2' by 'k' to 'm1'.
 * The result is put into 'm'.
 *
 * 'zMat3DAddDRC()' adds 'm2' directly to the 'm1'.
 *
 * 'zMat3DSubDRC()' subtracts 'm2' directly from 'm1'.
 *
 * 'zMat3DRevDRC()' reverses 'm' directly.
 *
 * 'zMat3DMulDRC()' multiplies 'm' directly by 'k'.
 *
 * 'zMat3DDivDRC()' divides 'm' directly  by 'k'.
 *
 * 'zMat3DCatDRC()' concatenates multiplied 'm2' by 'k'
 * directly to 'm1'.
 * [RETURN VALUE]
 * Each function returns a pointer the resultant matrix.
 *
 * And, 'zMat3DDiv()' and 'zMat3DDivDRC()' return the
 * NULL pointer if the given scalar value is 0.
 */
__EXPORT zMat3D *zMat3DAdd(zMat3D *m1, zMat3D *m2, zMat3D *m);
__EXPORT zMat3D *zMat3DSub(zMat3D *m1, zMat3D *m2, zMat3D *m);
__EXPORT zMat3D *zMat3DRev(zMat3D *m, zMat3D *rm);
__EXPORT zMat3D *zMat3DMul(zMat3D *m, double k, zMat3D *mm);
__EXPORT zMat3D *zMat3DDiv(zMat3D *m, double k, zMat3D *dm);
__EXPORT zMat3D *zMat3DCat(zMat3D *m1, double k, zMat3D *m2, zMat3D *m);

#define zMat3DAddDRC(m1,m2)   zMat3DAdd(m1,m2,m1)
#define zMat3DSubDRC(m1,m2)   zMat3DSub(m1,m2,m1)
#define zMat3DRevDRC(m)       zMat3DRev(m,m)
#define zMat3DMulDRC(m,k)     zMat3DMul(m,k,m)
#define zMat3DDivDRC(m,k)     zMat3DDiv(m,k,m)
#define zMat3DCatDRC(m1,k,m2) zMat3DCat(m1,k,m2,m1)

/* METHOD:
 * zMat3DDyad - dyadic product.
 *
 * 'zMat3DDyad()' calculates a dyadic product of 'v1' and 'v2'.
 * The result is put into 'dyad' ( i.e. 'dyad' = 'v1' 'v2'^T ).
 * [RETURN VALUE]
 * 'zMat3DDyad()' returns a pointer 'dyad'.
 */
__EXPORT zMat3D *zMat3DDyad(zVec3D *v1, zVec3D *v2, zMat3D *dyad);

/* METHOD:
 * zVec3DOuterProdMat3D, zVec3DOuterProd2Mat3D
 * - creation of a matrix equivalent to the outer product
 *   of 3D vector.
 * [SYNOPSIS]
 * zMat3D *zVec3DOuterProdMat3D(zVec3D *v, zMat3D *m);
 * zMat3D *zVec3DOuterProd2Mat3D(zVec3D *v1, zVec3D *v2, zMat3D *m);
 * [DESCRIPTION]
 * 'zVec3DOuterProdMat3D()' computes an outer-product
 * 3D skew-symmetric matrix of a 3D vector 'v' ( i.e.
 * 'm' 'a' is equivalent to 'v' x 'a' ).
 *
 * Suppose 'v'=[vx vy vz]^T, 'm' will be
 *  |   0 -vz  vy |
 *  |  vz   0 -vx |
 *  | -vy  vx   0 |
 *
 * 'zVec3DOuterProd2Mat3D()' computes a twice-outer-product
 * 3D matrix of a 3D vector 'v1' and 'v2' ( i.e. 'm'
 * 'a' is equivalent to 'v1' x ( 'v2' x 'a' ) ).
 *
 * 'm' will be 'v2''v1'^T - 'v1'^T'v2' E, where E is a
 * 3x3 identity matrix.
 * [RETURN VALUE]
 * 'zVec3DOuterProdMat3D()' and 'zVec3DOuterProd2Mat3D()'
 * return a pointer 'm'.
 */
__EXPORT zMat3D *zVec3DOuterProdMat3D(zVec3D *v, zMat3D *m);
__EXPORT zMat3D *zVec3DOuterProd2Mat3D(zVec3D *v1, zVec3D *v2, zMat3D *m);

/* ********************************************************** */
/* inverse of a 3x3 matrix
 * ********************************************************** */

/*! \brief determinant of a 3x3 matrix.
 *
 * zMat3DDet() computes the determinant of an arbitrary 3x3 matrix \a m.
 * \retval the determinant of 'm'
 */
__EXPORT double zMat3DDet(zMat3D *m);

/*! \brief inverse of a 3x3 matrix.
 *
 * zMat3DInv() computes the inverse matrix of an arbitrary 3x3 matrix
 * \a m, and puts it into \a im. It does not assume that \a m is an
 * orthonormal matrix.
 * \retval a pointer \a im, if succeeds.
 * \retval the null pointer, if \a m is singular.
 * \notes
 * \a im has to point to a different address from \a m.
 * When \a im is equal to \a m, anything might happen.
 */
__EXPORT zMat3D *zMat3DInv(zMat3D *m, zMat3D *im);

/* ********************************************************** */
/* multiplication of a 3D vector by a 3x3 matrix
 * ********************************************************** */

/* METHOD:
 * zMulMatVec3D, zMulMatTVec3D,
 * zMulMatVec3DDRC, zMulMatTVec3DDRC,
 * zMulInvMatVec3D
 * - multiplication of 3D vector and 3D matrix.
 *
 * 'zMulMatVec3D()' multiplies a 3D vector 'v' by a 3D
 * matrix 'm'. The result is put into 'mv'.
 *
 * 'zMulMatTVec3D()' multiplies 'v' by the transpose
 * matrix of 'm'. The result is put into 'mv'.
 *
 * 'zMulMatVec3DDRC()' directly multiplies 'v' by 'm'.
 *
 * 'zMulMatTVec3DDRC()' directly multiplies 'v' by the
 * transpose matrix of 'm'.
 *
 * 'zMulInvMatVec3D()' multiplies 'v' by the inverse
 * of 'm'. The result is put into 'miv'.
 * [RETURN VALUE]
 * Each function returns the pointer to the result.
 */
__EXPORT zVec3D *zMulMatVec3D(zMat3D *m, zVec3D *v, zVec3D *mv);
__EXPORT zVec3D *zMulMatTVec3D(zMat3D *m, zVec3D *v, zVec3D *mv);

#define zMulMatVec3DDRC(m,v)  zMulMatVec3D(m,v,v)
#define zMulMatTVec3DDRC(m,v) zMulMatTVec3D(m,v,v)

__EXPORT zVec3D *zMulInvMatVec3D(zMat3D *m, zVec3D *v, zVec3D *miv);

/* METHOD:
 * zVec3DCatRatio - concatenate ratio of vector.
 *
 * 'zVec3DCatRatio()' calculates the concatenate
 * ratio of vector 'v' for three bases 'v1', 'v2'
 * and 'v3'. Namely, an arbitrary vector 'v' is
 * represented by a concatenation of three bases
 * vectors 'v1', 'v2' and 'v3' as follows.
 *  'v' = 'r1'*'v1' + 'r2'*'v2' + 'r3'*'v3'
 * 'zVec3DCatRatio()' computes 'r1', 'r2' and 'r3'
 * in the above equation and puts them into 'ratio'.
 * [NOTES]
 * Consequently, the array 'ratio' must have three
 * components at least. If not, anything might
 * happen.
 *
 * This function fails if 'v1', 'v2' and 'v3' are
 * not independent with each other.
 * [RETURN VALUE]
 * 'zVec3DCatRatio()' returns a pointer 'ratio'.
 */
__EXPORT double *zVec3DCatRatio(zVec3D *v1, zVec3D *v2, zVec3D *v3, zVec3D *v, double ratio[]);

/* ********************************************************** */
/* multiplication of a 3x3 matrix by another 3x3 matrix
 * ********************************************************** */

/* METHOD:
 * zMulMatMat3D, zMulMatTMat3D, zMulMatMatT3D,
 * zMulMatMat3DDRC, zMulMatTMat3DDRC, zMulMatMatT3DDRC,
 * zMulInvMatMat3D
 * - multiplication of 3D matrices.
 *
 * 'zMulMatMat3D()' multiplies a 3D matrix 'm2' by the other
 * 'm1' from leftside. The result is put into 'm'.
 *
 * 'zMulMatTMat3D()' multiplies 'm2' by the transpose matrix
 * of 'm1' from leftside. The result is put into 'm'.
 *
 * 'zMulMatMatT3D()' multiplies 'm1' by the transpose of
 * 'm2' from rightside. The result is put into 'm'.
 *
 * 'zMulMatMat3DDRC()' directly multiplies 'm2' by the
 * other 'm1' from leftside.
 *
 * 'zMulMatTMat3DDRC()' directly multiplies 'm2' by
 * the transpose of the other 'm1' from leftside.
 *
 * 'zMulMatMatT3DDRC()' directly multiplies 'm1' by
 * the transpose of 'm2' from rightside.
 *
 * 'zMulInvMatMat3D()' multiplies 'm2' by the inverse
 * of 'm1' from left side. The result is put into 'm'.
 * [RETURN VALUE]
 * Each function returns the pointer to the result.
 */
__EXPORT zMat3D *zMulMatMat3D(zMat3D *m1, zMat3D *m2, zMat3D *m);
__EXPORT zMat3D *zMulMatTMat3D(zMat3D *m1, zMat3D *m2, zMat3D *m);
__EXPORT zMat3D *zMulMatMatT3D(zMat3D *m1, zMat3D *m2, zMat3D *m);

#define zMulMatMat3DDRC(m1,m2)  zMulMatMat3D(m1,m2,m2)
#define zMulMatTMat3DDRC(m1,m2) zMulMatTMat3D(m1,m2,m2)
#define zMulMatMatT3DDRC(m1,m2) zMulMatMatT3D(m1,m2,m1)

__EXPORT zMat3D *zMulInvMatMat3D(zMat3D *m1, zMat3D *m2, zMat3D *m);

/* ********************************************************** */
/* multiplication of a 6D spatial vector by a 3x3 matrix
 * ********************************************************** */

/* METHOD:
 * zMulMatVec6D, zMulMatTVec6D,
 * zMulMatVec6DDRC, zMulMatTVec6DDRC,
 * - multiplication of 6D vector and 3D matrix.
 *
 * 'zMulMatVec6D()' multiplies a 6D vector 'v' by a 3D
 * matrix 'm'. The result is put into 'mv'.
 *
 * 'zMulMatTVec6D()' multiplies 'v' by the transpose matrix
 * of 'm'. The result is put into 'mv'.
 *
 * 'zMulMatVec6DDRC()' directly multiplies 'v' by 'm'.
 *
 * 'zMulMatTVec6DDRC()' directly multiplies 'v' by the
 * transpose matrix of 'm'.
 * [RETURN VALUE]
 * Each function returns the pointer to the result 6D vector.
 */
__EXPORT zVec6D *zMulMatVec6D(zMat3D *m, zVec6D *v, zVec6D *mv);
__EXPORT zVec6D *zMulMatTVec6D(zMat3D *m, zVec6D *v, zVec6D *mv);

#define zMulMatVec6DDRC(m,v)  zMulMatVec6D(m,v,v)
#define zMulMatTVec6DDRC(m,v) zMulMatTVec6D(m,v,v)

/* ********************************************************** */
/* rotation
 * ********************************************************** */

/* METHOD:
 * zMat3DRotRoll, zMat3DRotPitch, zMat3DRotYaw,
 * zMat3DRotRollSC, zMat3DRotPitchSC, zMat3DRotYawSC,
 * zMat3DRotRollDRC, zMat3DRotPitchDRC, zMat3DRotYawDRC,
 * zMat3DRotRollSCDRC, zMat3DRotPitchSCDRC, zMat3DRotYawSCDRC
 * - rotate matrix about a base specified axis.
 *
 * 'zMat3DRotRoll()', 'zMat3DRotPitch()' and 'zMat3DRotYaw()'
 * rotate a matrix 'm' with the angle 'angle' about x-axis,
 * y-axis and z-axis, respectively.
 * 'angle' is in radian.
 * The result is put into 'rm'.
 *
 * 'zMat3DRotRollSC()', 'zMat3DRotPitchSC()' and 'zMat3DRotYawSC()'
 * rotate a matrix 'm' not with the angle but with its trigonometric
 * values. 's' and 'c' are for sine and cosine values, respectively.
 *
 * 'zMat3DRotRollDRC()', 'zMat3DRotPitchDRC()', 'zMat3DRotYawDRC()',
 * 'zMat3DRotRollSCDRC()', 'zMat3DRotPitchSCDRC()' and
 * 'zMat3DRotYawSCDRC()' are destructive versions of the above
 * functions, directly updating a given matrix 'm'.
 * [RETURN VALUE]
 * 'zMat3DRotRoll()', 'zMat3DRotPitch()', 'zMat3DRotYaw()'
 * 'zMat3DRotRollSC()', 'zMat3DRotPitchSC()' and 'zMat3DRotYawSC()'
 * return a pointer 'rm'.
 *
 * 'zMat3DRotRollDRC()', 'zMat3DRotPitchDRC()', 'zMat3DRotYawDRC()'
 * 'zMat3DRotRollSCDRC()', 'zMat3DRotPitchSCDRC()' and
 * 'zMat3DRotYawSCDRC()' return a pointer 'm'.
 */
__EXPORT zMat3D *zMat3DRotRoll(zMat3D *m, double angle, zMat3D *rm);
__EXPORT zMat3D *zMat3DRotPitch(zMat3D *m, double angle, zMat3D *rm);
__EXPORT zMat3D *zMat3DRotYaw(zMat3D *m, double angle, zMat3D *rm);
__EXPORT zMat3D *zMat3DRotRollSC(zMat3D *m, double s, double c, zMat3D *rm);
__EXPORT zMat3D *zMat3DRotPitchSC(zMat3D *m, double s, double c, zMat3D *rm);
__EXPORT zMat3D *zMat3DRotYawSC(zMat3D *m, double s, double c, zMat3D *rm);
__EXPORT zMat3D *zMat3DRotRollDRC(zMat3D *m, double theta);
__EXPORT zMat3D *zMat3DRotPitchDRC(zMat3D *m, double theta);
__EXPORT zMat3D *zMat3DRotYawDRC(zMat3D *m, double theta);
__EXPORT zMat3D *zMat3DRotRollSCDRC(zMat3D *m, double s, double c);
__EXPORT zMat3D *zMat3DRotPitchSCDRC(zMat3D *m, double s, double c);
__EXPORT zMat3D *zMat3DRotYawSCDRC(zMat3D *m, double s, double c);

/* METHOD:
 * zMat3DZYX, zMat3DZYXSC, zMat3DToZYX,
 * zMat3DZYZ, zMat3DZYZSC, zMat3DToZYZ,
 * zMat3DAA, zMat3DToAA
 * - 3D attitude alternation with matrix.
 *
 * 'zMat3DZYX()' calculates a 3D attitude matrix from z-y-x
 * Eulerian angle. The identity matrix is rotated firstly by
 * 'azim' about z-axis, secondly 'elev' about y-axis, and
 * finally 'tilt' about x-axis.
 * The result is put into 'm'.
 *
 * 'zMat3DZYXSC()' directly accepts trigonometric values of
 * z-y-x Eulerian angle. 'sa'/'ca', 'se'/'ce' and 'st'/'ct'
 * are for azimuth, elevation and tilt angles, respectively.
 *
 * 'zMat3DToZYX()' is the inverse transformation of
 * 'zMat3DZYX()' from a matrix to a quasi 3D vector which
 * respresents z-y-x Eulerian angle 'angle' in the order of
 * azimuth, elevation and tilt. The result is put into 'angle'.
 *
 * 'zMat3DZYZ()' calculates a 3D attitude matrix from z-y-z
 * Eulerian angle. The identity matrix is rotated firstly by
 * 'heading' about z-axis, secondly 'pitch' about y-axis,
 * and finally 'bank' about z-axis.
 * The result is put into 'm'.
 *
 * 'zMat3DZYZSC()' directly accepts trigonometric values of
 * z-y-x Eulerian angle. 'sh'/'ch', 'sp'/'cp' and 'sb'/'cb'
 * are for heading, pitch and bank angles, respectively.
 *
 * 'zMat3DToZYZ()' is the inverse transformation of
 * 'zMat3DZYZ()' from a matrix to a quasi 3D vector which
 * respresents z-y-z Eulerian angle 'angle' in the order of
 * heading, pitch and bank. The result is put into 'angle'.
 *
 * 'zMat3DAA()' calculates a 3D attitude matrix from the
 * equivalent angle-axis vector 'aa'. The identity matrix
 * is rotated about the direction of 'aa' with the angle
 * equal to the norm of 'aa'.
 * The result is put into 'm'.
 *
 * 'zMat3DToAA()' is the inverse transformation of
 * 'zMat3DAA()' from a matrix to the equivalent angle-axis
 * vector. The result is put into 'aa'.
 * [RETURN VALUE]
 * 'zMat3DZYX()', 'zMat3DZYXSC()', 'zMat3DZYZ()', 'zMat3DZYZSC()'
 * and 'zMat3DAA()' return a pointer 'm'.
 *
 * 'zMat3DToZYX()' and 'zMat3DToZYZ()' return a pointer 'angle'.
 *
 * 'zMat3DToAA()' returns a pointer 'aa'.
 */
__EXPORT zMat3D *zMat3DZYX(zMat3D *m, double azim, double elev, double tilt);
__EXPORT zMat3D *zMat3DZYXSC(zMat3D *m, double sa, double ca, double se, double ce, double st, double ct);
__EXPORT zVec3D *zMat3DToZYX(zMat3D *m, zVec3D *angle);
__EXPORT zMat3D *zMat3DZYZ(zMat3D *m, double heading, double pitch, double bank);
__EXPORT zMat3D *zMat3DZYZSC(zMat3D *m, double sh, double ch, double sp, double cp, double sb, double cb);
__EXPORT zVec3D *zMat3DToZYZ(zMat3D *m, zVec3D *angle);
__EXPORT zMat3D *zMat3DAA(zMat3D *m, zVec3D *aa);
__EXPORT zVec3D *zMat3DToAA(zMat3D *m, zVec3D *aa);

/* METHOD:
 * zRotMat3D, zRotMat3DInv, zRotMat3DDRC, zRotMat3DInvDRC,
 * - rotational multiplication of 3D matrices.
 *
 * 'zRotMat3D()' multiplies a 3D matrix 'm' by 'r' from
 * leftside, then does it by a transpose of 'r' from
 * rightside, and puts the result into 'rm'.
 * Namely, 'rm' = 'r' 'm' 'r'^T.
 *
 * 'zRotMat3DInv()' is an opposite computation of
 * 'zRotMat3D()'. Namely, 'rm' = 'r'^T 'm' 'r'.
 *
 * 'zRotMat3DDRC()' is the same computation with
 * 'zRotMat3D()' except that it puts the result directly
 * into 'm'.
 * 'zRotMat3DInvDRC()' is the same computation with
 * 'zRotMat3DInv()' except that it puts the result
 * directly into 'm'.
 * [RETURN VALUE]
 * Each function returns the pointer to the result.
 */
__EXPORT zMat3D *zRotMat3D(zMat3D *r, zMat3D *m, zMat3D *rm);
__EXPORT zMat3D *zRotMat3DInv(zMat3D *r, zMat3D *m, zMat3D *rm);

#define zRotMat3DDRC(r,m)    zRotMat3D(r,m,m)
#define zRotMat3DInvDRC(r,m) zRotMat3DInv(r,m,m)

/* METHOD:
 * zMulVecOPMat3D, zMulVecOPMat3DDRC
 * - multiply cross product of vector and matrix.
 *
 * 'zMulVecOPMat3D()' multiplies a matrix 'm' by the
 * outer product of a vector 'ohm'. The result is put
 * into 'mv'. Namely, 'mv = ohm x m'.
 *
 * 'zMulVecOPMat3DDRC()' directly multiplies 'm' by
 * the outer product of 'ohm'.
 * [RETURN VALUE]
 * 'zMulVecOPMat3D()' returns a pointer 'mv'.
 * 'zMulVecOPMat3DDRC()' returns a pointer 'm'.
 */
__EXPORT zMat3D *zMulVecOPMat3D(zVec3D *ohm, zMat3D *m, zMat3D *mv);
#define zMulVecOPMat3DDRC(o,m) zMulVecOPMat3D( o, m, m )

/* METHOD:
 * zMat3DRot - rotate matrix about arbitrary axis.
 *
 * 'zMat3DRot()' rotates a 3D attitude matrix 'm' by
 * angle-axis vector 'aa'. The axis of rotation is
 * in parallel to 'aa', and the norm of 'aa' is the
 * rotation angle in radian. The direction of rotation
 * is according to right-handed screw rule.
 * The result is put into 'rm'.
 * [RETURN VALUE]
 * 'zMat3DRot()' returns a pointer 'rm'.
 */
__EXPORT zMat3D *zMat3DRot(zMat3D *m, zVec3D *aa, zMat3D *rm);

__EXPORT zMat3D *zMat3DRotCat(zMat3D *m, zVec3D *omega, double dt, zMat3D *rm);

/*! \brief cascade an angle-axis vector to another.
 */
__EXPORT zVec3D *zAACascade(zVec3D *aa1, zVec3D *aa2, zVec3D *aa);

/*! \brief error vector between two attitude matrices.
 *
 * zMat3DError() calculates the error vector, namely,
 * the equivalent angle-axis vector from \a m2 to \a m1
 * (note the order). The result is put into \a err.
 * \retval \a err
 */
__EXPORT zVec3D *zMat3DError(zMat3D *m1, zMat3D *m2, zVec3D *err);

/*! \brief error between two angle-axis vectors.
 *
 * zAAError() calculates the error vector, namely, the angle-axis
 * vector from an attitude represented by \a a2 to another \a a1
 * (note the order). The result is put into \a err.
 * \retval \a err
 */
__EXPORT zVec3D *zAAError(zVec3D *a1, zVec3D *a2, zVec3D *err);

/* ********************************************************** */
/* differential kinematics
 * ********************************************************** */

__EXPORT zVec3D *zMat3DDif(zMat3D *m, zMat3D *mnew, double dt, zVec3D *omega);

/* ********************************************************** */
/* eigensystem
 * ********************************************************** */

/* METHOD:
 * zMat3DSymEig
 * - eigenvalues of a symmetric 3x3 matrix by Jacobi's method.
 *
 * 'zMat3DSymEig()' calculates eigenvalues and
 * eigenvectors of a symmetric 3x3 matrix 'm'
 * with Jacobi s method.
 * Each eigenvalue and eigenvector are stored
 * in 'eval' and 'evec' in a corresponding order.
 * [RETURN VALUE]
 * 'zMat3DSymEig()' returns no value.
 * [NOTES]
 * 'm' must be symmetric. Otherwise, the correct
 * result will not be expected.
 */
__EXPORT void zMat3DSymEig(zMat3D *m, double eval[], zVec3D evec[]);

/* ********************************************************** */
/* I/O
 * ********************************************************** */

/* METHOD:
 * zMat3DFRead, zMat3DRead, zMat3DFWrite, zMat3DWrite
 * - input/output of 3D matrix.
 *
 * 'zMat3DFRead()' reads nine values from the current
 * position of the file 'fp', and creates the 3D matrix
 * 'm' from them.
 * 'zMat3DRead()' simply reads values from the standard
 * input.
 *
 * 'zMat3DFWrite()' writes the 3D matrix 'm' to the current
 * position of the file 'fp' in the following style.
 *  {
 *   a11, a12, a13
 *   a21, a22, a23
 *   a31, a32, a33
 *  }
 * When the NULL pointer is given, it writes the following string.
 *  (null 3D matrix)
 * 'zMat3DWrite()' simply writes 'm' to the standard out.
 * [RETURN VALUE]
 * 'zMat3DFRead()' and 'zMat3DRead()' return a pointer 'm'.
 *
 * 'zMat3DFWrite()' and 'zMat3DWrite()' return no value.
 */
__EXPORT zMat3D *zMat3DFRead(FILE *fp, zMat3D *m);
#define zMat3DRead(m) zMat3DFRead( stdin, (m) )
__EXPORT void zMat3DFWrite(FILE *fp, zMat3D *m);
#define zMat3DWrite(m) zMat3DFWrite( stdout, (m) )

/* METHOD:
 * zMat3DFWriteXML - xml output.
 * ... yet testing.
 */
__EXPORT void zMat3DFWriteXML(FILE *fp, zMat3D *m);

__END_DECLS

#include <zeo/zeo_vec3d_pca.h> /* principal component analysis */

#endif /* __ZEO_MAT3D_H__ */
