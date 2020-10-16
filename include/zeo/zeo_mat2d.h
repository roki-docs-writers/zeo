/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_mat2d.h - 2D matrix.
 */

#ifndef __ZEO_MAT2D_H__
#define __ZEO_MAT2D_H__

#include <zeo/zeo_vec2d.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zMat2D
 * 2D matrix class
 * ********************************************************** */

typedef double zMat2D[4];

/* OBJECT:
 * zmat2Dzero, zmat2Dident
 * - 2D zero matrix and identity matrix.
 */
extern const zMat2D zmat2Dzero;
extern const zMat2D zmat2Dident;
#define ZMAT2DZERO  ( (zMat2D)zmat2Dzero )
#define ZMAT2DIDENT ( (zMat2D)zmat2Dident )

/* METHOD:
 * zMat2DCreate, zMat2DCopy, zMat2DClear, zMat2DIdent
 * - creation, copy, cleanup and identification of 2D matrix.
 *
 * 'zMat2DCreate()' creates a 2D matrix like the following.
 *  | 'a11' 'a12' |
 *  | 'a21' 'a22' |
 * #
 * 'zMat2DCopy()' copies 2D matrix 'src' to 'dest'.
 * #
 * 'zMat2DClear()' sets all factors of a 2D matrix 'm' as 0.
 * #
 * 'zMat2DIdent()' sets a 2D matrix 'm' for the identity matrix.
 * [RETURN VALUE]
 * 'zMat2DCreate()', 'zMat2DClear()' and 'zMat2DIdent()'
 * return a pointer 'm'.
 * #
 * 'zMat2DCopy()' returns no value.
 * [NOTES]
 * It is also possible to write simply *dest = *src instead of
 * 'zMat2DCopy()'. Actually, 'zMat2DCopy()' is defined
 * as macro(see "zeo_mat2d.h").
 * #
 * 'zMat2DClear()' is a macro for zMat2DCreate( v, 0, 0, 0, 0 ).
 * #
 * 'zMat2DIdent()' is a macro for zMat2DCreate( v, 1, 0, 0, 1 ).
 */
__EXPORT double *zMat2DCreate(zMat2D m, double a11, double a12, double a21, double a22);
#define zMat2DCopy(src,dest) memcpy( (dest), (src), sizeof(double)*4 )
#define zMat2DClear(m)       zMat2DCreate( (m), 0, 0, 0, 0 )
#define zMat2DIdent(m)       zMat2DCreate( (m), 1, 0, 0, 1 )

/*! \brief transpose a 2D matrix.
 *
 * zMat2DT() transposes a 2D matrix \a m and puts it into \a tm.
 * \retval \a tm
 */
__EXPORT double *zMat2DT(zMat2D m, zMat2D tm);

/* METHOD:
 * zMat2DRow, zMat2DCol
 * - abstraction of row/column vector from 2D matrix.
 *
 * 'zMat2DRow()' abstract the row vectors of the 2D matrix
 * 'm'. The first and second roww of 'm' is abstracted
 * to 'r1' and 'r2', respectively. When 'r1' or 'r2'
 * is the null pointer, it is ignored.
 * #
 * 'zMat2DCol()' abstract the column vectors of the 2D
 * matrix 'm'. The first and second columns of 'm' is
 * abstracted to 'c1' and 'c2', respectively. When 'c1'
 * or 'r2' is the null pointer, it is ignored.
 * [RETURN VALUE]
 * Neither 'zMat2DRow()' nor 'zMat2DCol()' return any
 * values.
 */
__EXPORT void zMat2DRow(zMat2D m, zVec2D r1, zVec2D r2);
__EXPORT void zMat2DCol(zMat2D m, zVec2D c1, zVec2D c2);

/* ********************************************************** */
/* arithmetics
 * ********************************************************** */

/* METHOD:
 * zMat2DAdd, zMat2DSub, zMat2DRev, zMat2DMul, zMat2DDiv, zMat2DCat,
 * zMat2DAddDRC, zMat2DSubDRC, zMat2DRevDRC,
 * zMat2DMulDRC, zMat2DDivDRC, zMat2DCatDRC
 * - four rules of the arithmetics for 2D matrix.
 *
 * 'zMat2DAdd()' adds the two 2D matrices, 'm1' and 'm2'.
 * The result is put into 'm'.
 * #
 * 'zMat2DSub()' subtracts 'm2' from 'm1'.
 * The result is put into 'm'.
 * #
 * 'zMat2DRev()' reverses 'm'. The result is put into 'rm'.
 * #
 * 'zMat2DMul()' multiplies 'm' by a scalar value 'k'.
 * The result is put into 'mm'.
 * #
 * 'zMat2DDiv()' divides 'm' by 'k'.
 * The result is put into 'dm'.
 * #
 * 'zMat2DCat()' concatenates multiplied 'm2' by 'k' to 'm1'.
 * The result is put into 'm'.
 * #
 * 'zMat2DAddDRC()' adds 'm2' directly to the 'm1'.
 * #
 * 'zMat2DSubDRC()' subtracts 'm2' directly from 'm1'.
 * #
 * 'zMat2DRevDRC()' reverses 'm' directly.
 * #
 * 'zMat2DMulDRC()' multiplies 'm' directly by 'k'.
 * #
 * 'zMat2DDivDRC()' divides 'm' directly  by 'k'.
 * #
 * 'zMat2DCatDRC()' concatenates multiplied 'm2' by 'k'
 * directly to 'm1'.
 * [RETURN VALUE]
 * Each function returns a pointer the resultant matrix.
 * #
 * And, 'zMat2DDiv()' and 'zMat2DDivDRC()' return the
 * NULL pointer if the given scalar value is 0.
 */
__EXPORT double *zMat2DAdd(zMat2D m1, zMat2D m2, zMat2D m);
__EXPORT double *zMat2DSub(zMat2D m1, zMat2D m2, zMat2D m);
__EXPORT double *zMat2DRev(zMat2D m, zMat2D rm);
__EXPORT double *zMat2DMul(zMat2D m, double k, zMat2D mm);
__EXPORT double *zMat2DDiv(zMat2D m, double k, zMat2D dm);
__EXPORT double *zMat2DCat(zMat2D m1, double k, zMat2D m2, zMat2D m);

#define zMat2DAddDRC(m1,m2)   zMat2DAdd(m1,m2,m1)
#define zMat2DSubDRC(m1,m2)   zMat2DSub(m1,m2,m1)
#define zMat2DRevDRC(m)       zMat2DRev(m,m)
#define zMat2DMulDRC(m,k)     zMat2DMul(m,k,m)
#define zMat2DDivDRC(m,k)     zMat2DDiv(m,k,m)
#define zMat2DCatDRC(m1,k,m2) zMat2DCat(m1,k,m2,m1)

/* METHOD:
 * zMat2DDyad - dyadic product.
 *
 * 'zMat2DDyad()' calculates a dyadic product of 'v1' and 'v2'.
 * The result is put into 'dyad' ( i.e. 'dyad' = 'v1' 'v2'^T ).
 * [RETURN VALUE]
 * 'zMat2DDyad()' returns a pointer 'dyad'.
 */
__EXPORT double *zMat2DDyad(zVec2D v1, zVec2D v2, zMat2D dyad);

/* ********************************************************** */
/* inverse of a 2x2 matrix
 * ********************************************************** */

/* METHOD:
 * zMat2DDet
 * - determinant of a 2x2 matrix.
 *
 * 'zMat2DDet()' computes the determinant of a 2x2 matrix 'm'.
 * [RETURN VALUE]
 * 'zMat2DDet()' returns the determinant of 'm'.
 */
__EXPORT double zMat2DDet(zMat2D m);

/*! \brief inverse of a 2D matrix.
 *
 * zMat2DInv() computes the inverse of an arbitrary 2x2 matrix
 * \a m. The result is put into \a im.
 * \retval \a im
 * \notes
 * \a im has to be different from \a m. If \a tm is equal to \a m,
 * anything might happen.
 */
__EXPORT double *zMat2DInv(zMat2D m, zMat2D im);

/* ********************************************************** */
/* multiplication of a 2D vector by a 2x2 matrix
 * ********************************************************** */

/* METHOD:
 * zMulMatVec2D, zMulMatTVec2D, zMulInvMatVec2D
 * - multiplication of 2D vector and 2D matrix.
 *
 * 'zMulMatVec2D()' multiplies a 2D vector 'v' by a 2D
 * matrix 'm'. The result is put into 'mv'.
 * #
 * 'zMulMatTVec2D()' multiplies 'v' by the transpose
 * matrix of 'm'. The result is put into 'mv'.
 * #
 * 'zMulInvMatVec2D()' multiplies 'v' by the inverse
 * matrix of 'm'. The result is put into 'mv'.
 * [RETURN VALUE]
 * Each function returns the pointer to the result.
 */
__EXPORT double *zMulMatVec2D(zMat2D m, zVec2D v, zVec2D mv);
__EXPORT double *zMulMatTVec2D(zMat2D m, zVec2D v, zVec2D mv);
__EXPORT double *zMulInvMatVec2D(zMat2D m, zVec2D v, zVec2D mv);

/* ********************************************************** */
/* multiplication of a 2x2 matrix by another 2x2 matrix
 * ********************************************************** */

/* METHOD:
 * zMulMatMat2D, zMulMatTMat2D, zMulMatMatT2D, zMulInvMatMat2D
 * - multiplication of 2D matrices.
 *
 * 'zMulMatMat2D()' multiplies a 2D matrix 'm2' by the other
 * 'm1' from leftside. The result is put into 'm'.
 * #
 * 'zMulMatTMat2D()' multiplies 'm2' by the transpose matrix
 * of 'm1' from leftside. The result is put into 'm'.
 * #
 * 'zMulMatMatT2D()' multiplies 'm1' by the transpose of
 * 'm2' from rightside. The result is put into 'm'.
 * #
 * 'zMulMatMatT2D()' multiplies 'm2' by the inverse of
 * 'm1' from leftside. The result is put into 'm'.
 * [RETURN VALUE]
 * Each function returns the pointer to the result.
 */
__EXPORT double *zMulMatMat2D(zMat2D m1, zMat2D m2, zMat2D m);
__EXPORT double *zMulMatTMat2D(zMat2D m1, zMat2D m2, zMat2D m);
__EXPORT double *zMulMatMatT2D(zMat2D m1, zMat2D m2, zMat2D m);
__EXPORT double *zMulInvMatMat2D(zMat2D m1, zMat2D m2, zMat2D m);

/* ********************************************************** */
/* rotation
 * ********************************************************** */

/* METHOD:
 * zMat2DRot, zMat2DRotSC
 * - rotate matrix.
 *
 * 'zMat2DRot()' rotates a matrix 'm' with the angle
 * 'angle'. 'angle' is in radian. The result is put into 'rm'.
 * #
 * 'zMat2DRotSC()' rotates 'm' not with the angle but
 * with its trigonometric values. 's' and 'c' are for
 * sine and cosine values, respectively.
 * [RETURN VALUE]
 * 'zMat2DRot()' and 'zMat2DRotSC()' return a pointer 'rm'.
 */
__EXPORT double *zMat2DRot(zMat2D m, double angle, zMat2D rm);
__EXPORT double *zMat2DRotSC(zMat2D m, double s, double c, zMat2D rm);

/*! \brief error vector between two attitude matrices.
 *
 * zMat2DError() calculates the error angle from \a m2 to \a m1
 * (note the order).
 * \retval the error angle computed.
 */
__EXPORT double zMat2DError(zMat2D m1, zMat2D m2);

/* ********************************************************** */
/* I/O
 * ********************************************************** */

/* METHOD:
 * zMat2DFRead, zMat2DRead, zMat2DFWrite, zMat2DWrite
 * - input/output of 2D matrix.
 *
 * 'zMat2DFRead()' reads four values from the current
 * position of the file 'fp', and creates the 2D matrix
 * 'm' from them. 'zMat2DRead()' simply reads values
 * from the standard input.
 * #
 * 'zMat2DFWrite()' writes the 2D matrix 'm' to the current
 * position of the file 'fp' in the following style.
 *  {
 *   a11, a12
 *   a21, a22
 *  }
 * When the NULL pointer is given, it writes the following string.
 *  (null 2D matrix)
 * 'zMat2DWrite()' simply writes 'm' to the standard out.
 * [RETURN VALUE]
 * 'zMat2DFRead()' and 'zMat2DRead()' return a pointer 'm'.
 * #
 * 'zMat2DFWrite()' and 'zMat2DWrite()' return no value.
 */
__EXPORT double *zMat2DFRead(FILE *fp, zMat2D m);
#define zMat2DRead(m) zMat2DFRead( stdin, (m) )
__EXPORT void zMat2DFWrite(FILE *fp, zMat2D m);
#define zMat2DWrite(m) zMat2DFWrite( stdout, (m) )

__END_DECLS

#endif /* __ZEO_MAT2D_H__ */
