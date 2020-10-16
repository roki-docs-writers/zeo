/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_frame - 3D frame.
 */

#ifndef __ZEO_FRAME_H__
#define __ZEO_FRAME_H__

#include <zeo/zeo_mat3d.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zFrame3D
 * 3D frame class
 * ********************************************************** */

typedef struct{
  zVec3D pos;
  zMat3D att;
} zFrame3D;

#define zFrame3DPos(f) ( &(f)->pos )
#define zFrame3DAtt(f) ( &(f)->att )

#define zFrame3DSetPos(f,p) zVec3DCopy( p, zFrame3DPos(f) )
#define zFrame3DSetAtt(f,r) zMat3DCopy( r, zFrame3DAtt(f) )

/* OBJECT:
 * zframe3Dident
 * - the identity frame.
 */
extern const zFrame3D zframe3Dident;
#define ZFRAME3DIDENT ( (zFrame3D *)&zframe3Dident )

/* METHOD:
 * zFrame3DCreate, zFrame3DCopy, zFrame3DIdent
 * - create, copy and initialize 3D frame.
 *
 * 'zFrame3DCreate()' creates a 3D frame 'f' with a
 * position of the original point 'p' and the attitude
 * matrix of 'm'.
 *
 * 'zFrame3DCopy()' copies a 3D frame 'src' to 'dest'.
 *
 * 'zFrame3DIdent()' initializes 'f', setting the
 * original vector for { 0, 0, 0 } and the attitude
 * for the identity matrix.
 * [RETURN VALUE]
 * 'zFrame3DCreate()' returns a pointer to 'f'.
 *
 * 'zFrame3DCopy()' and 'zFrame3DIdent()' return no value.
 * [NOTES]
 * It is also possible to write simply *'dest' = *'src'
 * instead of 'zFrame3DCopy()'. Actually, 'zFrame3DCopy()'
 * is defined as macro(see "zeo_frame.h").
 *
 * 'zFrame3DIdent()' is a macro zFrame3DCopy( ZIDENTFRAME3D, 'f' ).
 */
__EXPORT zFrame3D *zFrame3DCreate(zFrame3D *f, zVec3D *p, zMat3D *m);
#define zFrame3DCopy(src,dest) ( *(dest) = *(src) )
#define zFrame3DIdent(f) zFrame3DCopy( ZIDENTFRAME3D, f )

/* METHOD:
 * zFrame3DInv, zFrame3DCascade, zFrame3DXfer
 * - inverse and cascaded frame.
 *
 * 'zFrame3DInv()' calculates the inverse transformation
 * frame of the given 3D frame 'f'.
 * The result is put into 'fi'.
 *
 * 'zFrame3DCascade()' cascades the 3D frame 'c1' to
 * the other 'f2'. The result is put into 'f'.
 * That is, suppose 'c1' is the transformation from a
 * frame named S1 to another frame S0, and 'f2'
 * for that from a frame S2 to S1, 'f' is the
 * transformation from S2 to S0.
 *
 * 'zFrame3DXfer()' calculates the transformation
 * frame from the frame 'f1' to the other 'f2'.
 * The result is put into 'f'.
 * Suppose 'f1' is the transformation from a frame named S1
 * to another frame S0, and 'f2' is that from a frame S2
 * to S0, 'f' is the transformation from S2 to S1, namely,
 * the cascade of the inverse of 'f1' and 'f2'.
 * [RETURN VALUE]
 * Each of these functions returns a pointer to the resultant
 * frame.
 * [NOTES]
 * 'zFrame3DInv()' expects that the attitude matrix of
 * 'f' is a homogeneous unitary matrix.
 * For these three functions, it is not permitted to let
 * any argument point to the same address with the other
 * arguments.
 * When some of them are equal, anything might happen.
 */
__EXPORT zFrame3D *zFrame3DInv(zFrame3D *f, zFrame3D *fi);
__EXPORT zFrame3D *zFrame3DCascade(zFrame3D *f1, zFrame3D *f2, zFrame3D *f);
__EXPORT zFrame3D *zFrame3DXfer(zFrame3D *f1, zFrame3D *f2, zFrame3D *f);

/* METHOD:
 * zXfer3D, zXfer3DInv, zXfer3DDRC, zXfer3DInvDRC
 * - transformation of 3D vector.
 *
 * 'zXfer3D()' transforms the given 3D vector 'v' by
 * the frame 'f'. The result is put into 'tv'.
 * In other word, suppose 'f' is for the transformation from
 * a frame S1 to S0, 'tv' points where 'v' points in S1 with
 * respect to S0.
 *
 * On the same assumption for 'f', 'zXfer3DInv()'
 * transforms the given 3D vector 'v' in S0 to what is with
 * respect to S1 by the inverse frame of 'f'.
 * The result is put into 'tv'.
 *
 * 'zXfer3DDRC()' directly transforms the 3D vector 'v'
 * by the frame 'f'.
 *
 * 'zXfer3DInvDRC()' directly transforms the 3D
 * vector 'v' by the inverse frame of 'f'.
 * [RETURN VALUE]
 * 'zXfer3D()', 'zXfer3DInv()', 'zXfer3DDRC()' and
 * 'zXfer3DInvDRC()' return a pointer to the result.
 * [NOTES]
 * 'zXfer3DInv()' expects that the attitude matrix of
 * 'f' is a homogeneous unitary matrix.
 */
__EXPORT zVec3D *zXfer3D(zFrame3D *f, zVec3D *v, zVec3D *tv);
__EXPORT zVec3D *zXfer3DInv(zFrame3D *f, zVec3D *v, zVec3D *tv);

#define zXfer3DDRC(f,v)    zXfer3D(f,v,v)
#define zXfer3DInvDRC(f,v) zXfer3DInv(f,v,v)

/*! \brief transfer 6D vector
 */
__EXPORT zVec6D *zXfer6DLin(zFrame3D *f, zVec6D *v, zVec6D *vc);
__EXPORT zVec6D *zXfer6DAng(zFrame3D *f, zVec6D *v, zVec6D *vc);

/* zFrame3DTwist
 * - twist a frame by a torsion vector
 *   (position offset & angle-axis rotation).
 */
__EXPORT zFrame3D *zFrame3DTwist(zFrame3D *f1, zVec6D *t, zFrame3D *f2);

/*! \brief error between two frames.
 *
 * zFrame3DError() calculates the error vector between two frames
 * from \a f2 to \a f1 (note the order).
 * The result is put into \a err. Intuitively, \a err equals to
 * \a f1 minus \a f2.
 * \retval \a err
 * \sa zMat3DError
 */
__EXPORT zVec6D *zFrame3DError(zFrame3D *f1, zFrame3D *f2, zVec6D *err);

/* METHOD:
 * zFrame3DZYX, zFrame3DZYZ, zFrame3DDH
 * - creation of frame from handy expression.
 *
 * 'zFrame3DZYX()' creates a frame whose original point
 * is ( 'x', 'y', 'z' ) and attitude is expressed by
 * z-y-x Eulerian angle (refer zMat3DZYX).
 *
 * 'zFrame3DZYZ()' creates a frame whose original point
 * is ( 'x', 'y', 'z' ) and attitude is expressed by
 * z-y-z Eulerian angle (refer zMat3DZYZ).
 *
 * 'zFrame3DDH()' creates a frame from Denaviet-Hartenberg(D-H)
 * parameters ('a', 'alpha', 'd', 'theta').
 *
 * For any of these functions, the result is put into 'f'.
 * [RETURN VALUE]
 * All these functions return a pointer to 'f'.
 * [SEE ALSO]
 * zMat3DZYX, zMat3DZYZ
 */
__EXPORT zFrame3D *zFrame3DZYX(zFrame3D *f, double x, double y, double z, double azim, double elev, double tilt);
__EXPORT zFrame3D *zFrame3DZYZ(zFrame3D *f, double x, double y, double z, double heading, double pitch, double bank);
__EXPORT zFrame3D *zFrame3DAA(zFrame3D *f, double x, double y, double z, double xa, double ya, double za);
__EXPORT zFrame3D *zFrame3DDH(zFrame3D *f, double a, double alpha, double d, double theta);

/* METHOD:
 * zArrayToFrame3DZYX, zFrame3DZYXToArray,
 * zArrayToFrame3DZYZ, zFrame3DZYZToArray,
 * zVec6DToFrame3D, zFrame3DToVec6D,
 * zVec6DToFrame3DZYX, zFrame3DZYXToVec6D,
 * zVec6DToFrame3DZYZ, zFrame3DZYZToVec6D
 * - conversion from/to 3D frame to/from array or 6D vector.
 *
 * 'zArrayToFrame3DZYX()' converts the given array of
 * double-precision floating point values with
 * [ x, y, z, azimuth, elevation, tilt ] pointed by 'array'
 * to a 3D frame. The contents of 'array' is for the position
 * of the original point and z-y-x Eulerian angle.
 * The result is put into 'f'.
 *
 * 'zFrame3DZYXToArray()' converts the frame 'f' to an array
 * of double-precision floating point value
 * [ x, y, z, azimuth, elevation, tilt ], of which the
 * first three values are for the position of the original
 * point and the last three are for z-y-x Eulerian angle.
 * The result is put into the buffer pointed by 'array'.
 *
 * 'zVec6DToFrame3DZYX()' converts the given 6D vector
 * 'v' to the frame. 'v' is a quasi vector for
 * [ x, y, z, azimuth, elevation, tilt ], which is for
 * the position of the original point and z-y-x Eulerian angle.
 * The result is put into 'f'.
 *
 * 'zFrame3DZYXToVec6D()' converts the given frame 'f'
 * to a 6D quasi vector which is for
 * [ x, y, z, azimuth, elevation, tilt ], of which the first
 * three values are for the position of the original point
 * and the last three are for z-y-x Eulerian angle.
 * The result is put into 'v'.
 *
 * 'zArrayToFrame3DZYZ()' converts the given array of
 * double-precision floating point values with
 * [ x, y, z, heading, pitch, bank ] pointed by 'array' to a 3D
 * frame. The components of 'array' is for the position
 * of the original point and z-y-z Eulerian angle.
 * The result is put into 'f'.
 *
 * 'zFrame3DZYZToArray()' converts the frame 'f'
 * to an array of double-precision floating point values
 * [ x, y, z, heading, pitch, bank ], of which the first three
 * values are for the position of the original point
 * and the last three are for z-y-z Eulerian angle.
 * The result is put into the buffer pointed by 'array'.
 *
 * 'zVec6DToFrame3DZYZ()' converts the given 6D vector
 * 'v' to the frame 'f'.
 * 'v' is a quasi vector for [ x, y, z, heading, pitch, bank ],
 * which is for the position of the original point and
 * z-y-z Eulerian angle.
 * The result is put into 'f'.
 *
 * 'zFrame3DZYZToVec6D()' converts the given frame 'f'
 * to a 6D quasi vector 'v', which is for
 * [ x, y, z, heading, pitch, bank ], of which the first
 * three values are for the position of the original point
 * and the last three are for z-y-z Eulerian angle.
 * The result is put into 'v'.
 * [RETURN VALUE]
 * 'zArrayToFrame3DZYX()', 'zVec6DToFrame3DZYX()',
 * 'zArrayToFrame3DZYZ()' and 'zVec6DToFrame3DZYZ()'
 * return a pointer to 'f'.
 *
 * 'zFrame3DZYXToArray()' and 'zFrame3DZYZToArray()'
 * return a pointer to 'array'.
 *
 * 'zFrame3DZYXToVec6D()' and 'zFrame3DZYZToVec6D()'
 * return a pointer to 'v'.
 * [NOTES]
 * For 'zFrame3DZYXToArray()' and 'zFrame3DZYZToArray()',
 * the buffer pointed by 'array' must have enough size
 * - array of more than six values.
 * If not, anything may happen.
 */
__EXPORT zFrame3D *zArrayToFrame3DZYX(double *array, zFrame3D *f);
__EXPORT double *zFrame3DToArrayZYX(zFrame3D *f, double *array);
__EXPORT zFrame3D *zVec6DToFrame3DZYX(zVec6D *v, zFrame3D *f);
__EXPORT zVec6D *zFrame3DToVec6DZYX(zFrame3D *f, zVec6D *v);

__EXPORT zFrame3D *zArrayToFrame3DZYZ(double *array, zFrame3D *f);
__EXPORT double *zFrame3DToArrayZYZ(zFrame3D *f, double *array);
__EXPORT zFrame3D *zVec6DToFrame3DZYZ(zVec6D *v, zFrame3D *f);
__EXPORT zVec6D *zFrame3DToVec6DZYZ(zFrame3D *f, zVec6D *v);

__EXPORT zFrame3D *zArrayToFrame3DAA(double *array, zFrame3D *f);
__EXPORT double *zFrame3DToArrayAA(zFrame3D *f, double *array);
__EXPORT zFrame3D *zVec6DToFrame3DAA(zVec6D *v, zFrame3D *f);
__EXPORT zVec6D *zFrame3DToVec6DAA(zFrame3D *f, zVec6D *v);

/* METHOD:
 * zFrame3DFRead, zFrame3DRead, zFrame3DDHFRead, zFrame3DDHRead,
 * zFrame3DFWrite, zFrame3DWrite
 * - input/output of 3D frame.
 *
 * 'zFrame3DFRead()' reads 12 values from the current
 * position of the file 'fp', and create the 3D frame
 * from them. The meaning of the sequencial value is for
 *  a11, a12, a13, x,
 *  a21, a22, a23, y,
 *  a31, a32, a33, z
 * The left 3x3 values are for the 3D attitude matrix, or
 * each { a11, a21, a31 }, { a12, a22, a32 } and { a13, a23, a33 }
 * means bases of the frame, while the right 3x1 values
 * are for the position of the original point { x, y, z }.
 * The result is put into 'f'.
 * 'zFrame3DRead()' simply reads values from the standard
 * input. The result is put into 'f'.
 *
 * 'zFrame3DDHFRead()' reads four values from the current
 * position of the file 'fp', and create 3D frame from them.
 * The meaning of the sequencial value is for a DH parameters
 * ( a, alpha, d, theta ).
 * The result is put into 'f'.
 * The unit of both 'alpha' and 'theta' is degree.
 * 'zFrame3DDHRead()' simply reads four values from the
 * standard input. The result is put into 'f'.
 *
 * 'zFrame3DFWrite()' writes the given 3D frame 'f' to
 * the current position of the file 'fp' in the following style.
 *  {
 *   a11, a12, a13, x,
 *   a21, a22, a23, y
 *   a31, a32, a33, z
 *  }
 * When the NULL pointer is given, it writes the following string.
 *  (null 3D frame)
 * 'zFrame3DWrite()' simply writes the 3D frame 'f' to
 * the standard output.
 * [RETURN VALUE]
 * 'zFrame3DFRead()', 'zFrame3DRead()',
 * 'zFrame3DDHFRead()' and 'zFrame3DDHRead()' return
 * a pointer to 'f'.
 *
 * 'zFrame3DFWrite()' and 'zFrame3DWrite()' return
 * no value.
 */
__EXPORT zFrame3D *zFrame3DFRead(FILE *fp, zFrame3D *f);
#define zFrame3DRead(f) zFrame3DFRead( stdin, (f) )
__EXPORT zFrame3D *zFrame3DDHFRead(FILE *fp, zFrame3D *f);
#define zFrame3DDHRead(f) zFrame3DDHFRead( stdin, (f) )
__EXPORT void zFrame3DFWrite(FILE *fp, zFrame3D *f);
#define zFrame3DWrite(f) zFrame3DFWrite( stdout, (f) )

/* METHOD:
 * zFrame3DFWriteXML - xml output.
 * ... yet testing.
 */
__EXPORT void zFrame3DFWriteXML(FILE *fp, zFrame3D *f);

__END_DECLS

#endif /* __ZEO_FRAME_H__ */
