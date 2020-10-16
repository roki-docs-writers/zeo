/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_optic - optical properties
 */

#ifndef __ZEO_OPTIC_H__
#define __ZEO_OPTIC_H__

#include <cure/cure.h>
#include <zeo/zeo_color.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zOpticalInfo
 * class for the information about optical characterization
 * parameters
 * ********************************************************** */

typedef struct{
  Z_NAMED_CLASS
  zRGB amb;     /* coefficients of reflection for ambient */
  zRGB dif;     /* coefficients of diffuse reflection */
  zRGB spc;     /* coefficients of specular reflection */
  double ns;    /* Phong's exponential for specular reflection */
  double sns;   /* shininess */
  double alpha; /* alpha value */
} zOpticalInfo;

#define zOpticalInfoAmb(o)   ( &(o)->amb )
#define zOpticalInfoDif(o)   ( &(o)->dif )
#define zOpticalInfoSpc(o)   ( &(o)->spc )
#define zOpticalInfoExp(o)   (o)->ns
#define zOpticalInfoSns(o)   (o)->sns
#define zOpticalInfoAlpha(o) (o)->alpha

#define zOpticalInfoSetAmb(o,a)   zRGBCopy(a,zOpticalInfoAmb(o))
#define zOpticalInfoSetDif(o,d)   zRGBCopy(d,zOpticalInfoDif(o))
#define zOpticalInfoSetSpc(o,s)   zRGBCopy(s,zOpticalInfoSpc(o))
#define zOpticalInfoSetExp(o,n)   ( zOpticalInfoExp(o) = (n) )
#define zOpticalInfoSetSns(o,s)   ( zOpticalInfoSns(o) = (s) )
#define zOpticalInfoSetAlpha(o,a) ( zOpticalInfoAlpha(o) = (a) )

/* METHOD:
 * zOpticalInfoCreate, zOpticalInfoInit,
 * zOpticalInfoCopy, zOpticalInfoDestroy
 * - creation, initialization, copy and destruction of optical
 * parameter set.
 *
 * 'zOpticalInfoCreate()' creates an optical parameter set 'oi'
 * from the arguments. The meanings of each argument are as follows -
 *  'amb'  : coefficients of reflection for ambient.
 *  'dif'  : coefficients of diffuse reflection.
 *  'spc'  : coefficients of specular reflection.
 *  'ns'   : Phong's specular exponential value.
 *  'sns'  : shininess.
 *  'alpha': alpha value.
 *  'name' : a name of the optical set.
 * #
 * 'zOpticalInfoInit()' initializes 'oi'. All parameters will
 * be set for 1.0. Actually, 'zOpticalInfoInit()' is defined as
 * a macro(see zoptic.h).
 * #
 * 'zOpticalInfoCopy()' copies parameters of 'src' to 'dest'.
 * #
 * 'zOpticalInfoDestroy()' destroys the instance 'oi'.
 * [RETURN VALUE]
 * 'zOpticalInfoCreate()' and 'zOpticalInfoInit()' return a
 * pointer to 'oi'.
 * #
 * 'zOpticalInfoCopy()' returns a pointer to 'dest'.
 * #
 * 'zOpticalInfoDestroy()' returns no value.
 */
__EXPORT zOpticalInfo *zOpticalInfoCreate(zOpticalInfo *oi, float ar, float ag, float ab, float dr, float dg, float db, float sr, float sg, float sb, double ns, double sns, double alpha, char *name);
#define zOpticalInfoCreateSimple(o,r,g,b,n) \
  zOpticalInfoCreate( (o), 0.5*r, 0.5*g, 0.5*b, r, g, b, 0, 0, 0, 0, 0, 1, (n) )
#define zOpticalInfoInit(o) \
  zOpticalInfoCreate( (o), 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, NULL )
__EXPORT zOpticalInfo *zOpticalInfoCopy(zOpticalInfo *src, zOpticalInfo *dest);
__EXPORT zOpticalInfo *zOpticalInfoClone(zOpticalInfo *src, zOpticalInfo *dest);
#define zOpticalInfoDestroy(o) do{\
  zNameDestroy(o);\
  zOpticalInfoInit(o);\
} while(0)

__EXPORT zOpticalInfo *zOpticalInfoMul(zOpticalInfo *oi1, zOpticalInfo *oi2, zOpticalInfo *oi);

/* tag to identify optical info. */
#define ZOPTIC_TAG "optic"

/* METHOD:
 * zOpticalInfoFRead, zOpticalInfoRead,
 * zOpticalInfoFWrite, zOpticalInfoWrite
 * - input/output of the information of the optical parameter set.
 *
 * 'zOpticalInfoFRead()' reads an information about the optical
 * parameter set from the current position of the file 'fp',
 * and copies them to 'oi'.
 * An acceptable data file format is as follows.
 * #
 *  name     <string> <- the name of optical parameter set
 *  ambient  <value> <value> <value>
 *    ^ coefficients of reflection for ambient
 *  diffuse  <value> <value> <value>
 *    ^ coefficients of diffuse reflection
 *  specular <value> <value> <value>
 *    ^ coefficients of specular reflection
 *  exp      <value>
 *    ^ Phong s specular exponential value
 *  alpha <value>
 *    ^ alpha value
 * #
 * 'zOpticalInfoRead()' reads an information about the optical
 * parameter set simply from the standard in and copies them
 * to 'oi'.
 * #
 * 'zOpticalInfoFWrite()' writes an information of 'oi' to the
 * current position of the file 'fp'.
 * #
 * 'zOpticalInfoWrite()' writes an information of 'oi' simply
 * to the standard out.
 * [RETURN VALUE]
 * 'zOpticalInfoFRead()' and 'zOpticalInfoRead()' return a
 * pointer to 'oi'.
 * #
 * 'zOpticalInfoFWrite()' and 'zOpticalInfoWrite()' return
 * no value.
 */
__EXPORT zOpticalInfo *zOpticalInfoFRead(FILE *fp, zOpticalInfo *oi);
#define zOpticalInfoRead(i) zOpticalInfoFRead( stdin, (i) )
__EXPORT void zOpticalInfoFWrite(FILE *fp, zOpticalInfo *oi);
#define zOpticalInfoWrite(i) zOpticalInfoFWrite( stdout, (i) )

/* METHOD:
 * zOpticalInfoFWriteXML - xml output.
 * ... yet testing.
 */
__EXPORT void zOpticalInfoFWriteXML(FILE *fp, zOpticalInfo *oi, int indent);

__END_DECLS

#endif /* __ZEO_OPTIC_H__ */
