/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_color - color elements
 */

#ifndef __ZEO_COLOR_H__
#define __ZEO_COLOR_H__

#include <cure/cure.h>
#include <zeo/zeo_misc.h>

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zColor
 * color class - color space expression.
 * ********************************************************** */

typedef struct{
  float v1, v2, v3;
} zColor;

/* ********************************************************** */
/* CLASS: zRGB
 * RGB class - expression of color with RGB intensity set.
 * ********************************************************** */

typedef zColor zRGB;

#define zRed(c)   (c)->v1
#define zGreen(c) (c)->v2
#define zBlue(c)  (c)->v3

#define zRedSet(c,r)   ( zRed(c) = (r) )
#define zGreenSet(c,g) ( zGreen(c) = (g) )
#define zBlueSet(c,b)  ( zBlue(c) = (b) )

/* OBJECT:
 * zblackrgb, zwhitergb
 * - black & white RGB set.
 */
extern const zRGB zblackrgb, zwhitergb;
#define ZBLACKRGB ( (zRGB *)&zblackrgb )
#define ZWHITERGB ( (zRGB *)&zwhitergb )

/* METHOD:
 * zRGBSet, zRGBCopy, zRGBGet
 * - creation, initialization and copy of RGB set.
 *
 * 'zRGBSet()' sets 'red', 'green' and 'blue' to the RGB
 * parameter set 'rgb'.
 * #
 * 'zRGBCopy()' copies the RGB parameter set pointed by 'src'
 * to the other set pointed by 'dest'. Actually, 'zRGBInit()'
 * is defined as a macro(see zoptic.h).
 * #
 * 'zRGBGet()' abstracts red, green and blue values of 'rgb'
 * and puts them into 'red', 'green' and 'blue', respectively.
 * [RETURN VALUE]
 * 'zRGBSet()' and 'zRGBInit()' returns a pointer to 'rgb'.
 * #
 * 'zRGBCopy()' returns no value.
 * #
 * 'zRGBGet()' returns no value.
 */
__EXPORT zRGB *zRGBSet(zRGB *rgb, float red, float green, float blue);
#define zRGBCopy(s,d) ( *(d) = *(s) )
#define zRGBGet(c,r,g,b) do{\
  *(r) = zRed(c);\
  *(g) = zGreen(c);\
  *(b) = zBlue(c);\
} while(0)

/* grayscale mode. */
#define zGS(c)       ( ( zRed(c) + zGreen(c) + zBlue(c) ) / 3.0 )
#define zGSSet(c,i)  ( zRed(c) = zGreen(c) = zBlue(c) = (i) )

__EXPORT zRGB *zRGBMul(zRGB *rgb1, zRGB *rgb2, zRGB *rgb);

/* zRGBDec - decode a string to RGB.
 *
 * 'zRGBDec()' decodes a string to RGB. If the string begins
 * from '#', the string is decoded as hexadecimal expression.
 * Otherwise, it is decoded as floating-point value expression.
 * #
 * For hexadecimal expression, the string has to be a sequence
 * of 0-9 or a-f/A-F, and will be segmented into three values
 * with the same length. Consequently, the length of the string
 * should be a multiple of three.
 * Ex. "0f3df8" will be a vivid bluish color.
 * #
 * For floating-point value expression, the string has to be
 * a sequence of three floating-point values, each of which
 * should be in the range of 0-1.
 * Ex. "0.1:0.3:0.8" will be a smoky bluish color.
 * #
 * The resultant RGB is set where 'rgb' points.
 * [RETURN VALUE]
 * 'zRGBDec()' returns a pointer 'rgb'.
 */
__EXPORT zRGB *zRGBDec(zRGB *rgb, char *str);

/* METHOD:
 * zRGBFRead, zRGBRead, zRGBFWrite, zRGBWrite
 * - input/output of a RGB parameter set.
 *
 * 'zRGBFRead()' reads the sequence of 3 values as the
 * information about the RGB parameter set from the current
 * position of the file 'fp', and copies them to 'rgb'.
 * #
 * 'zRGBFRead()' reads an information about the RGB parameter
 * set simply from the standard in and copies them to 'rgb'.
 * #
 * 'zRGBFWrite()' writes the RGB parameter set of 'rgb' to the
 * current position of the file 'fp'.
 * #
 * 'zRGBWrite()' writes the RGB parameter set of 'rgb' simply
 * to the standard out.
 * [RETURN VALUE]
 * 'zRGBFRead()' and 'zRGBRead()' return a pointer to 'rgb'.
 * #
 * 'zRGBFWrite()' and 'zRGBWrite()' return no value.
 */
__EXPORT zRGB *zRGBFRead(FILE *fp, zRGB *rgb);
#define zRGBRead(c)  zRGBFRead( stdin, (c) )
__EXPORT void zRGBFWrite(FILE *fp, zRGB *rgb);
#define zRGBWrite(c) zRGBFWrite( stdout, (c) )

/* METHOD:
 * zRGBFWriteXML - xml output.
 * ... yet testing.
 */
__EXPORT void zRGBFWriteXML(FILE *fp, zRGB *rgb);

/* ********************************************************** */
/* CLASS: zHSV
 * HSV class - expression of color by hue, saturation and value.
 * ********************************************************** */

typedef zColor zHSV;

#define zHSVHue(c) (c)->v1 /* hue */
#define zHSVSat(c) (c)->v2 /* saturation */
#define zHSVVal(c) (c)->v3 /* value */

__EXPORT zHSV *zRGB2HSV(zRGB *rgb, zHSV *hsv);
__EXPORT zRGB *zHSV2RGB(zHSV *hsv, zRGB *rgb);
__EXPORT void zHSVFWrite(FILE *fp, zHSV *hsv);
#define zHSVWrite(hsv) zHSVFWrite( stdout, hsv )

__END_DECLS

#endif /* __ZEO_COLOR_H__ */
