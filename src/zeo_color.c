/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_color - color elements
 */

#include <zeo/zeo_color.h>
#include <math.h>

/* ********************************************************** */
/* CLASS: zRGB
 * RGB class - expression of color with RGB intensity set.
 * ********************************************************** */

static zRGB *_zRGBRatio(zRGB *rgb, char *ratio);
static zRGB *_zRGBHex(zRGB *rgb, char *hex);

/* OBJECT:
 * zblackrgb, zwhitergb
 * - black & white RGB set.
 */
const zRGB zblackrgb = { 0, 0, 0 }, zwhitergb = { 1, 1, 1 };

/* zRGBSet
 * - set RGB parameters.
 */
zRGB *zRGBSet(zRGB *rgb, float red, float green, float blue)
{
  zRedSet(   rgb, zLimit( red,   0, 1 ) );
  zGreenSet( rgb, zLimit( green, 0, 1 ) );
  zBlueSet(  rgb, zLimit( blue,  0, 1 ) );
  return rgb;
}

/* zRGBMul
 * - multiply a set of RGB parameters by another.
 */
zRGB *zRGBMul(zRGB *rgb1, zRGB *rgb2, zRGB *rgb)
{
  return zRGBSet( rgb, zRed(rgb1)*zRed(rgb2), zGreen(rgb1)*zGreen(rgb2), zBlue(rgb1)*zBlue(rgb2) );
}

/* (static)
 * _zRGBRatio
 * - decode a string of sequential floating-point values to RGB.
 */
zRGB *_zRGBRatio(zRGB *rgb, char *ratio)
{
  float r, g, b;

  sscanf( ratio, "%f:%f:%f", &r, &g, &b );
  return zRGBSet( rgb, r, g, b );
}

/* (static)
 * _zRGBHex
 * - decode a hexadecimal string to RGB.
 */
zRGB *_zRGBHex(zRGB *rgb, char *hex)
{
  int len, i;
  uint r, g, b, d;

  len = strlen( hex );
  if( len % 3 != 0 )
    ZRUNWARN( ZEO_ERR_RGB, hex, len );
  len /= 3;
  for( r=g=b=d=0, i=0; i<len; i++ ){
    r <<= 4; r |= atox_c( hex[i] );
    g <<= 4; g |= atox_c( hex[i+len] );
    b <<= 4; b |= atox_c( hex[i+len*2] );
    d <<= 4; d |= 0xf;
  }
  return zRGBSet( rgb, (double)r/d, (double)g/d, (double)b/d );
}

/* zRGBDec
 * - decode a string to RGB.
 */
zRGB *zRGBDec(zRGB *rgb, char *str)
{
  return str[0] == '#' ? _zRGBHex( rgb, str+1 ) : _zRGBRatio( rgb, str );
}

/* zRGBFRead
 * - input RGB parameter set.
 */
zRGB *zRGBFRead(FILE *fp, zRGB *rgb)
{
  double r, g, b;

  r = zFDouble(fp);
  g = zFDouble(fp);
  b = zFDouble(fp);
  return zRGBSet( rgb, r, g, b );
}

/* zRGBFWrite
 * - output RGB parameter set.
 */
void zRGBFWrite(FILE *fp, zRGB *rgb)
{
  fprintf( fp, "%.10g:%.10g:%.10g\n", zRed(rgb), zGreen(rgb), zBlue(rgb) );
}

/* zRGBFWriteXML
 * - xml output.
 */
void zRGBFWriteXML(FILE *fp, zRGB *rgb)
{
  fprintf( fp, "\"%.10g %.10g %.10g\"", zRed(rgb), zGreen(rgb), zBlue(rgb) );
}

/* ********************************************************** */
/* CLASS: zHSV
 * HSV class - expression of color by hue, saturation and value.
 * ********************************************************** */

/* zRGB2HSV
 * - convert RGB to HSV.
 */
zHSV *zRGB2HSV(zRGB *rgb, zHSV *hsv)
{
  float max, min, d, phase, phase0;

  /* find principal primary color */
  if( zRed(rgb) > zGreen(rgb) ){
    if( zRed(rgb) > zBlue(rgb) ){
      max = zRed(rgb);
      min = zMin( zGreen(rgb), zBlue(rgb) );
      phase = zGreen(rgb) - zBlue(rgb);
      phase0 = 0;
    } else{
      max = zBlue(rgb);
      min = zGreen(rgb);
      phase = zRed(rgb) - zGreen(rgb);
      phase0 = 240;
    }
  } else{
    if( zGreen(rgb) > zBlue(rgb) ){
      max = zGreen(rgb);
      min = zMin( zRed(rgb), zBlue(rgb) );
      phase = zBlue(rgb) - zRed(rgb);
      phase0 = 120;
    } else{
      max = zBlue(rgb);
      min = zRed(rgb);
      phase = zRed(rgb) - zGreen(rgb);
      phase0 = 240;
    }
  }
  d = max - min;
  /* conversion */
  zHSVHue(hsv) =   d == 0 ? 0 : 60.0 * phase / d + phase0;
  if( zHSVHue(hsv) < 0 ) zHSVHue(hsv) += 360;
  zHSVSat(hsv) = max == 0 ? 0 : d / max;
  zHSVVal(hsv) = max;
  return hsv;
}

/* zHSV2RGB
 * - convert HSV to RGB.
 */
zRGB *zHSV2RGB(zHSV *hsv, zRGB *rgb)
{
  float h, f, p, q, t;
  int i;

  h = zHSVHue(hsv) / 60;
  i = (int)floor( h ) % 6;
  f = h - i;
  p = zHSVVal(hsv) * ( 1 - zHSVSat(hsv) );
  q = zHSVVal(hsv) * ( 1 - f * zHSVSat(hsv) );
  t = zHSVVal(hsv) * ( 1 - (1-f) * zHSVSat(hsv) );
  switch( i ){
  case 0: zRGBSet( rgb, zHSVVal(hsv), t, p ); break;
  case 1: zRGBSet( rgb, q, zHSVVal(hsv), p ); break;
  case 2: zRGBSet( rgb, p, zHSVVal(hsv), t ); break;
  case 3: zRGBSet( rgb, p, q, zHSVVal(hsv) ); break;
  case 4: zRGBSet( rgb, t, p, zHSVVal(hsv) ); break;
  case 5: zRGBSet( rgb, zHSVVal(hsv), p, q ); break;
  default: ;
  }
  return rgb;
}

/* zHSVFWrite
 * - output HSV parameter set.
 */
void zHSVFWrite(FILE *fp, zHSV *hsv)
{
  fprintf( fp, "%.10g:%.10g:%.10g\n", zHSVHue(hsv), zHSVSat(hsv), zHSVVal(hsv) );
}
