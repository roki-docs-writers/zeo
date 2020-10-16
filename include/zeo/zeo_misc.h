/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_misc - miscellanies.
 */

#ifndef __ZEO_MISC_H__
#define __ZEO_MISC_H__

#include <cure/cure.h>
#include <zm/zm.h>

#include <zeo/zeo_errmsg.h>

__BEGIN_DECLS

/* TYPE: zAxis
 * axis direction name
 */
typedef byte zAxis;
enum{
  zX=0, zY, zZ, zXA, zYA, zZA,
};

/* METHOD:
 * zAxisExpr, zAxisByStr - conversion between string and axis.
 * [SYNOPSIS]
 * char *zAxisExpr(zAxis axis);
 * zAxis zAxisByStr(char str[]);
 * [DESCRIPTION]
 * 'zAxisExpr()' returns a string for the name of 'axis',
 * while 'zAxisByStr()' returns an axis identifier for
 * a string 'str'. The correspondency between strings and
 * axis identifiers are as follows:
 *  zX   <-> "x"
 *  zY   <-> "y"
 *  zZ   <-> "z"
 *  zXA  <-> "tilt"
 *  zYA  <-> "elev"
 *  zZA  <-> "azim"
 * [RETURN VALUE]
 * 'zAxisExpr()' returns a pointer to the string which expresses
 * the name of 'axis'.
 * 'zAxisByStr()' returns the corresponding axis identifier.
 */
__EXPORT char *zAxisExpr(zAxis axis);
__EXPORT zAxis zAxisByStr(char str[]);

/* TYPE: zDir
 * direction
 */
typedef byte zDir;
enum{
  zNONE=0, zRIGHT, zLEFT, zFORWARD, zBACKWARD, zUP, zDOWN
};

/* METHOD:
 * zDirExpr - expression for the name of direction.
 * [SYNOPSIS]
 * char *zDirExpr(zDir dir);
 * [DESCRIPTION]
 * 'zDirExpr()' returns a string for the name of 'dir',
 * which is a type of direction in the followings.
 *  zNONE     -> "none"
 *  zRIGHT    -> "right"
 *  zLEFT     -> "left"
 *  zFORWARD  -> "forward"
 *  zBACKWARD -> "backward"
 *  zUP       -> "up"
 *  zDOWN     -> "down"
 * [RETURN VALUE]
 * 'zDirExpr()' returns a pointer to the string which
 * expresses the name of 'dir'.
 */
__EXPORT char *zDirExpr(zDir dir);

/* METHOD:
 * zDirRev - reverse direction.
 * [SYNOPSIS]
 * zDir zDirRev(zDir dir);
 * [DESCRIPTION]
 * 'zDirRev()' returns the reverse direction of 'dir', namely:
 * zLEFT for zRIGHT, zRIGHT for zLEFT,
 * zBACKWARD for zFORWARD, zFORWARD for zBACKWARD,
 * zDOWN for zUP and zUP for zDOWN.
 * When 'dir' is zNONE, zNONE is returned.
 * [RETURN VALUE]
 * See DESCRIPTION.
 */
__EXPORT zDir zDirRev(zDir dir);

__END_DECLS

#endif /* __ZEO_MISC_H__ */
