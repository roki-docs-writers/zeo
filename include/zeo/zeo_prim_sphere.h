/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_prim_sphere - primitive 3D shapes: sphere.
 */

#ifndef __ZEO_PRIM_SPHERE_H__
#define __ZEO_PRIM_SPHERE_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zSphere3D
 * 3D sphere class
 * ********************************************************** */

typedef struct{
  zVec3D center;
  double radius;
  int div;
} zSphere3D;

#define zSphere3DCenter(s)      ( &(s)->center )
#define zSphere3DRadius(s)      (s)->radius
#define zSphere3DDiv(s)         (s)->div

#define zSphere3DSetCenter(s,c) zVec3DCopy( c, zSphere3DCenter(s) )
#define zSphere3DSetRadius(s,r) ( zSphere3DRadius(s) = (r) )
#define zSphere3DSetDiv(s,d)    ( zSphere3DDiv(s) = (d) )

/* METHOD:
 * zSphere3DInit, zSphere3DCreate, zSphere3DCopy
 * - initialization, creation and copy of 3D sphere.
 *
 * 'zSphere3DInit()' initializes a 3D sphere 'sphere',
 * setting its center for the original point and radius
 * for zero.
 * #
 * 'zSphere3DCreate()' creates a 3D sphere whose center
 * is 'c' and radius is 'r'.
 * \a div is the number of division for polyhedral approximation.
 * When zero is given for \a div, ZEO_PRIM_DEFAULT_DIV is
 * set instead.
 * #
 * 'zSphere3DCopy()' copies a 3D sphere 'src' to the
 * other 'dest'.
 * [RETURN VALUE]
 * Each of 'zSphere3DInit()' and 'zSphere3DCreate()'
 * returns a pointer 'sphere'.
 * 'zSphere3DCopy()' returns a pointer the copied 'dest'.
 */
__EXPORT zSphere3D *zSphere3DCreate(zSphere3D *sphere, zVec3D *c, double r, int div);
__EXPORT zSphere3D *zSphere3DInit(zSphere3D *sphere);
__EXPORT zSphere3D *zSphere3DCopy(zSphere3D *src, zSphere3D *dest);
__EXPORT zSphere3D *zSphere3DMirror(zSphere3D *src, zSphere3D *dest, zAxis axis);

/* METHOD:
 * zSphere3DXfer, zSphere3DXferInv
 * - transformation of 3D sphere.
 *
 * 'zSphere3DXfer()' transforms a 3D sphere 'src' by
 * a transformation frame 'f', put it into 'dest'.
 * #
 * 'zSphere3DXferInv()' transforms 'src' by the
 * inverse transformation frame of 'f', put it into 'dest'.
 * [RETURN VALUE]
 * Each of 'zSphere3DXfer()' and 'zSphere3DXferInv()'
 * returns a pointer 'dest'.
 */
__EXPORT zSphere3D *zSphere3DXfer(zSphere3D *src, zFrame3D *f, zSphere3D *dest);
__EXPORT zSphere3D *zSphere3DXferInv(zSphere3D *src, zFrame3D *f, zSphere3D *dest);

/* METHOD:
 * zSphere3DClosest, zSphere3DPointDist, zSphere3DPointIsInside
 * - distance from 3D point to 3D sphere.
 *
 * 'zSphere3DClosest()' calculates the closest point
 * from a 3D point 'p' to a 3D sphere 'sphere', and put it
 * into 'cp'. When 'p' is inside of 'sphere', it just
 * copies 'p' to 'cp'.
 * #
 * 'zSphere3DPointDist()' calculates the distance from
 * a 3D point 'p' to a 3D sphere 'sphere'.
 * #
 * 'zSphere3DPointIsInside()' checks if a 3D point'p' is
 * inside of a 3D sphere 'sphere'. The point on the
 * surface of 'sphere' is judged to be inside of 'sphere'
 * if the true value is given for 'rim'.
 * [RETURN VALUE]
 * Each of 'zSphere3DClosest()' and 'zSphere3DPointDist()'
 * returns the signed distance from 'p' to 'sphere'.
 * The result is
 *  - a positive value when 'p' is outside of 'sphere', or
 *  - a negative value when 'p' is inside of 'sphere'.
 * #
 * 'zSphere3DPointIsInside()' returns the true value if
 * 'p' is inside of 'sphere', or the false value otherwise.
 */
__EXPORT double zSphere3DClosest(zSphere3D *sphere, zVec3D *p, zVec3D *cp);
__EXPORT double zSphere3DPointDist(zSphere3D *sphere, zVec3D *p);
__EXPORT bool zSphere3DPointIsInside(zSphere3D *sphere, zVec3D *p, bool rim);

/* METHOD:
 * zSphere3DVolume, zSphere3DInertia
 * - volume and inertia of 3D sphere.
 *
 * 'zSphere3DVolume()' calculates the volume of a
 * 3D sphere 'sphere'.
 * #
 * 'zSphere3DInertia()' calculates the inertia tensor
 * of a 3D sphere 'sphere' around its center. It treats
 * 'sphere' as a solid model. The result is put into 'inertia'.
 * [RETURN VALUE]
 * 'zSphere3DVolume()' returns the volume calculated.
 * 'zSphere3DInertia()' returns a pointer 'inertia'.
 */
__EXPORT double zSphere3DVolume(zSphere3D *sphere);
__EXPORT zMat3D *zSphere3DInertia(zSphere3D *sphere, zMat3D *inertia);

/* METHOD:
 * zSphere3DToPH - convert sphere to polyhedron.
 *
 * 'zSphere3DToPH()' converts a sphere 'sphere' to
 * a polyhedron 'ph' as a polygon model.
 * 'ph' should be initialized in advance.
 * It approximately divides the face into rectangles
 * by the stored number of division.
 * [RETURN VALUE]
 * zSphere3DToPH() returns a pointer \a ph.
 * [SEE ALSO]
 * zCyl3DToPH, zCyl3DToPHDRC,
 * zCone3DToPH, zCone3DToPHDRC
 */
/* default longitudinal & latitudinal division number are the same. */
__EXPORT zPH3D *zSphere3DToPH(zSphere3D *sphere, zPH3D *ph);

/* METHOD:
 * zSphere3DFRead, zSphere3DRead, zSphere3DFWrite, zSphere3DWrite,
 * - input/output of 3D sphere.
 *
 * 'zSphere3DFRead()' reads the information of a 3D sphere
 * from the current position of the file 'fp', and creates
 * the new sphere 'sphere'.
 * An acceptable data file format is as follows.
 * #
 *  center: <x> <y> <z>
 *  radius: <r>
 *  div: <div>
 * #
 * Each bracketed value must be substituted for a real number.
 * 'zSphere3DRead()' reads the information for 'sphere'
 * simply from the standard input.
 * The field div is skippable.
 * #
 * 'zSphere3DFWrite()' writes the information of 'sphere' to
 * the current position of the file 'fp' in the same format
 * with the above. 'zSphere3DWrite()' writes the information
 * of 'sphere' simply to the standard out.
 * [RETURN VALUE]
 * Each of 'zSphere3DFRead()' and 'zSphere3DRead()' returns
 * a pointer 'sphere'.
 * #
 * Neither 'zSphere3DFWrite()' nor 'zSphere3DWrite()' returns
 * any values.
 */
__EXPORT zSphere3D *zSphere3DFRead(FILE *fp, zSphere3D *sphere);
#define zSphere3DRead(s) zSphere3DFRead( stdin, (s) )
__EXPORT void zSphere3DFWrite(FILE *fp, zSphere3D *sphere);
#define zSphere3DWrite(s) zSphere3DFWrite( stdout, (s) )

/* methods for abstraction */
extern zPrimCom zprim_sphere3d_com;

__END_DECLS

#endif /* __ZEO_PRIM_SPHERE_H__ */
