/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_col_ph - collision checking: polyhedron.
 */

#include <zeo/zeo_col.h>

/* zColChkPH3D
 * - check if two polyhedra intersect with each other.
 */
bool zColChkPH3D(zPH3D *ph1, zPH3D *ph2, zVec3D *p1, zVec3D *p2)
{
  zVec3D _p1, _p2;

  if( p1 == NULL ) p1 = &_p1;
  if( p2 == NULL ) p2 = &_p2;
  return zGJK( zPH3DVertBuf(ph1), zPH3DVertNum(ph1), zPH3DVertBuf(ph2), zPH3DVertNum(ph2), p1, p2 );
}

/* Muller-Preparata's algorithm */

static zVec3D *_zTri3DDualXfer(zTri3D *t, zVec3D *p, zVec3D *dx);
static zVec3D *_zTri3DDualXfer_a(zTri3D *t, zVec3D *c, zVec3D *dx);
static zVec3D *_zTri3DDualXfer_b(zTri3D *t, zVec3D *dx);
static zPH3D *_zIntersectPH3D(zPH3D *ph1, zPH3D *ph2, zPH3D *phcol, zAABox3D *ib);

/* (static)
 * _zTri3DDualXfer
 * - dual transfer from/to point to/from plane.
 */
zVec3D *_zTri3DDualXfer(zTri3D *t, zVec3D *p, zVec3D *dx)
{
  return zVec3DDiv( zTri3DNorm(t), zVec3DInnerProd(zTri3DNorm(t),p), dx );
}

/* (static)
 * _zTri3DDualXfer_a
 * - dual transfer from/to point to/from plane shifting origin.
 */
zVec3D *_zTri3DDualXfer_a(zTri3D *t, zVec3D *c, zVec3D *dx)
{
  zVec3D p;

  zVec3DSub( zTri3DVert(t,0), c, &p );
  return _zTri3DDualXfer( t, &p, dx );
}

/* (static)
 * _zTri3DDualXfer_b
 * - simple dual transfer from/to point to/from plane.
 */
zVec3D *_zTri3DDualXfer_b(zTri3D *t, zVec3D *dx)
{
  return _zTri3DDualXfer( t, zTri3DVert(t,0), dx );
}

/* (static)
 * _zIntersectPH3D
 * - intersection of convices by Muller-Preparata's algorithm.
 */
zPH3D *_zIntersectPH3D(zPH3D *ph1, zPH3D *ph2, zPH3D *phcol, zAABox3D *ib)
{
  zVec3D *v, p1, p2, err;
  zTri3D *tri;
  zPH3D ch;
  register int i, n;

  zPH3DInit( phcol );
  /* compute proximity pair */
  if( !zColChkPH3D( ph1, ph2, &p1, &p2 ) ||
      !zVec3DIsTiny( zVec3DSub( &p1, &p2, &err ) ) ) return NULL;
  /* transfer to dual space */
  n = zPH3DFaceNum(ph1) + zPH3DFaceNum(ph2);
  if( !( v = zAlloc( zVec3D, n ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  /* dual-transfer triangles intersecting with
     the roughly-estimated intersection volume */
  n = 0;
  if( ib ){
    for( i=0; i<zPH3DFaceNum(ph1); i++ )
      if( zColChkTriAABox3D( ( tri = zPH3DFace(ph1,i) ), ib ) )
        _zTri3DDualXfer_a( tri, &p1, &v[n++] );
    for( i=0; i<zPH3DFaceNum(ph2); i++ )
      if( zColChkTriAABox3D( ( tri = zPH3DFace(ph2,i) ), ib ) )
        _zTri3DDualXfer_a( tri, &p1, &v[n++] );
  } else{
    for( i=0; i<zPH3DFaceNum(ph1); i++ )
      _zTri3DDualXfer_a( zPH3DFace(ph1,i), &p1, &v[n++] );
    for( i=0; i<zPH3DFaceNum(ph2); i++ )
      _zTri3DDualXfer_a( zPH3DFace(ph2,i), &p1, &v[n++] );
  }
  /* convex hull in dual space */
  if( !zCH3D( &ch, v, n ) ){
    phcol = NULL;
    goto TERMINATE;
  }
  zFree( v );
  /* re-transfer to the original space */
  if( !( v = zAlloc( zVec3D, zPH3DFaceNum(&ch) ) ) ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  for( i=0; i<zPH3DFaceNum(&ch); i++ ){
    _zTri3DDualXfer_b( zPH3DFace(&ch,i), &v[i] );
    zVec3DAddDRC( &v[i], &p1 );
  }
  if( !zCH3D( phcol, v, zPH3DFaceNum(&ch) ) ) phcol = NULL;

 TERMINATE:
  zPH3DDestroy( &ch );
  zFree( v );
  return phcol;
}

/* zIntersectPH3D
 * - intersection of convices by Muller-Preparata's algorithm.
 */
zPH3D *zIntersectPH3D(zPH3D *ph1, zPH3D *ph2, zPH3D *phcol)
{
  return _zIntersectPH3D( ph1, ph2, phcol, NULL );
}

/* zIntersectPH3DFast
 * - intersection of convices by Muller-Preparata's algorithm
 *   with a focusing-acceleration.
 */
zPH3D *zIntersectPH3DFast(zPH3D *ph1, zPH3D *ph2, zPH3D *phcol)
{
  zAABox3D ib;

  /* rough estimation of intersection volume by an axis-aligned box */
  if( !zIntersectPH3DBox( ph1, ph2, &ib ) ) return NULL;
  return _zIntersectPH3D( ph1, ph2, phcol, &ib );
}
