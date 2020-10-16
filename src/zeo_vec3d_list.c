/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_vec3d_list - 3D vector list.
 */

#include <zeo/zeo_mat3d.h>

/* ********************************************************** */
/* CLASS: zVec3DList
 * 3D vector list
 * ********************************************************** */

/* zVec3DListInsert
 * - insert 3D vector list cell.
 */
zVec3DListCell *zVec3DListInsert(zVec3DList *list, zVec3D *v, bool flag)
{
  zVec3DListCell *cell;

  if( !( cell = zAlloc( zVec3DListCell, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( flag ){
    if( !( cell->data = zAlloc( zVec3D, 1 ) ) ){
      ZALLOCERROR();
      zFree( cell );
      return NULL;
    }
    zVec3DCopy( v, cell->data );
  } else
    cell->data = v;
  zListInsertHead( list, cell );
  return cell;
}

/* zVec3DListCreate
 * - create vector list.
 */
zVec3DList *zVec3DListCreate(zVec3DList *list, zVec3D v[], int num, bool flag)
{
  zListInit( list );
  while( num-- > 0 )
    if( !zVec3DListInsert( list, v++, flag ) ){
      ZALLOCERROR();
      break;
    }
  return list;
}

/* zVec3DListDestroy
 * - destroy vector list.
 */
void zVec3DListDestroy(zVec3DList *list, bool flag)
{
  zVec3DListCell *cell;

  while( !zListIsEmpty( list ) ){
    zListDeleteHead( list, &cell );
    if( flag ) zFree( cell->data );
    zFree( cell );
  }
}

/* zVec3DListQuickSort
 * - a quick sort routine for vector list class.
 */
zListQuickSortDef( zVec3DList, zVec3DListCell )

/* zVec3DListFWrite
 * - output 3D vector list.
 */
void zVec3DListFWrite(FILE *fp, zVec3DList *list)
{
  zVec3DListCell *cp;

  fprintf( fp, "%d\n", zListNum(list) );
  zListForEach( list, cp )
    zVec3DFWrite( fp, cp->data );
}

/* zVec3DListDataFWrite
 * - output 3D vector list.
 */
void zVec3DListDataFWrite(FILE *fp, zVec3DList *list)
{
  zVec3DListCell *cp;

  zListForEach( list, cp )
    zVec3DDataNLFWrite( fp, cp->data );
}

/* ********************************************************** */
/* utilities for point clouds
 * ********************************************************** */

/* zVec3DBarycenterPL
 * - barycenter of vector cloud given by a list.
 */
zVec3D *zVec3DBarycenterPL(zVec3DList *vl, zVec3D *c)
{
  zVec3DListCell *vc;

  zVec3DClear( c );
  zListForEach( vl, vc )
    zVec3DAddDRC( c, vc->data );
  return zVec3DDivDRC( c, zListNum(vl) );
}

/* zVec3DBarycenter
 * - barycenter of vector cloud given by an array.
 */
zVec3D *zVec3DBarycenter(zVec3D v[], int num, zVec3D *c)
{
  register int i;

  zVec3DClear( c );
  for( i=0; i<num; i++ )
    zVec3DAddDRC( c, &v[i] );
  return zVec3DDivDRC( c, num );
}

/* zVec3DPCA_PL
 * - PCA against vector cloud given by a list.
 */
zVec3D *zVec3DPCA_PL(zVec3DList *vl, zVec3D evec[])
{
  zMat3D pv, vm;
  double eval[3];
  zVec3DListCell *vc;

  zMat3DClear( &vm );
  zListForEach( vl, vc ){
    zMat3DDyad( vc->data, vc->data, &pv );
    zMat3DAddDRC( &vm, &pv );
  }
  zMat3DSymEig( &vm, eval, evec );
  return evec;
}

/* zVec3DPCA
 * - PCA against vector cloud given by an array.
 */
zVec3D *zVec3DPCA(zVec3D v[], int num, zVec3D evec[])
{
  zMat3D pv, vm;
  double eval[3];
  register int i;

  zMat3DClear( &vm );
  for( i=0; i<num; i++ ){
    zMat3DDyad( &v[i], &v[i], &pv );
    zMat3DAddDRC( &vm, &pv );
  }
  zMat3DSymEig( &vm, eval, evec );
  return evec;
}

/* zVec3DBaryPCA_PL
 * - barycenter of and PCA against vector cloud given by a list.
 */
zVec3D *zVec3DBaryPCA_PL(zVec3DList *vl, zVec3D *c, zVec3D evec[])
{
  zMat3D pv, vm;
  double eval[3];
  zVec3DListCell *vc;
  zVec3D dp;

  zVec3DBarycenterPL( vl, c );
  zMat3DClear( &vm );
  zListForEach( vl, vc ){
    zVec3DSub( vc->data, c, &dp );
    zMat3DDyad( &dp, &dp, &pv );
    zMat3DAddDRC( &vm, &pv );
  }
  zMat3DSymEig( &vm, eval, evec );
  return c;
}

/* zVec3DBaryPCA
 * - barycenter of and PCA against vector cloud given by an array.
 */
zVec3D *zVec3DBaryPCA(zVec3D v[], int num, zVec3D *c, zVec3D evec[])
{
  zMat3D pv, vm;
  double eval[3];
  zVec3D dp;
  register int i;

  zVec3DBarycenter( v, num, c );
  zMat3DClear( &vm );
  for( i=0; i<num; i++ ){
    zVec3DSub( &v[i], c, &dp );
    zMat3DDyad( &dp, &dp, &pv );
    zMat3DAddDRC( &vm, &pv );
  }
  zMat3DSymEig( &vm, eval, evec );
  return c;
}

/* zVec3DListNN
 * - a naive algorithm to find the nearest neighbor in a 3D vector list.
 */
zVec3D *zVec3DListNN(zVec3DList *list, zVec3D *v, double *dmin)
{
  zVec3DListCell *cell;
  double d;
  zVec3D *nn = NULL;

  *dmin = HUGE_VAL;
  zListForEach( list, cell )
    if( ( d = zVec3DSqrDist( cell->data, v ) ) < *dmin ){
      *dmin = d;
      nn = cell->data;
    }
  *dmin = sqrt( *dmin );
  return nn;
}
