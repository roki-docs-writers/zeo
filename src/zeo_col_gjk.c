/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_col_gjk - collision checking: Gilbert-Johnson-Keerthi's algorithm.
 */

#include <zeo/zeo_col.h>

typedef struct{
  bool sw_w; /* included in W (the smallest simplex) */
  bool sw_y; /* included in Y (the updated simplex) */
  zVec3D w;  /* support map of Minkowski's sum */
  zVec3D *p1; /* corresponding vertex on object 1 to the support map */
  zVec3D *p2; /* corresponding vertex on object 2 to the support map */
  double s;  /* linear sum coefficient */
} zGJKSlot;

typedef struct{
  int n;
  zGJKSlot slot[4];
} zGJKSimplex;

static void _zGJKSlotInit(zGJKSlot *slot);
#ifdef DEBUG
static void _zGJKSlotVertFWrite(FILE *fp, zGJKSlot *slot);
static void _zGJKSlotWrite(zGJKSlot *slot);
#endif

static void _zGJKSimplexInit(zGJKSimplex *s);
#ifdef DEBUG
static void _zGJKSimplexYWrite(zGJKSimplex *s);
static void _zGJKSimplexWWrite(zGJKSimplex *s);
static void _zGJKSimplexYVertFWrite(FILE *fp, zGJKSimplex *s);
static void _zGJKSimplexWVertFWrite(FILE *fp, zGJKSimplex *s);
#endif

static bool _zGJKSimplexCheckSlot(zGJKSimplex *s, zGJKSlot *slot);
static int _zGJKSimplexAddSlot(zGJKSimplex *s, zGJKSlot *slot);
static zVec3D *_zGJKSimplexClosest(zGJKSimplex *s, zVec3D *v);
static void _zGJKSimplexMinimize(zGJKSimplex *s);
static zVec3D *_zGJKSupportMapOne(zVec3D p[], int n, zVec3D *v);
static zVec3D *_zGJKSupportMap(zGJKSlot *s, zVec3D p1[], int n1, zVec3D p2[], int n2, zVec3D *v);
static zVec3D *_zGJKSupportMapOnePL(zVec3DList *pl, zVec3D *v);
static zVec3D *_zGJKSupportMapPL(zGJKSlot *s, zVec3DList *pl1, zVec3DList *pl2, zVec3D *v);
static void _zGJKPair(zGJKSimplex *s, zVec3D *c1, zVec3D *c2);

/* (static)
 * _zGJKSlotInit
 * - initialize GJK slot.
 */
void _zGJKSlotInit(zGJKSlot *slot)
{
  slot->sw_w = slot->sw_y = false;
  zVec3DClear( &slot->w );
  slot->p1 = slot->p2 = NULL;
  slot->s = 0;
}

#ifdef DEBUG
/* (static)
 * _zGJKSlotWrite
 * - output slot information.
 */
void _zGJKSlotWrite(zGJKSlot *slot)
{
  printf( " w: " ); zVec3DWrite( &slot->w );
  printf( " p1: " ); zVec3DWrite( slot->p1 );
  printf( " p2: " ); zVec3DWrite( slot->p2 );
  printf( " s = %g\n", slot->s );
}

/* (static)
 * _zGJKSlotVertFWrite
 * - output vertices of a slot.
 */
void _zGJKSlotVertFWrite(FILE *fp, zGJKSlot *slot)
{
  zVec3DDataFWrite( fp, slot->p1 );
  zVec3DDataFWrite( fp, slot->p2 );
}
#endif

/* (static)
 * _zGJKSimplexInit
 * - initialize GJK simplex set.
 */
void _zGJKSimplexInit(zGJKSimplex *s)
{
  s->n = 0;
  _zGJKSlotInit( &s->slot[0] );
  _zGJKSlotInit( &s->slot[1] );
  _zGJKSlotInit( &s->slot[2] );
  _zGJKSlotInit( &s->slot[3] );
}

#ifdef DEBUG
/* (static)
 * _zGJKSimplexYWrite
 * - output GJK simplex set Y (non-minimum simplex)
 */
void _zGJKSimplexYWrite(zGJKSimplex *s)
{
  register int i;

  printf( "[Y]\n" );
  for( i=0; i<4; i++ )
    if( s->slot[i].sw_y ){
      printf( "<slot %d>\n", i );
      _zGJKSlotWrite( &s->slot[i] );
    }
}

/* (static)
 * _zGJKSimplexWWrite
 * - output GJK simplex set W (minimum simplex)
 */
void _zGJKSimplexWWrite(zGJKSimplex *s)
{
  register int i;

  printf( "[W]\n" );
  for( i=0; i<4; i++ )
    if( s->slot[i].sw_w ){
      printf( "<slot %d>\n", i );
      _zGJKSlotWrite( &s->slot[i] );
    }
}

/* (static)
 * _zGJKSimplexYVertFWrite
 * - output vertices of GJK simplex set Y (non-minimum simplex)
 */
void _zGJKSimplexYVertFWrite(FILE *fp, zGJKSimplex *s)
{
  register int i;

  for( i=0; i<4; i++ )
    if( s->slot[i].sw_y )
      _zGJKSlotVertFWrite( fp, &s->slot[i] );
}

/* (static)
 * _zGJKSimplexWVertFWrite
 * - output vertices GJK simplex set W (minimum simplex)
 */
void _zGJKSimplexWVertFWrite(FILE *fp, zGJKSimplex *s)
{
  register int i;

  for( i=0; i<4; i++ )
    if( s->slot[i].sw_w )
      _zGJKSlotVertFWrite( fp, &s->slot[i] );
}
#endif

/* (static)
 * _zGJKSimplexCheckSlot
 * - check if the specified slot is already in the previous testing simplex.
 */
bool _zGJKSimplexCheckSlot(zGJKSimplex *s, zGJKSlot *slot)
{
  register int i;

  for( i=0; i<4; i++ )
    if( s->slot[i].sw_y && zVec3DMatch( &s->slot[i].w, &slot->w ) )
      return true;
  return false;
}

/* (static)
 * _zGJKSimplexAddSlot
 * - add a slot into the current testing simplex.
 */
int _zGJKSimplexAddSlot(zGJKSimplex *s, zGJKSlot *slot)
{
  register int i;

  for( i=0; i<4; i++ )
    s->slot[i].sw_y = s->slot[i].sw_w;
  for( i=0; i<4; i++ )
    if( !s->slot[i].sw_y ){
      s->slot[i].sw_y = true;
      zVec3DCopy( &slot->w, &s->slot[i].w );
      s->slot[i].p1 = slot->p1;
      s->slot[i].p2 = slot->p2;
      s->slot[i].s = slot->s;
      return i;
    }
  ZRUNERROR( ZEO_ERR_FATAL );
  return -1;
}

/* (static)
 * _zGJKSimplexClosest
 * - find the closest point in/on the testing simplex.
 */
zVec3D *_zGJKSimplexClosest(zGJKSimplex *s, zVec3D *v)
{
  double _q[9], _c[3];
  double _a[] = { -1.0, -1.0, -1.0 }, _b[] = { -1.0 }, _l[3];
  zMatStruct q, a;
  zVecStruct c, b, l;
  zVec3D dp[3];
  double s0, cost;
  register int i, j, n, n1;
  int index[4];

  /* create index */
  for( n=0, i=0; i<4; i++ )
    if( s->slot[i].sw_y ) index[n++] = i;
  /* setup & solve QP */
  n1 = n - 1;
  zMatBuf(&q) = _q; zMatSetSize(&q,n1,n1);
  zVecBuf(&c) = _c; zVecSetSize(&c,n1);
  zMatBuf(&a) = _a; zMatSetSize(&a,1,n1);
  zVecBuf(&b) = _b; zVecSetSize(&b,1);
  zVecBuf(&l) = _l; zVecSetSize(&l,n1);
  for( i=1; i<n; i++ )
    zVec3DSub( &s->slot[index[i]].w, &s->slot[index[0]].w, &dp[i-1] );
  for( i=0; i<n1; i++ ){
    for( j=0; j<n1; j++ )
      zMatElem(&q,i,j) = zVec3DInnerProd(&dp[i],&dp[j]);
    zVecElem(&c,i) = zVec3DInnerProd(&s->slot[index[0]].w,&dp[i]);
  }
  zQPSolveLemke( &q, &c, &a, &b, &l, &cost );
  /* solve QP */
  zVec3DCopy( &s->slot[index[0]].w, v );
  s0 = 0;
  for( i=0; i<n1; i++ ){
    zVec3DCatDRC( v, zVecElem(&l,i), &dp[i] );
    s0 += ( s->slot[index[i+1]].s = zVecElem(&l,i) );
  }
  s->slot[index[0]].s = 1.0 - s0;
  return v;
}

/* (static)
 * _zGJKSimplexMinimize
 * - minimize the testing simplex.
 */
void _zGJKSimplexMinimize(zGJKSimplex *s)
{
  register int i;
  zVec3D ws;

  for( s->n=0, i=0; i<4; i++ ){
    zVec3DMul( &s->slot[i].w, s->slot[i].s, &ws );
    if( zVec3DIsTiny( &ws ) ){
      s->slot[i].s = 0;
      s->slot[i].sw_w = false;
    } else{
      s->slot[i].sw_w = true;
      s->n++;
    }
  }
}

/* (static)
 * _zGJKSupportMapOne
 * - support map of a set of points with respect to a direction vector.
 */
zVec3D *_zGJKSupportMapOne(zVec3D p[], int n, zVec3D *v)
{
  register int i;
  double d, d_max;
  zVec3D *sp;

  if( n <= 0 ){
    ZRUNWARN( ZEO_ERR_EMPTYSET );
    return NULL;
  }
  sp = &p[0];
  d_max = zVec3DInnerProd( sp, v );
  for( i=1; i<n; i++ )
    if( ( d = zVec3DInnerProd( &p[i], v ) ) > d_max ){
      sp = &p[i];
      d_max = d;
    }
  return sp;
}

/* (static)
 * _zGJKSupportMap
 * - support map of Minkowski difference.
 */
zVec3D *_zGJKSupportMap(zGJKSlot *s, zVec3D p1[], int n1, zVec3D p2[], int n2, zVec3D *v)
{
  zVec3D nv;

  zVec3DRev( v, &nv );
  s->p1 = _zGJKSupportMapOne( p1, n1, &nv );
  s->p2 = _zGJKSupportMapOne( p2, n2,   v );
  zVec3DSub( s->p1, s->p2, &s->w );
  return &s->w;
}

/* (static)
 * _zGJKSupportMapOnePL
 * - support map of a set of points with respect to a direction vector.
 */
zVec3D *_zGJKSupportMapOnePL(zVec3DList *pl, zVec3D *v)
{
  double d, d_max;
  zVec3D *sp;
  zVec3DListCell *cp;

  if( zListIsEmpty(pl) ){
    ZRUNWARN( ZEO_ERR_EMPTYSET );
    return NULL;
  }
  sp = zListTail(pl)->data;
  d_max = zVec3DInnerProd( sp, v );
  zListForEach( pl, cp ){
    if( ( d = zVec3DInnerProd( cp->data, v ) ) > d_max ){
      sp = cp->data;
      d_max = d;
    }
  }
  return sp;
}

/* (static)
 * _zGJKSupportMapPL
 * - support map of Minkowski difference.
 */
zVec3D *_zGJKSupportMapPL(zGJKSlot *s, zVec3DList *pl1, zVec3DList *pl2, zVec3D *v)
{
  zVec3D nv;

  zVec3DRev( v, &nv );
  s->p1 = _zGJKSupportMapOnePL( pl1, &nv );
  s->p2 = _zGJKSupportMapOnePL( pl2,   v );
  zVec3DSub( s->p1, s->p2, &s->w );
  return &s->w;
}

/* (static)
 * _zGJKPair
 * - a pair of points on the original convex hulls.
 */
void _zGJKPair(zGJKSimplex *s, zVec3D *c1, zVec3D *c2)
{
  register int i;

  zVec3DClear( c1 );
  zVec3DClear( c2 );
  for( i=0; i<4; i++ )
    if( s->slot[i].sw_w ){
      zVec3DCatDRC( c1, s->slot[i].s, s->slot[i].p1 );
      zVec3DCatDRC( c2, s->slot[i].s, s->slot[i].p2 );
    }
}

/* zGJK
 * - Gilbert-Johnson-Keerthi algorithm.
 */
bool zGJK(zVec3D p1[], int n1, zVec3D p2[], int n2, zVec3D *c1, zVec3D *c2)
{
  zGJKSimplex s; /* simplex */
  zGJKSlot slot;
  zVec3D v; /* proximity */
  double dv2;

  zVec3DSub( &p1[0], &p2[0], &v );
  dv2 = zVec3DSqrNorm( &v );
  _zGJKSimplexInit( &s );
  do{
    _zGJKSupportMap( &slot, p1, n1, p2, n2, &v );
    if( _zGJKSimplexCheckSlot( &s, &slot ) ||
        dv2 - zVec3DInnerProd(&slot.w,&v) <= zTOL ){
      break; /* succeed */
    }
    _zGJKSimplexAddSlot( &s, &slot );
    _zGJKSimplexClosest( &s, &v );
    _zGJKSimplexMinimize( &s );
    dv2 = zVec3DSqrNorm( &v );
  } while( s.n < 4 );
  _zGJKPair( &s, c1, c2 );
  return s.n == 4 ? true : false;
}

/* zGJKPL
 * - Gilbert-Johnson-Keerthi algorithm.
 */
bool zGJKPL(zVec3DList *pl1, zVec3DList *pl2, zVec3D *c1, zVec3D *c2)
{
  zGJKSimplex s; /* simplex */
  zGJKSlot slot;
  zVec3D v; /* proximity */
  double dv2;

  zVec3DSub( zListTail(pl1)->data, zListTail(pl2)->data, &v );
  dv2 = zVec3DSqrNorm( &v );
  _zGJKSimplexInit( &s );
  do{
    _zGJKSupportMapPL( &slot, pl1, pl2, &v );
    if( _zGJKSimplexCheckSlot( &s, &slot ) ||
        dv2 - zVec3DInnerProd(&slot.w,&v) <= zTOL )
      break; /* succeed */
    _zGJKSimplexAddSlot( &s, &slot );
    _zGJKSimplexClosest( &s, &v );
    _zGJKSimplexMinimize( &s );
    dv2 = zVec3DSqrNorm( &v );
  } while( s.n < 4 );
  _zGJKPair( &s, c1, c2 );
  return s.n == 4 ? true : false;
}
