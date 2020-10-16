/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_mat3d - 3x3 matrix.
 */

#include <zeo/zeo_mat3d.h>

/* ********************************************************** */
/* CLASS: zMat3D
 * 3D matrix class
 * ********************************************************** */

/* OBJECT:
 * zmat3Dzero, zmat3Dident
 * - 3D zero matrix and identity matrix.
 */
const zMat3D zmat3Dzero  = { { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } };
const zMat3D zmat3Dident = { { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } };

/* zMat3DCreate
 * - create 3D matrix.
 */
zMat3D *zMat3DCreate(zMat3D *m,
  double a11, double a12, double a13,
  double a21, double a22, double a23,
  double a31, double a32, double a33)
{
  zVec3DCreate( &m->v[0], a11, a21, a31 );
  zVec3DCreate( &m->v[1], a12, a22, a32 );
  zVec3DCreate( &m->v[2], a13, a23, a33 );
  return m;
}

/* zMat3DT
 * - transpose a 3D matrix.
 */
zMat3D *zMat3DT(zMat3D *m, zMat3D *tm)
{
  return zMat3DCreate( tm,
    m->e[0][0], m->e[1][0], m->e[2][0],
    m->e[0][1], m->e[1][1], m->e[2][1],
    m->e[0][2], m->e[1][2], m->e[2][2] );
}

/* zMat3DRow
 * - abstract row vector from 3D matrix.
 */
zVec3D *zMat3DRow(zMat3D *m, int i, zVec3D *v)
{
  return zVec3DCreate( v, m->e[0][i], m->e[1][i], m->e[2][i] );
}

/* zMat3DCol
 * - abstract column vector from 3D matrix.
 */
zVec3D *zMat3DCol(zMat3D *m, int i, zVec3D *v)
{
  return zVec3DCopy( &m->v[i], v );
}

/* ********************************************************** */
/* arithmetics
 * ********************************************************** */

/* zMat3DAdd
 * - add two 3D matrices.
 */
zMat3D *zMat3DAdd(zMat3D *m1, zMat3D *m2, zMat3D *m)
{
  zVec3DAdd( &m1->v[0], &m2->v[0], &m->v[0] );
  zVec3DAdd( &m1->v[1], &m2->v[1], &m->v[1] );
  zVec3DAdd( &m1->v[2], &m2->v[2], &m->v[2] );
  return m;
}

/* zMat3DSub
 * - subtract a 3D matrix from another.
 */
zMat3D *zMat3DSub(zMat3D *m1, zMat3D *m2, zMat3D *m)
{
  zVec3DSub( &m1->v[0], &m2->v[0], &m->v[0] );
  zVec3DSub( &m1->v[1], &m2->v[1], &m->v[1] );
  zVec3DSub( &m1->v[2], &m2->v[2], &m->v[2] );
  return m;
}

/* zMat3DRev
 * - reverse 3D matrix.
 */
zMat3D *zMat3DRev(zMat3D *m, zMat3D *rm)
{
  zVec3DRev( &m->v[0], &rm->v[0] );
  zVec3DRev( &m->v[1], &rm->v[1] );
  zVec3DRev( &m->v[2], &rm->v[2] );
  return rm;
}

/* zMat3DMul
 * - multiply 3D matrix by a scalar.
 */
zMat3D *zMat3DMul(zMat3D *m, double k, zMat3D *mm)
{
  zVec3DMul( &m->v[0], k, &mm->v[0] );
  zVec3DMul( &m->v[1], k, &mm->v[1] );
  zVec3DMul( &m->v[2], k, &mm->v[2] );
  return mm;
}

/* zMat3DDiv
 * - divide 3D matrix by a scalar.
 */
zMat3D *zMat3DDiv(zMat3D *m, double k, zMat3D *dm)
{
  if( k == 0 ){
    ZRUNWARN( ZEO_ERR_ZERODIV );
    return NULL;
  }
  k = 1.0 / k;
  zVec3DMul( &m->v[0], k, &dm->v[0] );
  zVec3DMul( &m->v[1], k, &dm->v[1] );
  zVec3DMul( &m->v[2], k, &dm->v[2] );
  return dm;
}

/* zMat3DCat
 * - concatenate 3D matrix.
 */
zMat3D *zMat3DCat(zMat3D *m1, double k, zMat3D *m2, zMat3D *m)
{
  zVec3DCat( &m1->v[0], k, &m2->v[0], &m->v[0] );
  zVec3DCat( &m1->v[1], k, &m2->v[1], &m->v[1] );
  zVec3DCat( &m1->v[2], k, &m2->v[2], &m->v[2] );
  return m;
}

/* zMat3DDyad
 * - dyadic product.
 */
zMat3D *zMat3DDyad(zVec3D *v1, zVec3D *v2, zMat3D *dyad)
{
  zVec3DMul( v1, v2->e[0], &dyad->v[0] );
  zVec3DMul( v1, v2->e[1], &dyad->v[1] );
  zVec3DMul( v1, v2->e[2], &dyad->v[2] );
  return dyad;
}

/* zVec3DOuterProdMat3D
 * - create a skew-symmetric outer-product matrix.
 */
zMat3D *zVec3DOuterProdMat3D(zVec3D *v, zMat3D *m)
{
  return zMat3DCreate( m,
           0,-v->e[zZ], v->e[zY],
    v->e[zZ],        0,-v->e[zX],
   -v->e[zY], v->e[zX],        0 );
}

/* zVec3DOuterProd2Mat3D
 * - create a twice-outer-product matrix.
 */
zMat3D *zVec3DOuterProd2Mat3D(zVec3D *v1, zVec3D *v2, zMat3D *m)
{
  double val;

  zMat3DDyad( v2, v1, m );
  val = zVec3DInnerProd( v1, v2 );
  m->e[0][0] -= val;
  m->e[1][1] -= val;
  m->e[2][2] -= val;
  return m;
}

/* ********************************************************** */
/* inverse of a 3x3 matrix
 * ********************************************************** */

static void _zMat3DInvRow(zMat3D *m, zMat3D *im, int i, int j, int k, double det);

/* zMat3DDet
 * - determinant of 3D matrix.
 */
double zMat3DDet(zMat3D *m)
{
  return ( m->e[1][1]*m->e[2][2] - m->e[1][2]*m->e[2][1] ) * m->e[0][0]
       + ( m->e[0][2]*m->e[2][1] - m->e[0][1]*m->e[2][2] ) * m->e[1][0]
       + ( m->e[0][1]*m->e[1][2] - m->e[0][2]*m->e[1][1] ) * m->e[2][0];
}

/* (static)
 * _zMat3DInvRow
 * - misc for real inverse of 3D matrix.
 */
void _zMat3DInvRow(zMat3D *m, zMat3D *im, int i, int j, int k, double det)
{
  im->e[i][i] = ( m->e[j][j]*m->e[k][k] - m->e[k][j]*m->e[j][k] ) / det;
  im->e[j][i] = ( m->e[k][i]*m->e[j][k] - m->e[j][i]*m->e[k][k] ) / det;
  im->e[k][i] = ( m->e[j][i]*m->e[k][j] - m->e[k][i]*m->e[j][j] ) / det;
}

/* zMat3DInv
 * - inverse matrix of 3x3 matrix.
 */
zMat3D *zMat3DInv(zMat3D *m, zMat3D *im)
{
  double det;

  if( zIsTiny( ( det = zMat3DDet( m ) ) ) ){
    ZRUNERROR( ZEO_ERR_SINGULARMAT );
    return NULL;
  }
  _zMat3DInvRow( m, im, 0, 1, 2, det );
  _zMat3DInvRow( m, im, 1, 2, 0, det );
  _zMat3DInvRow( m, im, 2, 0, 1, det );
  return im;
}

/* ********************************************************** */
/* multiplication of a 3D vector by a 3x3 matrix
 * ********************************************************** */

static double _zMulInvMatVec3DRow(zMat3D *m, zVec3D *v, int i, int j, int k, double det);
static zVec3D *_zMulInvMatVec3D(zMat3D *m, zVec3D *v, zVec3D *miv, double det);

/* zMulMatVec3D
 * - multiply a 3D vector by 3D matrix.
 */
zVec3D *zMulMatVec3D(zMat3D *m, zVec3D *v, zVec3D *mv)
{
  return zVec3DCreate( mv,
    m->e[0][0]*v->e[0] + m->e[1][0]*v->e[1] + m->e[2][0]*v->e[2],
    m->e[0][1]*v->e[0] + m->e[1][1]*v->e[1] + m->e[2][1]*v->e[2],
    m->e[0][2]*v->e[0] + m->e[1][2]*v->e[1] + m->e[2][2]*v->e[2] );
}

/* zMulMatTVec3D
 * - multiply a 3D vector by a transpose matrix.
 */
zVec3D *zMulMatTVec3D(zMat3D *m, zVec3D *v, zVec3D *mv)
{
  return zVec3DCreate( mv,
    m->e[0][0]*v->e[0] + m->e[0][1]*v->e[1] + m->e[0][2]*v->e[2],
    m->e[1][0]*v->e[0] + m->e[1][1]*v->e[1] + m->e[1][2]*v->e[2],
    m->e[2][0]*v->e[0] + m->e[2][1]*v->e[1] + m->e[2][2]*v->e[2] );
}

/* (static)
 * _zMulInvMatVec3DRow
 * - misc for multiplication of a vector by real inverse of 3D matrix.
 */
double _zMulInvMatVec3DRow(zMat3D *m, zVec3D *v, int i, int j, int k, double det)
{
  return
    ( ( m->e[j][j]*m->e[k][k] - m->e[k][j]*m->e[j][k] ) * v->e[i]
    + ( m->e[k][i]*m->e[j][k] - m->e[j][i]*m->e[k][k] ) * v->e[j]
    + ( m->e[j][i]*m->e[k][j] - m->e[k][i]*m->e[j][j] ) * v->e[k] ) / det;
}

/* (static)
 * _zMulInvMatVec3D
 * - multiply a 3D vector by inverse matrix of 3x3 matrix.
 */
zVec3D *_zMulInvMatVec3D(zMat3D *m, zVec3D *v, zVec3D *miv, double det)
{
  return zVec3DCreate( miv,
    _zMulInvMatVec3DRow( m, v, 0, 1, 2, det ),
    _zMulInvMatVec3DRow( m, v, 1, 2, 0, det ),
    _zMulInvMatVec3DRow( m, v, 2, 0, 1, det ) );
}

/* zMulInvMatVec3D
 * - multiply a 3D vector by inverse matrix of 3x3 matrix.
 */
zVec3D *zMulInvMatVec3D(zMat3D *m, zVec3D *v, zVec3D *miv)
{
  double det;

  if( zIsTiny( ( det = zMat3DDet( m ) ) ) ){
    ZRUNERROR( ZEO_ERR_SINGULARMAT );
    return NULL;
  }
  return _zMulInvMatVec3D( m, v, miv, det );
}

/* zVec3DCatRatio
 * - concatenate ratio of vector.
 */
double *zVec3DCatRatio(zVec3D *v1, zVec3D *v2, zVec3D *v3, zVec3D *v, double ratio[])
{
  zMat3D m, im;

  zVec3DCopy( v1, &m.v[0] );
  zVec3DCopy( v2, &m.v[1] );
  zVec3DCopy( v3, &m.v[2] );
  zMat3DInv( &m, &im );
  zMulMatVec3D( &im, v, (zVec3D*)ratio );
  return ratio;
}

/* ********************************************************** */
/* multiplication of a 3x3 matrix by another 3x3 matrix
 * ********************************************************** */

/* zMulMatMat3D
 * - multiply two 3D matrices.
 */
zMat3D *zMulMatMat3D(zMat3D *m1, zMat3D *m2, zMat3D *m)
{
  zMat3D tmp;

  zMulMatVec3D( m1, &m2->v[0], &tmp.v[0] );
  zMulMatVec3D( m1, &m2->v[1], &tmp.v[1] );
  zMulMatVec3D( m1, &m2->v[2], &tmp.v[2] );
  zMat3DCopy( &tmp, m );
  return m;
}

/* zMulMatTMat3D
 * - multiply a matrix by transpose matrix from leftside.
 */
zMat3D *zMulMatTMat3D(zMat3D *m1, zMat3D *m2, zMat3D *m)
{
  zMat3D tmp;

  zMulMatTVec3D( m1, &m2->v[0], &tmp.v[0] );
  zMulMatTVec3D( m1, &m2->v[1], &tmp.v[1] );
  zMulMatTVec3D( m1, &m2->v[2], &tmp.v[2] );
  zMat3DCopy( &tmp, m );
  return m;
}

/* zMulMatMatT3D
 * - multiply a matrix by transpose matrix from rightside.
 */
zMat3D *zMulMatMatT3D(zMat3D *m1, zMat3D *m2, zMat3D *m)
{
  register int i;
  zMat3D tmp;

  for( i=0; i<3; i++ ){
    zVec3DMul( &m1->v[0], m2->e[0][i], &tmp.v[i] );
    zVec3DCatDRC( &tmp.v[i], m2->e[1][i], &m1->v[1] );
    zVec3DCatDRC( &tmp.v[i], m2->e[2][i], &m1->v[2] );
  }
  zMat3DCopy( &tmp, m );
  return m;
}

/* _zMulInvMatMat3D
 * - multiply 3x3 matrix by inverse matrix of another 3x3 matrix.
 */
zMat3D *zMulInvMatMat3D(zMat3D *m1, zMat3D *m2, zMat3D *m)
{
  double det;

  if( zIsTiny( ( det = zMat3DDet( m1 ) ) ) ){
    ZRUNERROR( ZEO_ERR_SINGULARMAT );
    return NULL;
  }
  zMulInvMatVec3D( m1, &m2->v[0], &m->v[0] );
  zMulInvMatVec3D( m1, &m2->v[1], &m->v[1] );
  zMulInvMatVec3D( m1, &m2->v[2], &m->v[2] );
  return m;
}

/* ********************************************************** */
/* multiplication of a 6D spatial vector by a 3x3 matrix
 * ********************************************************** */

/* zMulMatVec6D
 * - multiply a 6D vector by 3D matrix.
 */
zVec6D *zMulMatVec6D(zMat3D *m, zVec6D *v, zVec6D *mv)
{
  zMulMatVec3D( m, zVec6DLin(v), zVec6DLin(mv) );
  zMulMatVec3D( m, zVec6DAng(v), zVec6DAng(mv) );
  return mv;
}

/* zMulMatTVec6D
 * - multiply a 6D vector by transpose matrix.
 */
zVec6D *zMulMatTVec6D(zMat3D *m, zVec6D *v, zVec6D *mv)
{
  zMulMatTVec3D( m, zVec6DLin(v), zVec6DLin(mv) );
  zMulMatTVec3D( m, zVec6DAng(v), zVec6DAng(mv) );
  return mv;
}

/* ********************************************************** */
/* rotation
 * ********************************************************** */

/* (static)
 * _zMat3DRotRPYSC
 * - rotate matrix along a base axis.
 */
static zMat3D *_zMat3DRotRPYSC(zMat3D *m, int a0, int a1, int a2, double s, double c, zMat3D *rm);
zMat3D *_zMat3DRotRPYSC(zMat3D *m, int a0, int a1, int a2, double s, double c, zMat3D *rm)
{
  rm->e[a0][a0] = m->e[a0][a0]*c - m->e[a0][a1]*s;
  rm->e[a1][a0] = m->e[a1][a0]*c - m->e[a1][a1]*s;
  rm->e[a2][a0] = m->e[a2][a0]*c - m->e[a2][a1]*s;
  rm->e[a0][a1] = m->e[a0][a0]*s + m->e[a0][a1]*c;
  rm->e[a1][a1] = m->e[a1][a0]*s + m->e[a1][a1]*c;
  rm->e[a2][a1] = m->e[a2][a0]*s + m->e[a2][a1]*c;
  rm->e[a0][a2] = m->e[a0][a2];
  rm->e[a1][a2] = m->e[a1][a2];
  rm->e[a2][a2] = m->e[a2][a2];
  return rm;
}

/* (static)
 * _zMat3DRotRPYSCDRC
 * - rotate matrix directly along a base axis.
 */
static zMat3D *_zMat3DRotRPYSCDRC(zMat3D *m, int a0, int a1, int a2, double s, double c);
zMat3D *_zMat3DRotRPYSCDRC(zMat3D *m, int a0, int a1, int a2, double s, double c)
{
  double r00, r01, r02, r10, r11, r12;

  r00 = m->e[a0][a0];
  r01 = m->e[a1][a0];
  r02 = m->e[a2][a0];
  r10 = m->e[a0][a1];
  r11 = m->e[a1][a1];
  r12 = m->e[a2][a1];
  m->e[a0][a0] = r00*c-r10*s;
  m->e[a1][a0] = r01*c-r11*s;
  m->e[a2][a0] = r02*c-r12*s;
  m->e[a0][a1] = r00*s+r10*c;
  m->e[a1][a1] = r01*s+r11*c;
  m->e[a2][a1] = r02*s+r12*c;
  return m;
}

/* zMat3DRotRollSC, zMat3DRotRoll, zMat3DRotRollSCDRC, zMat3DRotRollDRC
 * - rotate matrix along x-axis.
 */
zMat3D *zMat3DRotRollSC(zMat3D *m, double s, double c, zMat3D *rm){
  return _zMat3DRotRPYSC( m, 1, 2, 0, s, c, rm );
}
zMat3D *zMat3DRotRoll(zMat3D *m, double theta, zMat3D *rm){
  double s, c;
  zSinCos( theta, &s, &c );
  return zMat3DRotRollSC( m, s, c, rm );
}
zMat3D *zMat3DRotRollSCDRC(zMat3D *m, double s, double c){
  return _zMat3DRotRPYSCDRC( m, 1, 2, 0, s, c );
}
zMat3D *zMat3DRotRollDRC(zMat3D *m, double theta){
  double s, c;
  zSinCos( theta, &s, &c );
  return zMat3DRotRollSCDRC( m, s, c );
}

/* zMat3DRotPitchSC, zMat3DRotPitch, zMat3DRotPitchSCDRC, zMat3DRotPitchDRC
 * - rotate matrix along y-axis.
 */
zMat3D *zMat3DRotPitchSC(zMat3D *m, double s, double c, zMat3D *rm){
  return _zMat3DRotRPYSC( m, 2, 0, 1, s, c, rm );
}
zMat3D *zMat3DRotPitch(zMat3D *m, double theta, zMat3D *rm){
  double s, c;
  zSinCos( theta, &s, &c );
  return zMat3DRotPitchSC( m, s, c, rm );
}
zMat3D *zMat3DRotPitchSCDRC(zMat3D *m, double s, double c){
  return _zMat3DRotRPYSCDRC( m, 2, 0, 1, s, c );
}
zMat3D *zMat3DRotPitchDRC(zMat3D *m, double theta){
  double s, c;
  zSinCos( theta, &s, &c );
  return zMat3DRotPitchSCDRC( m, s, c );
}

/* zMat3DRotYawSC, zMat3DRotYaw, zMat3DRotYawSCDRC, zMat3DRotYawDRC
 * - rotate matrix along z-axis.
 */
zMat3D *zMat3DRotYawSC(zMat3D *m, double s, double c, zMat3D *rm){
  return _zMat3DRotRPYSC( m, 0, 1, 2, s, c, rm );
}
zMat3D *zMat3DRotYaw(zMat3D *m, double theta, zMat3D *rm){
  double s, c;
  zSinCos( theta, &s, &c );
  return zMat3DRotYawSC( m, s, c, rm );
}
zMat3D *zMat3DRotYawSCDRC(zMat3D *m, double s, double c){
  return _zMat3DRotRPYSCDRC( m, 0, 1, 2, s, c );
}
zMat3D *zMat3DRotYawDRC(zMat3D *m, double theta){
  double s, c;
  zSinCos( theta, &s, &c );
  return zMat3DRotYawSCDRC( m, s, c );
}

/* zMat3DZYXSC, zMat3DZYX
 * - 3D attitude matrix expressed by z-y-x Eulerian angle.
 */
zMat3D *zMat3DZYXSC(zMat3D *m, double sa, double ca, double se, double ce, double st, double ct)
{
  return zMat3DCreate( m,
    ce*ca, st*se*ca-ct*sa, ct*se*ca+st*sa,
    ce*sa, st*se*sa+ct*ca, ct*se*sa-st*ca,
   -se,    st*ce,          ct*ce );
}
zMat3D *zMat3DZYX(zMat3D *m, double azim, double elev, double tilt)
{
  double sa, ca, se, ce, st, ct;

  zSinCos( azim, &sa, &ca );
  zSinCos( elev, &se, &ce );
  zSinCos( tilt, &st, &ct );
  return zMat3DZYXSC( m, sa, ca, se, ce, st, ct );
}

/* zMat3DToZYX
 * - quasi 3D vector for the expression of 3D attitude
 *   as z-y-x Eulerian angle.
 */
zVec3D *zMat3DToZYX(zMat3D *m, zVec3D *angle)
{
  double azim, ca, sa;

  azim = atan2( m->e[0][1], m->e[0][0] );
  zSinCos( azim, &sa, &ca );
  return zVec3DCreate( angle,
    azim,
    atan2(-m->e[0][2], m->e[0][0]*ca+m->e[0][1]*sa ),
    atan2( m->e[2][0]*sa-m->e[2][1]*ca, -m->e[1][0]*sa+m->e[1][1]*ca ) );
}

/* zMat3DZYZSC, zMat3DZYZ
 * - 3D attitude matrix expressed by z-y-z Eulerian angle.
 */
zMat3D *zMat3DZYZSC(zMat3D *m, double sh, double ch, double sp, double cp, double sb, double cb)
{
  return zMat3DCreate( m,
    ch*cp*cb-sh*sb,-ch*cp*sb-sh*cb, ch*sp,
    sh*cp*cb+ch*sb,-sh*cp*sb+ch*cb, sh*sp,
      -sp*cb,          sp*sb,          cp );
}
zMat3D *zMat3DZYZ(zMat3D *m, double heading, double pitch, double bank)
{
  double sh, ch, sp, cp, sb, cb;

  zSinCos( heading, &sh, &ch );
  zSinCos( pitch, &sp, &cp );
  zSinCos( bank, &sb, &cb );
  return zMat3DZYZSC( m, sh, ch, sp, cp, sb, cb );
}

/* zMat3DToZYZ
 * - quasi 3D vector for the expression of 3D attitude
 *   as z-y-z Eulerian angle.
 */
zVec3D *zMat3DToZYZ(zMat3D *m, zVec3D *angle)
{
  double heading, sh, ch;

  heading = atan2( m->e[2][1], m->e[2][0] );
  zSinCos( heading, &sh, &ch );
  return zVec3DCreate( angle,
    heading,
    atan2( m->e[2][0]*ch+m->e[2][1]*sh, m->e[2][2] ),
    atan2(-m->e[0][0]*sh+m->e[0][1]*ch,-m->e[1][0]*sh+m->e[1][1]*ch ) );
}

/* zMat3DAA
 * - 3D attitude matrix expressed by angle-axis vector.
 */
zMat3D *zMat3DAA(zMat3D *m, zVec3D *aa)
{
  double l, c, s, r1, r2, r3;
  zVec3D n, nl;

  if( zVec3DIsTiny( aa ) ) return zMat3DIdent( m );
  l = zVec3DNorm( aa );
  zSinCos( l, &s, &c );
  zVec3DMul( aa, 1.0/l, &n );
  zVec3DMul( &n, 1-c, &nl );
  zVec3DOuterProd2Mat3D( &n, &nl, m );
  r1 = s * n.e[zX];
  r2 = s * n.e[zY];
  r3 = s * n.e[zZ];
  m->e[0][0] += 1.0;
  m->e[1][0] -= r3;
  m->e[2][0] += r2;
  m->e[0][1] += r3;
  m->e[1][1] += 1.0;
  m->e[2][1] -= r1;
  m->e[0][2] -= r2;
  m->e[1][2] += r1;
  m->e[2][2] += 1.0;
  return m;
}

/* zMat3DToAA
 * - convert 3x3 matrix to equivalent angle-axis vector.
 */
zVec3D *zMat3DToAA(zMat3D *m, zVec3D *aa)
{
  register int i;
  double l, a;

  zVec3DCreate( aa,
    m->e[1][2]-m->e[2][1], m->e[2][0]-m->e[0][2], m->e[0][1]-m->e[1][0] );
  l = zVec3DNorm( aa );
  a = atan2( l, m->e[0][0]+m->e[1][1]+m->e[2][2]-1 );
  if( zIsTiny( l ) )
    for( i=0; i<3; i++ ){
      zMat3DCol( m, i, aa );
      aa->e[i] += 1.0;
      if( !zIsTiny( l = zVec3DNorm(aa) ) )
        break;
    }
  return zVec3DMulDRC( aa, a / l );
}

/* zRotMat3D
 * - rotational multiplication for 3D matrices.
 */
zMat3D *zRotMat3D(zMat3D *r, zMat3D *m, zMat3D *rm)
{
  zMulMatMat3D( r, m, rm );
  return zMulMatMatT3DDRC( rm, r );
}

/* zRotMat3DInv
 * - inverse rotational multiplication for 3D matrices.
 */
zMat3D *zRotMat3DInv(zMat3D *r, zMat3D *m, zMat3D *rm)
{
  zMulMatMat3D( m, r, rm );
  return zMulMatTMat3DDRC( r, rm );
}

/* zMulVecOPMat3D
 * - multiply cross product of vector and matrix.
 */
zMat3D *zMulVecOPMat3D(zVec3D *ohm, zMat3D *m, zMat3D *mv)
{
  zVec3DOuterProd( ohm, &m->v[0], &mv->v[0] );
  zVec3DOuterProd( ohm, &m->v[1], &mv->v[1] );
  zVec3DOuterProd( ohm, &m->v[2], &mv->v[2] );
  return mv;
}

/* zMat3DRot
 * - rotate a 3D matrix along an arbitrary axis.
 */
zMat3D *zMat3DRot(zMat3D *m, zVec3D *aa, zMat3D *rm)
{
  zMat3D ma;

  zMat3DAA( &ma, aa );
  return zMulMatMat3D( &ma, m, rm );
}

/* zMat3DRotCat
 * - concatenate a 3D rotational vector to a 3D matrix.
 */
zMat3D *zMat3DRotCat(zMat3D *m, zVec3D *omega, double dt, zMat3D *rm)
{
  zVec3D aa;

  zVec3DMul( omega, dt, &aa );
  return zMat3DRot( m, &aa, rm );
}

/* zAACascade
 * - cascade an angle-axis vector to another.
 */
zVec3D *zAACascade(zVec3D *aa1, zVec3D *aa2, zVec3D *aa)
{
  zMat3D m1, m2, m;

  zMat3DAA( &m1, aa1 );
  zMat3DAA( &m2, aa2 );
  zMulMatMat3D( &m2, &m1, &m );
  return zMat3DToAA( &m, aa );
}

/* zMat3DError
 * - error vector between two attitude matrices.
 */
zVec3D *zMat3DError(zMat3D *m1, zMat3D *m2, zVec3D *err)
{
  zMat3D em;

  zMulMatMatT3D( m1, m2, &em );
  return zMat3DToAA( &em, err );
}

/* zAAError
 * - error vector between two angle-axis vectors.
 */
zVec3D *zAAError(zVec3D *a1, zVec3D *a2, zVec3D *err)
{
  zMat3D m1, m2;

  zMat3DAA( &m1, a1 );
  zMat3DAA( &m2, a2 );
  return zMat3DError( &m1, &m2, err );
}

/* ********************************************************** */
/* differential kinematics
 * ********************************************************** */

/* zMat3DDif
 * - numerical differentiation of 3D matrix.
 */
zVec3D *zMat3DDif(zMat3D *m, zMat3D *mnew, double dt, zVec3D *omega)
{
  zMat3DError( mnew, m, omega );
  return zVec3DDivDRC( omega, dt );
}

/* ********************************************************** */
/* eigensystem
 * ********************************************************** */

/* (static)
 * _zMat3DSymEigRot
 * - transformation of a symmetric matrix by Jacobi's rotation.
 */
static bool _zMat3DSymEigRot(zMat3D *m, zMat3D *r, int i, int j);
bool _zMat3DSymEigRot(zMat3D *m, zMat3D *r, int i, int j)
{
  register int k;
  double v, t, ti, c, s;
  double tmp1, tmp2;

  if( zIsTiny( m->e[j][i] ) ) return true;
  v = 0.5 * ( m->e[j][j] - m->e[i][i] ) / m->e[j][i];
  ti = sqrt( v * v + 1 ); /* sqrt */
  t = -v + ( v < 0 ? -ti : ti );
  s = t * ( c = 1 / sqrt( t*t + 1 ) );
  tmp1 = t * m->e[j][i];
  m->e[j][i] = m->e[i][j] = 0;
  m->e[i][i] -= tmp1;
  m->e[j][j] += tmp1;
  for( k=0; k<3; k++ ){
    /* update of transformation matrix */
    tmp1 = r->e[i][k];
    tmp2 = r->e[j][k];
    r->e[i][k] = c * tmp1 - s * tmp2;
    r->e[j][k] = s * tmp1 + c * tmp2;
    /* update of eigenmatrix */
    if( k == i || k == j ) continue;
    tmp1 = m->e[k][i];
    tmp2 = m->e[k][j];
    m->e[k][i] = m->e[i][k] = c * ( tmp1 - t * tmp2 );
    m->e[k][j] = m->e[j][k] = c * ( tmp2 + t * tmp1 );
  }
  return false;
}

/* zMat3DSymEig
 * - eigenvalues of a symmetric 3x3 matrix by Jacobi's method.
 */
void zMat3DSymEig(zMat3D *m, double eval[], zVec3D evec[])
{
  int n = 0;
  zMat3D l, r;
  bool ok;

  zMat3DCopy( m, &l );
  zMat3DIdent( &r );
  /* iterative elimination of non-diagonal components */
  do{
    ok = true;
    if( !_zMat3DSymEigRot( &l, &r, 1, 0 ) ) ok = false;
    if( !_zMat3DSymEigRot( &l, &r, 2, 0 ) ) ok = false;
    if( !_zMat3DSymEigRot( &l, &r, 2, 1 ) ) ok = false;
    if( n++ > Z_MAX_ITER_NUM ){
      ZITERWARN( Z_MAX_ITER_NUM );
      break;
    }
  } while( !ok );
  for( n=0; n<3; n++ ){
    eval[n] = l.e[n][n];
    zVec3DCopy( &r.v[n], &evec[n] );
  }
}

/* ********************************************************** */
/* I/O
 * ********************************************************** */

/* zMat3DFRead
 * - input a 3D matrix from file.
 */
zMat3D *zMat3DFRead(FILE *fp, zMat3D *m)
{
  register int i, j;

  for( i=0; i<3; i++ )
    for( j=0; j<3; j++ )
      m->e[j][i] = zFDouble( fp );
  return m;
}

/* zMat3DFWrite
 * - output a 3D matrix to file.
 */
void zMat3DFWrite(FILE *fp, zMat3D *m)
{
  if( !m )
    fprintf( fp, "(null 3D matrix)\n" );
  else{
    fprintf( fp, "{\n" );
    fprintf( fp, " %.10g, %.10g, %.10g\n", m->e[0][0], m->e[1][0], m->e[2][0] );
    fprintf( fp, " %.10g, %.10g, %.10g\n", m->e[0][1], m->e[1][1], m->e[2][1] );
    fprintf( fp, " %.10g, %.10g, %.10g\n", m->e[0][2], m->e[1][2], m->e[2][2] );
    fprintf( fp, "}\n" );
  }
}

/* METHOD:
 * zMat3DFWriteXML - xml output.
 * ... yet testing.
 */
void zMat3DFWriteXML(FILE *fp, zMat3D *m)
{
  fprintf( fp, "\"%.10g %.10g %.10g,\n", m->e[0][0], m->e[1][0], m->e[2][0] );
  fprintf( fp, " %.10g %.10g %.10g,\n",  m->e[0][1], m->e[1][1], m->e[2][1] );
  fprintf( fp, " %.10g %.10g %.10g\"",   m->e[0][2], m->e[1][2], m->e[2][2] );
}
