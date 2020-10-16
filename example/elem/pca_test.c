#include <zeo/zeo_elem.h>

#define N 1000
zVec3D vert[N];

void test_vert(double x, double y, double z, double a, double b, double c, double wx, double wy, double wz)
{
  register int i;
  zMat3D ori;

  zMat3DZYX( &ori, zDeg2Rad(a), zDeg2Rad(b), zDeg2Rad(c) );
  for( i=0; i<N; i++ ){
    zVec3DCreate( &vert[i],
      zRandF(-wx,wx), zRandF(-wy,wy), zRandF(-wz,wz) );
    zMulMatVec3DDRC( &ori, &vert[i] );
    zVec3DElem(&vert[i],zX) += x;
    zVec3DElem(&vert[i],zY) += y;
    zVec3DElem(&vert[i],zZ) += z;
  }
}

void output_vert(void)
{
  FILE *fp;
  register int i;
  zPlane3D pl;
  zVec3D v, c, d[3];

  fp = fopen( "a", "w" );
  for( i=0; i<N; i++ )
    zVec3DDataFWrite( fp, &vert[i] );
  fclose( fp );

  zVec3DPCA( vert, N, d );
  fp = fopen( "b", "w" );
  zPlane3DCreate( &pl, Z_ZEROVEC3D, &d[0] );
  for( i=0; i<N; i++ ){
    zPlane3DProj( &pl, &vert[i], &v );
    zVec3DDataFWrite( fp, &v );
  }
  fclose( fp );
  fp = fopen( "c", "w" );
  zPlane3DCreate( &pl, Z_ZEROVEC3D, &d[1] );
  for( i=0; i<N; i++ ){
    zPlane3DProj( &pl, &vert[i], &v );
    zVec3DDataFWrite( fp, &v );
  }
  fclose( fp );
  fp = fopen( "d", "w" );
  zPlane3DCreate( &pl, Z_ZEROVEC3D, &d[2] );
  for( i=0; i<N; i++ ){
    zPlane3DProj( &pl, &vert[i], &v );
    zVec3DDataFWrite( fp, &v );
  }
  fclose( fp );

  zVec3DBaryPCA( vert, N, &c, d );
  fp = fopen( "e", "w" );
  zPlane3DCreate( &pl, &c, &d[0] );
  for( i=0; i<N; i++ ){
    zPlane3DProj( &pl, &vert[i], &v );
    zVec3DDataFWrite( fp, &v );
  }
  fclose( fp );
  fp = fopen( "f", "w" );
  zPlane3DCreate( &pl, &c, &d[1] );
  for( i=0; i<N; i++ ){
    zPlane3DProj( &pl, &vert[i], &v );
    zVec3DDataFWrite( fp, &v );
  }
  fclose( fp );
  fp = fopen( "g", "w" );
  zPlane3DCreate( &pl, &c, &d[2] );
  for( i=0; i<N; i++ ){
    zPlane3DProj( &pl, &vert[i], &v );
    zVec3DDataFWrite( fp, &v );
  }
  fclose( fp );
}

#define X 1
#define Y 0.5
#define Z 1.5
int main(void)
{
  zRandInit();
  test_vert( X, Y, Z, 10, 20, 30, 2, 1, 3 );
  output_vert();
  return 0;
}
