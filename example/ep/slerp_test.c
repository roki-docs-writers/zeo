#include <zeo/zeo_ep.h>

#define N 10000

int main(void)
{
  zMat3D m1, m2, m;
  zVec3D e1, e2, e, em;
  double t;
  register int i;
#if 0
  FILE *fp[2];
#endif

  zRandInit();
  zMat3DZYX( &m1, zRandF(-zPI,zPI), zRandF(-zPI,zPI), zRandF(-zPI,zPI) );
  zMat3DZYX( &m2, zRandF(-zPI,zPI), zRandF(-zPI,zPI), zRandF(-zPI,zPI) );
  zMat3DToAA( &m1, &e1 );
  zMat3DToAA( &m2, &e2 );
#if 0
  fp[0] = fopen( "slerp.zvs", "w" );
  fp[1] = fopen( "vecip.zvs", "w" );
#endif
  for( i=0; i<=N; i++ ){
    t = (double)i / N;
    zMat3DInnerDiv( &m1, &m2, t, &m );
    zMat3DToAA( &m, &e );
    zVec3DInnerDiv( &e1, &e2, t, &em );
    printf( "%f %f %f %f ", zVec3DElem(&e1,zX), zVec3DElem(&e2,zX), zVec3DElem(&e,zX), zVec3DElem(&em,zX) );
    printf( "%f %f %f %f ", zVec3DElem(&e1,zY), zVec3DElem(&e2,zY), zVec3DElem(&e,zY), zVec3DElem(&em,zY) );
    printf( "%f %f %f %f\n", zVec3DElem(&e1,zZ), zVec3DElem(&e2,zZ), zVec3DElem(&e,zZ), zVec3DElem(&em,zZ) );
#if 0
    fprintf( fp[0], "0.1 3 " ); zVec3DDataFWrite(fp[0],&e);
    fprintf( fp[1], "0.1 3 " ); zVec3DDataFWrite(fp[1],&em);
#endif
  }
#if 0
  fclose( fp[0] );
  fclose( fp[1] );
#endif
  return 0;
}
