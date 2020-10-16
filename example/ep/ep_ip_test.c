#include <zeo/zeo.h>

#define N 100

int main(void)
{
  zEP ep1, ep2, ep;
  zMat3D m1, m2, m;
  zVec3D e1, e2, e;
  double t;
  register int i;

  zRandInit();
  zMat3DZYX( &m1, zRandF(-zPI,zPI), zRandF(-zPI,zPI), zRandF(-zPI,zPI) );
  zMat3DZYX( &m2, zRandF(-zPI,zPI), zRandF(-zPI,zPI), zRandF(-zPI,zPI) );
  zMat3DToZYX( &m1, &e1 );
  zMat3DToZYX( &m2, &e2 );
  zMat3DToEP( &m1, &ep1 );
  zMat3DToEP( &m2, &ep2 );
  for( i=0; i<=N; i++ ){
    t = (double)i / N;
    zEPInnerDiv( &ep1, &ep2, t, &ep );
    zMat3DEP( &m, &ep );
    zMat3DToZYX( &m, &e );
    printf( "%f %f %f ", zVec3DElem(&e1,zX), zVec3DElem(&e2,zX), zVec3DElem(&e,zX) );
    printf( "%f %f %f ", zVec3DElem(&e1,zY), zVec3DElem(&e2,zY), zVec3DElem(&e,zY) );
    printf( "%f %f %f\n", zVec3DElem(&e1,zZ), zVec3DElem(&e2,zZ), zVec3DElem(&e,zZ) );
  }
  return 0;
}
