#include <zeo/zeo_mat3d.h>

int main(void)
{
  zMat3D m;

  zRandInit();
  zMat3DCreate( &m,
    zRandF(-10,10), zRandF(-10,10), zRandF(-10,10),
    zRandF(-10,10), zRandF(-10,10), zRandF(-10,10),
    zRandF(-10,10), zRandF(-10,10), zRandF(-10,10) );
  zMat3DWrite( &m );

  zMat3DT( &m, &m );
  zMat3DWrite( &m );

  zMat3DT( &m, &m );
  zMat3DWrite( &m );
  return 0;
}
