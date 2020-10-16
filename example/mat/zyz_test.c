#include <zeo/zeo_mat3d.h>

int main(void)
{
  zVec3D v;
  zMat3D m;
  double heading, pitch, bank;

  do{
    printf( "heading  = " ); scanf( "%lf", &heading );
    printf( "pitch = " ); scanf( "%lf", &pitch );
    printf( "bank   = " ); scanf( "%lf", &bank );

    zMat3DWrite( zMat3DZYZ( &m, zDeg2Rad(heading), zDeg2Rad(pitch), zDeg2Rad(bank) ) );

    zMat3DToZYZ( &m, &v );
    heading  = zRad2Deg( zVec3DElem( &v, zX ) );
    pitch = zRad2Deg( zVec3DElem( &v, zY ) );
    bank   = zRad2Deg( zVec3DElem( &v, zZ ) );

    printf( "%f %f %f\n", heading, pitch, bank );
  } while( heading!=0 || pitch!=0 || bank!=0 );
  return 0;
}
