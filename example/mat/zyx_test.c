#include <zeo/zeo_mat3d.h>

int main(void)
{
  zVec3D v;
  zMat3D m;
  double azim, elev, tilt;

  do{
    printf( "azimuth   = " ); scanf( "%lf", &azim );
    printf( "elevation = " ); scanf( "%lf", &elev );
    printf( "tilt      = " ); scanf( "%lf", &tilt );

    zMat3DWrite( zMat3DZYX( &m, zDeg2Rad(azim), zDeg2Rad(elev), zDeg2Rad(tilt) ) );

    zMat3DToZYX( &m, &v );
    azim = zRad2Deg( zVec3DElem( &v, zX ) );
    elev = zRad2Deg( zVec3DElem( &v, zY ) );
    tilt = zRad2Deg( zVec3DElem( &v, zZ ) );

    printf( "%f %f %f\n", azim, elev, tilt );
  } while( azim!=0 || elev!=0 || tilt!=0 );
  return 0;
}
