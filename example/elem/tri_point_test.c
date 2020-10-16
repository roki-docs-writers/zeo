#include <zeo/zeo.h>

#define N  1000

void write_arrow(FILE *fp, zVec3D *p1, zVec3D *p2)
{
  fprintf( fp, "set arrow from %.10g, %.10g, %.10g to %.10g, %.10g, %.10g\n",
    zVec3DElem(p1,zX), zVec3DElem(p1,zY), zVec3DElem(p1,zZ),
    zVec3DElem(p2,zX), zVec3DElem(p2,zY), zVec3DElem(p2,zZ) );
}

int main(void)
{
  zVec3D v1, v2, v3, p, v;
  zTri3D t;
  int i;
  FILE *fp[2];

  zRandInit();
  zVec3DCreate( &v1, zRandF(0,1), zRandF(0,1), 0 );
  zVec3DCreate( &v2, zRandF(0,1), zRandF(0,1), 0 );
  zVec3DCreate( &v3, zRandF(0,1), zRandF(0,1), 0 );
  zTri3DCreate( &t, &v1, &v2, &v3 );
  fp[0] = fopen( "t", "w" );
  zVec3DDataFWrite( fp[0], &v1 );
  zVec3DDataFWrite( fp[0], &v2 );
  zVec3DDataFWrite( fp[0], &v3 );
  zVec3DDataFWrite( fp[0], &v1 );
  fclose( fp[0] );

  fp[0] = fopen( "d", "w" );
  fp[1] = fopen( "a", "w" );
  for( i=0; i<N; i++ ){
    zVec3DCreate( &p,  zRandF(0,1), zRandF(0,1), 0 );
    zVec3DElem(&p,zZ) = zTri3DClosest( &t, &p, &v );
    zVec3DDataFWrite( fp[0], &p );
    write_arrow( fp[1], &p, &v );
  }
  fclose( fp[0] );
  fclose( fp[1] );
  return 0;
}
