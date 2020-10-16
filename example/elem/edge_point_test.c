#include <zeo/zeo_elem.h>

#define N  1000
#define NN  100

void write_arrow(FILE *fp, zVec3D *p1, zVec3D *p2)
{
  fprintf( fp, "set arrow from %.10g, %.10g, %.10g to %.10g, %.10g, %.10g\n",
    zVec3DElem(p1,zX), zVec3DElem(p1,zY), zVec3DElem(p1,zZ),
    zVec3DElem(p2,zX), zVec3DElem(p2,zY), zVec3DElem(p2,zZ) );
}

int main(void)
{
  zVec3D v1, v2, p, v;
  zEdge3D e;
  int i;
  FILE *fp[3];

  zRandInit();
  zVec3DCreate( &v1, zRandF(0,1), zRandF(0,1), 0 );
  zVec3DCreate( &v2, zRandF(0,1), zRandF(0,1), 0 );
  zEdge3DCreate( &e, &v1, &v2 );
  fp[0] = fopen( "e", "w" );
  zVec3DDataFWrite( fp[0], &v1 );
  zVec3DDataFWrite( fp[0], &v2 );
  fclose( fp[0] );

  fp[0] = fopen( "d", "w" );
  fp[1] = fopen( "l", "w" );
  fp[2] = fopen( "a", "w" );
  for( i=0; i<N; i++ ){
    zVec3DCreate( &p,  zRandF(0,1), zRandF(0,1), 0 );
    if( zEdge3DPointIsOn( &e, &p ) )
      zVec3DDataFWrite( fp[1], &p );
    zEdge3DClosest( &e, &p, &v );
    zVec3DElem(&p,zZ) = zEdge3DPointDist( &e, &p );
    zVec3DDataFWrite( fp[0], &p );
    write_arrow( fp[2], &p, &v );
  }
  for( i=0; i<NN; i++ ){
    zVec3DMul( zEdge3DVec(&e), zRandF(-1.5,1.5), &p );
    zVec3DAddDRC( &p, zEdge3DVert(&e,0) );
    if( zEdge3DPointIsOn( &e, &p ) )
      zVec3DDataFWrite( fp[1], &p );
    zEdge3DClosest( &e, &p, &v );
    zVec3DElem(&p,zZ) = zEdge3DPointDist( &e, &p );
    zVec3DDataFWrite( fp[0], &p );
    write_arrow( fp[2], &p, &v );
  }
  fclose( fp[0] );
  fclose( fp[1] );
  fclose( fp[2] );
  return 0;
}
