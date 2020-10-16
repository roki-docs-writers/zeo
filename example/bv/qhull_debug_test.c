#include <zeo/zeo_bv.h>

#define N 27
zVec3D v[N];
void vec_create(void)
{
  register int i;
  FILE *fp;

  fp = fopen( "src", "w" );
  for( i=0; i<N; i++ ){
    zVec3DCreate( &v[i], i%3-1, (i/3)%3-1, (i/9)%3-1 );
    zVec3DDataFWrite( fp, &v[i] );
  }
  fclose( fp );
}

void output(zPH3D *ph)
{
  register int i;
  FILE *fp;

  eprintf( "%d vertices, %d faces.\n", zPH3DVertNum(ph), zPH3DFaceNum(ph) );
  fp = fopen( "ch", "w" );
  for( i=0; i<zPH3DFaceNum(ph); i++ ){
    zVec3DDataFWrite( fp, zPH3DFaceVert(ph,i,0) );
    zVec3DDataFWrite( fp, zPH3DFaceVert(ph,i,1) );
    fprintf( fp, "\n" );
    zVec3DDataFWrite( fp, zPH3DFaceVert(ph,i,2) );
    zVec3DDataFWrite( fp, zPH3DFaceVert(ph,i,2) );
    fprintf( fp, "\n\n" );
  }
  fclose( fp );
  fp = fopen( "chv", "w" );
  for( i=0; i<zPH3DVertNum(ph); i++ )
    zVec3DDataFWrite( fp, zPH3DVert(ph,i) );
  fclose( fp );
}

int main(void)
{
  zPH3D ch;

  vec_create();
  zCH3D( &ch, v, N );
  output( &ch );
  zPH3DDestroy( &ch );
  return 0;
}
