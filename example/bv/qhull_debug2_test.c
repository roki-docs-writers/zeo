#include <zeo/zeo_bv.h>

#define N 16
zVec3D v[N];
void vec_create(void)
{
  register int i;
  FILE *fp;

  /* This irregular case was provided by Toshiya Nishi. */
  zVec3DCreate( &v[0], 0.1203493299599752, 0.6944457444688849, 0.0000000000005107 );
  zVec3DCreate( &v[1], 0.1854702517066613, 0.6198785563037007, 0.0000000000005678 );
  zVec3DCreate( &v[2], 0.2397009340086134, 0.6672392266649270, 0.0000000000005698 );
  zVec3DCreate( &v[3], 0.1745800122619273, 0.7418064148301111, 0.0000000000005127 );
  zVec3DCreate( &v[4], 0.1203493299599715, 0.6944457444688886, 0.0092000000005107 );
  zVec3DCreate( &v[5], 0.1854702517066577, 0.6198785563037045, 0.0092000000005678 );
  zVec3DCreate( &v[6], 0.2397009340086098, 0.6672392266649307, 0.0092000000005698 );
  zVec3DCreate( &v[7], 0.1745800122619236, 0.7418064148301149, 0.0092000000005127 );
  zVec3DCreate( &v[8], 0.3669255481454933, 0.8322184531202354, 0.0000000000008041 );
  zVec3DCreate( &v[9], 0.4083570424248257, 0.7423049876004276, 0.0000000000008037 );
  zVec3DCreate( &v[10],0.3429654311376928, 0.7121729917609132, 0.0000000000008029 );
  zVec3DCreate( &v[11],0.3015339368583604, 0.8020864572807209, 0.0000000000008034 );
  zVec3DCreate( &v[12],0.3669255481454932, 0.8322184531202352, 0.0092000000008041 );
  zVec3DCreate( &v[13],0.4083570424248256, 0.7423049876004276, 0.0092000000008037 );
  zVec3DCreate( &v[14],0.3429654311376927, 0.7121729917609131, 0.0092000000008029 );
  zVec3DCreate( &v[15],0.3015339368583604, 0.8020864572807208, 0.0092000000008034 );

  fp = fopen( "src", "w" );
  for( i=0; i<N; i++ )
    zVec3DDataFWrite( fp, &v[i] );
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
