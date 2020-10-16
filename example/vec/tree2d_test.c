#include <zeo/zeo_vec3d.h>

/* +++ nearest neighbor search: naive algorithm +++ */
zVec3D *zVec3DListNN(zVec3DList *list, zVec3D *v)
{
  zVec3DListCell *cell;
  double d, dmin = HUGE_VAL;
  zVec3D *nn = NULL;

  zListForEach( list, cell )
    if( ( d = zVec3DSqrDist( cell->data, v ) ) < dmin ){
      dmin = d;
      nn = cell->data;
    }
  return nn;
}

void output_list(FILE *fp, zVec3DList *list)
{
  zVec3DListCell *cell;

  zListForEach( list, cell )
    zVec3DDataFWrite( fp, cell->data );
}

/* +++ kd-tree search +++ */
void output_node(FILE *fp, zVecTree3D *node) /* dummy */
{
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmin,zZ) );
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmax,zY), zVec3DElem(&node->vmin,zZ) );
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmax,zX), zVec3DElem(&node->vmax,zY), zVec3DElem(&node->vmin,zZ) );
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmax,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmin,zZ) );
  fprintf( fp, "%g %g %g\n\n", zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmin,zZ) );
}

void output(FILE *fp, zVecTree3D *tree)
{
  output_node( fp, tree );
  if( tree->s[0] )
    output( fp, tree->s[0] );
  if( tree->s[1] )
    output( fp, tree->s[1] );
}


#define N 100

int main(int argc, char *argv[])
{
  zVec3DList list; /* for comparison */
  zVecTree3D tree, *node;
  zVec3D v, *nn;
  int i;
  FILE *fp;

  zRandInit();
  zListInit( &list );
  zVecTree3DInit( &tree );
  zVec3DCreate( &tree.vmin,-10,-10, 0 );
  zVec3DCreate( &tree.vmax, 10, 10, 0 );
  fp = fopen( "src", "w" );
  for( i=0; i<N; i++ ){
    zVec3DCreate( &v, zRandF(-10,10), zRandF(-10,10), 0 );
    zVec3DListInsert( &list, &v, true );
    zVecTree3DAdd( &tree, &v );
    zVec3DDataFWrite( fp, &v );
  }
  fclose( fp );

  fp = fopen( "par", "w" );
  output( fp, &tree );
  fclose( fp );

  zVec3DCreate( &v, zRandF(-5,5), zRandF(-5,5), 0 );
  fp = fopen( "tst", "w" );
  node = zVecTree3DPart( &tree, &v );
  zVec3DDataFWrite( fp, &v );
  fprintf( fp, "\n" );
  output_node( fp, node );
  fclose( fp );

  zVecTree3DNN( &tree, &v, &node );
  printf( "kd-tree: " ); zVec3DWrite( &node->v );
  fp = fopen( "nn", "w" );
  zVec3DDataFWrite( fp, &v );
  zVec3DDataFWrite( fp, &node->v );
  fclose( fp );
  /* for comparison */
  nn = zVec3DListNN( &list, &v );
  printf( "naive  : " ); zVec3DWrite( nn );
  fp = fopen( "nnn", "w" );
  zVec3DDataFWrite( fp, &v );
  zVec3DDataFWrite( fp, nn );
  fclose( fp );

  zVec3DListDestroy( &list, true );
  zVecTree3DDestroy( &tree );
  return 0;
}
