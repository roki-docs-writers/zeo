#include <zeo/zeo_vec3d.h>
#include <sys/time.h>

/* +++ nearest neighbor search: naive algorithm +++ */
zVec3D *zVec3DListNN(zVec3DList *list, zVec3D *v, double *dmin)
{
  zVec3DListCell *cell;
  double d;
  zVec3D *nn = NULL;

  *dmin = HUGE_VAL;
  zListForEach( list, cell )
    if( ( d = zVec3DSqrDist( cell->data, v ) ) < *dmin ){
      *dmin = d;
      nn = cell->data;
    }
  *dmin = sqrt( *dmin );
  return nn;
}

void output_list(FILE *fp, zVec3DList *list)
{
  zVec3DListCell *cell;

  zListForEach( list, cell )
    zVec3DDataFWrite( fp, cell->data );
}

/* +++ kd-tree search +++ */
void output_node(FILE *fp, zVecTree3D *node)
{
  /* face 1 */
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmin,zZ) );
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmax,zY), zVec3DElem(&node->vmin,zZ) );
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmax,zX), zVec3DElem(&node->vmax,zY), zVec3DElem(&node->vmin,zZ) );
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmax,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmin,zZ) );
  fprintf( fp, "%g %g %g\n\n", zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmin,zZ) );
  /* face 2 */
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmax,zZ) );
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmax,zY), zVec3DElem(&node->vmax,zZ) );
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmax,zX), zVec3DElem(&node->vmax,zY), zVec3DElem(&node->vmax,zZ) );
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmax,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmax,zZ) );
  fprintf( fp, "%g %g %g\n\n", zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmax,zZ) );
  /* pole x 4 */
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmin,zZ) );
  fprintf( fp, "%g %g %g\n\n", zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmax,zZ) );
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmax,zY), zVec3DElem(&node->vmin,zZ) );
  fprintf( fp, "%g %g %g\n\n", zVec3DElem(&node->vmin,zX), zVec3DElem(&node->vmax,zY), zVec3DElem(&node->vmax,zZ) );
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmax,zX), zVec3DElem(&node->vmax,zY), zVec3DElem(&node->vmin,zZ) );
  fprintf( fp, "%g %g %g\n\n", zVec3DElem(&node->vmax,zX), zVec3DElem(&node->vmax,zY), zVec3DElem(&node->vmax,zZ) );
  fprintf( fp, "%g %g %g\n",   zVec3DElem(&node->vmax,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmin,zZ) );
  fprintf( fp, "%g %g %g\n\n", zVec3DElem(&node->vmax,zX), zVec3DElem(&node->vmin,zY), zVec3DElem(&node->vmax,zZ) );
}

void output(FILE *fp, zVecTree3D *tree)
{
  output_node( fp, tree );
  if( tree->s[0] )
    output( fp, tree->s[0] );
  if( tree->s[1] )
    output( fp, tree->s[1] );
}


int deltatime(struct timeval *tv1, struct timeval *tv2)
{
  return (int)( (tv2->tv_sec-tv1->tv_sec)*1000000+tv2->tv_usec-tv1->tv_usec );
}

#define N 1000

int main(int argc, char *argv[])
{
  zVec3DList list; /* for comparison */
  zVecTree3D tree, *node;
  zVec3D v, *nn;
  int i, n;
  double dmin1, dmin2;
  FILE *fp;
  struct timeval tv1, tv2;

  n = argc > 1 ? atoi( argv[1] ) : N;
  zRandInit();
  zListInit( &list );
  zVecTree3DInit( &tree );
  zVec3DCreate( &tree.vmin,-10,-10,-10 );
  zVec3DCreate( &tree.vmax, 10, 10, 10 );

  fp = fopen( "src", "w" );
  for( i=0; i<n; i++ ){
    zVec3DCreate( &v, zRandF(-10,10), zRandF(-10,10), zRandF(-10,10) );
    zVec3DListInsert( &list, &v, true );
    zVecTree3DAdd( &tree, &v );
    zVec3DDataFWrite( fp, &v );
  }
  fclose( fp );

  fp = fopen( "par", "w" );
  output( fp, &tree );
  fclose( fp );

  zVec3DCreate( &v, zRandF(-10,10), zRandF(-10,10), zRandF(-10,10) );

  fp = fopen( "tst", "w" );
  node = zVecTree3DPart( &tree, &v );
  zVec3DDataFWrite( fp, &v );
  fprintf( fp, "\n" );
  output_node( fp, node );
  fclose( fp );

  gettimeofday( &tv1, NULL );
  dmin1 = zVecTree3DNN( &tree, &v, &node );
  gettimeofday( &tv2, NULL );
  eprintf( "kd-tree: %g - ", dmin1 ); zVec3DFWrite( stderr, &node->v );
  printf( "%d %d ", n, deltatime(&tv1,&tv2) );
  fp = fopen( "nn", "w" );
  zVec3DDataFWrite( fp, &v );
  zVec3DDataFWrite( fp, &node->v );
  fclose( fp );
  /* for comparison */
  gettimeofday( &tv1, NULL );
  nn = zVec3DListNN( &list, &v, &dmin2 );
  gettimeofday( &tv2, NULL );
  eprintf( "naive  : %g - ", dmin2 ); zVec3DFWrite( stderr, nn );
  printf( "%d\n", deltatime(&tv1,&tv2) );
  fp = fopen( "nnn", "w" );
  zVec3DDataFWrite( fp, &v );
  zVec3DDataFWrite( fp, nn );
  fclose( fp );

  if( !zIsTiny( dmin1-dmin2 ) )
    eprintf( "FAILED! %g/%g (err.=%g)\n", dmin1, dmin2, dmin1-dmin2 );

  zVec3DListDestroy( &list, true );
  zVecTree3DDestroy( &tree );
  return 0;
}
