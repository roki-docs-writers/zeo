#include <zeo/zeo_nurbs3d.h>
#include <zm/zm_mat.h>

void adaption(zMat m, int x, int y, int xnum, int ynum, double z, int weight, int *cnt, double *work)
{
  if( x >= 0 && x < xnum && y >= 0 && y < ynum ){
    *cnt += weight;
    *work += weight * (zMatElem(m,x,y) - z);
  }
}

zNURBS3D *create_terra(zNURBS3D *ns, int xnum, int ynum, double gwidth, double zs, double zfin, int loop, int pointn)
{
  zMat md, ms;
  double alpha;
  register int c, p, i, j, n;
  int cnt;
  double work, z;

  zNURBS3DCreateGridXY( ns, 3, xnum, ynum, gwidth );

  md = zMatAlloc( xnum ,ynum );
  ms = zMatAlloc( xnum, ynum );
  zMatClear( md );
  zMatClear( ms );

  if( loop == 1 )
    alpha = log( zfin/zs );
  else
    alpha = log( zfin/zs ) / ( loop-1 );

  for( c=0; c<loop; c++ ){
    for( p=0; p<pointn*(c+1); p++ ){
      zMatElem(md,zRandI(0,xnum-1),zRandI(0,ynum-1)) += zs * exp( alpha*c );
      zMatElem(md,zRandI(0,xnum-1),zRandI(0,ynum-1)) -= zs * exp( alpha*c );
    }
    for( n=0; n<loop-c; n++ ){
      for( i=0; i<xnum; i++){
        for( j=0; j<ynum; j++ ){
          work = 0;
          cnt = 0;
          z = zMatElem(md,i,j);
          adaption( md, i-1, j-1, xnum, ynum, z, 1, &cnt, &work);
          adaption( md, i-1, j,   xnum, ynum, z, 2, &cnt, &work);
          adaption( md, i-1, j+1, xnum, ynum, z, 1, &cnt, &work);
          adaption( md, i,   j-1, xnum, ynum, z, 2, &cnt, &work);
          adaption( md, i,   j,   xnum, ynum, z, 4, &cnt, &work);
          adaption( md, i,   j+1, xnum, ynum, z, 2, &cnt, &work);
          adaption( md, i+1, j-1, xnum, ynum, z, 1, &cnt, &work);
          adaption( md, i+1, j,   xnum, ynum, z, 2, &cnt, &work);
          adaption( md, i+1, j+1, xnum, ynum, z, 1, &cnt, &work);

          zMatSetElem(ms, i, j, zMatElem(md,i,j) + work/cnt );
        }
      }
      zMatCopyNC( ms ,md );
    }
  }

  for( i=0; i<zNURBS3DCPSize(ns,0); i++ )
    for( j=0; j<zNURBS3DCPSize(ns,1); j++ ){
      zVec3DElem(zNURBS3DCP(ns,i,j),2) = zMatElem(md,i,j);
    }

  return ns;
}

#define DIV 200
void write4gnuplot(zNURBS3D *ns, char filename[])
{
  FILE *fp;
  register int i, j;
  double d1, d2;
  zVec3D v;

  fp = fopen( filename, "w" );
  d1 = ( zNURBS3DKnotE(ns,0) - zNURBS3DKnot0(ns,0) ) / DIV;
  d2 = ( zNURBS3DKnotE(ns,1) - zNURBS3DKnot0(ns,1) ) / DIV;
  for( i=0; i<DIV; i++ ){
    for( j=0; j<DIV; j++ ){
      zNURBS3DVec( ns, i*d1, j*d2, &v );
      fprintf( fp, "%f %f %f\n", zVec3DElem(&v,0), zVec3DElem(&v,1), zVec3DElem(&v,2) );
    }
  }
  fclose( fp );
}

int main(int argc, char *argv[])
{
  zNURBS3D ns;

  zNURBS3DInit( &ns );
  create_terra( &ns, 41, 61, 0.005, 0.5, 0.002, 20, 20 );

  zNURBS3DWriteFile( &ns, "nurbs.zn3d" );
  write4gnuplot( &ns, "nurbs_grid.dat" );


  zNURBS3DDestroy( &ns );
  return 0;
}

