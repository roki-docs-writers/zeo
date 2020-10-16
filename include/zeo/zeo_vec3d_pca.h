/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_vec_pca - principal component analysis.
 */

#ifndef __ZEO_VEC_PCA_H__
#define __ZEO_VEC_PCA_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* utilities
 * ********************************************************** */

/* METHOD:
 * zVec3DBarycenterPL, zVec3DBarycenter, zVec3DPCA_PL, zVec3DPCA,
 * zVec3DBaryPCA_PL, zVec3DBaryPCA
 * - barycenter of and PCA against vector cloud.
 *
 * 'zVec3DBarycenterPL()' and 'zVec3DBarycenter()'
 * computes the barycenter of a set of vectors. For
 * 'zVec3DBarycenterPL()', vectors are given by a list
 * 'vl', while given by an array 'v' for 'zVec3DBarycenter()'.
 * 'num' is the number of vectors in 'v'. The result
 * vector is put where pointed by 'c'.
 * #
 * 'zVec3DPCA_PL()' and 'zVec3DPCA()' examines principal
 * component analysis (PCA) for a set of vectors.
 * Each of the three principal components passes through
 * the original point. The vectors are also given by
 * 'vl' for 'zVec3DPCA_PL()', and 'v' and 'num' for
 * 'zVec3DPCA()', respectively. The result PCs are stored
 * into the array 'evec'.
 * #
 * 'zVec3DBaryPCA_PL()' and 'zVec3DBaryPCA()' are
 * combinations of the above functions. Each of the
 * three PCs 'evec' passes through the barycenter 'c'
 * of the vectors. The vectors are given by 'vl' for
 * 'zVec3DBaryPCA_PL()', and 'v' and 'num' for
 * 'zVec3DBaryPCA()', respectively.
 * [RETURN VALUE]
 * 'zVec3DBarycenterPL()', 'zVec3DBarycenter()',
 * 'zVec3DBaryPCA_PL()' and 'zVec3DBaryPCA()' return a
 * pointer 'c'.
 * #
 * 'zVec3DPCA_PL()' and 'zVec3DPCA()' return a pointer
 * * to the head of 'evec'.
 */
__EXPORT zVec3D *zVec3DBarycenterPL(zVec3DList *vl, zVec3D *c);
__EXPORT zVec3D *zVec3DBarycenter(zVec3D v[], int num, zVec3D *c);
__EXPORT zVec3D *zVec3DPCA_PL(zVec3DList *vl, zVec3D evec[]);
__EXPORT zVec3D *zVec3DPCA(zVec3D v[], int num, zVec3D evec[]);
__EXPORT zVec3D *zVec3DBaryPCA_PL(zVec3DList *vl, zVec3D *c, zVec3D evec[]);
__EXPORT zVec3D *zVec3DBaryPCA(zVec3D v[], int num, zVec3D *c, zVec3D evec[]);

__END_DECLS

#endif /* __ZEO_VEC_PCA_H__ */
