/* Zeo - Z/Geometry and optics computation library.
 * Copyright (C) 2005 Tomomichi Sugihara (Zhidao)
 *
 * zeo_vec_list - 3D vector list.
 */

#ifndef __ZEO_VEC_LIST_H__
#define __ZEO_VEC_LIST_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zVec3DList
 * 3D vector list
 * ********************************************************** */

zListClass( zVec3DList, zVec3DListCell, zVec3D* );

/*! \brief insert, create and destroy 3D vector list.
 *
 * zVec3DListInsert() inserts a new 3D vector list cell \a v
 * at the head of a vector list \a list.
 *
 * zVec3DListCreate() creates a vector list \a list from an
 * array of vectors \a v. \a num is the number of vectors in \a v.
 *
 * For these two functions, each cell inserted will have a copy
 * of \a v which is newly allocated if \a flag is the true value.
 * Otherwise, the cell has a pointer to \a v.
 *
 * zVec3DListDestroy() destroys \a list, freeing all cells.
 * If \a flag is the true value, it also frees the data held by
 * each cell. Otherwise, only the cells are freed.
 * \return
 * zVec3DListInsert() returns a pointer to the cell newly
 * inserted.
 * zVec3DListCreate() returns a pointer \a list.
 * zVec3DListDestroy() returns no value.
 */
__EXPORT zVec3DListCell *zVec3DListInsert(zVec3DList *list, zVec3D *v, bool flag);
__EXPORT zVec3DList *zVec3DListCreate(zVec3DList *list, zVec3D v[], int num, bool flag);
__EXPORT void zVec3DListDestroy(zVec3DList *list, bool flag);

/*! \brief a quick sort routine for vector list class.
 *
 * zVec3DListQuickSort() is a quick sort routine for zVec3DList
 * class.
 *
 * The cells of \a list will be sorted in ascending order
 * according to the comparison function \a cmp.
 * (The factor a in \a list is put after another factor b when
 * \a cmp(a,b,p) > 0, where p is for programmer's utility, given
 * by \a priv.)
 * \return
 * zVec3DListQuickSort() returns a pointer \a list.
 * \sa
 * zListQuickSortDef
 */
__EXPORT zVec3DList *zVec3DListQuickSort(zVec3DList *list, int (*cmp)(void*,void*,void*), void *priv);

/*! \brief output of 3D vector list.
 *
 * zVec3DListFWrite() writes the 3D vector list \a list to the
 * current position of the file \a fp in the following style.
 * n
 *  ( x1, y1, z1 )
 *   ...
 *  ( xn, yn, zn )
 * zVec3DWrite() writes \a list to the standard output.
 * \return
 * zVec3DListFWrite() and zVec3DListWrite() return no value.
 */
__EXPORT void zVec3DListFWrite(FILE *fp, zVec3DList *list);
#define zVec3DListWrite(l) zVec3DListFWrite( stdout, (l) )

__END_DECLS

#endif /* __ZEO_VEC_LIST_H__ */
