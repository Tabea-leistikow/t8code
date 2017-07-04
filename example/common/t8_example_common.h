/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** file t8_vtk.h
 * This header file collects macros that are needed for
 * the forest and cmesh vtk routines.
 * \see t8_forest_vtk.h \see t8_cmesh_vtk.h
 */

#ifndef T8_EXAMPLE_COMMON_H
#define T8_EXAMPLE_COMMON_H

#include <t8.h>

/** A levelset function in 3 space dimensions. */
typedef double      (*t8_example_level_set_fn) (double, double, double,
                                                void *);

/** Struct to handle refinement around a level-set function. */
typedef struct
{
  t8_example_level_set_fn L;  /**< The level set function. */
  void               *udata; /**< Data pointer that is passed to L */
  double              band_width; /**< Width of max_level elements around the zero-level set */
  int                 min_level; /**< The minimal refinement level. Elements with this level will not be coarsened. */
  int                 max_level; /**< The maximum refinement level. Elements with this level will not be refined. */
} t8_example_level_set_struct_t;

T8_EXTERN_C_BEGIN ();

/* function declarations */

/** Adapt a forest such that always the second child of the first
 * tree is refined and no other elements. This results in a highly
 * imbalanced forest.
 */
int                 t8_common_adapt_balance (t8_forest_t forest,
                                             t8_forest_t forest_from,
                                             t8_locidx_t which_tree,
                                             t8_eclass_scheme_c * ts,
                                             int num_elements,
                                             t8_element_t * elements[]);

/** Adapt a forest along a given level-set function.
 * The user data of forest must be a pointer to a \a t8_example_level_set_struct_t.
 * An element in the forest is refined, if it is in a band of \a band_with many
 * \a max_level elements around the zero level-set Gamma = { x | L(x) = 0}
 */
/* TODO: Currently the band_width control is not working yet. */
int                 t8_common_adapt_level_set (t8_forest_t forest,
                                               t8_forest_t forest_from,
                                               t8_locidx_t which_tree,
                                               t8_eclass_scheme_c * ts,
                                               int num_elements,
                                               t8_element_t * elements[]);

/** Compute the coordinates of the midpoint of an element.
 * \param [in]  forest  The forest in which the element is in (must be committed).
 * \param [in]  which_tree The local tree id of tree in which the element is in.
 * \param [in]  ts      The eclass scheme associated to the element.
 * \param [in]  element The element.
 * \param [in,out] elem_midpoint_f An array of 3 doubles. On output the coordinates
 *              of the midpoint of \a element are stored.
 * \note \a forest must be committed before calling this function.
 */
void                t8_common_midpoint (t8_forest_t forest,
                                        t8_locidx_t which_tree,
                                        t8_eclass_scheme_c * ts,
                                        t8_element_t * element,
                                        double elem_midpoint_f[3]);

T8_EXTERN_C_END ();

#endif /* !T8_EXAMPLE_COMMON_H */
