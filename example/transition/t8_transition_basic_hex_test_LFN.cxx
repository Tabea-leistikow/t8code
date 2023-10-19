/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/* This is step3 of the t8code tutorials.
 * After generating a coarse mesh (step1) and building a uniform forest
 * on it (step2), we will now adapt (= refine and coarsen) the forest
 * according to our own criterion.
 * 
 * The geometry (coarse mesh) is again a cube, this time modelled with
 * 6 tetrahedra, 6 prisms and 4 cubes.
 * We refine an element if its midpoint is whithin a sphere of given radius
 * around the point (0.5, 0.5, 1) and we coarsen outside of a given radius.
 * We will use non-recursive refinement, that means that the refinement level
 * of any element will change by at most +-1.
 * 
 * How you can experiment here:
 *   - Look at the paraview output files of the unifomr and the adapted forest.
 *     For the adapted forest you can apply a slice filter to look into the cube.
 *   - Run the program with different process numbers. You should see that refining is
 *     independent of the number of processes, but coarsening is not.
 *     This is due to the face that a family can only be coarsened if it is completely
 *     local to a single process and the distribution among the process may break this property.
 *   - Change the midpoint coordinates and the radii.
 *   - Change the adaptation criterion such that elements inside the sphere are coarsened
 *     and elements outside are refined.
 *   - Use t8_productionf to print the local number of elements on each process.
 *     Notice, that the uniform forest is evenly distributed, but that the adapted forest
 *     is not. This is due to the fact that we do not repartition our forest here.
 *   - Add a maximum refinement level to the adapt_data struct and use non-recursive refinement.
 *     Do not refine an element if it has reached the maximum level. (Hint: ts->t8_element_level)
 */

#include <t8.h>                 /* General t8code header, always include this. */
#include <t8_cmesh.h>           /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h> /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>        /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>     /* save forest */
#include <t8_forest/t8_forest_transition.h> 
#include <t8_forest/t8_forest_geometrical.h>    /* geometrical information of the forest */
#include <example/common/t8_example_common.h>
#include "t8_eclass.h"
#include <cmath>
#include <t8_vec.h>             /* Basic operations on 3D vectors. */
#include <t8_cmesh_vtk_writer.h> 
#include <t8_schemes/t8_default/t8_default_cxx.hxx> 
#define DO_TRANSITION_HEX_SCHEME 1

#if DO_TRANSITION_HEX_SCHEME 
#include <t8_schemes/t8_transition/t8_transition_conformal_hex/t8_transition_conformal_hex_cxx.hxx>
#include <t8_schemes/t8_transition/t8_transition_cxx.hxx>
#else
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#endif

T8_EXTERN_C_BEGIN ();

//------------------------adapt callback function -----------------
int
t8_adapt_callback (t8_forest_t forest,
                    t8_forest_t forest_from,
                    t8_locidx_t which_tree,
                    t8_locidx_t lelement_id,
                    t8_eclass_scheme_c *ts,
                    const int is_family,
                    const int num_elements, t8_element_t *elements[])
{
  if (((lelement_id) == 0) ) {
    /* Refine this element. */
    return 1;
  }
  /* Do not change this element. */
  return 0;
}


void
t8_print_vtk (t8_forest_t forest_adapt, char filename[BUFSIZ],
              int set_transition, int set_balance, int single_tree_mesh, int adaptation_count,
              t8_eclass_t eclass)
{
  if (set_transition) {
    if (single_tree_mesh)
      snprintf (filename, BUFSIZ, "forest_transitioned_LFN_hex_%i_%s",
                adaptation_count, t8_eclass_to_string[eclass]);
     
  if (set_balance) {
    if (single_tree_mesh)
      snprintf (filename, BUFSIZ, "forest_balanced_LFN_hex_%i_%s",
                adaptation_count, t8_eclass_to_string[eclass]);
  }
  t8_forest_write_vtk (forest_adapt, filename);
}
}

void t8_LFN_test(const t8_forest_t forest_adapt,
                 int adaptation_count, int num_adaptations){

    /* Collecting data of the adapted forest */
  t8_element_t       *current_element;
  t8_tree_t           current_tree;
  t8_locidx_t         forest_is_balanced = 1;
  t8_element_t      **neighbor_leafs;
  t8_locidx_t        *element_indices;
  t8_eclass_scheme_c *neigh_scheme;
  t8_eclass_t         eclass;
  t8_eclass_scheme_c *ts;

  int                *dual_faces;
  int                 num_neighbors;
  int                 face_id;
  int                 local_num_trees = t8_forest_get_num_local_trees (forest_adapt);   /* get the number of trees, this process knows about */
  int                 current_tree_num_elements;
  int                 subelement_count = 0;
  int                 LFN_call_count = 0;
  int                 tree_count;
  int                 elem_count;
  int                 neighbor_count;


  for (tree_count = 0; tree_count < local_num_trees; ++tree_count) {
    //eclass should be T8_ECLASS_HEX
    eclass = t8_forest_get_tree_class (forest_adapt, tree_count);
    //ts should be t8_scheme_new_transition_hex_cxx ()
    ts = t8_forest_get_eclass_scheme (forest_adapt, eclass);

    /* get the number of elements in the current tree */
    current_tree = t8_forest_get_tree (forest_adapt, tree_count);
    current_tree_num_elements =
      t8_forest_get_tree_element_count (current_tree);

    for (elem_count = 0; elem_count < current_tree_num_elements; ++elem_count) {

      /* determing the current element according to the given tree id and element id within the tree */
      current_element =
        t8_forest_get_element_in_tree (forest_adapt, tree_count, elem_count);

      if (ts->t8_element_is_subelement (current_element)) {
        subelement_count++;
      }



      for (face_id = 0; face_id < ts->t8_element_num_faces (current_element);
           ++face_id) {
      
        t8_forest_leaf_face_neighbors (forest_adapt, tree_count,
                                       current_element, &neighbor_leafs,
                                       face_id, &dual_faces, &num_neighbors,
                                       &element_indices, &neigh_scheme,forest_is_balanced);
  


        for( int neighbor_count = 0; neighbor_count < num_neighbors; neighbor_count++){
          t8_debugf ("\n_________"
                              "\nNeighbor: %i of %i at face %i: (dual face: %i | local index %i of %i (with ghosts)  | ghost, if >= %i):\n",
                              neighbor_count + 1, num_neighbors, face_id,
                              dual_faces[neighbor_count],
                              element_indices[neighbor_count],
                              t8_forest_get_local_num_elements (forest_adapt)
                              + t8_forest_get_num_ghosts (forest_adapt),
                              t8_forest_get_local_num_elements (forest_adapt)
                              - 1);                   
              ts->t8_element_debug_print (neighbor_leafs[neighbor_count]);                          
                                       }                                                                
        T8_FREE (element_indices);
        T8_FREE (neighbor_leafs);
        T8_FREE (dual_faces);

 
    } 
  }
}
                 }

void
t8_transition(void){

  t8_eclass_t         eclass = T8_ECLASS_HEX;

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];

  int level = 2;

  int num_adaptations = 3;

  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0,0,0);

  /* initialize a forest */
  t8_forest_init (&forest);

  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);

  t8_forest_set_level (forest, level);

 #if DO_TRANSITION_HEX_SCHEME
  t8_forest_set_scheme (forest, t8_scheme_new_transition_hex_cxx ());
#else
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
#endif


  t8_forest_commit (forest);


  int  adaptation_count;
  for (adaptation_count = 1; adaptation_count <= num_adaptations;
       ++adaptation_count) {

    t8_forest_init (&forest_adapt);

    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_callback, 0);

    t8_forest_set_transition (forest_adapt, forest, 0);

    t8_forest_commit (forest_adapt);    

    t8_LFN_test (forest_adapt, adaptation_count,
                  num_adaptations);
  
  /* Set forest to forest_adapt for the next step */
    forest = forest_adapt;
    }
  t8_print_vtk (forest, filename, 1, 1,1, 0,eclass);
  t8_forest_write_vtk(forest, "forest_adapt_LFN");
  t8_forest_unref (&forest);
}



int
main (int argc, char **argv)
{

  int                 mpiret;
  sc_MPI_Comm         comm;

  t8_locidx_t         local_num_trees;
  t8_gloidx_t         global_num_trees;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEBUG);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_DEBUG);
  comm = sc_MPI_COMM_WORLD;
  t8_transition();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();