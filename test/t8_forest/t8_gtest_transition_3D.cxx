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

/* In this test we test the t8_forest_set/get_user_data and
 * t8_forest_set/get_user_function functions.
 */

#include <gtest/gtest.h>
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
int transition_type;
    for (elem_count = 0; elem_count < current_tree_num_elements; ++elem_count) {

      /* determing the current element according to the given tree id and element id within the tree */
      current_element =
        t8_forest_get_element_in_tree (forest_adapt, tree_count, elem_count);
        transition_type = ts->t8_element_get_transition_type(current_element);

    if (elem_count >= 8 && elem_count < 17 ){
      EXPECT_EQ(transition_type, 32);
    }
    else if (elem_count >= 17 && elem_count < 26){
      EXPECT_EQ(transition_type, 8);
    }
    else if (elem_count >= 27 && elem_count < 36){
      EXPECT_EQ(transition_type, 2);
    }
    else{
      EXPECT_EQ(transition_type, 0);
    }
      if (ts->t8_element_is_subelement (current_element)) {

        subelement_count++;
      }



      for (face_id = 0; face_id < ts->t8_element_num_faces (current_element);
           ++face_id) {
      
        t8_forest_leaf_face_neighbors (forest_adapt, tree_count,
                                       current_element, &neighbor_leafs,
                                       face_id, &dual_faces, &num_neighbors,
                                       &element_indices, &neigh_scheme,forest_is_balanced);
  
        EXPECT_LE(num_neighbors, 2);

// Test for transition type 32
if ( (face_id == 0) && (elem_count == 8)){
EXPECT_EQ(dual_faces[0] , 0);
EXPECT_EQ(num_neighbors , 1);
EXPECT_EQ(ts->t8_element_get_subelement_id(neighbor_leafs[0]), 5);
}
if ( (face_id == 1) && (elem_count == 8)){
EXPECT_EQ(dual_faces[0] , 0);
EXPECT_EQ(num_neighbors , 1);
EXPECT_EQ(ts->t8_element_get_subelement_id(neighbor_leafs[0]), 1);
}
if ( (face_id == 2) && (elem_count == 8)){
EXPECT_EQ(dual_faces[0] , 0);
EXPECT_EQ(num_neighbors , 1);
EXPECT_EQ(ts->t8_element_get_subelement_id(neighbor_leafs[0]), 7);
}
if ( (face_id == 3) && (elem_count == 8)){
EXPECT_EQ(dual_faces[0] , 2);
EXPECT_EQ(num_neighbors , 1);
EXPECT_EQ(ts->t8_element_get_subelement_id(neighbor_leafs[0]), 2);
}

if ( (face_id == 0) && (elem_count == 13)){
  EXPECT_EQ(dual_faces[0] , 0);
  EXPECT_EQ(dual_faces[1] , 0);
  EXPECT_EQ(num_neighbors , 2);
  EXPECT_EQ(ts->t8_element_get_subelement_id(neighbor_leafs[0]), 0);
  EXPECT_EQ(ts->t8_element_get_subelement_id(neighbor_leafs[1]), 2);
}
                                                               
        T8_FREE (element_indices);
        T8_FREE (neighbor_leafs);
        T8_FREE (dual_faces);

 
    } 
  }
  EXPECT_EQ( subelement_count , 27);
}
                 }

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
t8_transition(void){

  t8_eclass_t         eclass = T8_ECLASS_HEX;

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];

  int level = 2;

  int num_adaptations = 1;

int mpi_rank; 
sc_MPI_Comm local_comm;

sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpi_rank);
sc_MPI_Comm_split (sc_MPI_COMM_WORLD, mpi_rank, mpi_rank, &local_comm);


  cmesh = t8_cmesh_new_hypercube (eclass, local_comm, 0,0,0);

  /* initialize a forest */
  t8_forest_init (&forest);

  t8_forest_set_cmesh (forest, cmesh, local_comm);

  t8_forest_set_level (forest, level);

 #if DO_TRANSITION_HEX_SCHEME
  t8_forest_set_scheme (forest, t8_scheme_new_transition_hex_cxx ());
#else
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
#endif



  t8_forest_commit (forest);
  // ASSERT_EQ(forest->mpisize, 2 );
  // if(forest->mpisize > 1){
  //   ASSERT_EQ(forest->mpisize, 2 );
  //    t8_forest_unref (&forest);
  //   return;
  // }


  int  adaptation_count;
  for (adaptation_count = 1; adaptation_count <= num_adaptations;
       ++adaptation_count) {

    t8_forest_init (&forest_adapt);

    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_callback, 0);

    t8_forest_set_ghost_ext (forest_adapt, 1 , T8_GHOST_FACES,
                               1);

    t8_forest_set_transition (forest_adapt, forest, 0);

    t8_forest_commit (forest_adapt);    

    t8_LFN_test (forest_adapt, adaptation_count,
                  num_adaptations);
  
  /* Set forest to forest_adapt for the next step */
    forest = forest_adapt;
    }
  t8_forest_unref (&forest);
  sc_MPI_Comm_free(&local_comm);
}


/* Test t8_forest_set/get_user_data.
 * We build a forest and set user data for it.
 * We then retrieve the data and check whether it is the same.
 */
TEST (transition_3D, test_transition_3D_LFN)
{
  t8_transition();


}



