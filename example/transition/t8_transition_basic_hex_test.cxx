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



/* TODO (JM): Copied this function from `t8_transition_local.cxx`. Adapt to your needs. */
//#ifdef T8_ENABLE_DEBUG
static int
t8_check_coordinates (double *coords)
{
  /* The initial hex_element is the unit hex with vertices (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1) ,(1,0,1), (0,1,1)  and (1,1,1).
   * We know that therefore, all children (even our subelements) will have vertices with coordinates 0, 0.5 or 1. */
  double              eps = 1e-126;     /* testing up to float precision */
  if ((fabs (coords[0] - 0.0) < eps || fabs (coords[0] - 0.5) < eps || fabs (coords[0] - 1.0) < eps) &&
      (fabs (coords[1] - 0.0) < eps || fabs (coords[1] - 0.5) < eps || fabs (coords[1] - 1.0) < eps) && 
      (fabs (coords[2] - 0.0) < eps || fabs (coords[2] - 0.5) < eps || fabs (coords[2] - 1.0) < eps)){
    return true;
  }
  return false;
}
int
t8_forest_transition_conformal_hex (t8_forest_t forest,
                                     t8_forest_t forest_from,
                                     t8_locidx_t ltree_id,
                                     t8_locidx_t lelement_id,
                                     t8_eclass_scheme_c *ts,
                                     const int is_family,
                                     int num_elements,
                                     t8_element_t *elements[])
{
  int                 iface, num_faces, neigh_face, transition_type = 0;
  t8_gloidx_t         neighbor_tree;
  t8_eclass_t         neigh_class;
  t8_eclass_scheme_c *neigh_scheme;
  t8_element_t       *element = elements[0], **face_neighbor;

  /* Hanging faces can only exist at non-maxlevel elements */
  if (forest_from->maxlevel_existing <= 0 ||
      ts->t8_element_level (element) < forest_from->maxlevel) {

    num_faces = ts->t8_element_num_faces (element);

    /* TODO: Update this comment to HEX. */
    /* We use a binary encoding (depending on the face enumeration), to determine which subelement type to use. 
     * Every face has a flag parameter, wich is set to 1, if there is a neighbor with a higher level 
     * and to 0, if the level of the neighbor is at most the level of the element.   
     *             
     *              f0                         1
     *        x - - x - - x              x - - x - - x       
     *        |           |              | \   |   / |
     *        |           |              |   \ | /   |                                                            | f3 | f2 | f1 | f0 |
     *    f3  x           | f2   -->   1 x - - x     | 0   -->   binary code (according to the face enumeration): |  1 |  0 |  0 |  1 | = 9 in base 10  
     *        |           |              |   /   \   |
     *        | elem      |              | /       \ |
     *        x - - - - - x              x - - - - - x
     *              f1                         0 
     *                      
     * Note, that this procedure is independent of the eclass (we only show an example for the quad scheme). 
     * Each neighbor-structure will lead to a unique binary code. 
     * Within the element scheme of the given eclass, this binary code is used to construct the right subelement type,
     * in order to remove hanging nodes from the mesh. */

    for (iface = 0; iface < num_faces; iface++) {
      /* Get the element class and scheme of the face neighbor */
      neigh_class = t8_forest_element_neighbor_eclass (forest_from,
                                                       ltree_id, element,
                                                       iface);

      neigh_scheme = t8_forest_get_eclass_scheme (forest_from, neigh_class);

      /* Allocate memory for the virtual face neighbor */
      // t8_element_t Array mit einem Element
      face_neighbor = T8_ALLOC (t8_element_t *, 1);

      neigh_scheme->t8_element_new (1, face_neighbor);

      /* Compute the virtual face neighbor of element at this face */
      neighbor_tree = t8_forest_element_face_neighbor (forest_from, ltree_id,
                                                       element,
                                                       face_neighbor[0],
                                                       neigh_scheme,
                                                       iface, &neigh_face);

      /* TODO: Update this code block for hex / sub-pyramids. */
      // if (neighbor_tree >= 0) {
      //   if (t8_forest_element_has_leaf_desc (forest_from, neighbor_tree,
      //                                        face_neighbor[0],
      //                                        neigh_scheme)) {
      //     /* Compute transition type as the decimal represenation of the binary concatenation */
      //     transition_type += 1 << ((num_faces - 1) - iface);
      //   }
      // }
      /* clean-up */
      neigh_scheme->t8_element_destroy (1, face_neighbor);
      T8_FREE (face_neighbor);
    }

    /* returning the right subelement types */
    if (transition_type == 0) { /* no hanging faces in this case */
      return 0;
    }
    else if (transition_type == 63) {   /* Six hanging faces in this case */
      return 1;
    }
    else {                      /* use a transition cell of subelements and add 1 to every type, to avoid refine = 1 */
      return transition_type + 1;
    }
  }
  return 0;                     /* if elem has maxlevel then keep it unchanged since there will never be hanging faces */
}
static t8_cmesh_t
t8_build_hex_coarse_mesh (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh; 
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_HEX, comm, 0, 0, 0);

  return cmesh;
}




//----------------vtk functions------------------------------------------- 
  static void
t8_write_cmesh_vtk (t8_cmesh_t cmesh, const char *prefix)
{
  t8_cmesh_vtk_write_file (cmesh, prefix, 1.0);
}

static void
t8_destroy_cmesh (t8_cmesh_t cmesh)
{
  t8_cmesh_destroy (&cmesh);
}

void
t8_print_forest_information (t8_forest_t forest)
{
  t8_locidx_t         local_num_elements;
  t8_gloidx_t         global_num_elements;

  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the local number of elements. */
  local_num_elements = t8_forest_get_local_num_elements (forest);
  /* Get the global number of elements. */
  global_num_elements = t8_forest_get_global_num_elements (forest);
  t8_global_productionf ("Local number of elements:\t\t%i\n",
                         local_num_elements);
  t8_global_productionf ("Global number of elements:\t%li\n",
                         global_num_elements);
}
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
  int test = t8_forest_transition_conformal_hex ( forest,
                                      forest_from,
                                      which_tree,
                                      lelement_id,
                                     ts, is_family,
                                      num_elements,
                                     elements);

  if ((which_tree == 0 )&& (lelement_id == 0)) {
    t8_debugf("ich werde verfeinert \n");
    /* Refine this element. */
    return 1;
  }
  /* Do not change this element. */
  return 0;
}

//------------------------adapt function -----------------------
t8_forest_t t8_adapt_forest (t8_forest_t forest)
{
  t8_forest_t         forest_adapt;

  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  forest_adapt =
    t8_forest_new_adapt (forest, t8_adapt_callback, 0, 0, NULL);

  return forest_adapt;
}
//-------------------------transition function -----------------
//-------------------------transition------------------------------------
void
t8_transition(void)
{
  t8_debugf ("~~~~~~~~~~ Into the t8_transition function ~~~~~~~~~~\n");
  t8_eclass_t         eclass = T8_ECLASS_HEX;  /* depending on the include file, this will be the transitioned or default hex implementation */

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];
  t8_eclass_scheme_c *ts; 

  int level = 1;
  int subelement_counter = 0;
  int elem_count;
  int tree_count = 0;
  int current_tree_num_elements;
  t8_tree_t           current_tree;


  int                 set_balance = 1;
  int                 set_transition = 1;

  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0,0,0);

    /* initialize a forest */
  t8_forest_init (&forest);

  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_level (forest, level);
  t8_forest_set_scheme (forest, t8_scheme_new_transition_hex_cxx ());
  t8_forest_commit (forest);

  t8_forest_init (&forest_adapt);

  t8_forest_set_adapt (forest_adapt, forest, t8_adapt_callback, 0);

    if (set_balance && !set_transition) {
      t8_forest_set_balance (forest_adapt, forest, 0);
    }
    if (set_transition) {
      t8_forest_set_transition (forest_adapt, forest, set_balance);
    }

  t8_forest_commit (forest_adapt);
  t8_write_cmesh_vtk (cmesh, "TEST");
  t8_forest_write_vtk(forest_adapt, "test_adapt");
  ts = t8_forest_get_eclass_scheme (forest_adapt, eclass);  

  current_tree = t8_forest_get_tree (forest_adapt, tree_count);

  current_tree_num_elements = t8_forest_get_tree_element_count (current_tree);

  for (elem_count = 0; elem_count < current_tree_num_elements; ++elem_count) {
    t8_element_t *current_element = t8_forest_get_element_in_tree (forest_adapt, tree_count, elem_count);
//       //Subelements counter 
       if (ts->t8_element_is_subelement (current_element)) {
         subelement_counter++;
         int sub_id = ts->t8_element_get_subelement_id(current_element);
         t8_productionf("sub ID %i\n", sub_id);
         }
  }
//       } t8_forest_get_element_in_tree (forest_adapt, 0, 0);

  int test1 = t8_forest_get_global_num_elements(forest_adapt);

  t8_productionf("subelement counter  %i \n", subelement_counter);
  t8_forest_unref (&forest_adapt);
}




int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  /* The prefix for our output files. */
  const char          prefix[BUFSIZ] = "t8_cmesh";
  const char         *prefix_uniform = "t8_uniform_forest";
   const char         *prefix_adapt= "t8_forest_adapt";
  t8_locidx_t         local_num_trees;
  t8_gloidx_t         global_num_trees;
  int level = 1;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);
  comm = sc_MPI_COMM_WORLD;
  t8_transition();

//   /* Build the coarse mesh */
//   cmesh = t8_build_hex_coarse_mesh (sc_MPI_COMM_WORLD);
//   /* Compute local and global number of trees. */
//   local_num_trees = t8_cmesh_get_num_local_trees (cmesh);
//   global_num_trees = t8_cmesh_get_num_trees (cmesh);

//   t8_global_productionf (" Local number of trees:\t%i\n",
//                          local_num_trees);
//   t8_global_productionf (" Global number of trees:\t%li\n",
//                          global_num_trees);
//   t8_write_cmesh_vtk (cmesh, prefix);
//   t8_global_productionf (" Wrote coarse mesh to vtu files: %s*\n",
//                          prefix);
// //---------------------create forest-------------------------
//   forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx(), level, 0,
//                            comm);
//   int maxlevel = t8_forest_get_maxlevel(forest);
//   t8_global_productionf (" Created uniform forest.\n");
//     t8_global_productionf (" Max level forest %i.\n", maxlevel);

//   /* Write forest to vtu files. */
//   t8_forest_write_vtk (forest, prefix_uniform);
//   t8_global_productionf ("Wrote uniform forest to vtu files: %s*\n",
//                          prefix_uniform);
// //-------------adapt the forest -----------------
// forest = t8_adapt_forest (forest);

//  /* Print information of our new forest. */
//   t8_global_productionf ("Adapted forest.\n");

//   /* Write forest to vtu files. */
//   t8_forest_write_vtk (forest, prefix_adapt);
//   t8_global_productionf (" Wrote adapted forest to vtu files: %s*\n",
//                          prefix_adapt);
                

// //-------------transition---------------------------------



// //-----------------------Destroy the forest ----------------------
//   t8_forest_unref (&forest);
//   t8_global_productionf ("Destroyed forest.\n");



  // t8_destroy_cmesh (cmesh);
  // t8_global_productionf (" Destroyed coarse mesh.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}





// int
// main (int argc, char **argv)
// {
//   int                 mpiret;
//   sc_MPI_Comm         comm;
//   t8_cmesh_t          cmesh;
//   t8_forest_t         forest;


//   /* Initialize MPI. This has to happen before we initialize sc or t8code. */
//   mpiret = sc_MPI_Init (&argc, &argv);
//   /* Error check the MPI return value. */
//   SC_CHECK_MPI (mpiret);

//   /* Initialize the sc library, has to happen before we initialize t8code. */
//   sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
//   /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
//   t8_init (SC_LP_PRODUCTION);


//   /* We will use MPI_COMM_WORLD as a communicator. */
//   comm = sc_MPI_COMM_WORLD;

//   /*
//    * Setup.
//    * Build cmesh and uniform forest.
//    */

//   /* Build a cube cmesh with tet, hex, and prism trees. */
//   cmesh = t8_cmesh_new_from_class (T8_ECLASS_HEX, comm);
//   const int           level = 3;
//   forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0, comm);


//   /*
//    *  Adapt the forest.
//    */
//   forest = t8_forest_new_adapt (forest, t8_adapt_callback, 0, 0, NULL);
//   t8_forest_write_vtk (forest, "t8_transition_basic_hex");

//   /*
//    * clean-up
//    */

//   /* Destroy the forest. */
//   t8_forest_unref (&forest);
//   sc_finalize ();

//   mpiret = sc_MPI_Finalize ();
//   SC_CHECK_MPI (mpiret);

//   return 0;
// }

// T8_EXTERN_C_END ();
