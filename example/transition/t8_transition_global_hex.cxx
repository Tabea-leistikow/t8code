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

/* Description:
 * This is the example file for refinement with transitioning. In this testcase, we are able to
 *     (i)   refine a mesh according to some refinement criterion and use transiton cells to make the mesh conformal
 *     (ii)  use multiple adaptation steps in which the refinement criterion changes (e.g. the geometry)
 *     (iii) decide, whether we want to check the LFN function for each mesh
 *     (iv)  decide, whether we want to get statistics printed out, regarding # of elements in the meshes and runtime infos of the several functions or other debugging information
 */

/* to switch between the default quad scheme and the transition implementation */
#include "t8_eclass.h"
#define DO_TRANSITION_HEX_SCHEME 1

#include "t8_forest/t8_forest_general.h"
#include <cstring>
#if DO_TRANSITION_HEX_SCHEME 
#include <t8_schemes/t8_transition/t8_transition_conformal_hex/t8_transition_conformal_hex_cxx.hxx>
#include <t8_schemes/t8_transition/t8_transition_cxx.hxx>
#else
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#endif
#include <t8_forest/t8_forest_io.h>     // to write vtk
#include <t8_vec.h>
#include <example/common/t8_example_common.h>
#include <t8_cmesh/t8_cmesh_examples.h> /* for cmesh initialization via for example t8_cmesh_new_hypercube */
#include <t8_cmesh_vtk_writer.h> 
/* In this example, a simple refinement criteria is used to construct an adapted and transitioned forest. 
 * Afterwards, we iterate through all elements and all faces of the this forest in order to test the leaf_face_neighbor function that will determine all neighbor elements. */

// typedef struct
// {
//   double              mid_point[3];
//   double              radius;
// } t8_basic_sphere_data_t;

/* Compute the distance to a sphere around a mid_point with given radius. */
// static double
// t8_distance_to_sphere (const double x[3], double t, void *data)
// {
//   t8_basic_sphere_data_t *sdata = (t8_basic_sphere_data_t *) data;
//   double             *M = sdata->mid_point;

//   return t8_vec_dist (M, x) - sdata->radius;
// }

// void
// t8_print_general_stats (double commit_time_total, int num_adaptations,
//                         int global_num_elements_accum, double total_time,
//                         double LFN_time_accum)
// {
//   t8_productionf
//     ("\n|++++++++++++++++++++++++ Commit statistics | total +++++++++++++++++++++++++++|\n"
//      "|    Average #elements:       %i\n"
//      "|    Time total [s]:          %.3f (%.2f %%)\n"
//      "|    Step time average [s]:   %.3f\n"
//      "|    LFN time total [s]:      %.3f (%.2f %%)\n"
//      "|    LFN time average [s]:    %.3f\n"
//      "|    Commit time total [s]:   %.3f (%.2f %%)\n"
//      "|    Commit time average [s]: %.3f\n"
//      "|    Rest [s]:                %.3f (%.2f %%)\n"
//      "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|\n\n",
//      global_num_elements_accum / num_adaptations, total_time, 100.,
//      total_time / (double) num_adaptations, LFN_time_accum,
//      100. * LFN_time_accum / total_time,
//      LFN_time_accum / (double) num_adaptations, commit_time_total,
//      100. * commit_time_total / total_time,
//      commit_time_total / (double) num_adaptations,
//      total_time - LFN_time_accum - commit_time_total,
//      100. * (total_time - LFN_time_accum - commit_time_total) / total_time);
// }

// void
// t8_print_commit_stats (double commit_time, int num_adaptations,
//                        int adaptation_count)
// {
//   t8_productionf
//     ("\n|++++++++++++++++++++ Commit statistics | adaptation %i of %i +++++++++++++++++++|\n"
//      "|    Commit time total [s]: %.9f\n"
//      "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|\n\n",
//      adaptation_count, num_adaptations, commit_time);
// }

// void
// t8_print_LFN_stats (int global_num_elements, int global_num_subelements,
//                     int local_num_elements, int local_num_subelements,
//                     int LFN_call_count, double time_LFN,
//                     double time_LFN_per_call, int adaptation_count,
//                     int num_adaptations)
// {
//   t8_productionf
//     ("\n|+++++++++++++++++++++ LFN statistics | adaptation %i of %i +++++++++++++++++++++|\n"
//      "|    Global #elements:         %i (#quads: %i, #subelements: %i)\n"
//      "|    Local #elements:          %i (#quads: %i, #subelements: %i)\n"
//      "|    #LFN calls:               %i\n"
//      "|    LFN runtime total [s]:    %f\n"
//      "|    LFN runtime per call [s]: %.9f\n"
//      "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|\n\n",
//      adaptation_count, num_adaptations, global_num_elements,
//      global_num_elements - global_num_subelements, global_num_subelements,
//      local_num_elements, local_num_elements - local_num_subelements,
//      local_num_subelements, LFN_call_count, time_LFN, time_LFN_per_call);
// }

void
t8_print_vtk (t8_forest_t forest_adapt, char filename[BUFSIZ],
              int set_transition, int set_balance, int single_tree_mesh,
              int multiple_tree_mesh, int adaptation_count,
              t8_eclass_t eclass)
{
  // if (set_transition) {
  //   if (single_tree_mesh)
  //     snprintf (filename, BUFSIZ, "forest_transitioned_hex_%i_%s",
  //               adaptation_count, t8_eclass_to_string[eclass]);
  //   else if (multiple_tree_mesh)
  //     snprintf (filename, BUFSIZ, "forest_transitioned_hex_%i_%s",
  //               adaptation_count, t8_eclass_to_string[T8_ECLASS_HEX]);
  //   else
  //     snprintf (filename, BUFSIZ, "forest_transitioned_hex_%i_hybrid",
  //               adaptation_count);
  // }
  // else if (set_balance) {
  //   if (single_tree_mesh)
  //     snprintf (filename, BUFSIZ, "forest_balanced_hex_%i_%s",
  //               adaptation_count, t8_eclass_to_string[eclass]);
  //   else if (multiple_tree_mesh)
  //     snprintf (filename, BUFSIZ, "forest_balanced_hex_%i_%s",
  //               adaptation_count, t8_eclass_to_string[T8_ECLASS_HEX]);
  //   else
  //     snprintf (filename, BUFSIZ, "forest_balanced_hex_%i_hybrid",
  //               adaptation_count);
  // }
  // else {
  //   if (single_tree_mesh)
  //     snprintf (filename, BUFSIZ, "forest_adapted_hex_%i_%s",
  //               adaptation_count, t8_eclass_to_string[eclass]);
  //   else if (multiple_tree_mesh)
  //     snprintf (filename, BUFSIZ, "forest_adapted_hex_%i_%s",
  //               adaptation_count, t8_eclass_to_string[T8_ECLASS_HEX]);
  //   else
  //     snprintf (filename, BUFSIZ, "forest_adapted_hex_%i_hybrid",
  //               adaptation_count);;
  // }
  t8_forest_write_vtk (forest_adapt, filename);
}

/* Compute neighbors of all elements in all trees at all faces */
void
t8_LFN_test (const t8_forest_t forest_adapt, int get_LFN_stats,
             int adaptation_count, int num_adaptations, int get_LFN_elem_info)
{
//   t8_debugf ("~~~~~~~~~~ Into the LFN test fucntion. ~~~~~~~~~~\n");

//   /* Collecting data of the adapted forest */
//   t8_element_t       *current_element;
//   t8_tree_t           current_tree;
//   t8_locidx_t         forest_is_balanced = 1;
//   t8_element_t      **neighbor_leafs;
//   t8_locidx_t        *element_indices;
//   t8_eclass_scheme_c *neigh_scheme;
//   t8_eclass_t         eclass;
//   t8_eclass_scheme_c *ts;

//   int                *dual_faces;
//   int                 num_neighbors;
//   int                 face_id;
//   int                 local_num_trees = t8_forest_get_num_local_trees (forest_adapt);   /* get the number of trees, this process knows about */
//   int                 current_tree_num_elements;
//   int                 subelement_count = 0;
//   int                 LFN_call_count = 0;
//   int                 tree_count;
//   int                 elem_count;
//   int                 neighbor_count;

//   double              time_LFN = 0;

//   for (tree_count = 0; tree_count < local_num_trees; ++tree_count) {
//     eclass = t8_forest_get_tree_class (forest_adapt, tree_count);
//     ts = t8_forest_get_eclass_scheme (forest_adapt, eclass);

//     /* get the number of elements in the current tree */
//     current_tree = t8_forest_get_tree (forest_adapt, tree_count);
//     current_tree_num_elements =
//       t8_forest_get_tree_element_count (current_tree);

//     for (elem_count = 0; elem_count < current_tree_num_elements; ++elem_count) {

//       /* determing the current element according to the given tree id and element id within the tree */
//       current_element =
//         t8_forest_get_element_in_tree (forest_adapt, tree_count, elem_count);
//       //Subelements counter 
//       if (ts->t8_element_is_subelement (current_element)) {
//         subelement_count++;
//       }

//       if (get_LFN_elem_info) {  /* print current element */
// #if T8_ENABLE_DEBUG
//         t8_productionf
//           ("\n\n________________"
//            "\nCurrent element: local elem index of this process: %i of %i (without ghosts)\n",
//            elem_count, t8_forest_get_local_num_elements (forest_adapt));
//         ts->t8_element_debug_print (current_element);
// #endif
//       }

//       for (face_id = 0; face_id < ts->t8_element_num_faces (current_element);
//            ++face_id) {
//         LFN_call_count++;
//         time_LFN -= sc_MPI_Wtime ();
//         t8_forest_leaf_face_neighbors (forest_adapt, tree_count,
//                                        current_element, &neighbor_leafs,
//                                        face_id, &dual_faces, &num_neighbors,
//                                        &element_indices, &neigh_scheme,
//                                        forest_is_balanced);
//         time_LFN += sc_MPI_Wtime ();
//         /* free memory if neighbors exist */
//         if (num_neighbors > 0) {
//           if (get_LFN_elem_info) {
//             /* print all neighbor elements */
//             for (neighbor_count = 0; neighbor_count < num_neighbors;
//                  neighbor_count++) {
// #if T8_ENABLE_DEBUG
//               t8_productionf ("\n_________"
//                               "\nNeighbor: %i of %i at face %i: (dual face: %i | local index %i of %i (with ghosts)  | ghost, if >= %i):\n",
//                               neighbor_count + 1, num_neighbors, face_id,
//                               dual_faces[neighbor_count],
//                               element_indices[neighbor_count],
//                               t8_forest_get_local_num_elements (forest_adapt)
//                               + t8_forest_get_num_ghosts (forest_adapt),
//                               t8_forest_get_local_num_elements (forest_adapt)
//                               - 1);
//               ts->t8_element_debug_print (neighbor_leafs[neighbor_count]);
// #endif
//             }
//           }

//           neigh_scheme->t8_element_destroy (num_neighbors, neighbor_leafs);

//           T8_FREE (element_indices);
//           T8_FREE (neighbor_leafs);
//           T8_FREE (dual_faces);
//         }
//         else {
//           if (get_LFN_elem_info) {
// #if T8_ENABLE_DEBUG
//             /* no neighbor in this case */
//             t8_productionf ("\n_________"
//                             "\nNeighbor: at face %i: There is no neighbor (domain boundary).\n",
//                             face_id);
// #endif
//           }
//         }
//       }                         /* end of face loop */
//     }                           /* end of element loop */
//   }                             /* end of tree loop */

//   T8_ASSERT (subelement_count ==
//              t8_forest_get_local_num_subelements (forest_adapt));

//   if (get_LFN_stats)
//     t8_print_LFN_stats (t8_forest_get_global_num_elements (forest_adapt),
//                         t8_forest_get_global_num_subelements (forest_adapt),
//                         t8_forest_get_local_num_elements (forest_adapt),
//                         t8_forest_get_local_num_subelements (forest_adapt),
//                         LFN_call_count, time_LFN,
//                         time_LFN / (double) LFN_call_count, adaptation_count,
//                         num_adaptations);

//   t8_debugf
//     ("~~~~~~~~~~ The LFN test function finshed successful ~~~~~~~~~~\n");
}                               /* end of t8_LFN_test */

/* Initializing and adapting a forest */
static void
t8_transition_global (void)
{
  t8_debugf
    ("~~~~~~~~~~ Into the t8_transition_global function ~~~~~~~~~~\n");

  /* At the moment, subelements are only implemented for T8_ECLASS_HEX and quads */
  t8_eclass_t         eclass = T8_ECLASS_HEX;  /* depending on the include file, this will be the transitioned or default hex implementation */

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];

  /* ************************************************* Case Settings ************************************************* */

  /* refinement setting */
  int                 initlevel = 2;    /* initial uniform refinement level */
  int                 adaptlevel = 3;
  int                 minlevel = initlevel;     /* lowest level allowed for coarsening (minlevel <= initlevel) */
  int                 maxlevel = initlevel + adaptlevel;        /* highest level allowed for refining */

  /* refinement/adaptation criteria settings */
  // double              circ_midpoint_x = 0.0;
  // double              circ_midpoint_y = 0.0;
  // double              circ_midpoint_z = 0.0;
  // double              start_radius = 0.0;
  // double              band_width = 2.0;

  int                 num_adaptations = 1;      /* 1 for a single adapted forest */
  // double              radius_increase = 0.2;

  /* adaptation setting */
  int                 set_balance = 1;
  int                 set_transition = 1;

  /* cmesh settings */
  int                 single_tree_mesh = 1;
  // int                 multiple_tree_mesh = 0, num_x_trees = 1, num_y_trees = 2;
  // int                 hybrid_tree_mesh = 0;

  int                 periodic_boundary = 0;    /* use periodic boundaries */

  /* partition setting */
  int                 do_partition = 1;

  /* ghost setting */
  int                 do_ghost = 1;     /* if do_LFN_test = 1, then do_ghost must be set to 1 as well when using multiple processes */
  int                 ghost_version = 1;        /* use v1 for transitioned forests */

  /* LFN settings */
  int                 do_LFN_test = 0;

  /* vtk setting */
  int                 do_vtk = 1;
  int                 do_vtk_cmesh = 1;
  int                 do_vtk_ghost = 0;

  /* Monitoring (only available in debug configuration) */
  // int                 get_LFN_stats = 0;
  // int                 get_LFN_elem_info = 0;
  // int                 get_commit_stats = 0;
  // int                 get_general_stats = 0;

  /* ************************************** Check settings ************************************** */

  SC_CHECK_ABORT (num_adaptations > 0,
                  "Setting-Check failed: Set num_adaptations > 0");
  // SC_CHECK_ABORT (single_tree_mesh + multiple_tree_mesh + hybrid_tree_mesh ==
  //                 1,
  //                 "Setting-check failed: choose only one of {single_tree, multiple_tree, hybrid_cmesh}");
  if (do_LFN_test == 1) {
    SC_CHECK_ABORT (set_balance == 1,
                    "LFN is not implemented for non-balanced forests.");
    SC_CHECK_ABORT (do_ghost == 1,
                    "Setting-check failed: set do_ghost to one when applying the LFN test");
    if (set_transition == 1) {
      SC_CHECK_ABORT (ghost_version == 1,
                      "Setting-check failed: use ghost version 1 when applying the LFN test for transitioned forests.");

    }
  }

  /* ************************************** Initializing refinement criterion ************************************** */

  /* user-data (minlevel, maxlevel) */
  // t8_example_level_set_struct_t ls_data;
  // t8_basic_sphere_data_t sdata;

  /* Midpoint and radius of a sphere */
  // sdata.mid_point[0] = circ_midpoint_x;
  // sdata.mid_point[1] = circ_midpoint_y;
  // sdata.mid_point[2] = circ_midpoint_z;
  // sdata.radius = start_radius;

  // /* refinement parameter */
  // ls_data.band_width = band_width;
  // ls_data.L = t8_distance_to_sphere;
  // ls_data.min_level = minlevel;
  // ls_data.max_level = maxlevel;
  // ls_data.udata = &sdata;

  /* ********************************************* Initializing cmesh ********************************************** */

  /* building the cmesh, using the initlevel */
  if (single_tree_mesh) {
    /* construct a single tree hex cmesh */
    cmesh =
      t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0,
                              periodic_boundary);
  }


  // else if (multiple_tree_mesh) {
  //   T8_ASSERT (eclass == T8_ECLASS_HEX);
  //   /* this is by default a 2D or 3D quad cmesh of multiple trees */
  //   p4est_connectivity_t *brick =
  //     p4est_connectivity_new_brick (num_x_trees, num_y_trees,
  //                                   periodic_boundary, periodic_boundary);
  //   cmesh = t8_cmesh_new_from_p4est (brick, sc_MPI_COMM_WORLD, 0);
  //   p4est_connectivity_destroy (brick);
  // }
  // else if (hybrid_tree_mesh) {
  //   T8_ASSERT (eclass == T8_ECLASS_HEX || eclass == T8_ECLASS_TRIANGLE);
  //   /* this is by default a hybrid 2D quad-triangle forest */
  //   cmesh = t8_cmesh_new_periodic_hybrid (sc_MPI_COMM_WORLD);
  // }
  // else {
  //   SC_ABORT ("Specify cmesh geometry.");
  // }

  /* initialize a forest */
//   t8_forest_init (&forest);
  
//   /* set forest parameter via cmesh */
//   t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
//   t8_forest_set_level (forest, initlevel);
// #if DO_TRANSITION_HEX_SCHEME
//   t8_forest_set_scheme (forest, t8_scheme_new_transition_hex_cxx ());
// #else
//   t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
// #endif

//   /* commit the forest */
//   t8_forest_commit (forest);
// t8_debugf ("~~~~~~~~~~ Bis hier hin( after commit)\n");
//   t8_debugf ("~~~~~~~~~~ cmesh has been build ~~~~~~~~~~\n");

  // if (do_vtk_cmesh) {
    // snprintf (filename, BUFSIZ, "forest_hex_cmesh");
    // t8_forest_write_vtk (cmesh, filename);
    // t8_debugf ("~~~~~~~~~~ vtk of cmesh has been constructed ~~~~~~~~~~\n");
 // }
   
  /* ********************************** Adaptation (possibly with multiple steps) ************************************ */

  // double              commit_time_accum = 0;
  // double              LFN_time_accum = 0;
  // double              total_time = 0;
  // int                 global_num_elements_accum = 0;

  // total_time -= sc_MPI_Wtime ();
  // int                 adaptation_count;
  // for (adaptation_count = 1; adaptation_count <= num_adaptations;
  //      ++adaptation_count) {

  //   t8_debugf ("~~~~~~~~~~ Into adaptation %i of %i ~~~~~~~~~~\n",
  //              adaptation_count, num_adaptations);

    // /* initialization */
    // t8_forest_init (&forest_adapt);

    // /* Adapt the mesh according to the user data */

    // // t8_forest_set_user_data (forest_adapt, &ls_data);

    // t8_forest_set_adapt (forest_adapt, forest, t8_common_adapt_level_set, 0);

    // if (set_balance && !set_transition) {
    //   t8_forest_set_balance (forest_adapt, forest, 0);
    // }
    // if (set_transition) {
    //   t8_forest_set_transition (forest_adapt, forest, set_balance);
    // }
    // if (do_ghost) {
    //   t8_forest_set_ghost_ext (forest_adapt, do_ghost, T8_GHOST_FACES,
    //                            ghost_version);
    // }
    // if (do_partition) {
    //   t8_forest_set_partition (forest_adapt, forest, 0);
    // }

    /* adapt the mesh and take runtime for monitoring */
    // double              commit_time = 0;
    // commit_time_accum -= sc_MPI_Wtime ();
    // commit_time -= sc_MPI_Wtime ();
    // t8_forest_commit (forest_adapt);    /* adapt the forest */
    // commit_time_accum += sc_MPI_Wtime ();
    // commit_time += sc_MPI_Wtime ();
    // if (get_commit_stats) {
    //   t8_print_commit_stats (commit_time, num_adaptations, adaptation_count);
    //   t8_debugf ("~~~~~~~~~~ forest has been adapted ~~~~~~~~~~\n");
    // }

    // if (do_vtk) {
    //   t8_print_vtk (forest_adapt, filename, set_transition, set_balance,
    //                 single_tree_mesh, 0, 0,
    //                 eclass);}
      t8_debugf
        ("~~~~~~~~~~ vtk of adapted forest has been constructed ~~~~~~~~~~\n");
    //   if (do_vtk_ghost) {
    //     snprintf (filename, BUFSIZ, "forest_ghost_%i_%s",
    //               adaptation_count, t8_eclass_to_string[T8_ECLASS_HEX]);
    //     t8_forest_write_vtk_ext (forest_adapt, filename, 1, 1, 1, 1, 1, 0, 0,
    //                              0, NULL);
    //     t8_debugf
    //       ("~~~~~~~~~~ vtk of ghost has been constructed ~~~~~~~~~~\n");
    //   }
    // }

    /* iterate through all elements of the adapted, transitioned forest and compute
     * their neighbors to all faces. */
    // if (do_LFN_test) {
    //   LFN_time_accum -= sc_MPI_Wtime ();
    //   t8_LFN_test (forest_adapt, get_LFN_stats, adaptation_count,
    //                num_adaptations, get_LFN_elem_info);
    //   LFN_time_accum += sc_MPI_Wtime ();
    //   t8_debugf
    //     ("~~~~~~~~~~ all neighbors have been identified ~~~~~~~~~~\n");
    // }
 //t8_forest_write_vtk (forest,"forest_adapt_hex");
    /* Set forest to forest_adapt for the next step */
  //  forest = forest_adapt;
   
    /* Increase the radius of the sphere for the next step */
    // sdata.radius += radius_increase;

    /* Monitoring the total number of elements in the forest */
    // global_num_elements_accum +=
    //   t8_forest_get_global_num_elements (forest_adapt);

  //}                             /* end of adaptation loop */
  // total_time += sc_MPI_Wtime ();

 // t8_forest_unref (&forest_adapt);


  // if (get_general_stats) {
  //   t8_print_general_stats (commit_time_accum, num_adaptations,
  //                           global_num_elements_accum, total_time,
  //                           LFN_time_accum);
  // }

  t8_debugf
    ("~~~~~~~~~~ The t8_transition_global function finshed successful ~~~~~~~~~~\n");
}                               /* end of t8_transition_global */
    static void
t8_write_cmesh_vtk (t8_cmesh_t cmesh, const char *prefix)
{
  t8_cmesh_vtk_write_file (cmesh, prefix, 1.0);
}
int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_DEFAULT);
  t8_init (SC_LP_DEFAULT);

  t8_transition_global ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
