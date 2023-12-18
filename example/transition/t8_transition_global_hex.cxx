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
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
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



typedef struct
{
  double              mid_point[3];
  double              radius;
} t8_basic_sphere_data_t;

/* Compute the distance to a sphere around a mid_point with given radius. */
static double
t8_distance_to_sphere (const double x[3], double t, void *data)
{
  t8_basic_sphere_data_t *sdata = (t8_basic_sphere_data_t *) data;
  double             *M = sdata->mid_point;

  return t8_vec_dist (M, x) - sdata->radius;
}


void
t8_print_vtk (t8_forest_t forest_adapt, char filename[BUFSIZ],
              int set_transition, int set_balance, int single_tree_mesh,
              t8_eclass_t eclass, int adaptation_count)
{
  snprintf (filename, BUFSIZ, "forest_transitioned_%i",
                adaptation_count);
  t8_forest_write_vtk (forest_adapt, filename);
}


void
t8_print_commit_stats (double commit_time, int num_adaptations,
                       int adaptation_count)
{
  t8_productionf
    ("\n|++++++++++++++++++++ Commit statistics | adaptation %i of %i +++++++++++++++++++|\n"
     "|    Commit time total [s]: %.9f\n"
     "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|\n\n",
     adaptation_count, num_adaptations, commit_time);
}




void
t8_print_general_stats (double commit_time_total, int num_adaptations,
                        int global_num_elements_accum, double total_time,
                        double LFN_time_accum)
{
  t8_productionf
    ("\n|++++++++++++++++++++++++ Commit statistics | total +++++++++++++++++++++++++++|\n"
     "|    Average #elements:       %i\n"
     "|    Time total [s]:          %.3f (%.2f %%)\n"
     "|    Step time average [s]:   %.3f\n"
     "|    LFN time total [s]:      %.3f (%.2f %%)\n"
     "|    LFN time average [s]:    %.3f\n"
     "|    Commit time total [s]:   %.3f (%.2f %%)\n"
     "|    Commit time average [s]: %.3f\n"
     "|    Rest [s]:                %.3f (%.2f %%)\n"
     "|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|\n\n",
     global_num_elements_accum / num_adaptations, total_time, 100.,
     total_time / (double) num_adaptations, LFN_time_accum,
     100. * LFN_time_accum / total_time,
     LFN_time_accum / (double) num_adaptations, commit_time_total,
     100. * commit_time_total / total_time,
     commit_time_total / (double) num_adaptations,
     total_time - LFN_time_accum - commit_time_total,
     100. * (total_time - LFN_time_accum - commit_time_total) / total_time);
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


  if ((lelement_id % 2 ) == 0 ) {
        /* Refine this element. */
    return 1;
  }
  /* Do not change this element. */
  return 0;
}






/* Initializing and adapting a forest */
static void
t8_transition_global1 (void)
{
  // t8_debugf
  //   ("~~~~~~~~~~ Into the t8_transition_global function ~~~~~~~~~~\n");

  /* At the moment, subelements are only implemented for T8_ECLASS_HEX and quads */
  t8_eclass_t         eclass = T8_ECLASS_HEX;  /* depending on the include file, this will be the transitioned or default hex implementation */

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_forest_t         forest_adapt2;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];

  /* ************************************************* Case Settings ************************************************* */

  /* refinement setting */
  int                 initlevel = 1;    /* initial uniform refinement level */
  int                 adaptlevel = 3;

  /* adaptation setting */
  int                 set_balance = 1;
  int                 set_transition = 1;

  /* cmesh settings */
  int                 single_tree_mesh = 1;

  int                 periodic_boundary = 0;    /* use periodic boundaries */

  /* partition setting */
  int                 do_partition = 1;

  /* vtk setting */
  int                 do_vtk = 1;
  int                 do_vtk_cmesh = 1;


  /* Monitoring (only available in debug configuration) */
  // int                 get_commit_stats = 0;
  // int                 get_general_stats = 0;



  /* ********************************************* Initializing cmesh ********************************************** */


    /* construct a single tree hex cmesh */
    cmesh =
      t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0,
                              periodic_boundary);

  /* initialize a forest */
  t8_forest_init (&forest);
  
  /* set forest parameter via cmesh */
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_level (forest, initlevel);
#if DO_TRANSITION_HEX_SCHEME
  t8_forest_set_scheme (forest, t8_scheme_new_transition_hex_cxx ());
#else
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
#endif
// t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  /* commit the forest */
  t8_forest_commit (forest);

  int num_adaptions = 2;
  int  adaptation_count;
  for (adaptation_count = 1; adaptation_count <= num_adaptions;
       ++adaptation_count) {

    t8_forest_init (&forest_adapt);

    // forest_adapt = t8_forest_new_adapt (forest, t8_adapt_callback, 0, 0, 0);
    t8_debugf("------------------------------------\n");     
    t8_forest_set_adapt (forest_adapt, forest, t8_adapt_callback, 0);
    // t8_forest_set_balance (forest_adapt, forest, 1);   
  //  t8_debugf("---------------Into transition ---------------------\n");    
   t8_forest_set_transition (forest_adapt, forest, 1);

    //t8_forest_set_partition (forest_adapt, forest, 0);

    /* adapt the mesh and take runtime for monitoring */
    double              commit_time = 0;
    // commit_time_accum -= sc_MPI_Wtime ();
    commit_time -= sc_MPI_Wtime ();
    t8_forest_commit (forest_adapt);    /* adapt the forest */
    // commit_time_accum += sc_MPI_Wtime ();
    commit_time += sc_MPI_Wtime ();


    t8_print_commit_stats (commit_time, num_adaptions, adaptation_count);
  /* Set forest to forest_adapt for the next step */
    forest = forest_adapt;
    }
    //  t8_forest_init (&forest_adapt);
    // //  t8_forest_set_adapt (forest_adapt, forest, t8_adapt_callback, 0);
    //   t8_forest_set_balance (forest_adapt, forest, 0);
    //   t8_forest_set_transition(forest_adapt, forest, 1);
    //   t8_forest_commit (forest_adapt); 
    // // intsnprintf (filename, BUFSIZ, "forest_hex_adapted_mesh");
    // // t8_pr_vtk(forest_adapt, filename, 0, 0, 1, eclass);

    // // snprintf (filename, BUFSIZ, "forest_hex_balanced_mesh");
    // // t8_print_vtk(forest_adapt, filename, 0, 1, 1, eclass);

    //snprintf (filename, BUFSIZ, "forest_hex_transition_mesh");
    t8_print_vtk(forest, filename, 0, 0, 1, eclass, adaptation_count);

    t8_forest_unref (&forest);

}                               /* end of t8_transition_global */







/* Initializing and adapting a forest */
static void
t8_transition_global (void)
{
  t8_debugf
    ("~~~~~~~~~~ Into the t8_transition_global function ~~~~~~~~~~\n");

  /* At the moment, subelements are only implemented for T8_ECLASS_QUADS */
  t8_eclass_t         eclass = T8_ECLASS_HEX;  /* depending on the include file, this will be the transitioned or default quad implementation */

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];

  /* ************************************************* Case Settings ************************************************* */

  /* refinement setting */
  int                 initlevel = 3;    /* initial uniform refinement level */
  int                 adaptlevel = 3;
  int                 minlevel = initlevel;     /* lowest level allowed for coarsening (minlevel <= initlevel) */
  int                 maxlevel = initlevel + adaptlevel;        /* highest level allowed for refining */

  /* refinement/adaptation criteria settings */
  double              circ_midpoint_x = 0.5;
  double              circ_midpoint_y = 0.5;
  double              circ_midpoint_z = 0.5;
  double              start_radius = 0.0;
  double              band_width = 2.0;

  int                 num_adaptations = 10;      /* 1 for a single adapted forest */
  double              radius_increase = 0.05;

  /* adaptation setting */
  int                 set_balance = 1;
  int                 set_transition = 1;

  /* cmesh settings */
  int                 single_tree_mesh = 1;
  int                 multiple_tree_mesh = 0, num_x_trees = 1, num_y_trees =  2;
  int                 hybrid_tree_mesh = 0;

  int                 periodic_boundary = 0;    /* use periodic boundaries */

  /* partition setting */
  int                 do_partition = 1;

  /* ghost setting */
  int                 do_ghost = 0;     /* if do_LFN_test = 1, then do_ghost must be set to 1 as well when using multiple processes */
  int                 ghost_version = 0;        /* use v1 for transitioned forests */

  /* LFN settings */
  int                 do_LFN_test = 0;

  /* vtk setting */
  int                 do_vtk = 1;
  int                 do_vtk_cmesh = 0;
  int                 do_vtk_ghost = 0;

  /* Monitoring (only available in debug configuration) */
  int                 get_LFN_stats = 0;
  int                 get_LFN_elem_info = 0;
  int                 get_commit_stats = 0;
  int                 get_general_stats = 0;

  /* ************************************** Check settings ************************************** */

  SC_CHECK_ABORT (num_adaptations > 0,
                  "Setting-Check failed: Set num_adaptations > 0");
  SC_CHECK_ABORT (single_tree_mesh + multiple_tree_mesh + hybrid_tree_mesh ==
                  1,
                  "Setting-check failed: choose only one of {single_tree, multiple_tree, hybrid_cmesh}");
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
  t8_example_level_set_struct_t ls_data;
  t8_basic_sphere_data_t sdata;

  /* Midpoint and radius of a sphere */
  sdata.mid_point[0] = circ_midpoint_x;
  sdata.mid_point[1] = circ_midpoint_y;
  sdata.mid_point[2] = circ_midpoint_z;
  sdata.radius = start_radius;

  /* refinement parameter */
  ls_data.band_width = band_width;
  ls_data.L = t8_distance_to_sphere;
  ls_data.min_level = minlevel;
  ls_data.max_level = maxlevel;
  ls_data.udata = &sdata;

  /* ********************************************* Initializing cmesh ********************************************** */

  /* building the cmesh, using the initlevel */
  if (single_tree_mesh) {
    /* construct a single tree quad cmesh */
    cmesh =
      t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0,
                              periodic_boundary);
  }
  
  /* initialize a forest */
  t8_forest_init (&forest);

  /* set forest parameter via cmesh */
  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  t8_forest_set_level (forest, initlevel);
#if DO_TRANSITION_HEX_SCHEME
  t8_forest_set_scheme (forest, t8_scheme_new_transition_hex_cxx ());
#else
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
#endif

  /* commit the forest */
  t8_forest_commit (forest);

  t8_debugf ("~~~~~~~~~~ cmesh has been build ~~~~~~~~~~\n");

  if (do_vtk_cmesh) {
    snprintf (filename, BUFSIZ, "forest_cmesh");
    t8_forest_write_vtk (forest, filename);
    t8_debugf ("~~~~~~~~~~ vtk of cmesh has been constructed ~~~~~~~~~~\n");
  }

  /* ********************************** Adaptation (possibly with multiple steps) ************************************ */

  double              commit_time_accum = 0;
  double              LFN_time_accum = 0;
  double              total_time = 0;
  int                 global_num_elements_accum = 0;

  total_time -= sc_MPI_Wtime ();
  int                 adaptation_count;
  for (adaptation_count = 1; adaptation_count <= num_adaptations;
       ++adaptation_count) {

    t8_debugf ("~~~~~~~~~~ Into adaptation %i of %i ~~~~~~~~~~\n",
               adaptation_count, num_adaptations);

    /* initialization */
    t8_forest_init (&forest_adapt);

    /* Adapt the mesh according to the user data */
    t8_forest_set_user_data (forest_adapt, &ls_data);
    t8_forest_set_adapt (forest_adapt, forest, t8_common_adapt_level_set, 1);

    if (set_balance && !set_transition) {
      t8_forest_set_balance (forest_adapt, forest, 0);
    }
    if (set_transition) {
      t8_forest_set_transition (forest_adapt, forest, set_balance);
    }
    if (do_ghost) {
      t8_forest_set_ghost_ext (forest_adapt, do_ghost, T8_GHOST_FACES,
                               ghost_version);
    }
    if (do_partition) {
      t8_forest_set_partition (forest_adapt, forest, 0);
    }

    /* adapt the mesh and take runtime for monitoring */
    double              commit_time = 0;
    commit_time_accum -= sc_MPI_Wtime ();
    commit_time -= sc_MPI_Wtime ();
    t8_forest_commit (forest_adapt);    /* adapt the forest */
    commit_time_accum += sc_MPI_Wtime ();
    commit_time += sc_MPI_Wtime ();

    t8_print_commit_stats(commit_time, num_adaptations, adaptation_count);
    t8_print_vtk(forest_adapt, filename, 0, 0, 1, eclass, adaptation_count);
    /* Set forest to forest_adapt for the next step */
    forest = forest_adapt;

    /* Increase the radius of the sphere for the next step */
    sdata.radius += radius_increase;

    /* Monitoring the total number of elements in the forest */
    global_num_elements_accum +=
      t8_forest_get_global_num_elements (forest_adapt);

  }                             /* end of adaptation loop */
  total_time += sc_MPI_Wtime ();

  t8_forest_unref (&forest);


  t8_print_general_stats (commit_time_accum, num_adaptations,
                          global_num_elements_accum, total_time,
                          LFN_time_accum);

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
  //p4est_init (NULL, SC_LP_DEFAULT);
  t8_init (SC_LP_DEFAULT);

//Decide whether to take t8_transition_global t8_transition_global1
  t8_transition_global ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
