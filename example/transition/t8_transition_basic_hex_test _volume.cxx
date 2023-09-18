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



#ifdef T8_ENABLE_DEBUG
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
#endif
void
t8_print_vtk (t8_forest_t forest_adapt, char filename[BUFSIZ],
              int set_transition, int set_balance, int single_tree_mesh, int adaptation_count,
              t8_eclass_t eclass)
{
  if (set_transition) {
    if (single_tree_mesh)
      snprintf (filename, BUFSIZ, "forest_transitioned_hex_%i_%s",
                adaptation_count, t8_eclass_to_string[eclass]);
     
  if (set_balance) {
    if (single_tree_mesh)
      snprintf (filename, BUFSIZ, "forest_balanced_hex_%i_%s",
                adaptation_count, t8_eclass_to_string[eclass]);
  }
  t8_forest_write_vtk (forest_adapt, filename);
}
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

/* If a subelement is given, we apply the callback function to its parent */
  if (ts->t8_element_is_subelement (elements[0])) {
    t8_element_t      **parent = T8_ALLOC (t8_element_t *, 1);
    ts->t8_element_new (1, parent);
    ts->t8_element_parent (elements[0], parent[0]);
    T8_FREE (parent);
  }
  if ((lelement_id == 0) ) {
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

//-------------------print general status function ------------------------
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


//-------------------- print commit status ---------------------
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
//------------refinement criteria --------------------------------
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

//-------------------------transition function -----------------

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

  int level = 3;
  int subelement_counter = 0;
  int elem_count;
  int tree_count = 0;
  int current_tree_num_elements;
  t8_tree_t           current_tree;


  int                 set_balance = 1;
  int                 set_transition = 1;

  /* refinement setting */
  int                 initlevel = 1;    /* initial uniform refinement level */
  int                 adaptlevel = 1;
  int                 minlevel = initlevel;     /* lowest level allowed for coarsening (minlevel <= initlevel) */
  int                 maxlevel = initlevel + adaptlevel;        /* highest level allowed for refining */

  /* refinement/adaptation criteria settings */
  double              circ_midpoint_x = 0.5;
  double              circ_midpoint_y = 0.5;
  double              circ_midpoint_z = 0.5;
  double              start_radius = 0.0;
  double              band_width = 2.0;

  int                 num_adaptations = 1;      /* 1 for a single adapted forest */
  double              radius_increase = 0.2;

  /* Monitoring (only available in debug configuration) */
  int                 get_LFN_stats = 0;
  int                 get_LFN_elem_info = 0;
  int                 get_commit_stats = 0;
  int                 get_general_stats = 0;

  //----------------- measurements-----------------

double              commit_time_accum = 0;
double              total_time = 0;
int                 global_num_elements_accum = 0;
total_time -= sc_MPI_Wtime ();

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

  cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0,0,0);

    /* initialize a forest */
  t8_forest_init (&forest);

  t8_forest_set_cmesh (forest, cmesh, sc_MPI_COMM_WORLD);
  //t8_forest_set_level (forest, level);
  t8_forest_set_level (forest, initlevel);

 #if DO_TRANSITION_HEX_SCHEME
  t8_forest_set_scheme (forest, t8_scheme_new_transition_hex_cxx ());
#else
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
#endif


  t8_forest_commit (forest);


  int                 adaptation_count;
  for (adaptation_count = 1; adaptation_count <= num_adaptations;
       ++adaptation_count) {
  t8_forest_init (&forest_adapt);
  t8_forest_set_user_data (forest_adapt, &ls_data);
  //t8_forest_set_adapt (forest_adapt, forest, t8_adapt_callback, 0);
  t8_forest_set_adapt (forest_adapt, forest, t8_common_adapt_level_set, 1);

    if (set_balance && !set_transition) {
      t8_forest_set_balance (forest_adapt, forest, 0);
    }
    if (set_transition) {
      t8_forest_set_transition (forest_adapt, forest, set_balance);
    }

   /* adapt the mesh and take runtime for monitoring */
    double              commit_time = 0;
    commit_time_accum -= sc_MPI_Wtime ();
    commit_time -= sc_MPI_Wtime ();
    t8_forest_commit (forest_adapt);    /* adapt the forest */
    commit_time_accum += sc_MPI_Wtime ();
    commit_time += sc_MPI_Wtime ();

    t8_print_commit_stats (commit_time, 1, 1);


 // t8_forest_commit (forest_adapt);
  t8_print_vtk (forest_adapt, filename, set_transition, 1,1, 0,eclass);
  t8_forest_write_vtk(forest_adapt, "forest_adapt_wo_transition");
  ts = t8_forest_get_eclass_scheme (forest_adapt, eclass);  

//   current_tree = t8_forest_get_tree (forest_adapt, tree_count);

//   current_tree_num_elements = t8_forest_get_tree_element_count (current_tree);
//   //t8_productionf("num elements %i\n", current_tree_num_elements);

//   for (elem_count = 0; elem_count < current_tree_num_elements; ++elem_count) {
//     t8_element_t *current_element = t8_forest_get_element_in_tree (forest_adapt, tree_count, elem_count);
// //       //Subelements counter 
//        if (ts->t8_element_is_subelement (current_element)) {
//          subelement_counter++;
//          int sub_id = ts->t8_element_get_subelement_id(current_element);
//          //t8_productionf("sub ID %i\n", sub_id);
//          }
//   }
     /* Monitoring the total number of elements in the forest */
  global_num_elements_accum +=
  t8_forest_get_global_num_elements (forest_adapt);

    }                             /* end of adaptation loop */

  total_time += sc_MPI_Wtime ();

  t8_forest_unref (&forest_adapt);

      t8_print_general_stats (commit_time_accum, 1,
                            global_num_elements_accum, total_time,
                            0);

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
  //int level = 1;

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

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}


