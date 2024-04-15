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
#include <cstdint>
#include <sc_options.h>
#include <sc_statistics.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_vtk.h>
#include <example/common/t8_example_common.h>
#include <t8_cmesh.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_vec.h>
#include <t8_cmesh/t8_cmesh_examples.h> /* this is for t8_cmesh functions */
#include <t8_forest/t8_forest_profiling.h>
#include <string.h>             // for strlen fct

#define DO_TRANSITION_HEX_SCHEME 1

#if DO_TRANSITION_HEX_SCHEME 
#include <t8_schemes/t8_transition/t8_transition_conformal_hex/t8_transition_conformal_hex_cxx.hxx>
#include <t8_schemes/t8_transition/t8_transition_cxx.hxx>
#else
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#endif
/* helpful options for debugging */
#define GET_DEBUG_OUTPUT 1      /* get print statements of LFN, etc. */
#define DO_DEBUGGING_EXAMPLE 0  /* simulate a small example (4 timesteps, two adaptations) */
#define ENABLE_GLOBAL_CONSERVATION_CHECKS 1     /* this throws warnings for unused variables when 0 */

#define MAX_FACES 8             /* The maximum number of faces of an element */
/* TODO: This is not memory efficient. If we run out of memory, we can optimize here. */
T8_EXTERN_C_BEGIN ();


// /* Enum for statistics. */
// typedef enum
// {
//   ADVECT_ADAPT = 0,             /* adapt runtime */
//   ADVECT_PARTITION,             /* partition runtime */
//   ADVECT_PARTITION_PROCS,       /* number of processes sent to in partition */
//   ADVECT_PARTITION_DATA,        /* data partitioning runtime */
//   ADVECT_BALANCE,               /* balance runtime */
//   ADVECT_BALANCE_ROUNDS,        /* number of rounds in balance */
//   ADVECT_GHOST,                 /* ghost runtime */
//   ADVECT_GHOST_SENT,            /* number of ghosts sent to other processes */
//   ADVECT_GHOST_EXCHANGE,        /* ghost exchange runtime */
//   ADVECT_GHOST_WAIT,            /* ghost exchange waittime */
//   ADVECT_REPLACE,               /* forest_iterate_replace runtime */
//   ADVECT_IO,                    /* vtk runtime */
//   ADVECT_ELEM_AVG,              /* average global number of elements (per time step) */
//   ADVECT_INIT,                  /* initialization runtime */
//   ADVECT_AMR,                   /* AMR runtime (adapt+partition+ghost+balance) including data exchange (partition/ghost) */
//   ADVECT_NEIGHS,                /* neighbor finding runtime */
//   ADVECT_FLUX,                  /* flux computation runtime */
//   ADVECT_DUMMY,                 /* dummy operations to increase load (see -s option) */
//   ADVECT_SOLVE,                 /* solver runtime */
//   ADVECT_TOTAL,                 /* overall runtime */
//   ADVECT_ERROR_INF,             /* l_infty error */
//   ADVECT_ERROR_2,               /* L_2 error */
//   ADVECT_VOL_LOSS,              /* The loss in volume (region with LS < 0) in percent */
//   ADVECT_NUM_STATS              /* The number of statistics that we measure */
// } advect_stats_t


// /* Names of statistics that we measure */
// const char         *advect_stat_names[ADVECT_NUM_STATS] = {
//   "adapt",
//   "partition",
//   "partition_procs_sent",
//   "partition_data",
//   "balance",
//   "balance_rounds",
//   "ghost",
//   "ghost_sent",
//   "ghost_exchange",
//   "ghost_exchange_wait",
//   "replace",
//   "vtk_print",
//   "number_elements",
//   "init",
//   "AMR",
//   "neighbor_finding",
//   "flux_computation",
//   "dummy_ops",
//   "solve",
//   "total",
//   "l_infty_error",
//   "L_2",
//   "volume_loss_[%]"
// };


// /** The description of the problem configuration.
//  *  We store all necessary parameters, such as the initial level-set function, the flow function,
//  *  data needed for adaptation etc.
//  *  We also store the current forest and element data here.
//  */
// typedef struct
// {
//   t8_flow_function_3d_fn u; /**< Fluid field */
//   t8_example_level_set_fn phi_0; /**< Initial condition for phi */
//   void               *udata_for_phi; /**< User data passed to phi */
//   t8_forest_t         forest; /**< The forest in use */
//   t8_forest_t         forest_adapt; /**< The forest after adaptation */
//   sc_array_t         *element_data; /**< Array of type t8_advect_element_data_t of length
//                               num_local_elements + num_ghosts */
//   /* TODO: shorten element_data by number of ghosts */
//   sc_array_t         *element_data_adapt; /**< element_data for the adapted forest, used during adaptation to interpolate values */
//   /* We store the phi values in an extra array, since this data must exist for ghost
//    * element as well and is communicated with other processes in ghost_exchange. */
//   sc_array_t         *phi_values; /**< For each element and ghost its phi value. */
//   sc_array_t         *phi_values_adapt; /**< phi values for the adapted forest, used during adaptaption to interpolate values. */
//   sc_MPI_Comm         comm; /**< MPI communicator used */
//   sc_statinfo_t       stats[ADVECT_NUM_STATS]; /**< Runtimes and other statistics. */
//   double              t; /**< Current simulation time */
//   double              T; /**< End time */
//   double              cfl; /**< CFL number */
//   double              delta_t; /**< Current time step */
//   double              min_grad, max_grad; /**< bounds for refinement */
//   double              min_vol; /**< minimum element volume at level 'level' */
//   double              band_width; /**< width of the refinement band */
//   int                 num_time_steps; /**< Number of time steps computed so far.
//                                         (If delta_t is constant then t = num_time_steps * delta_t) */
//   int                 vtk_count; /**< If vtk output is enabled, count the number of pvtu files written. */
//   int                 level; /**< Initial refinement level */
//   int                 maxlevel; /**< Maximum refinement level */
//   int                 volume_refine; /**< If >= refine elements only if their volume is greater
//                                        than the minimum volume at level 'level + volume_refine' */
//   int                 dim; /**< The dimension of the mesh */
//   int                 dummy_op; /**< If true, we carry out more (but useless) operations
//                                      per element, in order to simulate more computation load */
//   int                 transition;       /* Flag to decide whether the forest should be transitioned or not */
// } t8_advect_problem_t;



// /** The per element data */
// typedef struct
// {
//   double              midpoint[3]; /**< coordinates of element midpoint in R^3 */
//   double              vol; /**< Volume of this element */
//   double              phi_new; /**< Value of solution at midpoint in next time step */
//   double             *fluxes[MAX_FACES]; /**< The fluxes to each neeighbor at a given face */
//   int                 flux_valid[MAX_FACES];  /**< If > 0, this flux was computed, if 0 memory was allocated
//                                                    for this flux, but not computed. If < 0, no memory was allocated. */
//   int                 level; /**< The refinement level of the element. */
//   int                 num_faces; /**< The number of faces */
//   int                 num_neighbors[MAX_FACES]; /**< Number of neighbors for each face */
//   int                *dual_faces[MAX_FACES]; /**< The face indices of the neighbor elements */
//   t8_locidx_t        *neighs[MAX_FACES]; /**< Indices of the neighbor elements */
//   int8_t              neigh_level[MAX_FACES]; /**< The level of the face neighbors at this face. */
// } t8_advect_element_data_t;


// double              time_interpolation = 0, time_adapt =
//   0, time_leaf_face_neighbors = 0;



//   /* Return the phi value of a given local or ghost element.
//  * 0 <= ielement < num_elements + num_ghosts
//  */
// static double
// t8_advect_element_get_phi (const t8_advect_problem_t * problem,
//                            t8_locidx_t ielement)
// {
//   return *((double *)
//            t8_sc_array_index_locidx (problem->phi_values, ielement));
// }

// static double
// t8_advect_element_get_phi_adapt (const t8_advect_problem_t * problem,
//                                  t8_locidx_t ielement)
// {
//   return *((double *)
//            t8_sc_array_index_locidx (problem->phi_values_adapt, ielement));
// }

// static double
// t8_advect_get_global_phi (const t8_advect_problem_t * problem)
// {
//   t8_locidx_t         lelement, num_local_elem;
//   double              scaled_phi_global = 0;
//   t8_advect_element_data_t *elem_data;

//   num_local_elem = t8_forest_get_local_num_elements (problem->forest);
//   T8_ASSERT (num_local_elem <=
//              (t8_locidx_t) problem->element_data->elem_count);
//   /* iterate over all all elements */
//   for (lelement = 0; lelement < num_local_elem; lelement++) {

//     /* Get element data */
//     elem_data = (t8_advect_element_data_t *)
//       t8_sc_array_index_locidx (problem->element_data, lelement);

//     scaled_phi_global +=
//       t8_advect_element_get_phi (problem, lelement) * elem_data->vol;
//   }

//   return scaled_phi_global;
// }

// /* Set the phi value of an element to a given entry */
// static void
// t8_advect_element_set_phi (const t8_advect_problem_t * problem,
//                            t8_locidx_t ielement, double phi)
// {
//   *((double *) t8_sc_array_index_locidx (problem->phi_values, ielement)) =
//     phi;
// }

// /* Set the phi value of an element in the adapted forest to a given entry */
// static void
// t8_advect_element_set_phi_adapt (const t8_advect_problem_t * problem,
//                                  t8_locidx_t ielement, double phi)
// {
//   *((double *) t8_sc_array_index_locidx (problem->phi_values_adapt, ielement))
//     = phi;
// }

// /* Adapt the forest. Elements are refined and coarsened randomly. */
// static int
// t8_advect_adapt_random (t8_forest_t forest, t8_forest_t forest_from,
//                         t8_locidx_t ltree_id, t8_locidx_t lelement_id,
//                         t8_eclass_scheme_c *ts, int isfamily,
//                         int num_elements, t8_element_t *elements[])
// {
//   t8_advect_problem_t *problem;
//   int                 level;

//   /* Get a pointer to the problem from the user data pointer of forest */
//   problem = (t8_advect_problem_t *) t8_forest_get_user_data (forest);

//   /* Get the element's level */
//   level = ts->t8_element_level (elements[0]);

//   int                 r = rand () % 99; /* random number between 0 and 99 */

//   if (level > problem->level && num_elements > 1 && r < 100) {
//     /* It is not possible to refine this level */
//     return -1;
//   }
//   else if (level < problem->maxlevel && r < 20) {
//     return 1;
//   }
//   else {
//     return 0;
//   }
// }


// /* Adapt the forest. We refine if the level-set function is close to zero
//  * and coarsen if it is larger than a given threshhold. */
// static int
// t8_advect_adapt (t8_forest_t forest, t8_forest_t forest_from,
//                  t8_locidx_t ltree_id, t8_locidx_t lelement_id,
//                  t8_eclass_scheme_c *ts, int isfamily, int num_elements,
//                  t8_element_t *elements[])
// {
//   t8_advect_problem_t *problem;
//   t8_advect_element_data_t *elem_data;
//   double              band_width;
//   double              elem_diam;
//   int                 level;
//   t8_locidx_t         offset;
//   double              phi;
//   double              vol_thresh;

//   /* Get a pointer to the problem from the user data pointer of forest */
//   problem = (t8_advect_problem_t *) t8_forest_get_user_data (forest);
//   /* Get the element's level */
//   level = ts->t8_element_level (elements[0]);
//   if (level >= problem->maxlevel && num_elements == 1) {
//     /* It is not possible to refine this level */
//     return 0;
//   }
//   if (level < problem->level && num_elements == 1) {
//     /* refine up to min init refinement level */
//     return 1;
//   }

//   /* Compute the volume threshold. Elements larger than this and
//    * close to the 0 level-set are refined */
//   if (problem->volume_refine >= 0) {
//     vol_thresh =
//       problem->min_vol / (1 << (problem->dim * problem->volume_refine));
//   }
//   else {
//     vol_thresh = 0;
//   }
//   /* Get the value of phi at this element */
//   offset = t8_forest_get_tree_element_offset (forest_from, ltree_id);
//   phi = t8_advect_element_get_phi (problem, lelement_id + offset);

//   if (true) {
//     float               phi_threshhold = 0.3;
//     if (phi >= phi_threshhold && level < problem->maxlevel) {
//       return 1;
//     }
//     if (phi <= phi_threshhold && num_elements > 1 && level > problem->level) {
//       return -1;
//     }
//     /* if no rule applies, do nothing */
//     return 0;
//   }
//   else {
//     /* Get a pointer to the element data */
//     elem_data = (t8_advect_element_data_t *)
//       t8_sc_array_index_locidx (problem->element_data, lelement_id + offset);

//     /* Refine if close to levelset, coarsen if not */
//     band_width = problem->band_width;
//     elem_diam = t8_forest_element_diam (forest_from, ltree_id, elements[0]);
//     if (fabs (phi) > 2 * band_width * elem_diam) {
//       /* coarsen if this is a family and level is not too small */
//       return -(num_elements > 1 && level > problem->level);
//     }
//     else if (fabs (phi) < band_width * elem_diam
//              && elem_data->vol > vol_thresh) {
//       /* refine if level is not too large */
//       return level < problem->maxlevel;
//     }
//     /* if no rule applies, do nothing */
//     return 0;
//   }                             /* refine areas with high values */
// }                               /* end of t8_advect_adapt */


// /* Initial geometric adapt scheme */
// static int
// t8_advect_adapt_init (t8_forest_t forest, t8_forest_t forest_from,
//                       t8_locidx_t ltree_id, t8_locidx_t lelement_id,
//                       t8_eclass_scheme_c *ts, int isfamily, int num_elements,
//                       t8_element_t *elements[])
// {
// #if 1                           /* do nothing */
//   return 0;
// #endif

// #if 0                           /* refine all lower right elements */
//   int                 coord[3] = { };
//   ts->t8_element_anchor (elements[0], coord);

//   if (coord[0] > coord[1]) {
//     return 1;
//   }
//   return 0;
// #endif

// #if 0                           /* refinement diag */
//   int                 coord[3] = { };
//   ts->t8_element_anchor (elements[0], coord);

//   if (coord[0] == coord[1]) {
//     return 1;
//   }
//   return 0;
// #endif

// #if 0                           /* refinement every second element */
//   if (lelement_id % 2 == 0) {
//     return 1;
//   }
//   return 0;
// #endif

// #if 0                           /* refinement all left elements */
//   int                 coord[3] = { };
//   ts->t8_element_anchor (elements[0], coord);
//   int                 len = P8EST_ROOT_LEN;
//   if (coord[0] < len / 2) {
//     return 1;
//   }
//   return 0;
// #endif
// }


// /* Compute the total volume of the elements with negative phi value */
// static double
// t8_advect_level_set_volume (const t8_advect_problem_t * problem)
// {
//   t8_locidx_t         num_local_elements, ielem;
//   t8_advect_element_data_t *elem_data;
//   double              volume = 0, global_volume = 0;
//   double              phi;

//   num_local_elements = t8_forest_get_local_num_elements (problem->forest);

//   for (ielem = 0; ielem < num_local_elements; ielem++) {
//     elem_data = (t8_advect_element_data_t *)
//       t8_sc_array_index_locidx (problem->element_data, ielem);
//     phi = t8_advect_element_get_phi (problem, ielem);
//     if (phi < 0) {
//       volume += elem_data->vol;
//     }
//   }
//   sc_MPI_Allreduce (&volume, &global_volume, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
//                     problem->comm);
//   return global_volume;
// }


// /* Compute the mean l_1 error of the stored phi values compared to a
//  * given analytical function at time problem->t */
// static double
// t8_advect_l1_error_mean (const t8_advect_problem_t * problem,
//                          t8_example_level_set_fn analytical_sol)
// {
//   t8_locidx_t         num_local_elements, ielem;
//   t8_advect_element_data_t *elem_data;
//   double              phi;
//   double              diff, ana_sol;
//   double              error = 0;

//   num_local_elements = t8_forest_get_local_num_elements (problem->forest);
//   for (ielem = 0; ielem < num_local_elements; ielem++) {
//     elem_data = (t8_advect_element_data_t *)
//       t8_sc_array_index_locidx (problem->element_data, ielem);
//     /* Compute the analytical solution: assume the analytical solution at time T equals the initial condition */
//     double              ana_sol_transported[3] = { };
//     ana_sol_transported[0] = elem_data->midpoint[0];
//     ana_sol_transported[1] = elem_data->midpoint[1];
//     ana_sol_transported[2] = elem_data->midpoint[2];
//     ana_sol =
//       analytical_sol (ana_sol_transported, problem->t,
//                       problem->udata_for_phi);

//     /* Compute the error as the stored value at the midpoint of this element minus the solution at this midpoint */
//     phi = t8_advect_element_get_phi (problem, ielem);
//     diff = phi - ana_sol;
//     error += fabs (diff);       /* add absolute difference */
//   }
//   /* return mean l1 error */
//   return error / num_local_elements;
// }


// /* Compute the relative l_2 error of the stored phi values compared to a
//  * given analytical function at time problem->t */
// static double
// t8_advect_l_2_rel (const t8_advect_problem_t * problem,
//                    t8_example_level_set_fn analytical_sol, double distance)
// {
//   t8_locidx_t         num_local_elements, ielem, count = 0;
//   t8_advect_element_data_t *elem_data;
//   double              phi;
//   double              diff, ana_sol;
//   double              error[2] = {
//     0, 0
//   }, el_error, global_error[2];

//   num_local_elements = t8_forest_get_local_num_elements (problem->forest);
//   for (ielem = 0; ielem < num_local_elements; ielem++) {
//     elem_data = (t8_advect_element_data_t *)
//       t8_sc_array_index_locidx (problem->element_data, ielem);
//     /* Compute the analytical solution */
//     ana_sol =
//       analytical_sol (elem_data->midpoint, problem->t,
//                       problem->udata_for_phi);
// #if 1
//     if (fabs (ana_sol) < distance)
// #endif
//     {
//       count++;
//       /* Compute the error as the stored value at the midpoint of this element
//        * minus the solution at this midpoint */
//       phi = t8_advect_element_get_phi (problem, ielem);
//       diff = fabs (phi - ana_sol);
//       el_error = diff * diff * elem_data->vol;
//       error[0] += el_error;
//       error[1] += ana_sol * ana_sol * elem_data->vol;
//     }
//   }
//   t8_debugf ("[advect] L_2 %e  %e\n", error[0], error[1]);
//   t8_debugf ("[advect] L_2 %i elems\n", count);
//   /* Compute the maximum of the error among all processes */
//   sc_MPI_Allreduce (&error, &global_error, 2, sc_MPI_DOUBLE, sc_MPI_SUM,
//                     problem->comm);

//   /* Return the relative error, that is the l_infty error divided by
//    * the l_infty norm of the analytical solution */
//   return sqrt (global_error[0]) / sqrt (global_error[1]);
// }


// static double
// t8_advect_l_infty_rel (const t8_advect_problem_t * problem,
//                        t8_example_level_set_fn analytical_sol,
//                        double distance)
// {
//   t8_locidx_t         num_local_elements, ielem;
//   t8_advect_element_data_t *elem_data;
//   double              phi;
//   double              ana_sol;
//   double              error[2] = {
//     -1, 0
//   }, el_error, global_error[2];

//   num_local_elements = t8_forest_get_local_num_elements (problem->forest);
//   for (ielem = 0; ielem < num_local_elements; ielem++) {
//     elem_data = (t8_advect_element_data_t *)
//       t8_sc_array_index_locidx (problem->element_data, ielem);

//     /* Compute the analytical solution */
//     ana_sol =
//       analytical_sol (elem_data->midpoint, problem->t,
//                       problem->udata_for_phi);
// #if 1
//     if (fabs (ana_sol) < distance)
// #endif
//     {
//       /* Compute the error as the stored value at the midpoint of this element
//        * minus the solution at this midpoint */
//       phi = t8_advect_element_get_phi (problem, ielem);
//       el_error = fabs ((phi - ana_sol));
//       error[0] = SC_MAX (error[0], el_error);
//       /* Compute the l_infty norm of the analytical solution */
//       error[1] = SC_MAX (error[1], ana_sol);
//     }
//   }
//   /* Compute the maximum of the error among all processes */
//   sc_MPI_Allreduce (&error, &global_error, 2, sc_MPI_DOUBLE, sc_MPI_MAX,
//                     problem->comm);

//   /* Return the relative error, that is the l_infty error divided by
//    * the l_infty norm of the analytical solution */
//   return global_error[0] / global_error[1];
// }

// static double
// t8_advect_flux_upwind_1d (const t8_advect_problem_t * problem,
//                           const t8_locidx_t el_plus,
//                           const t8_locidx_t el_minus, int face)
// {
//   double              x_j_half[3];
//   int                 idim;
//   double              u_at_x_j_half[3];
//   double              phi;
//   int                 sign;
//   t8_advect_element_data_t *el_data_plus;

//   /*
//    *    | --x-- | --x-- |   Two elements, midpoints marked with 'x'
//    *       x_j     x_j+1
//    *          x_j_half
//    */
//   /* Compute x_j_half */
//   el_data_plus = (t8_advect_element_data_t *)
//     t8_sc_array_index_locidx (problem->element_data, el_plus);
//   for (idim = 0; idim < 3; idim++) {
//     x_j_half[idim] =
//       (el_data_plus->midpoint[idim] -
//        (idim == 0 ? el_data_plus->vol / 2 : 0));
//   }
//   /* Compute u at the interval boundary. */
//   problem->u (x_j_half, problem->t, u_at_x_j_half);
//   /* In 1D we are only interested in the firs coordinate of u */

//   sign = face == 0 ? -1 : 1;
//   if (sign * u_at_x_j_half[0] >= 0) {
//     /* we have outflow */
//     phi = -t8_advect_element_get_phi (problem, el_plus);
//   }
//   else {
//     /* we have inflow */
//     /* el_minus may be negative, in this case, el_plus is at the boundary
//      * and we use phi = 0. */
//     phi = el_minus >= 0 ? t8_advect_element_get_phi (problem, el_minus) : 0;
//   }
//   return u_at_x_j_half[0] * phi;
// }

// /* Compute the flux across a given face between two elements */
// /* face is the face number as seen from el_data_plus */
// /* This works also if element_plus hangs on element_minus.
//  * It does not work if it hangs the other way around. */
// static double
// t8_advect_flux_upwind (const t8_advect_problem_t * problem,
//                        double el_plus_phi,
//                        double el_minus_phi,
//                        t8_locidx_t ltreeid,
//                        const t8_element_t *element_plus,
//                        const double *tree_vertices, int face)
// {
//   double              face_center[3];
//   double              u_at_face_center[3];
//   double              normal[3], normal_times_u;
//   double              area;

//   /*
//    *    | --x-- | --x-- |   Two elements, midpoints marked with 'x'
//    *       x_j     x_j+1
//    *          face_center
//    */

//   /* Compute the center coordinate of the face */
//   t8_forest_element_face_centroid (problem->forest, ltreeid, element_plus,
//                                    face, face_center);
//   /* Compute u at the face center. */
//   problem->u (face_center, problem->t, u_at_face_center);
//   /* Compute the normal of the element at this face */
//   t8_forest_element_face_normal (problem->forest, ltreeid, element_plus,
//                                  face, normal);
//   /* Compute the area of the face */
//   area =
//     t8_forest_element_face_area (problem->forest, ltreeid, element_plus,
//                                  face);

//   /* Compute the dot-product of u and the normal vector */
//   normal_times_u = t8_vec_dot (normal, u_at_face_center);

// #if GET_DEBUG_OUTPUT
//   /* Output, mainly for debugging */
//   t8_productionf ("[advect] face %i\n", face);
//   t8_productionf ("[advect] face normal vec %f %f %f\n", normal[0], normal[1],
//                   normal[2]);
//   t8_productionf ("[advect] face center %f %f %f\n", face_center[0],
//                   face_center[1], face_center[2]);
//   t8_productionf ("[advect] u %f %f %f\n", u_at_face_center[0],
//                   u_at_face_center[1], u_at_face_center[2]);
//   t8_productionf ("[advect] normal_times_u: %f\n", normal_times_u);
//   t8_productionf ("[advect] area %f\n", area);
//   t8_productionf ("[advect] phi+ (phi at current element) %f\n", el_plus_phi);
//   t8_productionf ("[advect] phi- (phi at neighbor) %f\n", el_minus_phi);
// #endif

//   if (normal_times_u >= 0) {
// #if GET_DEBUG_OUTPUT
//     /* u flows out of the element_plus */
//     t8_debugf ("[advect] flux (out of current element): %f\n",
//                -el_plus_phi * normal_times_u * area);
// #endif
//     return -el_plus_phi * normal_times_u * area;
//   }
//   else {
//     /* u flows into the element_plus */
// #if GET_DEBUG_OUTPUT
//     t8_debugf ("[advect] flux (into current element): %f\n",
//                -el_minus_phi * normal_times_u * area);
// #endif
//     return -el_minus_phi * normal_times_u * area;
//   }
// }

// /* Compute the flux if the element has hanging neighbors:
//  *
//  *  x -- x ---- x
//  *  |    |      |
//  *  x -- x      |
//  *  |    |      |
//  *  x -- x ---- x
//  *
//  * neighs  el_hang
//  *
//  */
// int                 a = 0;
// static double
// t8_advect_flux_upwind_hanging (const t8_advect_problem_t * problem,
//                                t8_locidx_t iel_hang,
//                                t8_locidx_t ltreeid,
//                                t8_element_t *element_hang,
//                                int face, int adapted_or_partitioned)
// {
//   int                 i, num_face_children, child_face;
//   t8_eclass_scheme_c *ts;
//   t8_eclass           eclass;
//   t8_element_t      **face_children;
//   t8_advect_element_data_t *neigh_data;
//   double              flux = 0;
//   int                 dual_face;
//   t8_locidx_t         neigh_id;
//   int                 neigh_is_ghost;
//   t8_advect_element_data_t *el_hang;
//   double              phi_plus, phi_minus;

//   /* Get a pointer to the element */
//   el_hang = (t8_advect_element_data_t *)
//     t8_sc_array_index_locidx (problem->element_data, iel_hang);
//   /* Get the eclass and the scheme for the element */
//   eclass = t8_forest_get_tree_class (problem->forest, ltreeid);
//   ts = t8_forest_get_eclass_scheme (problem->forest, eclass);
//   /* Compute the children of the element at the face */
//   num_face_children = ts->t8_element_num_face_children (element_hang, face);
//   T8_ASSERT (num_face_children == el_hang->num_neighbors[face]);

//   face_children = T8_ALLOC (t8_element_t *, num_face_children);
//   ts->t8_element_new (num_face_children, face_children);
//   ts->t8_element_children_at_face (element_hang, face, face_children,
//                                    num_face_children, NULL);

//   /* Store the phi value of el_hang. We use it as the phi value of the
//    * children to compute the flux */
//   phi_plus = t8_advect_element_get_phi (problem, iel_hang);

//   for (i = 0; i < num_face_children; i++) {
//     child_face = ts->t8_element_face_child_face (element_hang, face, i);
//     /* Get a pointer to the neighbor's element data */
//     neigh_id = el_hang->neighs[face][i];
//     neigh_data = (t8_advect_element_data_t *)
//       t8_sc_array_index_locidx (problem->element_data, neigh_id);
//     neigh_is_ghost =
//       neigh_id >= t8_forest_get_local_num_elements (problem->forest);
//     phi_minus = t8_advect_element_get_phi (problem, neigh_id);
//     /* Compute the flux */
//       el_hang->fluxes[face][i] =
//       t8_advect_flux_upwind (problem, phi_plus, phi_minus, ltreeid,
//                              face_children[face_children_count],
//                              tree_vertices, child_face);
//     el_hang->fluxes[face][i] =
//       t8_advect_flux_upwind (problem, phi_plus, phi_minus, ltreeid,
//                              face_children[i], child_face);
//     // if (a == 1) printf  ("%i %i %f\n",face, i, el_hang->fluxes[face][i]);
//     /* Set the flux of the neighbor element */
//     dual_face = el_hang->dual_faces[face][i];
//     if (!adapted_or_partitioned && !neigh_is_ghost) {

//       if (neigh_data->flux_valid[dual_face] < 0) {
//         /* We need to allocate the fluxes */
//         neigh_data->fluxes[dual_face] = T8_ALLOC (double, 1);
//       }
//       // printf ("face %i neigh %i df %i\n", face, neigh_id, dual_face);
//       SC_CHECK_ABORT (dual_face < neigh_data->num_faces, "num\n");
//       // SC_CHECK_ABORT (neigh_data->num_neighbors[dual_face] == 1, "entry\n");
//       neigh_data->num_neighbors[dual_face] = 1;
//       neigh_data->fluxes[dual_face][0] = -el_hang->fluxes[face][i];
//       neigh_data->flux_valid[dual_face] = 1;
//     }

//     flux += el_hang->fluxes[face][i];
//   }

//   el_hang->flux_valid[face] = 1;
//   /* clean-up */
//   ts->t8_element_destroy (num_face_children, face_children);
//   T8_FREE (face_children);

//   a = 2;
//   return flux;
// }

// /* If an element is at the domain boundary, we encode boundary conditions
//  * by setting a phi value for an imaginative neighbor element.
//  * We currently set the phi value
//  * to the value of the element itself. */
// static void
// t8_advect_boundary_set_phi (const t8_advect_problem_t * problem,
//                             t8_locidx_t ielement, double *boundary_phi)
// {

//   *boundary_phi = t8_advect_element_get_phi (problem, ielement);
// }

// #if 0
// static double
// t8_advect_lax_friedrich_alpha (const t8_advect_problem_t * problem,
//                                const t8_advect_element_data_t *
//                                el_data_plus,
//                                const t8_advect_element_data_t * el_data_minus)
// {
//   double              alpha;
//   double              dist, u_plus[3], u_minus[3];

//   /* We compute alpha as the derivative of u at the midpoint between
//    * the cells */

//   /* The distance between the two cells is the sum of their length divided by two */

//   dist = (el_data_plus->vol + el_data_minus->vol) / 2.;
//   /* Approximate the derivative of u */

//   problem->u (el_data_plus->midpoint, problem->t, u_plus);
//   problem->u (el_data_minus->midpoint, problem->t, u_minus);
//   /* in 1D we are only interested in the first coordinate of u */
//   alpha = fabs ((u_plus[0] - u_minus[0]) / dist);

//   return alpha;
// }

// static double
// t8_advect_flux_lax_friedrich_1d (const t8_advect_problem_t * problem,
//                                  const t8_advect_element_data_t *
//                                  el_data_plus,
//                                  const t8_advect_element_data_t *
//                                  el_data_minus)
// {
//   double              alpha = 0;        /* TODO: Choose alpha according to a reasonable criterion */
//   double              x_j_half[3];
//   int                 idim;
//   double              u_at_x_j_half[3];
//   double              phi_sum, phi_diff;

//   /*
//    *    | --x-- | --x-- |   Two elements, midpoints marked with 'x'
//    *       x_j     x_j+1
//    *          x_j_half
//    */
//   /* Compute x_j_half */
//   for (idim = 0; idim < 3; idim++) {
//     x_j_half[idim] =
//       (el_data_plus->midpoint[idim] -
//        (idim == 0 ? el_data_plus->vol / 2 : 0));
//   }

//   /* Compute u at the interval boundary. */
//   problem->u (x_j_half, problem->t, u_at_x_j_half);

//   /* Compute the sum of both phi values */
//   phi_sum = el_data_minus->phi + el_data_plus->phi;
//   /* Compute the difference of both */
//   phi_diff = el_data_plus->phi - el_data_minus->phi;

//   /* Compute alpha */
//   alpha =
//     t8_advect_lax_friedrich_alpha (problem, el_data_plus, el_data_minus);
//   /* in 1D only the first coordinate of u is interesting */
//   return .5 * (u_at_x_j_half[0] * phi_sum - alpha * phi_diff);
// }
// #endif


// static void
// t8_advect_advance_element (t8_advect_problem_t * problem,
//                            t8_locidx_t lelement)
// {
//   int                 iface, ineigh;
//   double              flux_sum = 0;
//   double              phi;
//   t8_advect_element_data_t *elem_data;

//   /* Get a pointer to the element */
//   elem_data = (t8_advect_element_data_t *)
//     t8_sc_array_index_locidx (problem->element_data, lelement);
//   /* Get the phi value of the element */
//   phi = t8_advect_element_get_phi (problem, lelement);
//   T8_ASSERT (3 <= elem_data->num_faces && elem_data->num_faces <= 4);
//   /* Sum all the fluxes */
//   for (iface = 0; iface < elem_data->num_faces; iface++) {
//     for (ineigh = 0; ineigh < SC_MAX (1, elem_data->num_neighbors[iface]);
//          ineigh++) {
//       flux_sum += elem_data->fluxes[iface][ineigh];
//     }
//   }
//   /* Phi^t = dt/dx * (f_(j-1/2) - f_(j+1/2)) + Phi^(t-1) */
//   elem_data->phi_new = (problem->delta_t / elem_data->vol) * flux_sum + phi;
// #if GET_DEBUG_OUTPUT
//   t8_productionf
//     ("advance el with delta_t %f vol %f phi %f  flux %f to %f\n",
//      problem->delta_t, elem_data->vol, phi, flux_sum, elem_data->phi_new);
// #endif
// }

// /* Compute element midpoint and vol and store at element_data field.
//  * tree_vertices can be NULL, if not it should point to the vertex coordinates of the tree */
// static void
// t8_advect_compute_element_data (t8_advect_problem_t * problem,
//                                 t8_advect_element_data_t * elem_data,
//                                 t8_element_t *element,
//                                 t8_locidx_t ltreeid,
//                                 t8_eclass_scheme_c *ts,
//                                 const double *tree_vertices)
// {
//   if (tree_vertices == NULL) {
//     /* Get the vertices of the coarse tree */
//     tree_vertices =
//       t8_cmesh_get_tree_vertices (t8_forest_get_cmesh (problem->forest),
//                                   t8_forest_ltreeid_to_cmesh_ltreeid
//                                   (problem->forest, ltreeid));
//   }
//   /* Compute the midpoint coordinates of element */
//   t8_forest_element_centroid (problem->forest, ltreeid, element,
//                               elem_data->midpoint);
//   /* Compute the volume (length in case of line) of this element */
//   elem_data->vol =
//     t8_forest_element_volume (problem->forest, ltreeid, element);
// }

// bool
// t8_advect_global_conservation_check (double scaled_global_phi_beginning,
//                                      double scaled_global_phi_end)
// {
//   double              eps = 1e-6;
//   double              diff_phi_global =
//     scaled_global_phi_beginning - scaled_global_phi_end;
//   double              abs_diff_phi_global =
//     ((diff_phi_global < 0) ? -diff_phi_global : diff_phi_global);
// #if GET_DEBUG_OUTPUT
//   t8_productionf
//     ("\nglobal conservation check:\n"
//      "global_phi_beginning: %e global_phi_end: %e abs_diff_phi_global: %e, threshhold: %e, return %d\n",
//      scaled_global_phi_beginning, scaled_global_phi_end, abs_diff_phi_global,
//      eps, ((abs_diff_phi_global < eps) ? true : false));
// #endif
//   return ((abs_diff_phi_global < eps) ? true : false);
// }

// bool
// t8_advect_conservation_check_phi (double outgoing_phi, double incoming_phi)
// {
//   double              eps = 1e-6;
//   double              diff_phi = outgoing_phi - incoming_phi;
//   double              abs_diff_phi = ((diff_phi < 0) ? -diff_phi : diff_phi);
// #if GET_DEBUG_OUTPUT
//   t8_productionf
//     ("\nlocal conservation check:\n"
//      "outgoing_phi: %e incoming_phi: %e abs_diff_phi: %e, threshhold: %e, return %d\n",
//      outgoing_phi, incoming_phi, abs_diff_phi, eps,
//      ((abs_diff_phi < eps) ? true : false));
// #endif
//   return ((abs_diff_phi < eps) ? true : false);
// }

// bool
// t8_advect_conservation_check_volume (double outgoing_volume,
//                                      double incoming_volume)
// {
//   double              eps = 1e-6;
//   double              diff_volume = outgoing_volume - incoming_volume;
//   double              abs_diff_volume =
//     ((diff_volume < 0) ? -diff_volume : diff_volume);
// #if GET_DEBUG_OUTPUT
//   t8_productionf
//     ("outgoing_volume: %e incoming_volume: %e abs_diff_volume: %e, threshhold: %e, return %d\n",
//      outgoing_volume, incoming_volume, abs_diff_volume, eps,
//      ((abs_diff_volume < eps) ? true : false));
// #endif
//   return ((abs_diff_volume < eps) ? true : false);
// }

// /* Replace callback to interpolate a refined or coarsened element.
//  * If an element is refined, each child gets the phi value of its parent.
//  * If elements are coarsened, the parent gets the average phi value of the children.
//  */
// /* outgoing are the old elements and incoming the new ones */
// static void
// t8_advect_replace (t8_forest_t forest_old,
//                    t8_forest_t forest_new,
//                    t8_locidx_t which_tree,
//                    t8_eclass_scheme_c *ts,
//                    int isfamily,
//                    int num_outgoing,
//                    t8_locidx_t first_outgoing,
//                    int num_incoming, t8_locidx_t first_incoming)
// {
//   time_interpolation -= sc_MPI_Wtime ();
//   t8_advect_problem_t *problem;
//   t8_advect_element_data_t *elem_data_in, *elem_data_out, *elem_data_in_mem,
//     *elem_data_out_mem;
//   t8_locidx_t         first_incoming_data, first_outgoing_data;
//   t8_element_t       *element, *first_outgoing_elem, *first_incoming_elem,
//     *elem_out_iterate, *elem_in_iterate;;
//   int                 iface;
//   int                 incoming_count;
//   int                 outgoing_count;
//   double              phi_old;
//   const int           num_quater_points = 12;   /* this is only valid for 2D quad with sub scheme */

//   /* Get the problem description */
//   problem = (t8_advect_problem_t *) t8_forest_get_user_data (forest_new);
//   T8_ASSERT (forest_old == problem->forest);
//   T8_ASSERT (forest_new == problem->forest_adapt);

//   /* Get pointers to the element datas */
//   first_incoming_data =
//     first_incoming + t8_forest_get_tree_element_offset (forest_new,
//                                                         which_tree);
//   first_outgoing_data =
//     first_outgoing + t8_forest_get_tree_element_offset (forest_old,
//                                                         which_tree);
//   elem_data_out = (t8_advect_element_data_t *)
//     t8_sc_array_index_locidx (problem->element_data, first_outgoing_data);
//   elem_data_in = (t8_advect_element_data_t *)
//     t8_sc_array_index_locidx (problem->element_data_adapt,
//                               first_incoming_data);

//   /* get the first incoming and outgoing elements */
//   first_outgoing_elem =
//     t8_forest_get_element_in_tree (problem->forest, which_tree,
//                                    first_outgoing);
//   first_incoming_elem =
//     t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
//                                    first_incoming);

// #if GET_DEBUG_OUTPUT            /* for debugging */
// #if T8_ENABLE_DEBUG
//   t8_productionf ("first Element out:\n");
//   ts->t8_element_debug_print (first_outgoing_elem);
//   t8_productionf ("first Element in:\n");
//   ts->t8_element_debug_print (first_incoming_elem);
//   t8_productionf ("num_outgoing: %i  num_incoming %i\n", num_outgoing,
//                   num_incoming);
// #endif
// #endif

//   int                 elem_out_count;
//   double              outgoing_phi = 0, outgoing_volume = 0;
//   for (elem_out_count = 0; elem_out_count < num_outgoing; elem_out_count++) {
//     outgoing_phi +=
//       t8_advect_element_get_phi (problem,
//                                  first_outgoing_data + elem_out_count)
//       * elem_data_out[elem_out_count].vol;
//     outgoing_volume += elem_data_out[elem_out_count].vol;
//   }

//   /* In the following, we check the type of interpolation */
//   bool                elem_to_elem = false;
//   bool                transition_to_transition_same = false;
//   bool                transition_to_transition_diff = false;
//   bool                transition_refined = false;
//   bool                coarsened_to_transition = false;

//   /* TODO: use the maxlevel macro */
//   int                 maxlevel_allowed = 18;
//   if (ts->t8_element_level (first_outgoing_elem) <= maxlevel_allowed - 2) {
//     /* We check the following cases for the quad scheme with subelements in order to improve the element interpolation. 
//      * Otherwise we will just use the standard interpolation. */
//     if (t8_forest_get_tree_class (problem->forest, which_tree) ==
//         T8_ECLASS_HEX) {
//       /* check whether the old element stayed unchanged during the adapting process */
//       if (ts->t8_element_level (first_outgoing_elem) ==
//           ts->t8_element_level (first_incoming_elem)) {

//         if (ts->t8_element_is_subelement (first_outgoing_elem) &&
//             ts->t8_element_is_subelement (first_incoming_elem)) {

//           if (ts->t8_element_get_transition_type (first_outgoing_elem) ==
//               ts->t8_element_get_transition_type (first_incoming_elem)) {
//             transition_to_transition_same = true;
//           }
//           else {
//             transition_to_transition_diff = true;
//           }
//         }

//         if (!ts->t8_element_is_subelement (first_outgoing_elem) &&
//             !ts->t8_element_is_subelement (first_incoming_elem)) {
//           elem_to_elem = true;
//         }
//       }
//       /* check whether a transition cell is refined */
//       if (ts->t8_element_level (first_outgoing_elem) <
//           ts->t8_element_level (first_incoming_elem)) {

//         T8_ASSERT (ts->t8_element_level (first_outgoing_elem) + 1 ==
//                    ts->t8_element_level (first_incoming_elem));

//         if (ts->t8_element_is_subelement (first_outgoing_elem)) {
//           transition_refined = true;
//         }
//       }
//       /* check whether a set of elements is coarsened to a transition cell */
//       if (ts->t8_element_level (first_outgoing_elem) >
//           ts->t8_element_level (first_incoming_elem)) {

//         T8_ASSERT (ts->t8_element_level (first_outgoing_elem) ==
//                    ts->t8_element_level (first_incoming_elem) + 1);

//         if (ts->t8_element_is_subelement (first_incoming_elem)) {
//           coarsened_to_transition = true;
//         }
//       }
//     }
//   }
// /*
//    *
//    * START OF THE INTERPOLATION
//    *
//    */

//   /* Case 1/5: elem_old equals elem_new (either transition cell or standard element, same level) */
//   if (elem_to_elem || transition_to_transition_same) {
//     T8_ASSERT (num_outgoing == num_incoming);

//     for (incoming_count = 0; incoming_count < num_incoming; incoming_count++) {
//       phi_old =
//         t8_advect_element_get_phi (problem,
//                                    first_outgoing_data + incoming_count);
//       /* The element is not changed, copy phi and vol */
//       elem_data_out_mem = (t8_advect_element_data_t *)
//         t8_sc_array_index_locidx (problem->element_data,
//                                   first_outgoing_data + incoming_count);
//       elem_data_in_mem = (t8_advect_element_data_t *)
//         t8_sc_array_index_locidx (problem->element_data_adapt,
//                                   first_incoming_data + incoming_count);
//       memcpy (elem_data_in_mem, elem_data_out_mem,
//               sizeof (t8_advect_element_data_t));

//       t8_advect_element_set_phi_adapt (problem,
//                                        first_incoming_data + incoming_count,
//                                        phi_old);

//       /* Set the neighbor entries to uninitialized */
//       elem_data_in[incoming_count].num_faces =
//         ts->t8_element_num_faces (first_incoming_elem);

//       for (iface = 0; iface < elem_data_in[incoming_count].num_faces; iface++) {
//         elem_data_in[incoming_count].num_neighbors[iface] = 0;
//         elem_data_in[incoming_count].flux_valid[iface] = -1;
//         elem_data_in[incoming_count].dual_faces[iface] = NULL;
//         elem_data_in[incoming_count].fluxes[iface] = NULL;
//         elem_data_in[incoming_count].neighs[iface] = NULL;
//       }                         /* end of iface loop */

//       elem_data_in[incoming_count].level =
//         elem_data_out[incoming_count].level;

//     }                           /* end of num_incoming loop */

//   }                             /* end of if-case 1/5 */

//   /* Case 2/5: elem_old is a transition cell and elem_new is a transition element of an ohter type (same level) */
//   else if (transition_to_transition_diff) {

//     /* Iterate over the of subelements. Copy when they are equal and interpolate when they differ. */
//     int                 subelem_out_count = 0;
//     int                 subelem_in_count = 0;
//     for (subelem_in_count = 0; subelem_in_count < num_incoming;
//          subelem_in_count++) {

//       /* get the recent elements */
//       elem_out_iterate =
//         t8_forest_get_element_in_tree (problem->forest, which_tree,
//                                        first_outgoing + subelem_out_count);
//       elem_in_iterate =
//         t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
//                                        first_incoming + subelem_in_count);

// #if GET_DEBUG_OUTPUT            /* for debugging */
// #if T8_ENABLE_DEBUG
//       t8_productionf ("Out count: %i, In count: %i\n", subelem_out_count,
//                       subelem_in_count);
//       t8_productionf ("Element out:\n");
//       ts->t8_element_debug_print (elem_out_iterate);
//       t8_productionf ("Element in:\n");
//       ts->t8_element_debug_print (elem_in_iterate);
// #endif
// #endif

//       T8_ASSERT (ts->t8_element_is_subelement (elem_out_iterate));
//       T8_ASSERT (ts->t8_element_is_subelement (elem_in_iterate));

// int transition_type_out = ts->t8_element_get_transition_type(elem_out_iterate);
//       int location1[3];
//       int location2[3];
//       ts->t8_element_get_location_of_subelement(elem_out_iterate, location1);
//       ts->t8_element_get_location_of_subelement(elem_in_iterate, location2)
//       /* both subelements are identically */
//       if (location1[1] == location2[1]) {

//         phi_old =
//           t8_advect_element_get_phi (problem,
//                                      first_outgoing_data + subelem_out_count);

//         /* The element is not changed, copy phi and vol */
//         elem_data_out_mem = (t8_advect_element_data_t *)
//           t8_sc_array_index_locidx (problem->element_data,
//                                     first_outgoing_data + subelem_out_count);

//         elem_data_in_mem = (t8_advect_element_data_t *)
//           t8_sc_array_index_locidx (problem->element_data_adapt,
//                                     first_incoming_data + subelem_in_count);

//         memcpy (elem_data_in_mem, elem_data_out_mem,
//                 sizeof (t8_advect_element_data_t));

//         t8_advect_element_set_phi_adapt (problem,
//                                          first_incoming_data +
//                                          subelem_in_count, phi_old);

//         /* Set the neighbor entries to uninitialized */
//         elem_data_in[subelem_in_count].num_faces =
//           ts->t8_element_num_faces (first_incoming_elem);

//         for (iface = 0; iface < elem_data_in[subelem_in_count].num_faces;
//              iface++) {
//           elem_data_in[subelem_in_count].num_neighbors[iface] = 0;
//           elem_data_in[subelem_in_count].flux_valid[iface] = -1;
//           elem_data_in[subelem_in_count].dual_faces[iface] = NULL;
//           elem_data_in[subelem_in_count].fluxes[iface] = NULL;
//           elem_data_in[subelem_in_count].neighs[iface] = NULL;
//         }                       /* end of iface loop */

//         elem_data_in[subelem_in_count].level =
//           elem_data_out[subelem_in_count].level;

//         subelem_out_count++;

//       }
//       /* subelement outgoing is split and subelement in is not */
//       else if ((location1[1] == 1 ) && location2[1] == 0) {

//         // T8_ASSERT (ts->t8_element_get_face_number_of_hypotenuse
//         //            (elem_in_iterate)
//         //            == 1);

//         /* interpolate the recent outgoing subelement and the following to the incoming one */
//         double              phi = 0;
//         double              total_volume = 0;
//         int                 subelement_count;
//         int                 num_subelements = 4;
//         /* Compute average of phi (important in case that a transition cell goes out) */
//         for (subelement_count = 0; subelement_count < num_subelements;
//              subelement_count++) {
//           phi +=
//             t8_advect_element_get_phi (problem,
//                                        first_outgoing_data +
//                                        subelem_out_count +
//                                        subelement_count) *
//             elem_data_out[subelem_out_count + subelement_count].vol;
//           total_volume +=
//             elem_data_out[subelem_out_count + subelement_count].vol;
//         }
//         phi /= total_volume;

//         element =
//           t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
//                                          first_incoming + subelem_in_count);
//         /* Compute midpoint and vol of the new element */
//         t8_advect_compute_element_data (problem,
//                                         elem_data_in + subelem_in_count,
//                                         element, which_tree, ts, NULL);

//         t8_advect_element_set_phi_adapt (problem,
//                                          first_incoming_data +
//                                          subelem_in_count, phi);

//         /* Set the neighbor entries to uninitialized */
//         elem_data_in[subelem_in_count].num_faces =
//           ts->t8_element_num_faces (first_incoming_elem);

//         for (iface = 0; iface < elem_data_in[subelem_in_count].num_faces;
//              iface++) {
//           elem_data_in[subelem_in_count].num_neighbors[iface] = 0;
//           elem_data_in[subelem_in_count].flux_valid[iface] = -1;
//           elem_data_in[subelem_in_count].dual_faces[iface] = NULL;
//           elem_data_in[subelem_in_count].fluxes[iface] = NULL;
//           elem_data_in[subelem_in_count].neighs[iface] = NULL;
//         }                       /* end of iface loop */

//         elem_data_in[subelem_in_count].level =
//           elem_data_out[subelem_in_count].level;

//         subelem_out_count += 4;
//       }
//       /* subelement outgoing is not split and subelement in is */
//       else {

//         // T8_ASSERT (ts->t8_element_get_face_number_of_hypotenuse
//         //            (elem_out_iterate)
//         //            == 1);

//         // T8_ASSERT (ts->t8_element_get_face_number_of_hypotenuse
//         //            (elem_in_iterate)
//         //            == 0);

//         /* copy the value of the outgoing subelement to the recent incoming one and the following */
//         phi_old =
//           t8_advect_element_get_phi (problem,
//                                      first_outgoing_data + subelem_out_count);
//         int                 subelement_count;
//         int                 num_subelements_at_face = 4;
//         for (subelement_count = 0; subelement_count < num_subelements_at_face;
//              subelement_count++) {
//           element =
//             t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
//                                            first_incoming + subelem_in_count +
//                                            subelement_count);
//           /* Compute midpoint and vol of the new element */
//           t8_advect_compute_element_data (problem,
//                                           elem_data_in + subelem_in_count +
//                                           subelement_count, element,
//                                           which_tree, ts, NULL);
//           t8_advect_element_set_phi_adapt (problem,
//                                            first_incoming_data +
//                                            subelem_in_count +
//                                            subelement_count, phi_old);

//           /* Set the neighbor entries to uninitialized */
//           elem_data_in[subelem_in_count + subelement_count].num_faces =
//             ts->t8_element_num_faces (first_incoming_elem);
//           for (iface = 0;
//                iface <
//                elem_data_in[subelem_in_count + subelement_count].num_faces;
//                iface++) {
//             elem_data_in[subelem_in_count +
//                          subelement_count].num_neighbors[iface] = 0;
//             elem_data_in[subelem_in_count +
//                          subelement_count].flux_valid[iface] = -1;
//             elem_data_in[subelem_in_count +
//                          subelement_count].dual_faces[iface] = NULL;
//             elem_data_in[subelem_in_count + subelement_count].fluxes[iface] =
//               NULL;
//             elem_data_in[subelem_in_count + subelement_count].neighs[iface] =
//               NULL;
//           }
//           elem_data_in[subelem_in_count + subelement_count].level =
//             elem_data_out[subelem_in_count].level;

//         }                       /* end of subelement at face loop */

//         subelem_out_count = subelem_out_count + 3;
//         subelem_in_count = subelem_in_count + 3;

//       }

//     }                           /* end of subelement iteration */

//   }                             /* end of if-case 2/5 */

//   /* Case 3/5: elem_old is transition cell and is refined */
//   else if (transition_refined) {

//     /* iterate through all new elements and set their new phi values */
//     for (incoming_count = 0; incoming_count < num_incoming; incoming_count++) {
//       /* get the recent incoming element */
//       elem_in_iterate =
//         t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
//                                        first_incoming + incoming_count);

//       /* compute the vertices of the recent incoming element */
//       int                 elem_num_faces_in =
//         ts->t8_element_num_faces (elem_in_iterate);

//       /* initialization */
//       int                *corner_coords_in_x, *corner_coords_in_y,*corner_coords_in_z, *cap;
//       corner_coords_in_x = T8_ALLOC (int, elem_num_faces_in);
//       corner_coords_in_y = T8_ALLOC (int, elem_num_faces_in);
//       corner_coords_in_z = T8_ALLOC (int, elem_num_faces_in);
//       cap = T8_ALLOC (int, num_outgoing);

//       int                 corner_iterate_in;
//       for (corner_iterate_in = 0;
//            corner_iterate_in < ts->t8_element_num_faces (elem_in_iterate);
//            corner_iterate_in++) {
//         int                 corner_coords[3] = { };
//         ts->t8_element_vertex_coords (elem_in_iterate, corner_iterate_in,
//                                       corner_coords);
//         corner_coords_in_x[corner_iterate_in] = corner_coords[0];
//         corner_coords_in_y[corner_iterate_in] = corner_coords[1];
//         corner_coords_in_z[corner_iterate_in] = corner_coords[2];
//       }

//       for (outgoing_count = 0; outgoing_count < num_outgoing;
//            outgoing_count++) {
//         cap[outgoing_count] = 0;
//         /* get the recent incoming element */
//         elem_out_iterate =
//           t8_forest_get_element_in_tree (problem->forest, which_tree,
//                                          first_outgoing + outgoing_count);

//         T8_ASSERT (ts->t8_element_is_subelement (elem_out_iterate));

//         /* compute the vertices and quater face points of the recent outgoing element */
//         int                 corner_coords_out_x[num_quater_points] = { };
//         int                 corner_coords_out_y[num_quater_points] = { };
//         int                 corner_coords_out_z[num_quater_points] = { };

//         int                 corner_iterate_out;
//         for (corner_iterate_out = 0; corner_iterate_out < ts->t8_element_num_faces (elem_out_iterate); corner_iterate_out++) {  /* iterate over the corners */
//           int                 corner_coords[3] = { };
//           ts->t8_element_vertex_coords (elem_out_iterate, corner_iterate_out,
//                                         corner_coords);
//           corner_coords_out_x[corner_iterate_out] = corner_coords[0];
//           corner_coords_out_y[corner_iterate_out] = corner_coords[1];
//           corner_coords_out_z[corner_iterate_out] = corner_coords[2];
//         }

//         /* TODO think of a better way */
//         corner_coords_out_x[5] =
//           corner_coords_out_x[0] / 2 + corner_coords_out_x[1] / 2 + corner_coords_out_x[2] / 2;
//         corner_coords_out_y[5] =
//           corner_coords_out_y[0] / 2 + corner_coords_out_y[1] / 2;
//         corner_coords_out_z[5] =
//           corner_coords_out_y[0] / 2 + corner_coords_out_y[1] / 2;

//         corner_coords_out_x[4] =
//           corner_coords_out_x[1] / 2 + corner_coords_out_x[2] / 2;
//         corner_coords_out_y[4] =
//           corner_coords_out_y[1] / 2 + corner_coords_out_y[2] / 2;

//         corner_coords_out_x[5] =
//           corner_coords_out_x[2] / 2 + corner_coords_out_x[0] / 2;
//         corner_coords_out_y[5] =
//           corner_coords_out_y[2] / 2 + corner_coords_out_y[0] / 2;

//         corner_coords_out_x[6] =
//           corner_coords_out_x[0] / 2 + corner_coords_out_x[3] / 2;
//         corner_coords_out_y[6] =
//           corner_coords_out_y[0] / 2 + corner_coords_out_y[3] / 2;

//         corner_coords_out_x[7] =
//           corner_coords_out_x[3] / 2 + corner_coords_out_x[1] / 2;
//         corner_coords_out_y[7] =
//           corner_coords_out_y[3] / 2 + corner_coords_out_y[1] / 2;

//         corner_coords_out_x[8] =
//           corner_coords_out_x[1] / 2 + corner_coords_out_x[4] / 2;
//         corner_coords_out_y[8] =
//           corner_coords_out_y[1] / 2 + corner_coords_out_y[4] / 2;

//         corner_coords_out_x[9] =
//           corner_coords_out_x[4] / 2 + corner_coords_out_x[2] / 2;
//         corner_coords_out_y[9] =
//           corner_coords_out_y[4] / 2 + corner_coords_out_y[2] / 2;

//         corner_coords_out_x[10] =
//           corner_coords_out_x[2] / 2 + corner_coords_out_x[5] / 2;
//         corner_coords_out_y[10] =
//           corner_coords_out_y[2] / 2 + corner_coords_out_y[5] / 2;

//         corner_coords_out_x[11] =
//           corner_coords_out_x[5] / 2 + corner_coords_out_x[0] / 2;
//         corner_coords_out_y[11] =
//           corner_coords_out_y[5] / 2 + corner_coords_out_y[0] / 2;

//         /* compare the vertices and midpoints of the recent outgoing element and the recent incoming one */
//         int                 corner_check = 0;
//         int                 coord_iterate_in, coord_iterate_out;
//         for (coord_iterate_in = 0;
//              coord_iterate_in < ts->t8_element_num_faces (elem_in_iterate);
//              coord_iterate_in++) {

//           for (coord_iterate_out = 0; coord_iterate_out < num_quater_points;
//                coord_iterate_out++) {

//             T8_ASSERT (corner_coords_out_x[coord_iterate_out] >= 0
//                        && corner_coords_out_y[coord_iterate_out] >= 0
//                        && corner_coords_in_x[coord_iterate_in] >= 0
//                        && corner_coords_in_y[coord_iterate_in] >= 0);

//             if (corner_coords_out_x[coord_iterate_out] ==
//                 corner_coords_in_x[coord_iterate_in]
//                 && corner_coords_out_y[coord_iterate_out] ==
//                 corner_coords_in_y[coord_iterate_in]) {
//               /* in this case a matching coordinate is found */
//               corner_check++;
//             }
//           }
//         }

//         if (corner_check > 2) {
//           T8_ASSERT (corner_check == 3);
//           /* in this case we have three matching coordinates -> recent elem_out must intersect recent elem_in */
//           cap[outgoing_count] = 1;
//         }

//       }                         /* end of loop over outcoming subelements */

//       double              phi = 0;
//       int                 cap_sum = 0;
//       for (outgoing_count = 0; outgoing_count < num_outgoing;
//            ++outgoing_count) {

//         if (cap[outgoing_count] == 1) {
//           /* get the phi value of the k-th outgoing element */
//           phi +=
//             t8_advect_element_get_phi (problem,
//                                        first_outgoing_data + outgoing_count);
//           cap_sum += 1;
//         }

//       }

//       T8_ASSERT (cap_sum == 1 || cap_sum == 2);

//       phi /= cap_sum;

//       /* Compute midpoint and vol of the new element */
//       t8_advect_compute_element_data (problem, elem_data_in + incoming_count,
//                                       elem_in_iterate, which_tree, ts, NULL);

//       t8_advect_element_set_phi_adapt (problem,
//                                        first_incoming_data + incoming_count,
//                                        phi);

//       /* Set the neighbor entries to uninitialized */
//       elem_data_in[incoming_count].num_faces =
//         ts->t8_element_num_faces (elem_in_iterate);
//       for (iface = 0; iface < elem_data_in[incoming_count].num_faces; iface++) {
//         elem_data_in[incoming_count].num_neighbors[iface] = 0;
//         elem_data_in[incoming_count].flux_valid[iface] = -1;
//         elem_data_in[incoming_count].dual_faces[iface] = NULL;
//         elem_data_in[incoming_count].fluxes[iface] = NULL;
//         elem_data_in[incoming_count].neighs[iface] = NULL;
//       }                         /* end of face loop */

//       elem_data_in[incoming_count].level = elem_data_out[0].level + 1;

//       T8_FREE (corner_coords_in_x);
//       T8_FREE (corner_coords_in_y);
//       T8_FREE (cap);

//     }                           /* end of incoming loop */

//   }                             /* end of if-case 3/5 */

//   /* Case 4/5: elem_old is a set of elements that is coarsened into a transition cell elem_new */
//   else if (coarsened_to_transition) {
//     /* iterate through all new elements and set their phi values */
//     for (incoming_count = 0; incoming_count < num_incoming; incoming_count++) {
//       /* get the recent incoming element */
//       elem_in_iterate =
//         t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
//                                        first_incoming + incoming_count);

// #if T8_ENABLE_DEBUG
//       t8_debugf ("elem_in_itertate:\n");
//       ts->t8_element_debug_print (elem_in_iterate);
// #endif

//       T8_ASSERT (ts->t8_element_is_subelement (elem_in_iterate));

//       /* compute the vertices and quater face points of the recent outgoing element */
//       int                 corner_coords_in_x[num_quater_points] = { };
//       int                 corner_coords_in_y[num_quater_points] = { };

//       int                 corner_iterate_in;
//       for (corner_iterate_in = 0;
//            corner_iterate_in < ts->t8_element_num_faces (elem_in_iterate);
//            corner_iterate_in++) {
//         int                 corner_coords[2] = { };
//         ts->t8_element_vertex_coords (elem_in_iterate, corner_iterate_in,
//                                       corner_coords);
//         corner_coords_in_x[corner_iterate_in] = corner_coords[0];
//         corner_coords_in_y[corner_iterate_in] = corner_coords[1];
//       }

//       /* TODO think of a better way -> maybe just use a lookuptable */
//       corner_coords_in_x[3] =
//         corner_coords_in_x[0] / 2 + corner_coords_in_x[1] / 2;
//       corner_coords_in_y[3] =
//         corner_coords_in_y[0] / 2 + corner_coords_in_y[1] / 2;

//       corner_coords_in_x[4] =
//         corner_coords_in_x[1] / 2 + corner_coords_in_x[2] / 2;
//       corner_coords_in_y[4] =
//         corner_coords_in_y[1] / 2 + corner_coords_in_y[2] / 2;

//       corner_coords_in_x[5] =
//         corner_coords_in_x[2] / 2 + corner_coords_in_x[0] / 2;
//       corner_coords_in_y[5] =
//         corner_coords_in_y[2] / 2 + corner_coords_in_y[0] / 2;

//       corner_coords_in_x[6] =
//         corner_coords_in_x[0] / 2 + corner_coords_in_x[3] / 2;
//       corner_coords_in_y[6] =
//         corner_coords_in_y[0] / 2 + corner_coords_in_y[3] / 2;

//       corner_coords_in_x[7] =
//         corner_coords_in_x[3] / 2 + corner_coords_in_x[1] / 2;
//       corner_coords_in_y[7] =
//         corner_coords_in_y[3] / 2 + corner_coords_in_y[1] / 2;

//       corner_coords_in_x[8] =
//         corner_coords_in_x[1] / 2 + corner_coords_in_x[4] / 2;
//       corner_coords_in_y[8] =
//         corner_coords_in_y[1] / 2 + corner_coords_in_y[4] / 2;

//       corner_coords_in_x[9] =
//         corner_coords_in_x[4] / 2 + corner_coords_in_x[2] / 2;
//       corner_coords_in_y[9] =
//         corner_coords_in_y[4] / 2 + corner_coords_in_y[2] / 2;

//       corner_coords_in_x[10] =
//         corner_coords_in_x[2] / 2 + corner_coords_in_x[5] / 2;
//       corner_coords_in_y[10] =
//         corner_coords_in_y[2] / 2 + corner_coords_in_y[5] / 2;

//       corner_coords_in_x[11] =
//         corner_coords_in_x[5] / 2 + corner_coords_in_x[0] / 2;
//       corner_coords_in_y[11] =
//         corner_coords_in_y[5] / 2 + corner_coords_in_y[0] / 2;

//       /* initialization */
//       int                *cap;
//       cap = T8_ALLOC (int, num_outgoing);
//       for (outgoing_count = 0; outgoing_count < num_outgoing;
//            outgoing_count++) {
//         cap[outgoing_count] = 0;
//         /* get the recent outgoing element */
//         elem_out_iterate =
//           t8_forest_get_element_in_tree (problem->forest, which_tree,
//                                          first_outgoing + outgoing_count);

// #if T8_ENABLE_DEBUG
//         t8_debugf ("elem_out_itertate\n");
//         ts->t8_element_debug_print (elem_out_iterate);
// #endif

//         /* compute the vertices of the recent outgoing element */
//         int                 elem_num_faces_out =
//           ts->t8_element_num_faces (elem_out_iterate);

//         /* initialization */
//         int                *corner_coords_out_x, *corner_coords_out_y;
//         corner_coords_out_x = T8_ALLOC (int, elem_num_faces_out);
//         corner_coords_out_y = T8_ALLOC (int, elem_num_faces_out);

//         int                 corner_iterate_out;
//         for (corner_iterate_out = 0; corner_iterate_out < ts->t8_element_num_faces (elem_out_iterate); corner_iterate_out++) {  /* iterate over the corners */
//           int                 corner_coords[2] = { };
//           ts->t8_element_vertex_coords (elem_out_iterate, corner_iterate_out,
//                                         corner_coords);
//           corner_coords_out_x[corner_iterate_out] = corner_coords[0];
//           corner_coords_out_y[corner_iterate_out] = corner_coords[1];
//         }

//         /* compare the vertices of the recent outgoing element and the recent incoming one */
//         int                 corner_check = 0;
//         int                 coord_iterate_in, coord_iterate_out;
//         for (coord_iterate_in = 0; coord_iterate_in < num_quater_points;
//              coord_iterate_in++) {
//           for (coord_iterate_out = 0;
//                coord_iterate_out <
//                ts->t8_element_num_faces (elem_out_iterate);
//                coord_iterate_out++) {

//             T8_ASSERT (corner_coords_out_x[coord_iterate_out] >= 0
//                        && corner_coords_out_y[coord_iterate_out] >= 0
//                        && corner_coords_in_x[coord_iterate_in] >= 0
//                        && corner_coords_in_y[coord_iterate_in] >= 0);

//             if (corner_coords_out_x[coord_iterate_out] ==
//                 corner_coords_in_x[coord_iterate_in]
//                 && corner_coords_out_y[coord_iterate_out] ==
//                 corner_coords_in_y[coord_iterate_in]) {
//               /* in this case a matching coordinate is found */
//               corner_check++;
//             }
//           }
//         }
//         if (corner_check > 2) {
//           T8_ASSERT (corner_check == 3);
//           /* in this case we have three matching coordinates -> recent elem_out must intersect recent elem_in */
//           cap[outgoing_count] = 1;
//         }

//         T8_FREE (corner_coords_out_x);
//         T8_FREE (corner_coords_out_y);

//       }                         /* end of loop over outcoming subelements */

//       double              phi = 0;
//       int                 cap_sum = 0;
//       double              total_volume = 0;
//       for (outgoing_count = 0; outgoing_count < num_outgoing;
//            outgoing_count++) {
//         if (cap[outgoing_count] == 1) {
//           /* get the recent outgoing element */
//           elem_out_iterate =
//             t8_forest_get_element_in_tree (problem->forest, which_tree,
//                                            first_outgoing + outgoing_count);
//           if (ts->t8_element_is_subelement (elem_out_iterate)) {
//             /* get the phi value of the k-th outgoing element */
//             phi +=
//               t8_advect_element_get_phi (problem,
//                                          first_outgoing_data +
//                                          outgoing_count) *
//               elem_data_out[outgoing_count].vol;
//             total_volume += elem_data_out[outgoing_count].vol;
//           }
//           else {                /* if the recent outgoing element is a quad, then half of it will intersect the incoming triangle */
//             /* get the phi value of the k-th outgoing element */
//             phi +=
//               t8_advect_element_get_phi (problem,
//                                          first_outgoing_data +
//                                          outgoing_count) *
//               (elem_data_out[outgoing_count].vol / 2);;
//             total_volume += elem_data_out[outgoing_count].vol / 2;
//           }
//           cap_sum++;
//         }
//       }

//       T8_ASSERT (cap_sum > 0 && cap_sum <= 6);  /* this is specific for quad subelements */

//       phi /= total_volume;

//       /* Compute midpoint and vol of the new element */
//       t8_advect_compute_element_data (problem, elem_data_in + incoming_count,
//                                       elem_in_iterate, which_tree, ts, NULL);

//       t8_advect_element_set_phi_adapt (problem,
//                                        first_incoming_data + incoming_count,
//                                        phi);

//       /* Set the neighbor entries to uninitialized */
//       elem_data_in[incoming_count].num_faces =
//         ts->t8_element_num_faces (elem_in_iterate);
//       for (iface = 0; iface < elem_data_in[incoming_count].num_faces; iface++) {
//         elem_data_in[incoming_count].num_neighbors[iface] = 0;
//         elem_data_in[incoming_count].flux_valid[iface] = -1;
//         elem_data_in[incoming_count].dual_faces[iface] = NULL;
//         elem_data_in[incoming_count].fluxes[iface] = NULL;
//         elem_data_in[incoming_count].neighs[iface] = NULL;
//       }

//       elem_data_in[incoming_count].level = elem_data_out[0].level - 1;

//       T8_FREE (cap);

//     }                           /* end of incoming loop */

//   }                             /* end of if-case 4/5 */

//   /* Case 5/5: In every other case we will eihter copy the old value to the new elements or compute the mean of old values and copy it to the new element(s).
//    * The other cases are:
//    *   1) elem_old is a standard element and is refined into higher level elements
//    *   2) elem_old is standard element and is refined to a transition cell of the same level
//    *   3) elem_old is a set of elements (family of children or transition cell) that are coarsened to a standard element */
//   else {
//     /* get the mean value of all subelements of the transition cell */
//     double              phi = 0;
//     double              total_volume = 0;
//     /* Compute average of phi (important in case that a transition cell goes out) */
//     for (outgoing_count = 0; outgoing_count < num_outgoing; outgoing_count++) {
//       phi +=
//         t8_advect_element_get_phi (problem,
//                                    first_outgoing_data +
//                                    outgoing_count) *
//         elem_data_out[outgoing_count].vol;
//       total_volume += elem_data_out[outgoing_count].vol;
//     }                           /* end of outgoing loop */

//     phi /= total_volume;

//     /* iterate through all incoming elements and set the new phi value */
//     for (incoming_count = 0; incoming_count < num_incoming; incoming_count++) {
//       /* Get a pointer to the new element */
//       element =
//         t8_forest_get_element_in_tree (problem->forest_adapt, which_tree,
//                                        first_incoming + incoming_count);
//       /* Compute midpoint and vol of the new element */
//       t8_advect_compute_element_data (problem, elem_data_in + incoming_count,
//                                       element, which_tree, ts, NULL);
//       t8_advect_element_set_phi_adapt (problem,
//                                        first_incoming_data + incoming_count,
//                                        phi);
//       /* Set the neighbor entries to uninitialized */
//       elem_data_in[incoming_count].num_faces =
//         ts->t8_element_num_faces (element);
//       for (iface = 0; iface < elem_data_in[incoming_count].num_faces; iface++) {
//         elem_data_in[incoming_count].num_neighbors[iface] = 0;
//         elem_data_in[incoming_count].flux_valid[iface] = -1;
//         elem_data_in[incoming_count].dual_faces[iface] = NULL;
//         elem_data_in[incoming_count].fluxes[iface] = NULL;
//         elem_data_in[incoming_count].neighs[iface] = NULL;
//       }
//       if (ts->t8_element_level (first_outgoing_elem) ==
//           ts->t8_element_level (first_incoming_elem)) {
//         elem_data_in[incoming_count].level = elem_data_out->level;
//       }
//       else if (ts->t8_element_level (first_outgoing_elem) <
//                ts->t8_element_level (first_incoming_elem)) {
//         elem_data_in[incoming_count].level = elem_data_out->level + 1;
//       }
//       else {
//         elem_data_in[incoming_count].level = elem_data_out->level - 1;
//       }

//     }                           /* end of outgoing loop */

//   }                             /* end of if-case 5/5 */

//   double              incoming_phi = 0;
//   double              incoming_volume = 0;
//   for (incoming_count = 0; incoming_count < num_incoming; incoming_count++) {
//     incoming_phi +=
//       t8_advect_element_get_phi_adapt (problem,
//                                        first_incoming_data + incoming_count)
//       * elem_data_in[incoming_count].vol;
//     incoming_volume += elem_data_in[incoming_count].vol;
//   }

//   time_interpolation += sc_MPI_Wtime ();

//   /* Check that conservation is fulfilled for each interpolation step. 
//    * We require, that: 
//    *    1) the volume difference of the sum of all incoming and outgoing elements is at most 1e-12
//    *    2) the volume scaled phi difference of the sum of all incoming and outgoing elements is at most 1e-12
//    */
//   SC_CHECK_ABORT (t8_advect_conservation_check_volume
//                   (outgoing_volume, incoming_volume),
//                   "Conservation check volume failed at the end of t8_advect_replace");
//   SC_CHECK_ABORT (t8_advect_conservation_check_phi
//                   (outgoing_phi, incoming_phi),
//                   "Conservation check phi failed at the end of t8_advect_replace");

// }                               /* end of t8_advect_replace */





 



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
  

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();