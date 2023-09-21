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

/* The data that we want to store for each element.
 * In this example we want to store the element's level and volume. */
struct t8_data_per_element
{
  int                 level;
  double              volume;
  int                 height;
};



struct t8_adapt_data adapt_data = {
    double midpoint[3]                = {0.5, 0.5, 1},              /* Midpoints of the sphere. */
    double refine_if_inside_radius    = 0.2,                        /* Refine if inside this radius. */
    double coarsen_if_outside_radius  = 0.4                         /* Coarsen if outside this radius. */
  };



//--------------------vtu output ----------------
static void
t8_output_data_to_vtu (t8_forest_t forest,
                             struct t8_data_per_element *data,
                             const char *prefix)
{
  t8_locidx_t         num_elements =
    t8_forest_get_local_num_elements (forest);
  t8_locidx_t         ielem;
  /* We need to allocate a new array to store the volumes on their own.
   * This array has one entry per local element. */
  double             *element_volumes = T8_ALLOC (double, num_elements);
  double             *element_heights = T8_ALLOC (double, num_elements);
  /* The number of user defined data fields to write. */
  int                 num_data = 1;
  /* For each user defined data field we need one t8_vtk_data_field_t variable */
  t8_vtk_data_field_t vtk_data;
  /* Set the type of this variable. Since we have one value per element, we pick T8_VTK_SCALAR */
  vtk_data.type = T8_VTK_SCALAR;
  /* The name of the field as should be written to the file. */
  strcpy (vtk_data.description, "Element volume");
  vtk_data.data = element_volumes;
  strcpy (vtk_data.description, "Element height");
  vtk_data.data = element_heights;
  /* Copy the elment's volumes from our data array to the output array. */
  for (ielem = 0; ielem < num_elements; ++ielem) {
    element_volumes[ielem] = data[ielem].volume;
    element_heights[ielem] = data[ielem].height;
  }
  
  {
    /* To write user defined data, we need to extended output function t8_forest_vtk_write_file
     * from t8_forest_vtk.h. Despite writin user data, it also offers more control over which 
     * properties of the forest to write. */
    int                 write_treeid = 1;
    int                 write_mpirank = 1;
    int                 write_level = 1;
    int                 write_element_id = 1;
    int                 write_ghosts = 0;
    t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank,
                             write_level, write_element_id, write_ghosts,
                             0, 0, num_data, &vtk_data);
  }
  T8_FREE (element_volumes);
  T8_FREE (element_heights);
}

static struct t8_data_per_element *
t8_create_element_data (t8_forest_t forest)
{
  t8_locidx_t         num_local_elements;
  t8_locidx_t         num_ghost_elements;
 struct t8_data_per_element *element_data;

  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Get the number of local elements of forest. */
  num_local_elements = t8_forest_get_local_num_elements (forest);
  t8_productionf("get local num elements %i\n", num_local_elements);
  t8_productionf("global num subelements  %i\n", forest->global_num_subelements);
  /* Get the number of ghost elements of forest. */
  num_ghost_elements = t8_forest_get_num_ghosts (forest);
t8_productionf("get num ghost elements %i\n", num_ghost_elements);
  /* Now we need to build an array of our data that is as long as the number
   * of elements plus the number of ghosts. You can use any allocator such as
   * new, malloc or the t8code provide allocation macro T8_ALLOC. 
   * Note that in the latter case you need
   * to use T8_FREE in order to free the memory.
   */
  element_data =
    T8_ALLOC (struct t8_data_per_element,
              num_local_elements + num_ghost_elements);

  {
    t8_locidx_t         itree, num_local_trees;
    t8_locidx_t         current_index;
    t8_locidx_t         ielement, num_elements_in_tree;
    t8_eclass_t         tree_class;
    t8_eclass_scheme_c *eclass_scheme;
    const t8_element_t *element;

    /* Get the number of trees that have elements of this process. */
    num_local_trees = t8_forest_get_num_local_trees (forest);
    t8_productionf("get num local trees %i\n", num_local_trees);
    for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
      tree_class = t8_forest_get_tree_class (forest, itree);
      eclass_scheme = t8_forest_get_eclass_scheme (forest, tree_class);
      /* Get the number of elements of this tree. */
      num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);
      t8_productionf("get num of elements in  tree %i\n", num_elements_in_tree);
      for (ielement = 0; ielement < num_elements_in_tree;
           ++ielement, ++current_index) {
        /* This loop iterates through all the local elements of the forest in the current tree. */

        element = t8_forest_get_element_in_tree (forest, itree, ielement);

      const double        x = 0.5;
      const double        y = 0.5;
      const double        z = 0.5;
      const double        r = sqrt (x * x + y * y) * 20.0;      
       
        /* We want to store the elements level and its volume (and its height) as data. We compute these
         * via the eclass_scheme and the forest_element interface. */
        element_data[current_index].level =
          eclass_scheme->t8_element_level (element);
        element_data[current_index].volume =
          t8_forest_element_volume (forest, itree, element);
           t8_productionf("Das volumen ist %f \n",t8_forest_element_volume (forest, itree, element) );
        element_data[current_index].height = sin (2.0 * r) / r;
      }
    }
  }
  return element_data;
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

 double              centroid[3];  
 const struct t8_adapt_data *adapt_data =
    (const struct t8_adapt_data *) t8_forest_get_user_data (forest);

 double              dist;

  T8_ASSERT (adapt_data != NULL);


/* If a subelement is given, we apply the callback function to its parent */
  if (ts->t8_element_is_subelement (elements[0])) {
    t8_element_t      **parent = T8_ALLOC (t8_element_t *, 1);
    ts->t8_element_new (1, parent);
    ts->t8_element_parent (elements[0], parent[0]);
    T8_FREE (parent);
  }

  /* Compute the element's centroid coordinates. */
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);

  /* Compute the distance to our sphere midpoint. */
  dist = t8_vec_dist (centroid, adapt_data->midpoint);
  if (dist < adapt_data->refine_if_inside_radius) {
    /* Refine this element. */
    return 1;
  }
  else if (is_family && dist > adapt_data->coarsen_if_outside_radius) {
    /* Coarsen this family. Note that we check for is_family before, since returning < 0
     * if we do not have a family as input is illegal. */
    return -1;
  }
  /* Do not change this element. */
  return 0;
}
//   if ((lelement_id == 0) ) {
//     /* Refine this element. */
//     return 1;
//   }
//   /* Do not change this element. */
//   return 0;
// }


void t8_transition (void)
{

  t8_eclass_t         eclass = T8_ECLASS_HEX;  /* depending on the include file, this will be the transitioned or default hex implementation */

  t8_forest_t         forest;
  t8_forest_t         forest_adapt;
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];
  t8_eclass_scheme_c *ts; 
  int level = 3;

  t8_data_per_element *data;
  t8_adapt_data adapt_data;


  int set_balance = 1;
  int set_transition = 1;


  const char         *prefix_forest_with_data = "t8_forest_with_volume_data";
/* Creating cmesh as 3D hypercube with eclass = hex */
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

  t8_forest_init (&forest_adapt);
  //t8_forest_set_user_data (forest_adapt, &adapt_data);
  
  //t8_forest_set_adapt (forest_adapt, forest, t8_adapt_callback, 1);
  

if (set_balance && !set_transition) {
      t8_forest_set_balance (forest_adapt, forest, 0);
    }
    if (set_transition) {
      t8_forest_set_transition (forest_adapt, forest, set_balance);
    }

//t8_forest_commit (forest_adapt); 
forest_adapt = t8_forest_new_adapt (forest, t8_adapt_callback, 0, 0, &adapt_data);
data = t8_create_element_data (forest_adapt);
  t8_output_data_to_vtu (forest_adapt, data, prefix_forest_with_data);
T8_FREE (data);

  // t8_output_data_to_vtu (forest, data, prefix_forest_with_data);
  // t8_global_productionf (" Wrote forest and volume data to %s*.\n",
  //                        prefix_forest_with_data);

 //t8_productionf("das ist data %d\n", &data);
t8_forest_unref (&forest_adapt);
}




/* Write the forest as vtu and also write the element's volumes in the file.
 * 
 * t8code supports writing element based data to vtu as long as its stored
 * as doubles. Each of the data fields to write has to be provided in its own
 * array of length num_local_elements.
 * We support two types: T8_VTK_SCALAR - One double per element
 *                  and  T8_VTK_VECTOR - 3 doubles per element
 */


int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_forest_t         forest;
  /* The prefix for our output files. */
  const char         *prefix_forest = "t8_forest";
  /* The uniform refinement level of the forest. */
  const int           level = 1;
  /* The array that will hold our per element data. */
  t8_data_per_element *data;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;

  /*
   * Setup.
   * Build cmesh and uniform forest.
   * Adapt forest similar to step 3 & 4.
   */
  t8_transition();

  // t8_forest_write_vtk (forest, prefix_forest);
  // t8_global_productionf (" [step5] Wrote forest to vtu files: %s*\n",
  //                        prefix_forest);

  /*
   * Build data array and gather data for the local elements.
   */
  // data = t8_create_element_data (forest);

  // t8_global_productionf
  //   (" Computed level and volume data for local elements.\n");
  // if (t8_forest_get_local_num_elements (forest) > 0) {
  //   /* Output the stored data of the first local element (if it exists). */
  //   t8_global_productionf (" Element 0 has level %i and volume %e.\n",
  //                          data[0].level, data[0].volume);
  // }


  /*
   * Output the volume data to vtu.
   */

  // t8_global_productionf (" Wrote forest and volume data to %s*.\n",
  //                        prefix_forest_with_data);

  /*
   * clean-up
   */

  /* Free the data array. */
  T8_FREE (data);


  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();