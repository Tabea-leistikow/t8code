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

/* Description:
 * This is the low-level structure of 3D hexahedral elements with transition
 * cells of pyramidal subelements. */

#include <t8.h>
#include <p8est_bits.h>
#include <t8_schemes/t8_default/t8_default_line/t8_dline_bits.h>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>
#include "t8_transition_conformal_hex_cxx.hxx"
#include <cmath>

/* *INDENT-OFF* */
/* Connectivity of subelement faces depends on which type of subelement is considered (6 pyramids = 6 types): 
 *         |
 *    (Z)  | / (Y)
 *         |/_____ (X)
 *
 *     Type 0 = quadliteral face at face 0 from the hexahdron:    Type 3:
 *     f_0 <-> f_0                                                f_0 <-> f_1
 *     f_1 <-> f_0                                                f_1 <-> f_3
 *     f_2 <-> f_0                                                f_2 <-> f_3
 *     f_3 <-> f_0                                                f_3 <-> f_3 
 *     f_4 <-> f_4 (assuming a neighboring transition cell)       f_4 <-> f_4
 *     Type 1:                                                    Type 4:
 *     f_0 <-> f_1                                                f_0 <-> f_2
 *     f_1 <-> f_1                                                f_1 <-> f_2
 *     f_2 <-> f_1                                                f_2 <-> f_2
 *     f_3 <-> f_1                                                f_3 <-> f_2
 *     f_4 <-> f_4 (assuming a neighboring transition cell)       f_4 <-> f_4 
 *     Type 0:                                                    Type 5:
 *     f_0 <-> f_0                                                f_0 <-> f_3
 *     f_1 <-> f_2                                                f_1 <-> f_3
 *     f_2 <-> f_0                                                f_2 <-> f_3
 *     f_3 <-> f_2                                                f_3 <-> f_3
 *     f_4 <-> f_4 (assuming a neighboring transition cell)       f_4 <-> f_4
 */

const int           subelement_face_dual[6][5] = {
  {0, 0, 0, 0, 4},
  {1, 1, 1, 1, 4},
  {0, 2, 0, 2, 4},
  {1, 3, 3, 3, 4},
  {2, 2, 2, 2, 4},
  {3, 3, 3, 3, 4}
  };

/* Connectivity of a subelements location within a transition cell 
 * and the parent hexs faces:
 *     location[0] = 0 -> parents dual face = 1
 *     location[0] = 1 -> parents dual face = 0
 *     location[0] = 2 -> parents dual face = 3
 *     location[0] = 3 -> parents dual face = 2
 *     location[0] = 4 -> parents dual face = 5
 *     location[0] = 5 -> parents dual face = 4 */
const int           subelement_location_to_parent_dual_face[6] = { 1, 0, 3, 2, 5, 4 };

/* Connectivity of a subelements location within a transition cell (only if the subelements are not splitted)
 *      --> gets subelement_duals
 */
const int           subelement_face_to_dual_subelement[6][5] = {
  { 2, 3, 4, 5, -1 };
  { 2, 3, 4, 5, -1 };
  { 0, 1, 4, 5, -1 };
  { 0, 1, 4, 5, -1 };
  { 0, 1, 2, 3, -1 };
  { 0, 1, 2, 3, -1 };
}
/* Connectivity of a subelements location within a transition cell 
 * and the parent hexs faces: Wie in mit Koordinatensystem in Davids Masterarbeit
 * starte links, dann rechts, vorne, hinten, unten oben.
 *     location[0] = 0 (left)   -> parents face = 0
 *     location[0] = 1 (right)  -> parents face = 1
 *     location[0] = 2 (front)  -> parents face = 2
 *     location[0] = 3 (back)   -> parents face = 3
 *     location[0] = 4 (bottom) -> parents face = 4
 *     location[0] = 5 (up)     -> parents face = 5  
 */
const int           subelement_location_to_parent_face[6] = { 0, 1, 2, 3, 4, 5};
/* *INDENT-ON* */

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* This function is used by other element functions and we thus need to
 * declare it up here */
t8_linearidx_t      t8_element_get_linear_id (const t8_element_t *elem,
                                              int level);

int
t8_subelement_scheme_hex_c::t8_element_maxlevel (void) const
{
  return P8EST_QMAXLEVEL;
}

/* *INDENT-OFF* */
t8_eclass_t
t8_subelement_scheme_hex_c::t8_element_child_eclass (int childid) const
/* *INDENT-ON* */

{
  T8_ASSERT (0 <= childid && childid < P8EST_CHILDREN);

  return T8_ECLASS_HEX;
}

int
t8_subelement_scheme_hex_c::t8_element_level (const t8_element_t *elem) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  T8_ASSERT (t8_element_is_valid (elem));
  return (int) ((const t8_hex_with_subelements *) phex_w_sub)->p8q.level;
}

//Wieso ist nur q const aber nicht r? Ode bezieht sich das const auf beide 
/* JM: `r` ist die Ausgabevariable in die geschrieben bzw. rueberkopiert wird. */
static void
t8_element_copy_surround (const p8est_quadrant_t * q, p8est_quadrant_t * r)
{
  T8_HEX_SET_TDIM (r, T8_HEX_GET_TDIM (q));
  if (T8_HEX_GET_TDIM (q) == 3) {
    T8_HEX_SET_TNORMAL (r, T8_HEX_GET_TNORMAL (q));
    T8_HEX_SET_TCOORD (r, T8_HEX_GET_TCOORD (q));
  }
}

void
t8_subelement_scheme_hex_c::t8_element_copy (const t8_element_t *source,
                                              t8_element_t *dest) const
{
  const t8_hex_with_subelements *phex_w_sub_source =
    (const t8_hex_with_subelements *) source;
  t8_hex_with_subelements *phex_w_sub_dest =
    (t8_hex_with_subelements *) dest;

  const p8est_quadrant_t *q = &phex_w_sub_source->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_dest->p8q;

  T8_ASSERT (t8_element_is_valid (source));
  T8_ASSERT (t8_element_is_valid (dest));
  if (q == r &&
      phex_w_sub_source->transition_type ==
      phex_w_sub_dest->transition_type &&
      phex_w_sub_source->subelement_id == phex_w_sub_dest->subelement_id) {
    /* Do nothing if they are already the same hexahedra. */
    return;
  }
  *r = *q;

  t8_element_copy_subelement_values (source, dest);
  t8_element_copy_surround (q, r);
}

int
t8_subelement_scheme_hex_c::t8_element_compare (const t8_element_t *elem1,
                                                 const t8_element_t *elem2)
  const
{
  const t8_hex_with_subelements *phex_w_sub_elem1 =
    (const t8_hex_with_subelements *) elem1;
  const t8_hex_with_subelements *phex_w_sub_elem2 =
    (const t8_hex_with_subelements *) elem2;

  const p8est_quadrant_t *q = &phex_w_sub_elem1->p8q;
  const p8est_quadrant_t *r = &phex_w_sub_elem2->p8q;

  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));

  int                 compare = p8est_quadrant_compare (q, r);

  if (compare == 0 && (t8_element_is_subelement (elem1)
                       || t8_element_is_subelement (elem2))) {
    t8_debugf ("Caution, t8_element_compare is used with subelements.\n");
    if (t8_element_is_subelement (elem1)
        && t8_element_is_subelement (elem2)) {
      /* Caution: The compare function is used for two subelements. */

      if (phex_w_sub_elem1->transition_type ==
          phex_w_sub_elem2->transition_type
          && phex_w_sub_elem1->subelement_id ==
          phex_w_sub_elem2->subelement_id) {
        /* both subelements are identical */
        return 0;
      }
      /* return != 0 to avoid debug abortion in t8_ghost_add_remote */
      return 1;
    }
    else if (t8_element_is_subelement (elem1)) {
      return -1;                /* elem1 is subelement and therefore smaller */
    }
    else if (t8_element_is_subelement (elem2)) {
      return 1;                 /* elem2 is subelement and therefore smaller */
    }
  }

  /* Note that for subelements, their parent quadrant is compared at this point */
  return compare;

  // SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_parent (const t8_element_t *elem,
                                                t8_element_t *parent) const
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_parent =
    (t8_hex_with_subelements *) parent;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_parent->p8q;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (parent));

  if (t8_element_is_subelement (elem)) {
    phex_w_sub_parent->p8q = phex_w_sub_elem->p8q;
  }
  else {
    p8est_quadrant_parent (q, r);
  }

  /* the parent of any element will never be a subelement */
  t8_element_reset_subelement_values (parent);

  t8_element_copy_surround (q, r);

  // SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_sibling (const t8_element_t *elem,
                                                 int sibid,
                                                 t8_element_t *sibling) const
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_sibling =
    (t8_hex_with_subelements *) sibling;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_sibling->p8q;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (sibling));

  p8est_quadrant_sibling (q, r, sibid);
  t8_element_copy_surround (q, r);

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_num_faces (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  return (t8_element_is_subelement (elem) ? T8_HEX_SUBELEMENT_FACES :
          P8EST_FACES);
}

int
t8_subelement_scheme_hex_c::t8_element_max_num_faces (const t8_element_t
                                                       *elem) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  return P8EST_FACES;
}

int
t8_subelement_scheme_hex_c::t8_element_num_children (const t8_element_t
                                                      *elem) const
{
  /* Note that children of subelements equal the children of the parent hexahedra. 
   * Therefore, the number of children of a subelement equals P8EST_CHILDREN */
  T8_ASSERT (t8_element_is_valid (elem));
  return P8EST_CHILDREN;
}

/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */

int   
t8_subelement_scheme_hex_c::t8_element_num_siblings (const t8_element_t *
                                                      elem) const
/* *INDENT-ON* */

{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  //no hanging nodes 
  if (phex_w_sub->transition_type == 0){
    //eigentlich normale funktion t8_element_num_siblings
    //könnte auch 2hoch dim returnen 
    return P8EST_CHILDREN;
  }

  int                 num_hanging_faces = 0;
  int                 iface;
  for (iface = 0; iface < P8EST_FACES; iface++) {     
   /* Count the number of ones of the binary transition type.
    * This number equals the number of hanging faces. */

   // binary shift << 1 Left-shift, d.h. *2¹ 
   // Right shift >> 1 Right-shift, d.h. *2⁻¹
   // & (bitwise AND operator, num_hanging_faces wird nur erhöht, wenn mind. ein bit von transition_type != 0 ist. )
    num_hanging_faces +=
      (phex_w_sub->transition_type & (1 << iface)) >> iface;
  }

  return P8EST_CHILDREN + 3*num_hanging_faces;

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_num_face_children (const t8_element_t
                                                           *elem,
                                                           int face) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_is_valid (elem));

  /*  if we use hex scheme without set_transition, then we are only balanced 
  *   and four neighbors are possible */
  return 4;
}

/*
In header file steht, dass die Funktion anzahl von children zurück gibt,
hier aber Frage ob is sibling ja oder nein*/
  /** Check whether the neighbors of an element at a specic face are siblings
   *  \param [in] elem A valid element 
   *  \param [in] elem_face A valid face 
   *  \return true if the neighbor of elem at face elem_face is a sibling.
   */
//Wenn elements transition cell, dann immer alle seiten bis auf die Grundseite (f4)
int
t8_subelement_scheme_hex_c::t8_element_neighbor_is_sibling (const
                                                             t8_element_t
                                                             *elem,
                                                             const int face)
  const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_subelement (elem));

  if (face == 0 || face == 1 || face == 2 || face == 3) {
    return 1;
  }

  return 0;
  // SC_ABORT_NOT_REACHED();
}

/** Check whether the neighbors of an element at a specic face are siblings
   *  \param [in] elem A valid element 
   *  \param [in] elem_face A valid face 
   *  \return return the number of sibling neighbors at a given face.
   */
int
t8_subelement_scheme_hex_c::t8_element_get_num_sibling_neighbors_at_face (const t8_element_t *elem,
                                                                           const int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_subelement (elem));
  T8_ASSERT (face == 0 || face == 1 || face == 2 || face == 3);
  
  return 1;
  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_get_transition_refine_identifier () const
{
  return T8_TRANSITION_CONFORMAL_HEX_REFINE_FUNCTION;
}
/* *INDENT-ON* */

int
t8_subelement_scheme_hex_c::t8_element_get_face_corner (const t8_element_t
                                                         *elem, int face,
                                                         int corner) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  if (!t8_element_is_subelement (elem)) {
    /*
     *   2    f_2    3
     *     x -->-- x
     *     |       |
     *     ^       ^
     * f_0 |       | f_1
     *     x -->-- x
     *   0    f_3    1
     */
  //face has to be between 0 and 4
  //Corner has to be between 0 and 4
    T8_ASSERT (0 <= face && face < P8EST_FACES);
    T8_ASSERT (0 <= corner && corner < 8);

    return p8est_face_corners[face][corner];
  }
  else {
    int                 t8_face_corners_subelement[5][4] = {
      {0, 2, 4 ,-1},  //f0
      {1, 3, 4, -1},  //f1  
      {0, 1, 4, -1},  //f2
      {2, 3, 4, -1},  //f3
      {0, 1, 2, 3}    //f4
    };
    T8_ASSERT (0 <= face && face < T8_HEX_SUBELEMENT_FACES);
    T8_ASSERT (0 <= corner && corner < 3);

    return t8_face_corners_subelement[face][corner];
  }

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_get_corner_face (const t8_element_t
                                                         *elem, int corner,
                                                         int face) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= corner && corner < P8EST_CHILDREN);
  T8_ASSERT (0 <= face && face < 2);

  return p8est_corner_faces[corner][face];
}


void
t8_subelement_scheme_hex_c::t8_element_child (const t8_element_t *elem,
                                               int childid,
                                               t8_element_t *child) const
{
  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  /*
   *
   *         x - - - - - x        x - - x - - x 
   *         |           |        |     |     |
   *         |           |        |  2  |  3  |
   *         |   elem    |   =>   x - - x - - x
   *         |           |        |     |     |
   *         |           |        |  0  |  1  |
   *         x - - - - - x        x - - x - - x
   * 
   */
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_child =
    (t8_hex_with_subelements *) child;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_child->p8q;

  const p4est_qcoord_t shift = P8EST_QUADRANT_LEN (q->level + 1);

  /* it should not be possible to construct a child of a subelement */
  T8_ASSERT (!t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (child));
  T8_ASSERT (p8est_quadrant_is_extended (q));
  T8_ASSERT (q->level < P8EST_QMAXLEVEL);

  T8_ASSERT (childid >= 0 && childid < P8EST_CHILDREN);
// 0x01 is the least significant bit set. 
// | is the bitwise OR: The result of OR is 1 if any of the two bits is 1. 
  r->x = childid & 0x01 ? (q->x | shift) : q->x;
  r->y = childid & 0x02 ? (q->y | shift) : q->y;
  r->z = childid & 0x04 ? (q->z | shift) : q->z;
// Da r Kind von q muss Level von r um eins höher sein. 
  r->level = q->level + 1;

  T8_ASSERT (p8est_quadrant_is_parent (q, r));

  t8_element_reset_subelement_values (child);

  t8_element_copy_surround (q, r);

  // SC_ABORT_NOT_REACHED();
}

/* TODO: Implement this function. */
// void
// t8_subelement_scheme_hex_c::t8_element_get_sibling_neighbor_in_transition_cell (const t8_element_t
//                                                                                  *elem, const int face,
//                                                                                  const int num_neighbors,
//                                                                                  t8_element_t
//                                                                                  *neighbor_at_face[],
//                                                                                  int *neigh_face[])
// {
//   //SC_ABORT_NOT_REACHED ();
// }

/* FUNKTION NOCH NICHT RICHTIG !!! */
//Ein subelement kann mehr als einen face_neighbor haben! Wenn zum Bsp. ich selbst nicht gesplittet bin aber das subelement an dem mein face liegt schon. 
void
t8_subelement_scheme_hex_c::t8_element_get_sibling_neighbor_in_transition_cell (const t8_element_t
                                                                                 *elem, const int face,
                                                                                 const int num_neighbors,
                                                                                 t8_element_t
                                                                                 *neighbor_at_face[],
                                                                                 int *neigh_face[])
{
  T8_ASSERT (t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_neighbor_is_sibling (elem, face));
  T8_ASSERT (num_neighbors == 1);
  T8_ASSERT (t8_element_is_valid (neighbor_at_face[0]));

  /* source = elem, destination = neighbor at face. --> They have the same anchor node + Morton index + level.*/
  t8_element_copy (elem, neighbor_at_face[0]);
/* Expand neighbor_at_face to a subelement (subelement_id + transititon_type) + copy it into phex_w_sub_neighbor_at_face*/
  t8_hex_with_subelements *
    phex_w_sub_neighbor_at_face =
    (t8_hex_with_subelements *) neighbor_at_face[0];
  //iterator variable
  int iter;
  //Get informations about the location of the subelement.
  int
    num_siblings = t8_element_num_siblings (phex_w_sub_neighbor_at_face);
  int                 location[3] = { };
  //( location[0] = face_number of transition cell, location[1] = if splitted or not ( 1 = splitted ), location[2] = sub_id
    t8_element_get_location_of_subelement (phex_w_sub_neighbor_at_face, location);
  
  //Create a temporary variable to store the possible subelement_id of the neighbor
  int subelement_id_tmp = 0;
  int transition_type_tmp = 0;
  int amount_subelements;
  /* There are 4 cases that can happen:
   * 1. The subelement itself is not splitted, and its face neighbor is also not splitted.
   * 2. The subelement itself is not splitted, but its face neighbor is splitted. (Two face neighbors)
   * 3. The subelement itself is splitted, but its face neighbor is not splitted.
   * 4. The subelement itself is splitted, and its face neighbor is also splitted.
  */
  //First check if (the own) face is splitted.
  if(location[1] == 0){
    //get hex_face_number of the face_neighbored subelement
    *neigh_face[0] = subelement_face_dual[location[0]][face];
    subelement_id_tmp = subelement_face_to_dual_subelement[location[0]][*neigh_face[0]];
    /*Check if the dual subelement is splitted. If  it's splitted,the element has two neighbors 
    * as siblings. Then we always take the left, front or down (in this order) subelement.
    * when the transition type is = 1 at the hex_face, it's splitted.*/


//-----------------------------CASE 1--------------------------------------------------------
    if(phex_w_sub_neighbor_at_face->transition_type >> subelement_id_tmp == 0){ //neighbor not splitted
      for(iter = 0; iter <  *neigh_face[0]; iter++){
        //make rightshift until only the bits for the faces before our neighbors face are left.
          transition_type_tmp = phex_w_sub_neighbor_at_face->transition_type >> (5 - *neigh_face[0]);
        //Count the elemnts until our neighbored hex_face.
        if(transition_type_tmp >> iter == 1){
            amount_subelements += 4;
        }
        else{
          amount_subelements +=1;
        }
      }
      //The subelement_id of the neighbor is then the amount of subelements till then - 1
      subelement_id_tmp = amount_subelements -1 ;
    }


//-------------------------CASE 2--------------------------------------
    else{ 

    }
  }


  //-------------------------CASE 3--------------------------------------
    else{
    //get hex_face_number of the face_neighbored subelement
    *neigh_face[0] = subelement_face_dual[location[0]][face];
    //It's possible, that the neighbored subelement has the same face_hex number as the element itself
    //Now we need the subelement_id_type (location[2]) to determine the exact location of the 
    //subelement in the transition cell 

    /* We have to go through all hex_faces
    */
    if(location[0] == 0 || location[0] == 1){ // hex_face = 0,1
      if((location[2] & 2) == 1 ){//back  
          if(face == 0){ //then it's the element before.
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 1;
          } 
          if(location[2] & 1 == 1){//up
            if(face == 2){ //then it's the element below
              subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 2;
            }
          }
          else{ //down
            if(face == 3){// down
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 2;
            }
          } 
      }
      else{//front 
        if( face == 1){ //then it's the element after.
          subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 1;
        }
        if(location[2] & 1 == 1){//up
            if(face == 2 ){ //then it's the element below
              subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 2;
            }
          }
          else{ //down
            if(face == 3){//then its the element above
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 2;}
          }
      } 
    }
    if(location[0] == 2 || location[0] == 3){ // hex_face = 2,3
      if((location[2] & 4) == 1 ){//left 
          if(face == 1){ //then it's the next element.
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 1;
          } 
          if((location[2] & 1) == 1){//up
            if(face == 2){ //then it's the element below
              subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 2;
            }
          }
          else{ //down
            if(face == 3){// down
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 2;
            }
          } 
      }
      else{//right
        if( face == 0){ //then it's the element before.
          subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 1;
        }
        if((location[2] & 1) == 1){//up
            if(face == 2 ){ //then it's the element below
              subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 2;
            }
          }
          else{ //down
            if(face == 3){//then its the element above
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 2;}
          }
      } 
    }
    if(location[0] == 4 || location[0] == 5){ // hex_face = 4,5
      if((location[2] & 4) == 1 ){//left 
          if(face == 1){ //then it's the next element.
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 1;
          } 
          if((location[2] & 2) == 1){//back
            if(face == 2){ //then it's the element in front
              subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 2;
            }
          }
          else{ //front
            if(face == 3){
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 2;
            }
          } 
      }
      else{//right
        if( face == 0){ //then it's the element before.
          subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 1;
        }
        if((location[2] & 2) == 1){//back
            if(face == 2 ){ //then it's the element below
              subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id - 2;
            }
          }
          else{ //front
            if(face == 3){//then its the element above
            subelement_id_tmp = phex_w_sub_neighbor_at_face->subelement_id + 2;}
          }
      } 
    }
      
  } 
//-------------------------CASE 4--------------------------------------
//   else{
//     //First check hex_faces 0 and 1. 
//     //We need to know where exactly the subelements is. So we have to "decode" the sub_id (location[2])
//     if( location[0] == 0 || location[0] == 1){
//       if(location[2] == 0){
//         if( face % 2 == 0){
//           subelement_id_tmp = [location[0]][face];
//         }

//       }
//       if(location[2] == 0){
  
//       }
//       if(location[2] == 0){
  
//     }
//       if(location[2] == 0){
  
//     }
//     }
//   }

//   int                 num_subelements_faces = 0;
//   int                 iface = 0;
//   int                 hex_face = 0;
//   int                 ones_in_transition_type = 0;
//   int                 left_side, right_side, bottom_side, up_side, front_side, back_side = 0;
  
//   for (iface; iface < P8EST_FACES; iface++) {     
//    /* Check whether a subelement is splitted into four pyramids or not. If transition type at
//     * place 0 is equal to one, we know that there are four pyramids at face 0.
//     * If it's equal to 0, we know that there is one pyramid at this face. */
//     if(phex_w_sub_neighbor_at_face->transition_type & (1 << iface) == 1){
//         num_subelements += 4;
//         ones_in_transition_type += 1;
//     }
//     else{
//       num_subelements += 1;
//     }
//     /* Check whether the subelement ID is greater or equal than the number of pyramids till
//     *  then. If so, we know that the subelement is at face iface. We store this number in 
//     *  hex_face.  
//     *  number of zeros = abs(iface + ones_in_transition_type) 
//     */
//     if(num_subelements >= phex_w_sub_neighbor_at_face->subelement_id){
//       hex_face = iface;
//       break;
//     }
//   }

//   /* Check whether hex_face is splitted or not.
//   */
//   if(phex_w_sub_neighbor_at_face->transition_type & (1 << hex_face) == 1){
//   /* Check where exactly the subelement is (front, bottom etc. on that face (hex_face))
//   *  For hex_face \in {0,1} we have to decide whether its in the front or back, or up or down.
//   *  For hex_face \in {2,3} we have to decide whether its in the left or right, or up or down.
//   *  For hex_face \in {4,5} we have to decide whether its in the left or right, or front or back.
//   */
//   int subelement_id_is_even = phex_w_sub_neighbor_at_face->subelement_id % 2;
//   /*
//   First check left or right, then front or back and then bottom or up */
//   if( hex_face == 2){
//     //Check if hex_face = 2 is splited.
//     if(phex_w_sub_neighbor_at_face->transition_type & (1 << hex_face) == 1){
//       if( ones_in_transition_type % 2 == 0 ){
//         if(subelement_id_is_even = 0){ //even
//           right_side = 1;
//         }
//         else{
//           left_side = 1;
//         }
//       else{
//         if(subelement_id_is_even = 0){
//           left_side = 1;
//         }
//         else{
//           right_side = 1;
//         }}
//         //Check whether up or down
//         if( num_subelements - phex_w_sub_neighbor_at_face->subelement_id > 2){
//           bottom_side = 1;
//         }
//         else{ 
//           up_side = 1;
//         }
//       }
//     }
//     else {
//       //hier seiten nachbarn eintragen. Bsp. von seite 1 
//     } 
//   }


  


// //ich bin auf der linken Seite, linke vordere untere kleine Pyramide , d.h. subelement id = 0 , face 2 liegt nach unten.
// // Grenzt also an die Pyramide an die vom hex auf Seite f4 liegt.  
//   if(face == 2) {
//     if (phex_w_sub_neighbor_at_face->subelement_id == 0) {
//       phex_w_sub_neighbor_at_face->subelement_id += num_siblings - 1;
//     }
//   }

//   if (face == 0) {
//     /* adjust subelement id counter clockwise */
//     if (phex_w_sub_neighbor_at_face->subelement_id == 0) {
//       phex_w_sub_neighbor_at_face->subelement_id += num_siblings - 1;
//     }
//     else {
//       phex_w_sub_neighbor_at_face->subelement_id -= 1;
//     }
//   }
//   else {
//     /* adjust subelement id clockwise */
//     if (phex_w_sub_neighbor_at_face->subelement_id == num_siblings - 1) {
//       phex_w_sub_neighbor_at_face->subelement_id = 0;
//     }
//     else {
//       phex_w_sub_neighbor_at_face->subelement_id += 1;
//     }
//   }

//   /* return dual face with resprect to neighboring sibling subelement */
//   /* Compute the face number as seen from elem.
//    *  0 -> 2    2 -> 0
//    */
//   *neigh_face[0] = subelement_face_dual[face];

//   // SC_ABORT_NOT_REACHED();
// }
}
/* *INDENT-ON* */

void
t8_subelement_scheme_hex_c::t8_element_children (const t8_element_t *elem,
                                                  int length,
                                                  t8_element_t *c[]) const
{
  /* if elem is a subelement, then this function will construct the children of its parent p8est quadrant */
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements **phex_w_sub_children =
    (t8_hex_with_subelements **) c;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;

  int                 ichild;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (length == P8EST_CHILDREN);

#ifdef T8_ENABLE_DEBUG
  {
    int                 i;
    for (i = 0; i < P8EST_CHILDREN; i++) {
      T8_ASSERT (t8_element_is_valid (c[i]));
    }
  }
#endif

  /* set coordinates and levels of the children */
  p8est_quadrant_children (q, &phex_w_sub_children[0]->p8q,
                           &phex_w_sub_children[1]->p8q,
                           &phex_w_sub_children[2]->p8q,
                           &phex_w_sub_children[3]->p8q,
                           &phex_w_sub_children[4]->p8q,
                           &phex_w_sub_children[5]->p8q,
                           &phex_w_sub_children[6]->p8q,
                           &phex_w_sub_children[7]->p8q
);

  for (ichild = 0; ichild < P8EST_CHILDREN; ++ichild) {
    t8_element_reset_subelement_values (c[ichild]);
    t8_element_copy_surround (q, &phex_w_sub_children[ichild]->p8q);
  }

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_child_id (const t8_element_t *elem) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q = &phex_w_sub->p8q;

  T8_ASSERT (t8_element_is_valid (elem));

  return (t8_element_is_subelement (elem) ? phex_w_sub->subelement_id :
          p8est_quadrant_child_id (q));

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_ancestor_id (const t8_element_t *elem,
                                                     int level) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q = &phex_w_sub->p8q;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  return p8est_quadrant_ancestor_id (q, level);

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_is_family (t8_element_t **fam) const
{
  /* Note that this test is very rudimentary, especially when there subelements are in fam */
  t8_hex_with_subelements **phex_w_sub_family =
    (t8_hex_with_subelements **) fam;

#ifdef T8_ENABLE_DEBUG
  {
    int                 i;
    int                 num_siblings = t8_element_num_siblings (fam[0]);
    for (i = 0; i < num_siblings; i++) {
      T8_ASSERT (t8_element_is_valid (fam[i]));
    }
  }
#endif

  /* Subelements can not be refined into other elements of a higher level. 
   * So if the first element of fam is a subelement, we assume that the following num_siblings 
   * many elements are its siblings and therefore form a family. */
  if (phex_w_sub_family[0]->transition_type != 0) {
    return 1;
  }
  /* If the first element of fam is no subelement we check the following elements of fam */
  else {
    /* If any of the following elements is a subelement, then they can not form a family */
    if (phex_w_sub_family[1]->transition_type != 0 ||
        phex_w_sub_family[2]->transition_type != 0 ||
        phex_w_sub_family[3]->transition_type != 0 ||
        phex_w_sub_family[4]->transition_type != 0 ||
        phex_w_sub_family[5]->transition_type != 0 ||
        phex_w_sub_family[6]->transition_type != 0 ||
        phex_w_sub_family[7]->transition_type != 0) {
      return 0;
    }
    /* If all elements of fam are no subelements, then we can use the p8est check is_family */
    else {
      return p8est_quadrant_is_family (&phex_w_sub_family[0]->p8q,
                                       &phex_w_sub_family[1]->p8q,
                                       &phex_w_sub_family[2]->p8q,
                                       &phex_w_sub_family[3]->p8q,
                                       &phex_w_sub_family[4]->p8q,
                                       &phex_w_sub_family[5]->p8q,
                                       &phex_w_sub_family[6]->p8q,
                                       &phex_w_sub_family[7]->p8q);
    }
  }
}

void
t8_subelement_scheme_hex_c::t8_element_set_linear_id (t8_element_t *elem,
                                                       int level,
                                                       t8_linearidx_t id)
  const
{
  t8_hex_with_subelements *phex_w_sub = (t8_hex_with_subelements *) elem;
  p8est_quadrant_t   *q = &phex_w_sub->p8q;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);
  T8_ASSERT (0 <= id && id < ((t8_linearidx_t) 1) << P8EST_DIM * level);

  p8est_quadrant_set_morton (q, level, id);
  T8_HEX_SET_TDIM (q, 2);

  // SC_ABORT_NOT_REACHED();
}

t8_linearidx_t
  t8_subelement_scheme_hex_c::t8_element_get_linear_id (const t8_element_t
                                                         *elem,
                                                         int level) const
{
  t8_hex_with_subelements *phex_w_sub = (t8_hex_with_subelements *) elem;
  p8est_quadrant_t   *q = &phex_w_sub->p8q;

  /* Note that the id of a subelement equals the id of its parent quadrant.
   * Therefore, the binary search (for example used in the leaf_face_neighbor function) 
   * will find a random subelement of the transition cell which might not be the desired neighbor of a given element. */

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);

  return p8est_quadrant_linear_id (q, level);

  // SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_first_descendant (const t8_element_t
                                                          *elem,
                                                          t8_element_t *desc,
                                                          int level) const
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_desc =
    (t8_hex_with_subelements *) desc;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_desc->p8q;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);

  p8est_quadrant_first_descendant (q, r, level);
  T8_HEX_SET_TDIM (r, 3);

  /* We allow constructing a last descendant from a subelement. 
   * Keep in mind, that transforming a hex element to a subelement does not change the 
   * p8est quadrant. Therefore, we are constructing the last descendant of the parent 
   * hex element of the given subelement. Since the last descendant is not meant to be 
   * a subelement, we reset the corresponding subelement values. */
  t8_element_reset_subelement_values (desc);

  // SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_last_descendant (const t8_element_t
                                                         *elem,
                                                         t8_element_t *desc,
                                                         int level) const
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_desc =
    (t8_hex_with_subelements *) desc;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_desc->p8q;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (desc));
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);

  p8est_quadrant_last_descendant (q, r, level);
  T8_HEX_SET_TDIM (r, 2);

  /* We allow constructing a last descendant from a subelement. 
   * Keep in mind, that transforming a hex element to a subelement does not change the 
   * p8est quadrant. Therefore, we are constructing the last descendant of the parent 
   * hex element of the given subelement. Since the last descendant is not meant to be 
   * a subelement, we reset the corresponding subelement values. */
  t8_element_reset_subelement_values (desc);

  // SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_successor (const t8_element_t *elem1,
                                                   t8_element_t *elem2,
                                                   int level) const
{
  const t8_hex_with_subelements *phex_w_sub_elem1 =
    (const t8_hex_with_subelements *) elem1;
  t8_hex_with_subelements *phex_w_sub_elem2 =
    (t8_hex_with_subelements *) elem2;

  const p8est_quadrant_t *q = &phex_w_sub_elem1->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_elem2->p8q;

  t8_linearidx_t      id;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem1));

  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);

  id = p8est_quadrant_linear_id (q, level);
  T8_ASSERT (id + 1 < ((t8_linearidx_t) 1) << P8EST_DIM * level);
  t8_element_reset_subelement_values (elem2);
  p8est_quadrant_set_morton (r, level, id + 1);
  t8_element_copy_surround (q, r);

  // SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_nca (const t8_element_t *elem1,
                                             const t8_element_t *elem2,
                                             t8_element_t *nca) const
{
  const t8_hex_with_subelements *phex_w_sub_elem1 =
    (const t8_hex_with_subelements *) elem1;
  const t8_hex_with_subelements *phex_w_sub_elem2 =
    (const t8_hex_with_subelements *) elem2;
  t8_hex_with_subelements *phex_w_sub_nca =
    (t8_hex_with_subelements *) nca;

  const p8est_quadrant_t *q1 = &phex_w_sub_elem1->p8q;
  const p8est_quadrant_t *q2 = &phex_w_sub_elem2->p8q;
  p8est_quadrant_t   *r = &phex_w_sub_nca->p8q;

  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
#if 0
  /* TODO: This assertions throws an error since it expects a 3D hex.
   *       this does not make sense. investigate. */
  T8_ASSERT (t8_element_surround_matches (q1, q2));
#endif

  /* In case of subelements, we use the parent quadrant and construct nca of the parent quadrant */
  t8_element_reset_subelement_values (nca);
  p8est_nearest_common_ancestor (q1, q2, r);
  t8_element_copy_surround (q1, r);

  // SC_ABORT_NOT_REACHED();
}

//Nummerierung der Seiten(der Pyramiden) wie in Davids Masterarbeit
t8_element_shape_t
t8_subelement_scheme_hex_c::t8_element_face_shape (const t8_element_t *elem,
                                                    int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  if( face == 4){
    return T8_ECLASS_QUAD;
  }
  else{
    return T8_ECLASS_TRIANGLE;
  }
  
  
}

void
t8_subelement_scheme_hex_c::t8_element_children_at_face (const t8_element_t
                                                          *elem, int face,
                                                          t8_element_t
                                                          *children[],
                                                          int num_children,
                                                          int *child_indices)
  const
{
#ifdef T8_ENABLE_DEBUG
  {
    int                 i;
    for (i = 0; i < num_children; i++) {
      T8_ASSERT (t8_element_is_valid (children[i]));
    }
  }
#endif
  /* This function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (0 <= face && face < P8EST_FACES);
  T8_ASSERT (num_children == t8_element_num_face_children (elem, face));

  /*
   * Compute the child id of the first and second child at the face.
   *
   *            3
   *
   *      x - - x - - x           This picture shows a refined quadrant
   *      |     |     |           with child_ids and the label for the faces.
   *      | 2   | 3   |           For examle for face 2 (bottom face) we see
   * 0    x - - x - - x   1       first_child = 0 and second_child = 1.
   *      |     |     |
   *      | 0   | 1   |
   *      x - - x - - x
   *
   *            2
   */

  T8_ASSERT (num_children == 2);
  int                 first_child;
  int                 second_child;
  /* TODO: Think about a short and easy bitwise formula. */
  switch (face) {
  case 0:
    first_child = 0;
    second_child = 2;
    break;
  case 1:
    first_child = 1;
    second_child = 3;
    break;
  case 2:
    first_child = 0;
    second_child = 1;
    break;
  case 3:
    first_child = 2;
    second_child = 3;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  /* From the child ids we now construct the children at the faces. */
  /* We have to revert the order and compute second child first, since
   * the usage allows for elem == children[0].
   */
  this->t8_element_child (elem, second_child, children[1]);
  this->t8_element_child (elem, first_child, children[0]);
  if (child_indices != NULL) {
    child_indices[0] = first_child;
    child_indices[1] = second_child;
  }

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_face_child_face (const t8_element_t
                                                         *elem, int face,
                                                         int face_child) const
{
  // T8_ASSERT (t8_element_is_valid (elem));

  // if (t8_element_is_subelement (elem)) {
  //   T8_ASSERT (face == 1);
  //   return t8_element_face_parent_face (elem, face);
  // }
  // else {
  //   /* For quadrants the face enumeration of children is the same as for the parent. */
  //   return face;
  // }
  SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_face_parent_face (const t8_element_t
                                                          *elem,
                                                          int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (face >= -1 && face <= P8EST_FACES);

  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q = &phex_w_sub->p8q;

  int                 child_id;

  if (face == -1) {
    return -1;
  }

  /* For subelements we need to adjust the output of this function.
   * A subelements face is a subface of the parent quadrant (the transition cell) if and only if the face number is 1. */
  if (t8_element_is_subelement (elem)) {
    if (face == 1) {
      /* In this case the face is a subface of the parent. We use the location function in order
       * to determine which of the parents faces intersects the subelements face. */
      int                 location[3] = { };
      t8_element_get_location_of_subelement (elem, location);

      /* subelements in location are enumerated clockwise (not as quadrant faces) */
      return subelement_location_to_parent_face[location[0]];
    }
    else {
      return -1;
    }
  }

  if (q->level == 0) {
    return face;
  }
  /* Determine whether face is a subface of the parent.
   * This is the case if the child_id matches one of the faces corners */
  child_id = p8est_quadrant_child_id (q);
  if (child_id == p8est_face_corners[face][0]
      || child_id == p8est_face_corners[face][1]) {
    return face;
  }
  return -1;

  // SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_transform_face (const t8_element_t
                                                        *elem1,
                                                        t8_element_t *elem2,
                                                        int orientation,
                                                        int sign,
                                                        int is_smaller_face)
  const
{
  const t8_hex_with_subelements *phex_w_sub_elem1 =
    (const t8_hex_with_subelements *) elem1;
  t8_hex_with_subelements *phex_w_sub_elem2 =
    (t8_hex_with_subelements *) elem2;

  const p8est_quadrant_t *qin = &phex_w_sub_elem1->p8q;
  p8est_quadrant_t   *p = &phex_w_sub_elem2->p8q;

  const p8est_quadrant_t *q;
  p4est_qcoord_t      h = P8EST_QUADRANT_LEN (qin->level);
  p4est_qcoord_t      x = qin->x;       /* temp storage for x coordinate in case elem1 = elem 2 */

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem1));
  T8_ASSERT (t8_element_is_valid (elem1));
  T8_ASSERT (t8_element_is_valid (elem2));
  T8_ASSERT (0 <= orientation && orientation < P8EST_FACES);

  if (sign) {
    /* The tree faces have the same topological orientation, and
     * thus we have to perform a coordinate switch. */
    /* We use p as storage, since elem1 and elem2 are allowed to
     * point to the same hex */
    q = (const p8est_quadrant_t *) p;
    t8_element_copy_surround (qin, (p8est_quadrant_t *) q);
    ((p8est_quadrant_t *) q)->x = qin->y;
    ((p8est_quadrant_t *) q)->y = x;
    x = q->x;                   /* temp storage in case elem1 = elem 2 */
  }
  else {
    q = qin;
  }

  p->level = q->level;
  /*
   * The faces of the root quadrant are enumerated like this:
   *
   *   v_2      v_3
   *     x -->-- x
   *     |       |
   *     ^       ^
   *     |       |
   *     x -->-- x
   *   v_0      v_1
   *
   * Orientation is the corner number of the bigger face that coincides
   * with the corner v_0 of the smaller face.
   */
  /* If this face is not smaller, switch the orientation:
   *  sign = 0   sign = 1
   *  0 -> 0     0 -> 0
   *  1 -> 2     1 -> 1
   *  2 -> 1     2 -> 2
   *  3 -> 3     3 -> 3
   */
  if (!is_smaller_face && (orientation == 1 || orientation == 2) && !sign) {
    orientation = 3 - orientation;
  }

  switch (orientation) {
  case 0:                      /* Nothing to do */
    p->x = q->x;
    p->y = q->y;
    break;
  case 1:
    p->x = P8EST_ROOT_LEN - q->y - h;
    p->y = x;
    break;
  case 2:
    p->x = q->y;
    p->y = P8EST_ROOT_LEN - x - h;
    break;
  case 3:
    p->x = P8EST_ROOT_LEN - q->x - h;
    p->y = P8EST_ROOT_LEN - q->y - h;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  T8_HEX_SET_TDIM (p, 2);

  t8_element_reset_subelement_values (elem2);

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_extrude_face (const t8_element_t
                                                      *face,
                                                      const t8_eclass_scheme_c
                                                      *face_scheme,
                                                      t8_element_t *elem,
                                                      int root_face) const
{
  /* build (extrude) elem from a given face element */
  t8_hex_with_subelements *phex_w_sub = (t8_hex_with_subelements *) elem;
  p8est_quadrant_t   *q = &phex_w_sub->p8q;

  const t8_dline_t   *l = (const t8_dline_t *) face;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (T8_COMMON_IS_TYPE
             (face_scheme, const t8_default_scheme_line_c *));
  T8_ASSERT (face_scheme->eclass == T8_ECLASS_LINE);
  T8_ASSERT (face_scheme->t8_element_is_valid (elem));
  T8_ASSERT (0 <= root_face && root_face < P8EST_FACES);

  /*
   * The faces of the root quadrant are enumerated like this:
   *
   *        f_2
   *     x -->-- x
   *     |       |
   *     ^       ^
   * f_0 |       | f_1
   *     x -->-- x
   *        f_3
   *
   * The arrows >,^ denote the orientation of the faces.
   * We need to scale the coordinates since a root line may have a different
   * length than a root hex.
   */
  q->level = l->level;
  switch (root_face) {
  case 0:
    q->x = 0;
    q->y = ((int64_t) l->x * P8EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    break;
  case 1:
    q->x = P8EST_LAST_OFFSET (q->level);
    q->y = ((int64_t) l->x * P8EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    break;
  case 2:
    q->x = ((int64_t) l->x * P8EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    q->y = 0;
    break;
  case 3:
    q->x = ((int64_t) l->x * P8EST_ROOT_LEN) / T8_DLINE_ROOT_LEN;
    q->y = P8EST_LAST_OFFSET (q->level);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  t8_element_reset_subelement_values (elem);
  /* We return the face of q at which we extruded. This is the same number
   * as root_face. */
  return root_face;

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_tree_face (const t8_element_t *elem,
                                                   int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  /* If elem is a subelement, then this function should only be called together with 
   * face = 1 since other faces will never intersect a tree face. */
  if (t8_element_is_subelement (elem)) {
    T8_ASSERT (face == 1);

    return t8_element_face_parent_face (elem, face);
  }
  else {
    T8_ASSERT (0 <= face && face < P8EST_FACES);
    /* For quadrants the face and the tree face number are the same. */
    return face;
  }

  // SC_ABORT_NOT_REACHED();
}

/** Construct the first descendant of an element that touches a given face.   */
void
t8_subelement_scheme_hex_c::t8_element_first_descendant_face (const
                                                               t8_element_t
                                                               *elem,
                                                               int face,
                                                               t8_element_t
                                                               *first_desc,
                                                               int level)
  const
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_first_desc =
    (t8_hex_with_subelements *) first_desc;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *desc = &phex_w_sub_first_desc->p8q;

  int                 first_face_corner;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (0 <= face && face < P8EST_FACES);
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);

  /* Get the first corner of q that belongs to face */
  first_face_corner = p8est_face_corners[face][0];
  /* Construce the descendant in that corner */
  p8est_quadrant_corner_descendant (q, desc, first_face_corner, level);
  t8_element_reset_subelement_values (first_desc);

  // SC_ABORT_NOT_REACHED();
}

/** Construct the last descendant of an element that touches a given face.   */
void
t8_subelement_scheme_hex_c::t8_element_last_descendant_face (const
                                                              t8_element_t
                                                              *elem, int face,
                                                              t8_element_t
                                                              *last_desc,
                                                              int level) const
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_last_desc =
    (t8_hex_with_subelements *) last_desc;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *desc = &phex_w_sub_last_desc->p8q;

  int                 last_face_corner;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (!t8_element_is_subelement (last_desc));
  T8_ASSERT (0 <= face && face < P8EST_FACES);
  T8_ASSERT (0 <= level && level <= P8EST_QMAXLEVEL);

  /* Get the last corner of q that belongs to face */
  last_face_corner = p8est_face_corners[face][1];
  /* Construce the descendant in that corner */
  p8est_quadrant_corner_descendant (q, desc, last_face_corner, level);
  t8_element_reset_subelement_values (last_desc);

  // SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_boundary_face (const t8_element_t
                                                       *elem, int face,
                                                       t8_element_t *boundary,
                                                       const
                                                       t8_eclass_scheme_c
                                                       *boundary_scheme) const
{
  // const t8_hex_with_subelements *phex_w_sub =
  //   (const t8_hex_with_subelements *) elem;
  // const p8est_quadrant_t *q = &phex_w_sub->p8q;

  // t8_dline_t         *l = (t8_dline_t *) boundary;

  // T8_ASSERT (t8_element_is_valid (elem));
  // T8_ASSERT (T8_COMMON_IS_TYPE
  //            (boundary_scheme, const t8_default_scheme_line_c *));
  // T8_ASSERT (boundary_scheme->eclass == T8_ECLASS_LINE);
  // T8_ASSERT (boundary_scheme->t8_element_is_valid (boundary));

  // if (!t8_element_is_subelement (elem)) {
  //   T8_ASSERT (0 <= face && face < P8EST_FACES);
  //   /* The level of the boundary element is the same as the quadrant's level */
  //   l->level = q->level;
  //   /*
  //    * The faces of the quadrant are enumerated like this:
  //    *        f_2
  //    *     x ---- x
  //    *     |      |
  //    * f_0 |      | f_1
  //    *     x ---- x
  //    *        f_3
  //    *
  //    * If face = 0 or face = 1 then l->x = q->y
  //    * if face = 2 or face = 3 then l->x = q->x
  //    */
  //   l->x = ((face >> 1 ? q->x : q->y) *
  //           ((int64_t) T8_DLINE_ROOT_LEN) / P8EST_ROOT_LEN);
  // }
  // else {
  //   /* face number 1 is the only face of a subelement that points outward of the transition cell */
  //   T8_ASSERT (face == 1);
  //   /* boundary faces of subelements:
  //    *
  //    *         x - - - - - x
  //    *         | \       / |
  //    *         |   \   /   |
  //    *         x - - x  e2 | f1
  //    *         |e1 /   \   |
  //    *      f1 | /       \ |
  //    *         x - - x - - x
  //    *               
  //    * for a split subelement (e1), the boundary face has a higher level
  //    * for a non split element (e2), the boundary face has the same level. 
  //    */

  //   int                 location[3] = { };      /* location = {location of subelement (face number of transition cell), split, first or second element if split} */
  //   t8_element_get_location_of_subelement (elem, location);
  //   int                 split = location[1];
  //   int                 second = location[2];

  //   if (split) {                /* if the subelement lies at a split face */
  //     l->level = q->level + 1;
  //     int                 len =
  //       P8EST_QUADRANT_LEN (phex_w_sub->p8q.level + 1);
  //     if (second) {             /* second subelement */
  //       if (location[0] == 0) { /* left face */
  //         l->x = q->y + len;
  //       }
  //       else if (location[0] == 1) {    /* upper face */
  //         l->x = q->x + len;
  //       }
  //       else if (location[0] == 2) {    /* right face */
  //         l->x = q->y;
  //       }
  //       else {                  /* lower face */
  //         l->x = q->x;
  //       }
  //     }
  //     else {                    /* first subelement */
  //       if (location[0] == 0) { /* left face */
  //         l->x = q->y;
  //       }
  //       else if (location[0] == 1) {    /* upper face */
  //         l->x = q->x;
  //       }
  //       else if (location[0] == 2) {    /* right face */
  //         l->x = q->y + len;
  //       }
  //       else {                  /* lower face */
  //         l->x = q->x + len;
  //       }
  //     }
  //   }
  //   else {                      /* if the subelement is not split */
  //     l->level = q->level;
  //     if (location[0] == 0) {   /* left face */
  //       l->x = q->y;
  //     }
  //     else if (location[0] == 1) {      /* upper face */
  //       l->x = q->x;
  //     }
  //     else if (location[0] == 2) {      /* right face */
  //       l->x = q->y;
  //     }
  //     else {                    /* lower face */
  //       l->x = q->x;
  //     }
  //   }
  // }

  SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_boundary (const t8_element_t *elem,
                                                  int min_dim, int length,
                                                  t8_element_t **boundary)
  const
{
  SC_ABORT ("Not implemented\n");
#if 0
#ifdef T8_ENABLE_DEBUG
  int                 per_eclass[T8_ECLASS_COUNT];
#endif
  int                 iface;

  T8_ASSERT (length ==
             t8_eclass_count_boundary (T8_ECLASS_hex, min_dim, per_eclass));

  T8_ASSERT (length == P8EST_FACES);
  for (iface = 0; iface < P8EST_FACES; iface++) {
    t8_element_boundary_face (elem, iface, boundary[iface]);
  }
#endif
}

int
t8_subelement_scheme_hex_c::t8_element_is_root_boundary (const t8_element_t
                                                          *elem,
                                                          int face) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q = &phex_w_sub->p8q;

  p4est_qcoord_t      coord;

  /* In case of a subelement, we need to change its face number to the face number of the parent hex */
  if (t8_element_is_subelement (elem)) {
    if (face == 1) {
      /* adjust face of subelement to face of parent */
      face = t8_element_face_parent_face (elem, face);
    }
    else {                      /* in case of a subelement and face 0 or 2 the face is no subface of the root boundary */
      return 0;
    }
  }

  T8_ASSERT (0 <= face && face < P8EST_FACES);

  /* if face is 0 or 1 q->x
   *            2 or 3 q->y
   */
  coord = face >> 1 ? q->y : q->x;
  /* If face is 0 or 2 check against 0.
   * If face is 1 or 3  check against LAST_OFFSET */
  return coord == (face & 1 ? P8EST_LAST_OFFSET (q->level) : 0);
  
  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_face_neighbor_inside (const
                                                              t8_element_t
                                                              *elem,
                                                              t8_element_t
                                                              *neigh,
                                                              int face,
                                                              int *neigh_face)
  const
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (neigh));
  T8_ASSERT (0 <= face && face < P8EST_FACES);

  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements *phex_w_sub_neigh =
    (t8_hex_with_subelements *) neigh;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;
  p8est_quadrant_t   *n = &phex_w_sub_neigh->p8q;

  /* In case of a subelement one should construct the face neighbor of the face-corresponding child quadrant
   * of the subelements parent quadrant. Therefore we might want to adjust the level  and adapt the
   * anchor node. */
  if (t8_element_is_subelement (elem)) {        /* if elem is a subelement */

    T8_ASSERT (0 <= face && face < T8_HEX_SUBELEMENT_FACES);

    if (face == 0) {            /* in this case the face neighbor of the subelement is a sibling */
      /* level and anchor stay the same */
      n->x = q->x;
      n->y = q->y;
      n->level = q->level;
    }
    if (face == 2) {            /* in this case the face neighbor of the subelement is a sibling */
      /* level and anchor stay the same */
      n->x = q->x;
      n->y = q->y;
      n->level = q->level;
    }
    if (face == 1) {            /* in this case the face neighbor is no sibling */
      int                 location[3] = { };
      t8_element_get_location_of_subelement (elem, location);

      /* setting the anchor node of the neighbor element */
      n->x = q->x;
      n->y = q->y;

      /* half the side length of the transition cell of the subelement */
      const p4est_qcoord_t shift = P8EST_QUADRANT_LEN (q->level + 1);

      int                 split = location[1];
      int                 second = location[2];

      /* we need to take into account whether the subelement is split or not */
      if (split) {              /* split */

        /* increase the level by one */
        n->level = q->level + 1;

        /* adjust the anchor node of the neighbor of the subelement depending on its location */
        if (location[0] == 0) { /* left face */
          if (!second) {
            n->x = q->x - shift;
          }
          else {
            n->x = q->x - shift;
            n->y = q->y + shift;
          }
        }
        else if (location[0] == 2) {    /* right face */
          if (!second) {
            n->x = q->x + 2 * shift;
            n->y = q->y + shift;
          }
          else {
            n->x = q->x + 2 * shift;
          }
        }
        else if (location[0] == 3) {    /* lower face */
          if (!second) {
            n->x = q->x + shift;
            n->y = q->y - shift;
          }
          else {
            n->y = q->y - shift;
          }
        }
        else {                  /* upper face */
          if (!second) {
            n->y = q->y + 2 * shift;
          }
          else {
            n->x = q->x + shift;
            n->y = q->y + 2 * shift;
          }
        }
      }

      else {                    /* not split */
        /* level stays the same */
        n->level = q->level;

        /* adjust the anchor node of the neighbor of the subelement depending on its location */
        if (location[0] == 0) { /* left face */
          n->x = q->x - 2 * shift;
        }
        else if (location[0] == 2) {    /* right face */
          n->x = q->x + 2 * shift;
        }
        else if (location[0] == 3) {    /* lower face */
          n->y = q->y - 2 * shift;
        }
        else {                  /* upper face */
          n->y = q->y + 2 * shift;
        }
      }
    }
  }
  else {                        /* if elem is no subelement */
    /* Directly construct the face neighbor */
    p8est_quadrant_face_neighbor (q, face, n);
  }

  t8_element_reset_subelement_values (neigh);

  T8_HEX_SET_TDIM (n, 2);

  /* In the following we set the dual faces of our element at the given face. */
  if (t8_element_is_subelement (elem)) {
    if (face == 1) {
      /* return dual face with respect to neighboring hex element */
      int                 location[3] = { };
      t8_element_get_location_of_subelement (elem, location);
      /* if the face is pointing outwards, then we set the face equal to the transition cell face and determine its dual face.
       * Compute the face number as seen from q.
       *  0 -> 1    1 -> 2    2 -> 0    3 -> 3
       */
      *neigh_face = subelement_location_to_parent_dual_face[location[0]];
    }
    else {
      T8_ASSERT (face == 0 || face == 2);
      /* return dual face with resprect to neighboring sibling subelement (note that the constructed neigh is NOT a subelement but the parent hex) */
      /* Compute the face number as seen from q.
       *  0 -> 2    2 -> 0
       */

      /* TODO (JM): Changed to 2D array accessor simply set to zero. Please change! */
      *neigh_face = subelement_face_dual[face][0];
    }
  }
  else {
    /* Compute the face number as seen from q.
     *  0 -> 1    1 -> 0    2 -> 3    3 -> 2
     */
    T8_ASSERT (neigh_face != NULL);
    *neigh_face = p8est_face_dual[face];
  }

  /* return true if neigh is inside the root */
  return p8est_quadrant_is_inside_root (n);

  // SC_ABORT_NOT_REACHED();
}

void
t8_subelement_scheme_hex_c::t8_element_anchor (const t8_element_t *elem,
                                                int coord[3]) const
{
  t8_hex_with_subelements *phex_w_sub = (t8_hex_with_subelements *) elem;
  p8est_quadrant_t   *q = &phex_w_sub->p8q;

  /* this function is not implemented for subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));

  coord[0] = q->x;
  coord[1] = q->y;
  coord[2] = q->z;
  T8_HEX_SET_TDIM (q, 2);

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_root_len (const t8_element_t *elem) const
{
  return P8EST_ROOT_LEN;
}

int
t8_subelement_scheme_hex_c::t8_element_refines_irregular () const
{
  /* In general, subelements do not refine regularly */
  return 1;
}

void
t8_subelement_scheme_hex_c::t8_element_reference_coords (const t8_element_t
                                                          *elem,
                                                          const double
                                                          *ref_coords,
                                                          const void
                                                          *user_data,
                                                          double *out_coords)
  const
{
  SC_ABORT ("This function is not implemented for the given scheme.\n");
}

void
t8_subelement_scheme_hex_c::t8_element_vertex_reference_coords (const
                                                                 t8_element_t
                                                                 *t,
                                                                 int vertex,
                                                                 double
                                                                 coords[])
  const
{
  int                 coords_int[3] = { };
  t8_element_vertex_coords (t, vertex, coords_int);

  /* We divide the integer coordinates by the root length of the hex
   * to obtain the reference coordinates. */
  coords[0] = (double) coords_int[0] / (double) P8EST_ROOT_LEN;
  coords[1] = (double) coords_int[1] / (double) P8EST_ROOT_LEN;
  coords[2] = (double) coords_int[2] / (double) P8EST_ROOT_LEN;
}

void
t8_subelement_scheme_hex_c::t8_element_vertex_coords (const t8_element_t
                                                       *elem, int vertex,
                                                       int coords[]) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q1 = &phex_w_sub->p8q;

  T8_ASSERT (t8_element_is_valid (elem));

  if (!t8_element_is_subelement (elem)) {
    int                 len;

    T8_ASSERT (0 <= vertex && vertex < 8);
    /* Get the length of the quadrant */
    len = P8EST_QUADRANT_LEN (q1->level);

    /* Compute the x, y and z coordinates of the vertex depending on the
     * vertex number */
    coords[0] = q1->x + (vertex & 1 ? 1 : 0) * len;
    coords[1] = q1->y + (vertex & 2 ? 1 : 0) * len;
    coords[2] = q1->z + (vertex & 4 ? 1 : 0) * len;
  }
  else {
    t8_element_vertex_coords_of_subelement (elem, vertex, coords);
  }
}

void
t8_subelement_scheme_hex_c::t8_element_vertex_coords_of_subelement (const
                                                                     t8_element_t
                                                                     *elem,
                                                                     int
                                                                     vertex,
                                                                     int
                                                                     coords[])
  const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q1 = &phex_w_sub->p8q;
  
  int                 len;

  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_subelement (elem));
  T8_ASSERT (vertex >= 0 && vertex < T8_HEX_SUBELEMENT_FACES);      /* all subelements are pyramids so T8_HEX_SUBELEMENT_FACES = 5 */

  /* get the length of the current quadrant */
  len = P8EST_QUADRANT_LEN (q1->level);

  /* Compute the x,y and z coordinates of subelement vertices, depending on the transition type, id and vertex number 
   * (faces enumerated clockwise, starting at the center of the transition cell): 
   *
   *               f1                      V1
   *         x - - - - - x                 x
   *         | \   2   / |               / |
   *         | 1 \   / 3 |             / 3 |
   *      f0 x - - + - - x f2  -->   + - - x 
   *         | 0 / | \ 4 |           V0    V2
   *         | / 6 | 5 \ | 
   *         x - - x - - x
   *               f3
   * 
   * In this example, the below location array would contain the values [2, 1, 1] 
   * (second face, split, first subelement at this face) */

  /* get location information of the given subelement */
  int                 location[3] = { };
  t8_element_get_location_of_subelement (elem, location);

  /* the face number, the subelement is adjacent to */
  int                 face_number = location[0];
  /* = 1, if the adjacent face is split and = 0, if not */
  int                 split = location[1];
  /* subelement_id type. First bit (left) = 1 if right, = 0 if left, second bit( middle) = 1 if back , = 0 if front, third bit (right) = 0 if up and is = 1 if down 

   * second bit front = 0, back = 1, third bit 0 = bottom 1 = up. For example: 110 stands for right and back (so only hex face f_4 and f_5 are possible.)
   */
  int                 sub_face_id = location[2];

  // /* Check, whether the get_location function provides meaningful location data */
  // T8_ASSERT (face_number == 0 || face_number == 1 || face_number == 2
  //            || face_number == 3);
  // T8_ASSERT ((split == 0 && sub_face_id == 0)
  //            || (split == 1 && (sub_face_id == 0 || sub_face_id == 1)));

  coords[0] = q1->x;
  coords[1] = q1->y;
  coords[2] = q1->z;
//TEST --> Probiere anhand der subelement ID + transition type die Koordinaten zu berechnen 
             /* vertex 0 (the first vertex allways equals the center of the element) */

    switch(vertex){
    case 4: //vertex 4 always equals the center of the hexahedron
      coords[0] += len / 2;
      coords[1] += len / 2;
      coords[2] += len / 2; 
      break;
    case 0:
      if(split == 0){
        if(face_number == 0){
          coords[0] += len;
        }
        if(face_number == 3){
          coords[1] += len;
        }
      }
/* ----------- face 0 + 1 (splitted) --------------------*/
      else{
        if(face_number == 0 || face_number == 1){
          if(sub_face_id & 1 == 1){ // up
            coords[2] += len / 2;
            }
          if(sub_face_id & 2 != 0 ){ // back 
            coords[1] += len / 2;
          }
        }
      }
/* ----------- face 2 + 3 (splitted) --------------------*/
        if(face_number == 2 || face_number == 3){
          if(sub_face_id & 1 == 1){ // up
            coords[2] += len / 2;
            }
          if(sub_face_id & 4 == 0 ){ // right 
            coords[0] += len / 2;
          }
          if(face_number == 3){
            coords[1] += len;
          }
        }
/* ----------- face 4 + 5 (splitted) --------------------*/
        if(face_number == 4 || face_number == 5){
          if(sub_face_id & 2 != 0){ // back
            coords[1] += len / 2;
            }
          if(sub_face_id & 4 == 0 ){ // right 
            coords[0] += len / 2;
          }
          if(face_number == 5){
            coords[2] += len;
          }
        }
      break;

    case 1:
    if(split == 0){
        if(face_number == 0 || face_number == 1){
          coords[1] += len;
        }
        if(face_number > 1){
          coords[0] += len;
        }
        if(face_number == 5){
          coords[2] += len;
        }
      }
/* ----------- face 0 + 1 (splitted) --------------------*/
      else{
        if(face_number == 0 || face_number == 1){
          if(sub_face_id & 1 != 0){ // up
            coords[2] += len / 2;
            }
          if(sub_face_id & 2 == 0 ){ // front 
            coords[1] += len / 2;
          }
          else { //back
            coords[1] += len;
          }
          if( face_number == 1){
            coords[0] += len;
          }
        }
/* ----------- face 2 + 3 (splitted) --------------------*/
        if(face_number == 2 || face_number == 3){
          if(sub_face_id & 1 != 0){ // up
            coords[2] += len / 2;
            }
          if(sub_face_id & 4 == 0 ){ // right 
            coords[0] += len;
          }
          else{ //left
            coords[0] += len / 2;
          }
          if(face_number == 3){
            coords[1] += len;
          }        
          }

/* ----------- face 4 + 5 (splitted) --------------------*/
        if(face_number == 4 || face_number == 5){
          if(sub_face_id & 2 != 0){ // back
            coords[1] += len / 2;
            }
          if(sub_face_id & 4 == 0 ){ // right 
            coords[0] += len;
          }
          else{ //left
            coords[0] += len / 2;
          }
          if(face_number == 5){
            coords[2] += len;
          }
        }
     break;
    case 2:
    if(split == 0){
        if(face_number != 4){
          coords[2] += len;
        }
        if(face_number == 1){
          coords[0] += len;
        }
        if(face_number > 2){
          coords[1] += len;
        }
      }
/* ----------- face 0 + 1 (splitted) --------------------*/
      else{
        if(face_number == 0 || face_number == 1){
          if(sub_face_id & 2 != 0){ // back
            coords[1] += len / 2;
            }
          if(sub_face_id & 1 != 0 ){ // up
            coords[2] += len;
          }
          else { //bottom
            coords[2] += len / 2;
          }
          if(face_number == 1){
            coords[0] += len;
          }
        }

/* ----------- face 2 + 3 (splitted) --------------------*/
        if(face_number == 2 || face_number == 3){
          if(sub_face_id & 1 != 0){ // up
            coords[2] += len;
            }
          else{
            coords[2] += len / 2;
          }
          if(sub_face_id & 4 == 0 ){ // right 
            coords[0] += len / 2;
          }
          if(face_number == 3){
            coords[1] += len;
          }
        }  
/* ----------- face 4 + 5 (splitted) --------------------*/

        if(face_number == 4 || face_number == 5){
          if(sub_face_id & 2 != 0){ // back
            coords[1] += len;
            }
          else{ //front
            coords[1] += len / 2;
          }
          if(sub_face_id & 4 == 0 ){ // right 
            coords[0] += len / 2;
          }
          if(face_number == 5){
            coords[2] += len;
          }
        }
      }
     break;
    case 3:
    if(split == 0){
        if(face_number != 3){
          coords[1] += len;
        }
        if(face_number != 4){
          coords[2] += len;
        }
        if(face_number > 1){
          coords[0] += len;
        }
      }
      /* ----------- face 0 + 1 (splitted) --------------------*/
      else{
        if(face_number == 0 || face_number == 1){
          if(sub_face_id & 2 != 0){ // back
            coords[1] += len;
            }
          else { //front
            coords[1] += len / 2;
          }
          if(sub_face_id & 1 != 0 ){ // up
            coords[2] += len;
          }
          else{ // bottom
            coords[2] += len / 2;
          }
          if( face_number == 1){
            coords[0] += len;
          }
        }
/* ----------- face 2 + 3 (splitted) --------------------*/
        if(face_number == 2 || face_number == 3){
          if(sub_face_id & 1 != 0){ // up
            coords[2] += len;
            }
          else{
            coords[2] += len / 2;
          }
          if(sub_face_id & 4 == 0 ){ // right 
            coords[0] += len;
          }
          else{ // left 
            coords[0] += len / 2;
          }
          if(face_number == 3){
            coords[1] += len;
          }
        }

/* ----------- face 4 + 5 (splitted) --------------------*/
        if(face_number == 4 || face_number == 5){
          if(sub_face_id & 2 != 0){ // back
            coords[1] += len;
            }
          else{ //front
            coords[1] += len / 2;
          }
          if(sub_face_id & 4 == 0 ){ // right 
            coords[0] += len;
          }
          else{ // left
            coords[0] += len / 2;
          }
          if( face_number == 5){
            coords[2] += len;
          }
        }
      }
     break;
      }    
    }   

    
  /* using the location data to determine vertex coordinates */
  // if (vertex == 0) {            /* vertex 0 (the first vertex allways equals the center of the element) */
  //   coords[0] += len / 2;
  //   coords[1] += len / 2;
  //   coords[2] += len / 2;
  // }                             /* end of vertex == 0 */
  // else if (vertex == 1) {       /* vertex 1 */
  //   if (face_number == 0) {
  //     if (split && sub_face_id) {
  //       coords[1] += len / 2;
  //     }
  //   }
  //   else if (face_number == 1) {
  //     coords[1] += len;
  //     if (split && sub_face_id) {
  //       coords[0] += len / 2;
  //     }
  //   }
  //   else if (face_number == 2) {
  //     coords[0] += len;
  //     coords[1] += len;
  //     if (split && sub_face_id) {
  //       coords[1] -= len / 2;
  //     }
  //   }
  //   else {
  //     coords[0] += len;
  //     if (split && sub_face_id) {
  //       coords[0] -= len / 2;
  //     }
  //   }
  // }                             /* end of vertex == 1 */
  // else {                        /* vertex 2 */
  //   if (face_number == 0) {
  //     coords[1] += len;
  //     if (split && !sub_face_id) {
  //       coords[1] -= len / 2;
  //     }
  //   }
  //   else if (face_number == 1) {
  //     coords[0] += len;
  //     coords[1] += len;
  //     if (split && !sub_face_id) {
  //       coords[0] -= len / 2;
  //     }
  //   }
  //   else if (face_number == 2) {
  //     coords[0] += len;
  //     if (split && !sub_face_id) {
  //       coords[1] += len / 2;
  //     }
  //   }
  //   else {
  //     if (split && !sub_face_id) {
  //       coords[0] += len / 2;
  //     }
  //   }
  // }                             /* end of vertex == 2 */
}

void
t8_subelement_scheme_hex_c::t8_element_to_transition_cell (const t8_element_t
                                                            *elem, int transition_type,
                                                            t8_element_t *c[])
{
  const t8_hex_with_subelements *phex_w_sub_elem =
    (const t8_hex_with_subelements *) elem;
  t8_hex_with_subelements **phex_w_sub_subelement =
    (t8_hex_with_subelements **) c;

  const p8est_quadrant_t *q = &phex_w_sub_elem->p8q;

  /* this function should not be callable by subelements */
  T8_ASSERT (!t8_element_is_subelement (elem));
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (transition_type >= 0 && transition_type <= T8_SUB_HEX_MAX_TRANSITION_TYPE);

  int                 num_subelements =
    t8_element_get_number_of_subelements (transition_type);

#ifdef T8_ENABLE_DEBUG
  {
    int                 j;
    for (j = 0; j < num_subelements; j++) {
      T8_ASSERT (t8_element_is_valid (c[j]));
    }
  }
#endif

  /* get the length of a children-quadrant */
  const int8_t        level = (int8_t) (q->level);

  T8_ASSERT (p8est_quadrant_is_extended (q));
  T8_ASSERT (q->level < P8EST_QMAXLEVEL);

  /* Setting the parameter values for different subelements. 
   * The different subelement types (up to rotation) are:
   *                               
   *      x - - - - - - x         x - - - - - x        x - - - - - x        x - - - - - x        x - - x - - x        x - - x - - x
   *      |             |         | \   2   / |        | \       / |        | \       / |        | \   |   / |        | \   |   / |
   *      |             |         | 1 \   /   |        |   \   /   |        |   \   /   |        |   \ | /   |        |   \ | /   |
   *      |             |   -->   x - - X   3 |   or   x - - x     |   or   x - - x - - x   or   x - - x - - x   or   x - - x - - x
   *      |             |         | 0 /   \   |        |   / | \   |        |   /   \   |        |   /   \   |        |   / | \   |
   *      | elem        |         | /   4   \ |        | /   |   \ |        | /       \ |        | /       \ |        | /   |   \ |
   *      + - - - - - - x         x - - - - - x        x - - x - - x        x - - - - - x        x - - - - - x        x - - x - - x
   *           
   * Sub_ids are counted clockwise, starting with the (lower) left subelement with id 0.                    
   * Note, that we do not change the p8est quadrant. */

  int                 sub_id_counter = 0;
  for (sub_id_counter = 0; sub_id_counter < num_subelements; sub_id_counter++) {
    phex_w_sub_subelement[sub_id_counter]->p8q.x = q->x;
    phex_w_sub_subelement[sub_id_counter]->p8q.y = q->y;
    phex_w_sub_subelement[sub_id_counter]->p8q.z = q->z;
    phex_w_sub_subelement[sub_id_counter]->p8q.level = level;

    phex_w_sub_subelement[sub_id_counter]->transition_type = transition_type;

    //Überlegung  transition type umwandeln in subelement id--> nicht einfach counter
    phex_w_sub_subelement[sub_id_counter]->subelement_id = sub_id_counter;

    T8_ASSERT (t8_element_is_valid (c[sub_id_counter]));
    t8_element_copy_surround (q,
                              &phex_w_sub_subelement[sub_id_counter]->p8q);
  }
}

int
t8_subelement_scheme_hex_c::t8_element_get_number_of_subelements (int
                                                                   transition_type)
  const
{
  /* we could return 0 for transition type 0 but we will assert this case for safety reasons */
  T8_ASSERT (transition_type != 0);

  /* consider transition_type 16 = 010000 in base two -> there are 6 + (1)*3 = 9 subelements */
  int                 num_hanging_faces = 0;
  int                 ichild;
  for (ichild = 0; ichild < P8EST_FACES; ichild++) {    /* Count the number of ones of the binary transition type. This number equals the number of hanging faces. */
    num_hanging_faces += (transition_type & (1 << ichild)) >> ichild;
  }

  /* The number of subelements equals the number of neighbours: */
  return P8EST_FACES + num_hanging_faces*3;
}

void
t8_subelement_scheme_hex_c::t8_element_get_location_of_subelement (const
                                                                    t8_element_t
                                                                    *elem,
                                                                    int
                                                                    location
                                                                    []) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  /* this function only works for subelements */
  T8_ASSERT (t8_element_is_subelement (elem));

  T8_ASSERT (t8_element_is_valid (elem));

  /* Consider the following transition cell of type 13:
   *            
   *              f0                         1
   *        x - - x - - x              x - - x - - x           
   *        |           |              | \ 2 | 3 / |           faces:                                                      f3   f2   f1   f0
   *        |           |              | 1 \ | / 4 |           binary code:                                                 1    1    0    1   (=13)
   *     f3 x           x f2   -->   1 x - - x - - x 1   -->   rearrange binaries s.t. the faces are enumerated clockwise:  1    1    1    0
   *        |           |              | 0 /   \ 5 |           number subelements at face:                                  2    2    1    2
   *        | elem      |              | /   6   \ |           consider sub_id 3:                                                x -> second subelement on the upper face
   *        + - - - - - x              x - - - - - x
   *              f1                         0
   *           
   * We will use the binary representation to determine the location of the given subelement. 
   * 
   * We need to know: 
   *     i)   the face number of the first vertex (values: {0,1,2,3}).
   *     ii)  whether this face is split in half (values: {0,1}).
   *     iii) if the subelement is the first or second subelement at the face (values: {0,1}).
   * 
   * These informations are then saved in the location array which will be used by the element_vertex function, 
   * to automatically determine the vertex coordinates of the given subelement. 
   * 
   * The location array for the above example would be {1,1,1} (upper face, split = true, second subelement at the upper face). */

  /* 1) convert the transition type from a decimal to a binary representation */
  int                 type = phex_w_sub->transition_type;
  int                 binary_array[P8EST_FACES] = { };

  int                 iface;

   /* We need an array with 6 elements to store all subelement types of the hex scheme from 1 to 63 ({0, 0, 0, 0, 0, 1} to {1, 1, 1, 1, 1, 1}) */
  for (iface = 0; iface < P8EST_FACES; iface++) {      
    binary_array[(P8EST_FACES - 1) - iface] = (type & (1 << iface)) >> iface;
  }                             /* we now got a binary represenation of the transition type, bitwise stored in an array */

  /* 2) rearrange the binary representation to be in clockwise order */
  int                 binary_array_temp[P8EST_FACES] = { };

  for (iface = 0; iface < P8EST_FACES; iface++) {       /* copying the binary array */
    binary_array_temp[iface] = binary_array[iface];
  }

  // for (iface = 0; iface < P8EST_FACES; iface++) {       /* bringing the entries of binary array into clockwise order */
  //   binary_array[iface] =
  //     binary_array_temp[subelement_location_to_parent_face[iface]];
  // }

  /* 3) use the rearranged binary representation, and the sub_id to determine the location of the subelement and store these information in an array */
  /*     3.1) location[0] -> the face_number, the subelement is adjacent to */
  /*     3.2) location[1] -> if the face is split or not */
  /*     3.3) location[2] -> if the subelement is the left/right, front/back or bottom/up Same idea as with the transition type: first bit 0 = left, 1 = right,
  *                          second bit front = 0, back = 1, third bit 0 = bottom 1 = up. For example: 110 stands for right and back (so only hex face f_4 and f_5 are possible.)  */
  T8_ASSERT (phex_w_sub->subelement_id <
             t8_element_get_number_of_subelements
             (phex_w_sub->transition_type));

  int                 sub_id = phex_w_sub->subelement_id;
  int                 sub_face_id_array[3];
  int                 sub_face_id = 0;
  int                 face_number = -1;
  int                 split;

  int                 cum_neigh_array[P8EST_FACES] = { };

  /* construct a cumulative array of the number of neighbors from face 0 to face 5 */
  cum_neigh_array[0] = binary_array[0]*3 + 1;
  cum_neigh_array[1] = cum_neigh_array[0] + binary_array[1]*3 + 1;
  cum_neigh_array[2] = cum_neigh_array[1] + binary_array[2]*3 + 1;
  cum_neigh_array[3] = cum_neigh_array[2] + binary_array[3]*3 + 1;
  cum_neigh_array[4] = cum_neigh_array[3] + binary_array[4]*3 + 1;
  cum_neigh_array[5] = cum_neigh_array[4] + binary_array[5]*3 + 1;

  /* 3.1) we can use the cumulative array to determine the face number of the given subelement */
  if (sub_id < cum_neigh_array[0]) {
    face_number = 0;
  }
  else {
    for (iface = 0; iface < P8EST_FACES - 1; ++iface) {
      if (sub_id >= cum_neigh_array[iface]
          && sub_id < cum_neigh_array[iface + 1]) {
        face_number = iface + 1;
        break;
      }
    }
  }

  /* make sure that a face_number has been found */
  T8_ASSERT (face_number >= 0);

  /* 3.2) determine, whether the face is split or not */
  if (binary_array[face_number] == 0) {
    split = 0;                  /* the face is not split */
  }
  else {
    split = 1;                  /* the face is split */
  }
  t8_debugf("Is splitted? 1 = yes, 0 = not: %i \n", split);
  if(split == 1){

    /* 3.3) determine, whether the subelement is the left/right, front/back or bottom/up subelement at the face */
  //First left/ right (only for face number 2, 3, 4, 5)
  if( face_number >1){
    if ((sub_id + 1 == cum_neigh_array[face_number] ) || (sub_id + 3 == cum_neigh_array[face_number])) {
      sub_face_id_array[0] = 1;            /* right*/
    }
    else {
      sub_face_id_array[0] = 0;            /* left */
    } 
  }
  //Second check front or back (only for face numbers 0, 1, 4, 5)
  if( face_number <= 1 ){
    if ((sub_id + 1 == cum_neigh_array[face_number]) || (sub_id + 3 == cum_neigh_array[face_number])) {
      sub_face_id_array[1] = 1;            /* back subelement */
    }
    else {
      sub_face_id_array[1] = 0;            /* front subelement */
    } 
  }
    if( face_number >= 4 ){
    if ((sub_id + 1 == cum_neigh_array[face_number]) || (sub_id + 2 == cum_neigh_array[face_number])) {
      sub_face_id_array[1] = 1;            /* back subelement */
    }
    else {
      sub_face_id_array[1] = 0;            /* front subelement */
    } 
  }
  //Third check up or down (only for face numbers 0, 1, 2, 3)
  if( face_number < 4 ){
    if ((sub_id + 2 == cum_neigh_array[face_number] ) || (sub_id + 1 == cum_neigh_array[face_number] )) {
      sub_face_id_array[2] = 1;            /* up subelement */
    }
    else {
      sub_face_id_array[2] = 0;            /* bottom subelement */
    } 
  }
  //Calculate the sub_face_id out of the sub_face_id_array
  for(int i = 0; i < 3; i++ ){
    if( sub_face_id_array[i] == 1){
      sub_face_id += std::pow(2, 2-i);
    }
  }

  }

//Important: sub_face_id 0 7 is not valid, because one subelement can not be right, bottom, and front at the same time. 
//faces that can have subelements on the left/ right: 2,3,4,5
//faces that can have subelements on the bottom/ up: 0,1,2,3
//faces that can have subelements at the front/ back: 0,1,4,5
// --> if face_number \in {0,1} sub_face_id \in {0,1,2,3} (case 1) 
// --> if face_number \in {2,3} sub_face_id \in {0,2,4,6} (case 2)
// --> if face_number \in {4,5} sub_face_id \in {0,1,4,5} (case 3)
//T8_ASSERT(sub_face_id == 7);
//case 1
t8_debugf(" Das ist die face number %i\n", face_number);
t8_debugf(" Das ist die subface_id %i\n", sub_face_id);
// T8_ASSERT(face_number < 2 && sub_face_id > 3 );
// //case 2
// T8_ASSERT((face_number > 1 && face_number < 4) && (sub_face_id == 1 || sub_face_id == 3));
// //case 3
// T8_ASSERT((face_number > 3) && (sub_face_id == 2 || sub_face_id == 3));


  location[0] = face_number;
  location[1] = split;
  location[2] = sub_face_id;
}

void
t8_subelement_scheme_hex_c::t8_element_reset_subelement_values (t8_element *
                                                                 elem) const
{
  t8_hex_with_subelements *phex_w_sub = (t8_hex_with_subelements *) elem;

  phex_w_sub->transition_type = 0;
  phex_w_sub->subelement_id = 0;
}

void
t8_subelement_scheme_hex_c::t8_element_copy_subelement_values (const
                                                                t8_element *
                                                                source,
                                                                t8_element *
                                                                dest) const
{
  const t8_hex_with_subelements *phex_w_sub_source =
    (const t8_hex_with_subelements *) source;
  t8_hex_with_subelements *phex_w_sub_dest =
    (t8_hex_with_subelements *) dest;

  phex_w_sub_dest->transition_type = phex_w_sub_source->transition_type;
  phex_w_sub_dest->transition_type = phex_w_sub_source->transition_type;
  phex_w_sub_dest->subelement_id = phex_w_sub_source->subelement_id;
}

int
t8_subelement_scheme_hex_c::t8_element_is_subelement (const
                                                       t8_element *
                                                       elem) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  T8_ASSERT (phex_w_sub->transition_type >= 0);

  /* transition_type == 0 => elem is no subelement.
   * transition_type != 0 => elem is subelement 
   */
  return (phex_w_sub->transition_type == 0 ? false : true);
}

int
t8_subelement_scheme_hex_c::t8_element_get_transition_type (const
                                                             t8_element *
                                                             elem)
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  return phex_w_sub->transition_type;
}

int
t8_subelement_scheme_hex_c::t8_element_get_subelement_id (const
                                                           t8_element * elem)
{
  T8_ASSERT (t8_element_is_valid (elem));

  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  return phex_w_sub->subelement_id;
}

t8_element_shape_t
t8_subelement_scheme_hex_c::t8_element_shape (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  return (t8_element_is_subelement (elem) ? T8_ECLASS_PYRAMID :
          T8_ECLASS_HEX);
}

int
t8_subelement_scheme_hex_c::t8_element_num_corners (const t8_element_t *elem) const
{
  T8_ASSERT (t8_element_is_valid (elem));

  return (t8_element_is_subelement (elem) ? T8_HEX_SUBELEMENT_FACES :
          P8EST_FACES);
}

int
t8_subelement_scheme_hex_c::t8_element_find_neighbor_in_transition_cell
  (const t8_element_t *elem, const t8_element_t *pseudo_neigh, int elem_face)
{
  /* In this function, we assume pseudo_neigh to be a random subelement of a transition cell that includes
   * the real neighbor of elem at face elem_face. This function will output the subelement_id of the real neighbor of elem. */
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_valid (pseudo_neigh));

  /* we expect neigh to be a element in a transition cell, thus to be a subelement */
  T8_ASSERT (t8_element_is_subelement (pseudo_neigh));

  const t8_hex_with_subelements *
    phex_w_sub_elem = (const t8_hex_with_subelements *) elem;
  const t8_hex_with_subelements *
    phex_w_sub_pseudo_neigh =
    (const t8_hex_with_subelements *) pseudo_neigh;

  /* In the following, all possible neighbor configurations are defined, such that subelement neighbors can be
   * identified in LFN_transitioned. */
  if (phex_w_sub_elem->transition_type != 0
      && (elem_face == 0 || elem_face == 2)) {
    /* In this case, we have the following situation:
                  x - - - - - - - -x
                 / \              /|
               /    \           /  |
     *       /       \        /    |
     *      x - - - - - - - x      |
     *      |          \   /|      |
     *      |           +   |      |
     *      |         / \   |      |
     *      x       /    \  x
     *      |     /       \ |     /
     *      |   /          \|   /
     *      | /             | /
     *      x - - - - - - - x
     *
     * Elem and face 0 or face 2 is given and a random sibling subelement neigh is given, too. 
     * We are searching for the subelement id of the real neighbor N_f0 or N_f2, depending on the face number. */
    int
      shift;
    if (elem_face == 0) {
      shift = -1;
    }
    if (elem_face == 2) {
      shift = 1;
    }
    int
      num_subelements =
      t8_element_get_number_of_subelements
      (phex_w_sub_elem->transition_type);
    return ((phex_w_sub_elem->subelement_id + shift) + num_subelements) % num_subelements;     /* the neighbor is directly before or after elem modulo the number of subelements in the transition cell */
  }
  /* Below are the cases in which the neighbor can not be identified as simple as above. 
   * The idea is to fill a location array with the desired properties of the real neighbor. 
   * Togehter with the type of the transition cell of pseudo_neigh, we can then identify the sub_id of the right neighbor. */

  if (phex_w_sub_elem->transition_type != 0 && elem_face == 1) {
    /* In this case, we have the following situation:
     * 
     *      x - - - - - - - x - - - - - - - x
     *      | \           / | \           / |
     *      |   \       /   |   \       /   |
     *      |     \   /     |     \   /     |
     *      x - - - x neigh | elem  x       |
     *      |     /   \     |     / | \     |
     *      |   /pseudo \   |   /   |   \   |
     *      | /   neigh   \ | /     |     \ |
     *      x - - - - - - - x - - - x - - - x
     *
     * A subelement elem is given as well as a random subelement pseudo_neigh from a neighboring transition cell. 
     * We are searching for the subelement id of the real neighbor neigh. 
     * Note that both transition cells can have different levels. */

    /* get the location of elem */
    int
    location_elem[3] = { };     /* {face, is_split, number of subelement at face} */
    t8_element_get_location_of_subelement (elem, location_elem);

    /* Initialize the location array of the real neighbor. */
    int
    location_neigh[3] = { -1, -1, -1 };

    /* the pseudo_neigh tranaition cell has a lower level than the elem transition cell */
    if (phex_w_sub_pseudo_neigh->p8q.level < phex_w_sub_elem->p8q.level) {
      if (location_elem[0] == 0) {      /* left face of transition cell */
        if (phex_w_sub_pseudo_neigh->p8q.y == phex_w_sub_elem->p8q.y) {
          location_neigh[0] = 2;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 1;        /* second subelement at face */
        }
        if (phex_w_sub_pseudo_neigh->p8q.y != phex_w_sub_elem->p8q.y) {
          location_neigh[0] = 2;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 0;        /* first subelement at face */
        }
      }
      if (location_elem[0] == 1) {      /* upper face of transition cell */
        if (phex_w_sub_pseudo_neigh->p8q.x == phex_w_sub_elem->p8q.x) {
          location_neigh[0] = 3;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 1;        /* first or second subelement at face */
        }
        if (phex_w_sub_pseudo_neigh->p8q.x != phex_w_sub_elem->p8q.x) {
          location_neigh[0] = 3;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 0;        /* first subelement at face */
        }
      }
      if (location_elem[0] == 2) {      /* right face of transition cell */
        if (phex_w_sub_pseudo_neigh->p8q.y == phex_w_sub_elem->p8q.y) {
          location_neigh[0] = 0;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 0;        /* first subelement at face */
        }
        if (phex_w_sub_pseudo_neigh->p8q.y != phex_w_sub_elem->p8q.y) {
          location_neigh[0] = 0;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 1;        /* first or second subelement at face */
        }
      }
      if (location_elem[0] == 3) {      /* lower face of transition cell */
        if (phex_w_sub_pseudo_neigh->p8q.x == phex_w_sub_elem->p8q.x) {
          location_neigh[0] = 1;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 0;        /* first subelement at face */
        }
        if (phex_w_sub_pseudo_neigh->p8q.x != phex_w_sub_elem->p8q.x) {
          location_neigh[0] = 1;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 1;        /* second subelement at face */
        }
      }
    }
    /* the pseudo_neigh tranaition cell has not a lower level than the elem transition cell */
    if (phex_w_sub_pseudo_neigh->p8q.level >= phex_w_sub_elem->p8q.level) {
      if (location_elem[0] == 0) {      /* left face of transition cell */
        location_neigh[0] = 2;  /* face */
        location_neigh[1] = 0;  /* not split */
        location_neigh[2] = 0;  /* first (only) subelement at face */
      }
      if (location_elem[0] == 1) {      /* upper face of transition cell */
        location_neigh[0] = 3;  /* face */
        location_neigh[1] = 0;  /* not split */
        location_neigh[2] = 0;  /* first (only) subelement at face */
      }
      if (location_elem[0] == 2) {      /* right face of transition cell */
        location_neigh[0] = 0;  /* face */
        location_neigh[1] = 0;  /* not split */
        location_neigh[2] = 0;  /* first (only) subelement at face */
      }
      if (location_elem[0] == 3) {      /* lower face of transition cell */
        location_neigh[0] = 1;  /* face */
        location_neigh[1] = 0;  /* not split */
        location_neigh[2] = 0;  /* first (only) subelement at face */
      }
    }

    /* check, that a neighbor is found and the location array is adjusted */
    T8_ASSERT (location_neigh[0] >= 0 && location_neigh[1] >= 0
               && location_neigh[2] >= 0);

    /* Depending on the location of elem, we have filled location_neigh with the data of the real neighbor.
     * This data will be used to determine the sub_id of the neighbor within the transition cell of pseudo_neigh. */
    return
      t8_element_get_id_from_location (t8_element_get_transition_type
                                       (pseudo_neigh), location_neigh);
  }
  if (!t8_element_is_subelement (elem)) {
    /* In this case, we have the following situation:
     * 
     *      x - - - - - - - x - - - - - - - x
     *      | \           / |               |
     *      |   \       /   |               |
     *      |     \   /     |               |
     *      x - - - x neigh |     elem      |
     *      |     /   \     |               |
     *      |   /pseudo \   |               |
     *      | /   neigh   \ |               |
     *      x - - - - - - - x - - - - - - - x
     *
     * Subelement elem is given as well as a random subelement neigh from a neighboring transition cell. 
     * We are searching for the subelement id of the real neighbor neigh.
     * Note that the transition cell of pseudo_neigh and elem can have different levels. */

    /* Initialize the location array of the real neighbor. */
    int
    location_neigh[3] = { -1, -1, -1 };

    /* the pseudo_neigh tranaition cell has a lower level than elem */
    if (phex_w_sub_pseudo_neigh->p8q.level < phex_w_sub_elem->p8q.level) {
      if (elem_face == 0) {     /* left face */
        if (phex_w_sub_pseudo_neigh->p8q.y == phex_w_sub_elem->p8q.y) {
          location_neigh[0] = 2;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 1;        /* second subelement at face */
        }
        if (phex_w_sub_pseudo_neigh->p8q.y != phex_w_sub_elem->p8q.y) {
          location_neigh[0] = 2;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 0;        /* first subelement at face */
        }
      }
      if (elem_face == 1) {     /* right face */
        if (phex_w_sub_pseudo_neigh->p8q.y == phex_w_sub_elem->p8q.y) {
          location_neigh[0] = 0;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 0;        /* first subelement at face */
        }
        if (phex_w_sub_pseudo_neigh->p8q.y != phex_w_sub_elem->p8q.y) {
          location_neigh[0] = 0;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 1;        /* first or second subelement at face */
        }
      }
      if (elem_face == 2) {     /* lower face */
        if (phex_w_sub_pseudo_neigh->p8q.x == phex_w_sub_elem->p8q.x) {
          location_neigh[0] = 1;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 0;        /* first subelement at face */
        }
        if (phex_w_sub_pseudo_neigh->p8q.x != phex_w_sub_elem->p8q.x) {
          location_neigh[0] = 1;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 1;        /* second subelement at face */
        }
      }
      if (elem_face == 3) {     /* upper face */
        if (phex_w_sub_pseudo_neigh->p8q.x == phex_w_sub_elem->p8q.x) {
          location_neigh[0] = 3;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 1;        /* first or second subelement at face */
        }
        if (phex_w_sub_pseudo_neigh->p8q.x != phex_w_sub_elem->p8q.x) {
          location_neigh[0] = 3;        /* face */
          location_neigh[1] = 1;        /* split */
          location_neigh[2] = 0;        /* first subelement at face */
        }
      }
    }
    /* the pseudo_neigh tranaition cell has the same level as elem 
     * Note that the level of the trnasition cell can not be higher as the level of elem in this case, 
     * since elem would then be a subelement in a transtion cell. */
    if (phex_w_sub_pseudo_neigh->p8q.level == phex_w_sub_elem->p8q.level) {
      if (elem_face == 0) {     /* left face */
        location_neigh[0] = 2;  /* face */
        location_neigh[1] = 0;  /* not split */
        location_neigh[2] = 0;  /* first (only) subelement at face */
      }
      if (elem_face == 1) {     /* right face */
        location_neigh[0] = 0;  /* face */
        location_neigh[1] = 0;  /* not split */
        location_neigh[2] = 0;  /* first (only) subelement at face */
      }
      if (elem_face == 2) {     /* lower face */
        location_neigh[0] = 1;  /* face */
        location_neigh[1] = 0;  /* not split */
        location_neigh[2] = 0;  /* first (only) subelement at face */
      }
      if (elem_face == 3) {     /* upper face */
        location_neigh[0] = 3;  /* face */
        location_neigh[1] = 0;  /* not split */
        location_neigh[2] = 0;  /* first (only) subelement at face */
      }
    }

    /* check, that a neighbor is found and the location array is adjusted */
    T8_ASSERT (location_neigh[0] >= 0 && location_neigh[1] >= 0
               && location_neigh[2] >= 0);

    /* Depending on the location of elem, we have filled location_neigh with the data of the real neighbor.
     * This data will be used to determine the sub_id of the neighbor within the transition cell of pseudo_neigh. */
    return
      t8_element_get_id_from_location (t8_element_get_transition_type
                                       (pseudo_neigh), location_neigh);
  }
  return -1;                    /* return negative if no neighbor element could be found */

  // SC_ABORT_NOT_REACHED();
}

int
t8_subelement_scheme_hex_c::t8_element_get_id_from_location (int type,
                                                              int location[])
{
  T8_ASSERT (type >= 0 && type <= T8_SUB_HEX_MAX_TRANSITION_TYPE);

  int                 sub_id, subelements_count = 0;
  double              type_temp = double (type);        // would work for ints but we use libc pow(double, double)
  int                 binary_type[P8EST_FACES] = { };
  int                 binary_type_clockwise[P8EST_FACES] = { };

  /* get the type as a binary array */
  int                 iface;
  for (iface = 0; iface < P8EST_FACES; iface++) {
    if (type_temp >= pow (2.0, 4 - (iface + 1))) {
      binary_type[iface] = 1;
      type_temp -= pow (2.0, 4 - (iface + 1));
    }
    else {
      binary_type[iface] = 0;
    }
  }

  for (iface = 0; iface < P8EST_FACES; iface++) {       /* rearrange the binary type to be in clockwise order of the faces, starting with the left face */
    binary_type_clockwise[iface] =
      binary_type[subelement_location_to_parent_face[iface]];
  }

  /* count the number of elements up to the given location */
  int                 element_count;
  for (element_count = 0; element_count <= location[0]; element_count++) {
    if (element_count == location[0]) {
      if (location[1] == 0) {
        subelements_count += 1;
      }
      else {
        if (location[2] == 0) {
          subelements_count += 1;
        }
        else {
          subelements_count += 2;
        }
      }
    }
    else {
      subelements_count += binary_type_clockwise[element_count] + 1;
    }
  }

  /* get the sub_id */
  sub_id = subelements_count - 1;

  return sub_id;
}

int
t8_subelement_scheme_hex_c::t8_element_get_face_number_of_hypotenuse (const
                                                                       t8_element_t
                                                                       *elem)
{
  T8_ASSERT (t8_element_is_valid (elem));
  T8_ASSERT (t8_element_is_subelement (elem));

  int                 location[3] = { };
  t8_element_get_location_of_subelement (elem, location);

  int                 split = location[1];
  int                 second = location[2];

  if (!split) {                 /* if the face is not split, then the hypotenuse is always face number one */
    return 1;
  }
  else {
    if (!second) {              /* otherwise, the subelment is mirrored, depending on the value of 'second' */
      return 0;
    }
    else {
      return 2;
    }
  }
}

void
t8_subelement_scheme_hex_c::t8_element_new (int length, t8_element_t **elem) const
{
  /* allocate memory */
  t8_default_scheme_common_c::t8_element_new (length, elem);

  int                 elem_count;
  for (elem_count = 0; elem_count < length; elem_count++) {
    t8_hex_with_subelements *phex_w_sub =
      (t8_hex_with_subelements *) elem[elem_count];
    t8_element_init (1, elem[elem_count], 0);
    /* set dimension of hex to 2 */
    T8_HEX_SET_TDIM ((p8est_quadrant_t *) & phex_w_sub->p8q, 2);
  }
}

void
t8_subelement_scheme_hex_c::t8_element_init (int length, t8_element_t *elem,
                                              int new_called) const
{
  t8_hex_with_subelements *phex_w_sub = (t8_hex_with_subelements *) elem;

  int                 elem_count;

  for (elem_count = 0; elem_count < length; elem_count++) {
    /* initalize subelement parameters */
    phex_w_sub[elem_count].transition_type = 0;
    phex_w_sub[elem_count].subelement_id = 0;

#ifdef T8_ENABLE_DEBUG
    /* In debugging mode we iterate over all length many elements and 
     * set their hex to the level 0 hex with ID 0. */
    if (!new_called) {
      p8est_quadrant_t   *hex = &phex_w_sub[elem_count].p8q;
      p8est_quadrant_set_morton (hex, 0, 0);
      T8_HEX_SET_TDIM (hex, 2);
      T8_ASSERT (p8est_quadrant_is_extended (hex));
    }
#endif
  }
}

int
t8_subelement_scheme_hex_c::t8_element_scheme_supports_transitioning (void)
{
  return T8_HEX_TRANSITION_IS_IMPLEMENTED;
}

int
t8_subelement_scheme_hex_c::t8_element_transition_scheme_is_conformal (void)
{
  return T8_HEX_TRANSITION_SCHEME_IS_CONFORMAL;
}

#ifdef T8_ENABLE_DEBUG
void
t8_subelement_scheme_hex_c::t8_element_debug_print (const t8_element_t *elem) const
{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  t8_productionf ("\n|------------ t8_element_debug_print: ------------|"
                  "\n|    Transition Type:     %i"
                  "\n|    Subelement ID:       %i"
                  "\n|    Anchor (Morton):     (%i,%i)"
                  "\n|    Anchor (ref coords): (%lf,%lf)"
                  "\n|    Level:               %i"
                  "\n|-------------------------------------------------|\n",
                  phex_w_sub->transition_type, phex_w_sub->subelement_id,
                  phex_w_sub->p8q.x, phex_w_sub->p8q.y,
                  (double) phex_w_sub->p8q.x / (double) P8EST_ROOT_LEN,
                  (double) phex_w_sub->p8q.y / (double) P8EST_ROOT_LEN,
                  phex_w_sub->p8q.level);

  /* if the element is not valid, abort, but after printing */
  T8_ASSERT (t8_element_is_valid (elem));
}

/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_subelement_scheme_hex_c::t8_element_is_valid (const t8_element_t * elem) const 
/* *INDENT-ON* */

{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;
  const p8est_quadrant_t *q = &phex_w_sub->p8q;

  /* the 4pest quadrant AND the subelement values must be valid such that the whole element is valid */
  return (p8est_quadrant_is_extended (q)
          && t8_element_subelement_values_are_valid (elem));
}

/* *INDENT-OFF* */
/* indent bug, indent adds a second "const" modifier */
int
t8_subelement_scheme_hex_c::t8_element_subelement_values_are_valid (const
                                                                 t8_element_t *
                                                                 elem) const
/* *INDENT-ON* */

{
  const t8_hex_with_subelements *phex_w_sub =
    (const t8_hex_with_subelements *) elem;

  return ((phex_w_sub->transition_type >= 0 &&
           phex_w_sub->transition_type <= T8_SUB_HEX_MAX_TRANSITION_TYPE)
          || phex_w_sub->transition_type == 0) &&
    ((phex_w_sub->subelement_id >= 0 &&
      phex_w_sub->subelement_id <= T8_SUB_HEX_MAX_SUBELEMENT_ID)
     || phex_w_sub->subelement_id == 0);
}
#endif

/* Constructor */
t8_subelement_scheme_hex_c::t8_subelement_scheme_hex_c (void)
{
  eclass = T8_ECLASS_HEX;
  element_size = sizeof (t8_phex_sub_t);
  ts_context = sc_mempool_new (element_size);
}

t8_subelement_scheme_hex_c::~t8_subelement_scheme_hex_c ()
{
  /* This destructor is empty since the destructor of the
   * default_common scheme is called automatically and it
   * suffices to destroy the hex_scheme.
   * However we need to provide an implementation of the destructor
   * and hence this empty function. */
}

T8_EXTERN_C_END ();
