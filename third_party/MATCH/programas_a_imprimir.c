/************************************************************************
 * <AUTO EXTRACT>
 *
 * ROUTINE: atApplyTrans
 *
 * DESCRIPTION:
 * Given a list of s_star structures, apply the given TRANS to each item in
 * the list, modifying the "x" and "y" values.
 *
 * The TRANS structure has 6 coefficients, which are used as follows:
 *
 *       x' = A + Bx + Cx
 *       y' = D + Ex + Fy
 *
 * Actually, this is a top-level "driver" routine that calls smaller
 * functions to perform actual tasks.  It mostly creates the proper
 * inputs and outputs for the smaller routines.
 *
 *
 * RETURN:
 *    SH_SUCCESS         if all goes well
 *    SH_GENERIC_ERROR   if an error occurs
 *
 * </AUTO>
 */

int
atApplyTrans
   (
   int num,              /* I: number of stars in the linked list */
   s_star *star_list,    /* I/O: modify x,y coords of objects in this list */
   TRANS *trans          /* I: use this TRANS to transform the coords of */
                         /*       items in the list */
   )
{
   int i;
   struct s_star *ptr, *star_array;

   shAssert(star_list != NULL);

   /* 
    * convert the linked list to an array
    */
   star_array = list_to_array(num, star_list);

#ifdef DEBUG
   printf("before applying TRANS \n");
   print_star_array(star_array, num);
#endif

   /*
    * next, apply the transformation to each element of the array
    */
   apply_trans(star_array, num, trans);

#ifdef DEBUG
   printf("after applying TRANS \n");
   print_star_array(star_array, num);
#endif

   /*
    * transfer the coord values from the array back into the list
    */
   for (ptr = star_list, i = 0; i < num; i++, ptr = ptr->next) {
      shAssert(ptr != NULL);
      ptr->x = star_array[i].x;
      ptr->y = star_array[i].y;
   }

   /*
    * delete the array
    */
   free_star_array(star_array);

   /*
    * all done!
    */

   return(SH_SUCCESS);
}

/************************************************************************
 * <AUTO EXTRACT>
 *
 * ROUTINE: atMatchLists
 *
 * DESCRIPTION:
 * Given 2 lists of s_star structures, 
 * which have ALREADY been transformed so that the "x"
 * and "y" coordinates of each list are close to each other
 * (i.e. matching items from each list have very similar "x" and "y")
 * this routine attempts to find all instances of matching items
 * from the 2 lists.
 *
 * We consider a "match" to be the closest coincidence of centers 
 * which are within "radius" pixels of each other.  
 *
 * Use a slow, but sure, algorithm.
 *
 * We will match objects from A --> B.  It is possible to have several
 * As that match to the same B:
 *
 *           A1 -> B5   and A2 -> B5
 *
 * This function finds such multiple-match items and deletes all but
 * the closest of the matches.
 *
 * place the elems of A that are matches into output list J
 *                    B that are matches into output list K
 *                    A that are not matches into output list L
 *                    B that are not matches into output list M
 *
 * Place a count of the number of matching pairs into the final
 * argument, 'num_matches'.
 *
 *
 * RETURN:
 *    SH_SUCCESS         if all goes well
 *    SH_GENERIC_ERROR   if an error occurs
 *
 * </AUTO>
 */

int
atMatchLists
   (
   int numA,                /* I: number of stars in list A */
   s_star *listA,           /* I: first input list of items to be matched */
   int numB,                /* I: number of stars in list B */
   s_star *listB,           /* I: second list of items to be matched */
   double radius,           /* I: maximum radius for items to be a match */
   char *basename,          /* I: base of filenames used to store the */
                            /*      output; extension indicates which */
                            /*      .mtA    items from A that matched */
                            /*      .mtB    items from B that matched */
                            /*      .unA    items from A that didn't match */
                            /*      .unB    items from A that didn't match */
	int *num_matches         /* O: number of matching pairs we find */
   )
{
   s_star *star_array_A;
   int num_stars_A;
   s_star *star_array_B;
   int num_stars_B;
   s_star *star_array_J, *star_array_K, *star_array_L, *star_array_M;
   int num_stars_J, num_stars_K, num_stars_L, num_stars_M;
   char filename[100];

   shAssert(listA != NULL);
   shAssert(listB != NULL);

   /* convert from linked lists to arrays */
   num_stars_A = numA;
   num_stars_B = numB;
   star_array_A = list_to_array(numA, listA);
   star_array_B = list_to_array(numB, listB);
   shAssert(star_array_A != NULL);
   shAssert(star_array_B != NULL);

   /* reset the 'id' fields in the arrays to match those in the lists */
   reset_array_ids(listA, numA, star_array_A);
   reset_array_ids(listB, numB, star_array_B);


   /* do the matching process */
   if (match_arrays_slow(star_array_A, num_stars_A,
                         star_array_B, num_stars_B,
                         radius,
                         &star_array_J, &num_stars_J,
                         &star_array_K, &num_stars_K,
                         &star_array_L, &num_stars_L,
                         &star_array_M, &num_stars_M) != SH_SUCCESS) {
      shError("atMatchLists: match_arrays_slow fails");
      return(SH_GENERIC_ERROR);
   }

	/* 
	 * Set the 'num_matches' value to the number of matching pairs
	 *   (we could just as well use num_stars_K) 
	 */
	*num_matches = num_stars_J;

   /*
    * now write the output into ASCII text files, each of which starts
    * with 'basename', but has a different extension.  
    *
    *    basename.mtA    stars from list A that did match         array J
    *    basename.mtB    stars from list A that did match         array K
    *    basename.unA    stars from list A that did NOT match     array L
    *    basename.unB    stars from list A that did NOT match     array M
    */
   sprintf(filename, "%s.mtA", basename);
   write_array(num_stars_J, star_array_J, filename);
   sprintf(filename, "%s.mtB", basename);
   write_array(num_stars_K, star_array_K, filename);
   sprintf(filename, "%s.unA", basename);
   write_array(num_stars_L, star_array_L, filename);
   sprintf(filename, "%s.unB", basename);
   write_array(num_stars_M, star_array_M, filename);

   /*
    * all done!
    */
   free_star_array(star_array_J);
   free_star_array(star_array_K);
   free_star_array(star_array_L);
   free_star_array(star_array_M);

   return(SH_SUCCESS);
}

/************************************************************************
 * <AUTO EXTRACT>
 *
 * ROUTINE: atFindTrans
 *
 * DESCRIPTION:
 * This function is based on the algorithm described in Valdes et al.,
 * PASP 107, 1119 (1995).  It tries to 
 *         a. match up objects in the two chains
 *         a. find a coordinate transformation that takes coords in
 *               objects in chainA and changes to those in chainB.
 *
 *
 * Actually, this is a top-level "driver" routine that calls smaller
 * functions to perform actual tasks.  It mostly creates the proper
 * inputs and outputs for the smaller routines.
 *
 * RETURN:
 *    SH_SUCCESS         if all goes well
 *    SH_GENERIC_ERROR   if an error occurs
 *
 * </AUTO>
 */

int
atFindTrans
   (
   int numA,             /* I: number of stars in list A */
   struct s_star *listA, /* I: match this set of objects with list B */
   int numB,             /* I: number of stars in list B */
   struct s_star *listB, /* I: match this set of objects with list A */
   double radius,        /* I: max radius in triangle-space allowed for */
                         /*       a pair of triangles to match */
   int nobj,             /* I: max number of bright stars to use in creating */
                         /*       triangles for matching from each list */
   double min_scale,     /* I: minimum permitted relative scale factor */
                         /*       if -1, any scale factor is allowed */
   double max_scale,     /* I: maximum permitted relative scale factor */
                         /*       if -1, any scale factor is allowed */
   double rotation_deg,  /* I: desired relative angle of coord systems (deg) */
                         /*       if AT_MATCH_NOANGLE, any orientation is allowed */
   double tolerance_deg, /* I: allowed range of orientation angles (deg) */
                         /*       if AT_MATCH_NOANGLE, any orientation is allowed */
   int max_iter,         /* I: go through at most this many iterations */
                         /*       in the iter_trans() loop. */
   double halt_sigma,    /* I: halt the fitting procedure if the mean */
                         /*       residual becomes this small */
   TRANS *trans          /* O: place into this TRANS structure's fields */
                         /*       the coeffs which convert coords of chainA */
                         /*       into coords of chainB system. */
   )
{
   int i, nbright, min;
   int num_stars_A;          /* number of stars in chain A */
   int num_stars_B;          /* number of stars in chain B */
   int num_triangles_A;      /* number of triangles formed from chain A */
   int num_triangles_B;      /* number of triangles formed from chain B */
   int **vote_matrix; 
   int *winner_votes;        /* # votes gotten by top pairs of matched stars */
   int *winner_index_A;      /* elem i in this array is index in star array A */
                             /*    which matches ... */
   int *winner_index_B;      /* elem i in this array, index in star array B */
   int start_pairs;
   s_star *star_array_A = NULL;
   s_star *star_array_B = NULL;
   s_triangle *triangle_array_A = NULL;
   s_triangle *triangle_array_B = NULL;

   num_stars_A = numA;
   num_stars_B = numB;
   star_array_A = list_to_array(numA, listA);
   star_array_B = list_to_array(numB, listB);

#ifdef DEBUG3
   test_routine();
#endif

   shAssert(star_array_A != NULL);
   shAssert(star_array_B != NULL);

   switch (trans->order) {
   case AT_TRANS_LINEAR:
      start_pairs = AT_MATCH_STARTN_LINEAR;
      break;
   case AT_TRANS_QUADRATIC:
      start_pairs = AT_MATCH_STARTN_QUADRATIC;
      break;
   case AT_TRANS_CUBIC:
      start_pairs = AT_MATCH_STARTN_CUBIC;
      break;
   default:
      shError("atFindTrans: invalid trans->order %d ", trans->order);
      break;
   }

   

   /*
    * here we check to see if each list of stars contains a
    * required minimum number of stars.  If not, we return with
    * an error message, and SH_GENERIC_ERROR.  
    *
    * In addition, we check to see that each list has at least 'nobj'
    * items.  If not, we set 'nbright' to the minimum of the two
    * list lengths, and print a warning message so the user knows 
    * that we're using fewer stars than he asked.
    *
    * On the other hand, if the user specifies a value of "nobj" which
    * is too SMALL, then we ignore it and use the smallest valid
    * value (which is start_pairs).
    */
   min = (num_stars_A < num_stars_B ? num_stars_A : num_stars_B);
   if (min < start_pairs) {
      shError("atFindTrans: only %d stars in list(s), require at least %d",
	       min, start_pairs);
      free_star_array(star_array_A);
      free_star_array(star_array_B);
      return(SH_GENERIC_ERROR);
   }
   if (nobj > min) {
      shDebug(AT_MATCH_ERRLEVEL,
	       "atFindTrans: using only %d stars, fewer than requested %d",
	       min, nobj);
      nbright = min;
   }
   else {
      nbright = nobj;
   }
   if (nbright < start_pairs) {
      shDebug(AT_MATCH_ERRLEVEL,
	       "atFindTrans: must use %d stars, more than requested %d",
	       start_pairs, nobj);
      nbright = start_pairs;
   }

   /* this is a sanity check on the above checks */
   shAssert((nbright >= start_pairs) && (nbright <= min));


#ifdef DEBUG
   printf("here comes star array A\n");
   print_star_array(star_array_A, num_stars_A);
   printf("here comes star array B\n");
   print_star_array(star_array_B, num_stars_B);
#endif
   
   /*
    * we now convert each list of stars into a list of triangles, 
    * using only a subset of the "nbright" brightest items in each list.
    */
   triangle_array_A = stars_to_triangles(star_array_A, num_stars_A, 
	    nbright, &num_triangles_A);
   shAssert(triangle_array_A != NULL);
   triangle_array_B = stars_to_triangles(star_array_B, num_stars_B, 
	    nbright, &num_triangles_B);
   shAssert(triangle_array_B != NULL);


   /*
    * Now we prune the triangle arrays to eliminate those with 
    * ratios (b/a) > AT_MATCH_RATIO,
    * since Valdes et al. say that this speeds things up and eliminates
    * lots of closely-packed triangles.
    */
   prune_triangle_array(triangle_array_A, &num_triangles_A);
   prune_triangle_array(triangle_array_B, &num_triangles_B);
#ifdef DEBUG2
   printf("after pruning, here comes triangle array A\n");
   print_triangle_array(triangle_array_A, num_triangles_A,
	       star_array_A, num_stars_A);
   printf("after pruning, here comes triangle array B\n");
   print_triangle_array(triangle_array_B, num_triangles_B,
	       star_array_B, num_stars_B);
#endif


   /*
    * Next, we want to try to match triangles in the two arrays.
    * What we do is to create a "vote matrix", which is a 2-D array
    * with "nbright"-by-"nbright" cells.  The cell with
    * coords [i][j] holds the number of matched triangles in which
    * 
    *        item [i] in star_array_A matches item [j] in star_array_B
    *
    * We'll use this "vote_matrix" to figure out a first guess
    * at the transformation between coord systems.
    *
    * Note that if there are fewer than "nbright" stars
    * in either list, we'll still make the vote_matrix 
    * contain "nbright"-by-"nbright" cells ...
    * there will just be a lot of cells filled with zero.
    */
   vote_matrix = make_vote_matrix(star_array_A, num_stars_A,
                                  star_array_B, num_stars_B,
                                  triangle_array_A, num_triangles_A,
                                  triangle_array_B, num_triangles_B,
                                  nbright, radius, min_scale, max_scale,
                                  rotation_deg, tolerance_deg);


   /*
    * having made the vote_matrix, we next need to pick the 
    * top 'nbright' vote-getters.  We call 'top_vote_getters'
    * and are given, in its output arguments, pointers to three
    * arrays, each of which has 'nbright' elements pertaining
    * to a matched pair of STARS:
    * 
    *       winner_votes[]    number of votes of winners, in descending order
    *       winner_index_A[]  index of star in star_array_A 
    *       winner_index_B[]  index of star in star_array_B
    *
    * Thus, the pair of stars which matched in the largest number
    * of triangles will be 
    *
    *       star_array_A[winner_index_A[0]]    from array A
    *       star_array_B[winner_index_A[0]]    from array B
    *
    * and the pair of stars which matched in the second-largest number
    * of triangles will be 
    *
    *       star_array_A[winner_index_A[1]]    from array A
    *       star_array_B[winner_index_A[1]]    from array B
    * 
    * and so on.
    */
   top_vote_getters(vote_matrix, nbright, 
	       &winner_votes, &winner_index_A, &winner_index_B);

   /*
    * here, we disqualify any of the top vote-getters which have
    * fewer than AT_MATCH_MINVOTES votes.  This may decrease the
    * number of valid matched pairs below 'nbright', so we
    * re-set nbright if necessary.
    */
   for (i = 0; i < nbright; i++) {
      if (winner_votes[i] < AT_MATCH_MINVOTES) {
#ifdef DEBUG
         printf("disqualifying all winners after number %d, nbright now %d\n",
               i, i);
#endif
         nbright = i;
         break;
      }
   }


   /*
    * next, we take the "top" matched pairs of coodinates, and
    * figure out a transformation of the form
    *
    *       x' = A + Bx + Cx
    *       y' = D + Ex + Fy
    *
    * (i.e. a TRANS structure) which converts the coordinates
    * of objects in chainA to those in chainB.
    */
   if (iter_trans(nbright, star_array_A, num_stars_A, 
                       star_array_B, num_stars_B,
                       winner_votes, winner_index_A, winner_index_B, 
                       RECALC_NO, max_iter, halt_sigma, trans) != SH_SUCCESS) {

      shError("atFindTrans: iter_trans unable to create a valid TRANS");
      free_star_array(star_array_A);
      free_star_array(star_array_B);
      return(SH_GENERIC_ERROR);
   }

#ifdef DEBUG
   printf("  after calculating new TRANS structure, here it is\n");
   print_trans(trans);
#endif

   /*
    * clean up memory we allocated during the matching process 
    */
   free_star_array(star_array_A);
   free_star_array(star_array_B);

   return(SH_SUCCESS);
}

/************************************************************************
 * <AUTO EXTRACT>
 *
 * ROUTINE: atRecalcTrans
 *
 * DESCRIPTION:
 * Given two lists of stars which ALREADY have been matched,
 * this routine finds a coord transformation which takes coords
 * of stars in list A to those in list B.
 * 
 * We can skip all the matching-triangles business, which makes this
 * _much_ faster than atFindTrans.
 *
 * RETURN:
 *    SH_SUCCESS         if all goes well
 *    SH_GENERIC_ERROR   if an error occurs
 *
 * </AUTO>
 */

int
atRecalcTrans
   (
   int numA,             /* I: number of stars in list A */
   struct s_star *listA, /* I: match this set of objects with list B */
   int numB,             /* I: number of stars in list B */
   struct s_star *listB, /* I: match this set of objects with list A */
   int max_iter,         /* I: go through at most this many iterations */
                         /*       in the iter_trans() loop. */
   double halt_sigma,    /* I: halt the fitting procedure if the mean */
                         /*       residual becomes this small */
   TRANS *trans          /* O: place into this TRANS structure's fields */
                         /*       the coeffs which convert coords of chainA */
                         /*       into coords of chainB system. */
   )
{
   int i, nbright, min;
   int num_stars_A;          /* number of stars in chain A */
   int num_stars_B;          /* number of stars in chain B */
   int *winner_votes;        /* # votes gotten by top pairs of matched stars */
   int *winner_index_A;      /* elem i in this array is index in star array A */
                             /*    which matches ... */
   int *winner_index_B;      /* elem i in this array, index in star array B */
   int start_pairs;
   s_star *star_array_A = NULL;
   s_star *star_array_B = NULL;

   num_stars_A = numA;
   num_stars_B = numB;
   star_array_A = list_to_array(numA, listA);
   star_array_B = list_to_array(numB, listB);

   shAssert(star_array_A != NULL);
   shAssert(star_array_B != NULL);

   switch (trans->order) {
   case AT_TRANS_LINEAR:
      start_pairs = AT_MATCH_STARTN_LINEAR;
      break;
   case AT_TRANS_QUADRATIC:
      start_pairs = AT_MATCH_STARTN_QUADRATIC;
      break;
   case AT_TRANS_CUBIC:
      start_pairs = AT_MATCH_STARTN_CUBIC;
      break;
   default:
      shError("atRecalcTrans: invalid trans->order %d ", trans->order);
      break;
   }

   /*
    * here we check to see if each list of stars contains a
    * required minimum number of stars.  If not, we return with
    * an error message, and SH_GENERIC_ERROR.  
    *
    * We set 'nbright' to the minimum of the two list lengths
    */
   min = (num_stars_A < num_stars_B ? num_stars_A : num_stars_B);
   if (min < start_pairs) {
      shError("atRecalcTrans: only %d stars in list(s), require at least %d",
	       min, start_pairs);
      free_star_array(star_array_A);
      free_star_array(star_array_B);
      return(SH_GENERIC_ERROR);
   }
   nbright = min;

   /* this is a sanity check on the above checks */
   shAssert((nbright >= start_pairs) && (nbright <= min));


#ifdef DEBUG
   printf("here comes star array A\n");
   print_star_array(star_array_A, num_stars_A);
   printf("here comes star array B\n");
   print_star_array(star_array_B, num_stars_B);
#endif
   

   /*
    * We need to create dummy arrays for 'winner_votes', and the
    * 'winner_index' arrays.  We already know that all these stars
    * are good matches, and so we can just create some arrays
    * and fill them with identical numbers.  They aren't used by
    * iter_trans(), anyway.
    */
   winner_votes = (int *) shMalloc(nbright*sizeof(int));
   winner_index_A = (int *) shMalloc(nbright*sizeof(int));
   winner_index_B = (int *) shMalloc(nbright*sizeof(int));
   for (i = 0; i < nbright; i++) {
      winner_votes[i] = 100;
      winner_index_A[i] = i;
      winner_index_B[i] = i;
   }

   /*
    * next, we take ALL the matched pairs of coodinates, and
    * figure out a transformation of the form
    *
    *       x' = A + Bx + Cx
    *       y' = D + Ex + Fy
    *
    * (i.e. a TRANS structure) which converts the coordinates
    * of objects in list A to those in list B
    */
   if (iter_trans(nbright, star_array_A, num_stars_A, 
                       star_array_B, num_stars_B,
                       winner_votes, winner_index_A, winner_index_B, 
                       RECALC_YES, max_iter, halt_sigma, trans) != SH_SUCCESS) {

      shError("atRecalcTrans: iter_trans unable to create a valid TRANS");
      free_star_array(star_array_A);
      free_star_array(star_array_B);
      return(SH_GENERIC_ERROR);
   }

#ifdef DEBUG
   printf("  after calculating new TRANS structure, here it is\n");
   print_trans(trans);
#endif

   /*
    * clean up memory we allocated during the matching process 
    */
   shFree(winner_votes);
   shFree(winner_index_A);
   shFree(winner_index_B);
   free_star_array(star_array_A);
   free_star_array(star_array_B);

   return(SH_SUCCESS);
}

	/**********************************************************************
	 * PROCEDURE: prepare_to_recalc
	 *
	 * DESCRIPTION: This function sets us up to call "atRecalcTrans".
	 *              We have already found (or been given) a TRANS, and
	 *              used it to match up items from list A and list B.
	 *              Those matched items are in a pair of files
    *              with names based on 'outfile', but with
    *              extensions "mtA" and "mtB".  
	 *              Ex: if 'outfile' is "matched",
    *              then the two sets of matched items are in 
    *
    *                    matched.mtA      matched.mtB
    *
    *              The format of each file is one star per line, 
	 *              with 4 fields:
    *  
    *                    internal_ID     xval   yval     mag
    *
    *              where the coords (xval, yval) are in system of list B.
	 *
	 *              We are about to use these good, matched items to 
	 *              find an improved TRANS -- which should take objects
	 *              from coord system A to coord system B.
	 *
	 *              In order to do that, we need to 
	 *
	 *                  a. Read in the good, matched items.  The coordinates
	 *                     of objects in list A will have been transformed
	 *                     to their corresponding values in coord system 
	 *                     of list B, so ...
	 *
	 *                  b. We must then re-set the coords of the items in
	 *                     list A to their original values, so that we can
	 *                     re-calculate a TRANS which takes the coords
	 *                     from system A to system B.
	 *
	 *              We also take this opportunity to compare the transformed
	 *              positions of items in list A against the positions of
	 *              the matching objects in list B.  We calculate the
	 *              RMS of the differences in both "x" and "y" directions,
	 *              and place them into the "sx" and "sy" members of
	 *              the current TRANS.
	 *  
	 *
	 * RETURNS:
	 *    0             if all goes well
	 *    1             if there's an error
	 */

static int
prepare_to_recalc
	(
	char *outfile,                   /* I: stem of files with matched items */
	int *num_matched_A,              /* O: number of stars in matched set */
	                                 /*      from list A */
	struct s_star **matched_list_A,  /* O: fill this with matched items from */
	                                 /*      list A, in coord system B */
	int *num_matched_B,              /* O: number of stars in matched set */
	                                 /*      from list B */
	struct s_star **matched_list_B,  /* O: fill this with matched items from */
	                                 /*      list B, in coord system B */
	struct s_star *star_list_A_copy, /* O: fill this with matched items from */
	                                 /*      list A, with their orig coords  */
	TRANS *trans                     /* O: we calc herein the sx, sy fields  */
	                                 /*      so put them into this TRANS */
	)
{
	char matched_file_A[CMDBUFLEN + 4];
	char matched_file_B[CMDBUFLEN + 4];
	double Xrms,Yrms;

	shAssert(outfile != NULL);

   sprintf(matched_file_A, "%s.mtA", outfile);
   if (read_matched_file(matched_file_A, num_matched_A, matched_list_A)
                                                            != SH_SUCCESS) {
       shError("read_matched_file can't read data from file %s", 
                                 matched_file_A);
       return(1);
   }
   sprintf(matched_file_B, "%s.mtB", outfile);
   if (read_matched_file(matched_file_B, num_matched_B, matched_list_B)
                                                             != SH_SUCCESS) {
       shError("read_matched_file can't read data from file %s", 
                                 matched_file_B);
       return(1);
   }

   /* here we find the rms of those stars we read in -- JPB 17/Jan/02 */
   if (atCalcRMS(*num_matched_A, *matched_list_A, 
		 *num_matched_B, *matched_list_B, 
		 &Xrms, &Yrms) != SH_SUCCESS) {
      shFatal("atCalcRMS fails on matched pairs");
   }
   trans->sx = Xrms;
   trans->sy = Yrms;
   /************************************************/

   if (reset_A_coords(*num_matched_A, *matched_list_A, star_list_A_copy) != 0) {
      shError("prepare_to_recalc: reset_A_coords returns with error");
      return(1);
   }

   return(0);
}
