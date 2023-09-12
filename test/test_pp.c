/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (c) 2009 RELIC Authors
 *
 * This file is part of RELIC. RELIC is legal property of its developers,
 * whose names are not listed here. Please refer to the COPYRIGHT file
 * for contact information.
 *
 * RELIC is free software; you can redistribute it and/or modify it under the
 * terms of the version 2.1 (or later) of the GNU Lesser General Public License
 * as published by the Free Software Foundation; or version 2.0 of the Apache
 * License as published by the Apache Software Foundation. See the LICENSE files
 * for more details.
 *
 * RELIC is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the LICENSE files for more details.
 *
 * You should have received a copy of the GNU Lesser General Public or the
 * Apache License along with RELIC. If not, see <https://www.gnu.org/licenses/>
 * or <https://www.apache.org/licenses/>.
 */

/**
 * @file
 *
 * Tests for pairings defined over prime elliptic curves.
 *
 * @ingroup test
 */

#include <stdio.h>

#include "relic.h"
#include "relic_test.h"
#include "relic_bench.h"

static int doubling12(void) {
	printf("----> in doubling 12\n");
	int code = RLC_ERR;
	bn_t k, n;
	ep_t p;
	ep2_t q, r, s;
	fp12_t e1, e2;

	bn_null(k);
	bn_null(n);
	ep_null(p);
	ep2_null(q);
	ep2_null(r);
	ep2_null(s);
	fp12_null(e1);
	fp12_null(e2);

	/* RLC_TRY { */
		bn_new(n);
		bn_new(k);
		ep_new(p);
		ep2_new(q);
		ep2_new(r);
		ep2_new(s);
		fp12_new(e1);
		fp12_new(e2);

		ep_curve_get_ord(n);

		TEST_CASE("miller doubling is correct") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			pp_dbl_k12(e1, r, q, p);
			pp_norm_k12(r, r);
			ep2_dbl(s, q);
			ep2_norm(s, s);
			TEST_ASSERT(ep2_cmp(r, s) == RLC_EQ, end);
		} TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_CASE("miller doubling in affine coordinates is correct") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			fp12_zero(e1);
			fp12_zero(e2);
			fp_neg(p->y, p->y);
			pp_dbl_k12_basic(e2, r, q, p);
			pp_exp_k12(e2, e2);
#if EP_ADD == PROJC || EP_ADD == JACOB
			/* Precompute. */
			fp_dbl(p->z, p->x);
			fp_add(p->x, p->z, p->x);
#endif
			pp_dbl_k12(e1, r, q, p);
			pp_exp_k12(e1, e1);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

/* #if EP_ADD == PROJC || EP_ADD == JACOB || !defined(STRIP)              */
/* 		TEST_CASE("miller doubling in projective coordinates is correct") {  */
/* 			ep_rand(p);                                                         */
/* 			ep2_rand(q);                                                        */
/* 			ep2_rand(r);                                                        */
/* 			fp12_zero(e1);                                                      */
/* 			fp12_zero(e2);                                                      */
/* 			// Precompute.                                                  */
/* 			fp_neg(p->y, p->y);                                                 */
/* 			fp_dbl(p->z, p->x);                                                 */
/* 			fp_add(p->x, p->z, p->x);                                           */
/* 			pp_dbl_k12_projc(e2, r, q, p);                                      */
/* 			pp_exp_k12(e2, e2);                                                 */
/* #if EP_ADD == BASIC                                                    */
/* 			// Revert precomputing.                                          */
/* 			fp_hlv(p->x, p->z);                                                 */
/* #endif                                                                 */
/* 			pp_dbl_k12(e1, r, q, p);                                            */
/* 			pp_exp_k12(e1, e1);                                                 */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);                       */
/* 		} TEST_END;                                                          */

/* #if PP_EXT == BASIC || !defined(STRIP)                                 */
/* 		TEST_CASE("basic projective miller doubling is correct") {           */
/* 			ep_rand(p);                                                         */
/* 			ep2_rand(q);                                                        */
/* 			ep2_rand(r);                                                        */
/* 			fp12_zero(e1);                                                      */
/* 			fp12_zero(e2);                                                      */
/* 			pp_dbl_k12_projc(e1, r, q, p);                                      */
/* 			pp_dbl_k12_projc_basic(e2, r, q, p);                                */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);                       */
/* 		} TEST_END;                                                          */
/* #endif                                                                 */

/* #if PP_EXT == LAZYR || !defined(STRIP)                                 */
/* 		TEST_CASE("lazy-reduced projective miller doubling is consistent") { */
/* 			ep_rand(p);                                                         */
/* 			ep2_rand(q);                                                        */
/* 			ep2_rand(r);                                                        */
/* 			fp12_zero(e1);                                                      */
/* 			fp12_zero(e2);                                                      */
/* 			pp_dbl_k12_projc(e1, r, q, p);                                      */
/* 			pp_dbl_k12_projc_lazyr(e2, r, q, p);                                */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);                       */
/* 		} TEST_END;                                                          */
/* #endif                                                                 */
// #endif /* EP_ADD = PROJC */
	/* }                              */
	/* RLC_CATCH_ANY {                */
	/* 	util_print("FATAL ERROR!\n"); */
	/* 	RLC_ERROR(end);               */
	/* }                              */
	code = RLC_OK;
  end:
	bn_free(n);
	bn_free(k);
	ep_free(p);
	ep2_free(q);
	ep2_free(r);
	ep2_free(s);
	fp12_free(e1);
	fp12_free(e2);
	return code;
}

static int addition12(void) {
	int code = RLC_ERR;
	bn_t k, n;
	ep_t p;
	ep2_t q, r, s;
	fp12_t e1, e2;

	bn_null(k);
	bn_null(n);
	ep_null(p);
	ep2_null(q);
	ep2_null(r);
	ep2_null(s);
	fp12_null(e1);
	fp12_null(e2);

	/* RLC_TRY { */
		bn_new(n);
		bn_new(k);
		ep_new(p);
		ep2_new(q);
		ep2_new(r);
		ep2_new(s);
		fp12_new(e1);
		fp12_new(e2);

		ep_curve_get_ord(n);

		TEST_CASE("miller addition is correct") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			ep2_copy(s, r);
			pp_add_k12(e1, r, q, p);
			pp_norm_k12(r, r);
			ep2_add(s, s, q);
			ep2_norm(s, s);
			TEST_ASSERT(ep2_cmp(r, s) == RLC_EQ, end);
		} TEST_END;

#if EP_ADD == BASIC || !defined(STRIP)
		TEST_CASE("miller addition in affine coordinates is correct") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			ep2_copy(s, r);
			fp12_zero(e1);
			fp12_zero(e2);
			pp_add_k12(e1, r, q, p);
			pp_exp_k12(e1, e1);
			pp_add_k12_basic(e2, s, q, p);
			pp_exp_k12(e2, e2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if EP_ADD == PROJC || EP_ADD == JACOB || !defined(STRIP)
		TEST_CASE("miller addition in projective coordinates is correct") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			ep2_copy(s, r);
			fp12_zero(e1);
			fp12_zero(e2);
			pp_add_k12(e1, r, q, p);
			pp_exp_k12(e1, e1);
			pp_add_k12_projc(e2, s, q, p);
			pp_exp_k12(e2, e2);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;

#if PP_EXT == BASIC || !defined(STRIP)
		TEST_CASE("basic projective miller addition is consistent") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			ep2_copy(s, r);
			fp12_zero(e1);
			fp12_zero(e2);
			pp_add_k12_projc(e1, r, q, p);
			pp_add_k12_projc_basic(e2, s, q, p);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif

#if PP_EXT == LAZYR || !defined(STRIP)
		TEST_CASE("lazy-reduced projective miller addition is consistent") {
			ep_rand(p);
			ep2_rand(q);
			ep2_rand(r);
			ep2_copy(s, r);
			fp12_zero(e1);
			fp12_zero(e2);
			pp_add_k12_projc(e1, r, q, p);
			pp_add_k12_projc_lazyr(e2, s, q, p);
			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);
		} TEST_END;
#endif
#endif /* EP_ADD = PROJC */
	/* }                              */
	/* RLC_CATCH_ANY {                */
	/* 	util_print("FATAL ERROR!\n"); */
	/* 	RLC_ERROR(end);               */
	/* }                              */
	code = RLC_OK;
  end:
	bn_free(n);
	bn_free(k);
	ep_free(p);
	ep2_free(q);
	ep2_free(r);
	ep2_free(s);
	fp12_free(e1);
	fp12_free(e2);
	return code;
}

static int pairing12(void) {
	int j, code = RLC_ERR;
	bn_t k, n;
	ep_t p[2];
	ep2_t q[2], r;
	fp12_t e1, e2;

	bn_null(k);
	bn_null(n);
	fp12_null(e1);
	fp12_null(e2);
	ep2_null(r);

	/* RLC_TRY { */
		bn_new(n);
		bn_new(k);
		fp12_new(e1);
		fp12_new(e2);
		ep2_new(r);

		for (j = 0; j < 2; j++) {
			ep_null(p[j]);
			ep2_null(q[j]);
			ep_new(p[j]);
			ep2_new(q[j]);
		}

		ep_curve_get_ord(n);

		/* TEST_CASE("pairing non-degeneracy is correct") { */
			ep_rand(p[0]);
			ep2_rand(q[0]);
			pp_map_k12(e1, p[0], q[0]);
		/* //	TEST_ASSERT(fp12_cmp_dig(e1, 1) != RLC_EQ, end); */
		/* 	ep_set_infty(p[0]);                                */
		/* 	pp_map_k12(e1, p[0], q[0]);                        */
		/* //	TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end); */
		/* 	ep_rand(p[0]);                                     */
		/* 	ep2_set_infty(q[0]);                               */
		/* 	pp_map_k12(e1, p[0], q[0]);                        */
		//	TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);
		/* } TEST_END; */

		/* TEST_CASE("pairing is bilinear") {                */
		/* 	ep_rand(p[0]);                                   */
		/* 	ep2_rand(q[0]);                                  */
		/* 	bn_rand_mod(k, n);                               */
		/* 	ep2_mul(r, q[0], k);                             */
		/* 	pp_map_k12(e1, p[0], r);                         */
		/* 	pp_map_k12(e2, p[0], q[0]);                      */
		/* 	fp12_exp(e2, e2, k);                             */
		/* 	TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);    */
		/* 	ep_mul(p[0], p[0], k);                           */
		/* 	pp_map_k12(e2, p[0], q[0]);                      */
		/* 	TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);    */
		/* 	ep_dbl(p[0], p[0]);                              */
		/* 	pp_map_k12(e2, p[0], q[0]);                      */
		/* 	fp12_sqr(e1, e1);                                */
		/* 	TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);    */
		/* 	ep2_dbl(q[0], q[0]);                             */
		/* 	pp_map_k12(e2, p[0], q[0]);                      */
		/* 	fp12_sqr(e1, e1);                                */
		/* 	TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);    */
		/* } TEST_END;                                       */

		/* TEST_CASE("multi-pairing is correct") {           */
		/* 	ep_rand(p[i % 2]);                               */
		/* 	ep2_rand(q[i % 2]);                              */
		/* 	pp_map_k12(e1, p[i % 2], q[i % 2]);              */
		/* 	ep_rand(p[1 - (i % 2)]);                         */
		/* 	ep2_set_infty(q[1 - (i % 2)]);                   */
		/* 	pp_map_sim_k12(e2, p, q, 2);                     */
		/* 	TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);    */
		/* 	ep_set_infty(p[1 - (i % 2)]);                    */
		/* 	ep2_rand(q[1 - (i % 2)]);                        */
		/* 	pp_map_sim_k12(e2, p, q, 2);                     */
		/* 	TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);    */
		/* 	ep2_set_infty(q[i % 2]);                         */
		/* 	pp_map_sim_k12(e2, p, q, 2);                     */
		/* 	TEST_ASSERT(fp12_cmp_dig(e2, 1) == RLC_EQ, end); */
		/* 	ep_rand(p[0]);                                   */
		/* 	ep2_rand(q[0]);                                  */
		/* 	pp_map_k12(e1, p[0], q[0]);                      */
		/* 	ep_rand(p[1]);                                   */
		/* 	ep2_rand(q[1]);                                  */
		/* 	pp_map_k12(e2, p[1], q[1]);                      */
		/* 	fp12_mul(e1, e1, e2);                            */
		/* 	pp_map_sim_k12(e2, p, q, 2);                     */
		/* 	TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);    */
		/* 	ep_neg(p[1], p[0]);                              */
		/* 	ep2_copy(q[1], q[0]);                            */
		/* 	pp_map_sim_k12(e1, p, q, 2);                     */
		/* 	TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end); */
		/* } TEST_END;                                       */

/* #if PP_MAP == TATEP || !defined(STRIP)                         */
/* 		TEST_CASE("tate pairing non-degeneracy is correct") {        */
/* 			ep_rand(p[0]);                                              */
/* 			ep2_rand(q[0]);                                             */
/* 			pp_map_tatep_k12(e1, p[0], q[0]);                           */
/* 			TEST_ASSERT(fp12_cmp_dig(e1, 1) != RLC_EQ, end);            */
/* 			ep_set_infty(p[0]);                                         */
/* 			pp_map_tatep_k12(e1, p[0], q[0]);                           */
/* 			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);            */
/* 			ep_rand(p[0]);                                              */
/* 			ep2_set_infty(q[0]);                                        */
/* 			pp_map_tatep_k12(e1, p[0], q[0]);                           */
/* 			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);            */
/* 		} TEST_END;                                                  */

/* 		TEST_CASE("tate pairing is bilinear") {                      */
/* 			ep_rand(p[0]);                                              */
/* 			ep2_rand(q[0]);                                             */
/* 			bn_rand_mod(k, n);                                          */
/* 			ep2_mul(r, q[0], k);                                        */
/* 			pp_map_tatep_k12(e1, p[0], r);                              */
/* 			pp_map_tatep_k12(e2, p[0], q[0]);                           */
/* 			fp12_exp(e2, e2, k);                                        */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep_mul(p[0], p[0], k);                                      */
/* 			pp_map_tatep_k12(e2, p[0], q[0]);                           */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep_dbl(p[0], p[0]);                                         */
/* 			pp_map_tatep_k12(e2, p[0], q[0]);                           */
/* 			fp12_sqr(e1, e1);                                           */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep2_dbl(q[0], q[0]);                                        */
/* 			pp_map_tatep_k12(e2, p[0], q[0]);                           */
/* 			fp12_sqr(e1, e1);                                           */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 		} TEST_END;                                                  */

/* 		TEST_CASE("tate multi-pairing is correct") {                 */
/* 			ep_rand(p[i % 2]);                                          */
/* 			ep2_rand(q[i % 2]);                                         */
/* 			pp_map_tatep_k12(e1, p[i % 2], q[i % 2]);                   */
/* 			ep_rand(p[1 - (i % 2)]);                                    */
/* 			ep2_set_infty(q[1 - (i % 2)]);                              */
/* 			pp_map_sim_tatep_k12(e2, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep_set_infty(p[1 - (i % 2)]);                               */
/* 			ep2_rand(q[1 - (i % 2)]);                                   */
/* 			pp_map_sim_tatep_k12(e2, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep2_set_infty(q[i % 2]);                                    */
/* 			pp_map_sim_tatep_k12(e2, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp_dig(e2, 1) == RLC_EQ, end);            */
/* 			ep_rand(p[0]);                                              */
/* 			ep2_rand(q[0]);                                             */
/* 			pp_map_tatep_k12(e1, p[0], q[0]);                           */
/* 			ep_rand(p[1]);                                              */
/* 			ep2_rand(q[1]);                                             */
/* 			pp_map_tatep_k12(e2, p[1], q[1]);                           */
/* 			fp12_mul(e1, e1, e2);                                       */
/* 			pp_map_sim_tatep_k12(e2, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep_neg(p[1], p[0]);                                         */
/* 			ep2_copy(q[1], q[0]);                                       */
/* 			pp_map_sim_tatep_k12(e1, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);            */
/* 		} TEST_END;                                                  */
/* #endif                                                         */

/* #if PP_MAP == WEIL || !defined(STRIP)                          */
/* 		TEST_CASE("weil pairing non-degeneracy is correct") {        */
/* 			ep_rand(p[0]);                                              */
/* 			ep2_rand(q[0]);                                             */
/* 			pp_map_weilp_k12(e1, p[0], q[0]);                           */
/* 			TEST_ASSERT(fp12_cmp_dig(e1, 1) != RLC_EQ, end);            */
/* 			ep_set_infty(p[0]);                                         */
/* 			pp_map_weilp_k12(e1, p[0], q[0]);                           */
/* 			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);            */
/* 			ep_rand(p[0]);                                              */
/* 			ep2_set_infty(q[0]);                                        */
/* 			pp_map_weilp_k12(e1, p[0], q[0]);                           */
/* 			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);            */
/* 		} TEST_END;                                                  */

/* 		TEST_CASE("weil pairing is bilinear") {                      */
/* 			ep_rand(p[0]);                                              */
/* 			ep2_rand(q[0]);                                             */
/* 			bn_rand_mod(k, n);                                          */
/* 			ep2_mul(r, q[0], k);                                        */
/* 			pp_map_weilp_k12(e1, p[0], r);                              */
/* 			pp_map_weilp_k12(e2, p[0], q[0]);                           */
/* 			fp12_exp(e2, e2, k);                                        */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep_mul(p[0], p[0], k);                                      */
/* 			pp_map_weilp_k12(e2, p[0], q[0]);                           */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep_dbl(p[0], p[0]);                                         */
/* 			pp_map_weilp_k12(e2, p[0], q[0]);                           */
/* 			fp12_sqr(e1, e1);                                           */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep2_dbl(q[0], q[0]);                                        */
/* 			pp_map_weilp_k12(e2, p[0], q[0]);                           */
/* 			fp12_sqr(e1, e1);                                           */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 		} TEST_END;                                                  */

/* 		TEST_CASE("weil multi-pairing is correct") {                 */
/* 			ep_rand(p[i % 2]);                                          */
/* 			ep2_rand(q[i % 2]);                                         */
/* 			pp_map_weilp_k12(e1, p[i % 2], q[i % 2]);                   */
/* 			ep_rand(p[1 - (i % 2)]);                                    */
/* 			ep2_set_infty(q[1 - (i % 2)]);                              */
/* 			pp_map_sim_weilp_k12(e2, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep_set_infty(p[1 - (i % 2)]);                               */
/* 			ep2_rand(q[1 - (i % 2)]);                                   */
/* 			pp_map_sim_weilp_k12(e2, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep2_set_infty(q[i % 2]);                                    */
/* 			pp_map_sim_weilp_k12(e2, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp_dig(e2, 1) == RLC_EQ, end);            */
/* 			ep_rand(p[0]);                                              */
/* 			ep2_rand(q[0]);                                             */
/* 			pp_map_weilp_k12(e1, p[0], q[0]);                           */
/* 			ep_rand(p[1]);                                              */
/* 			ep2_rand(q[1]);                                             */
/* 			pp_map_weilp_k12(e2, p[1], q[1]);                           */
/* 			fp12_mul(e1, e1, e2);                                       */
/* 			pp_map_sim_weilp_k12(e2, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep_neg(p[1], p[0]);                                         */
/* 			ep2_copy(q[1], q[0]);                                       */
/* 			pp_map_sim_weilp_k12(e1, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);            */
/* 		} TEST_END;                                                  */
/* #endif                                                         */

/* #if PP_MAP == OATEP || !defined(STRIP)                         */
/* 		TEST_CASE("optimal ate pairing non-degeneracy is correct") { */
/* 			ep_rand(p[0]);                                              */
/* 			ep2_rand(q[0]);                                             */
/* 			pp_map_oatep_k12(e1, p[0], q[0]);                           */
/* 			TEST_ASSERT(fp12_cmp_dig(e1, 1) != RLC_EQ, end);            */
/* 			ep_set_infty(p[0]);                                         */
/* 			pp_map_oatep_k12(e1, p[0], q[0]);                           */
/* 			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);            */
/* 			ep_rand(p[0]);                                              */
/* 			ep2_set_infty(q[0]);                                        */
/* 			pp_map_oatep_k12(e1, p[0], q[0]);                           */
/* 			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);            */
/* 		} TEST_END;                                                  */

/* 		TEST_CASE("optimal ate pairing is bilinear") {               */
/* 			ep_rand(p[0]);                                              */
/* 			ep2_rand(q[0]);                                             */
/* 			bn_rand_mod(k, n);                                          */
/* 			ep2_mul(r, q[0], k);                                        */
/* 			pp_map_oatep_k12(e1, p[0], r);                              */
/* 			ep_mul(p[0], p[0], k);                                      */
/* 			pp_map_oatep_k12(e2, p[0], q[0]);                           */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep_dbl(p[0], p[0]);                                         */
/* 			pp_map_oatep_k12(e2, p[0], q[0]);                           */
/* 			fp12_sqr(e1, e1);                                           */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep2_dbl(q[0], q[0]);                                        */
/* 			pp_map_oatep_k12(e2, p[0], q[0]);                           */
/* 			fp12_sqr(e1, e1);                                           */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 		} TEST_END;                                                  */

/* 		TEST_CASE("optimal ate multi-pairing is correct") {          */
/* 			ep_rand(p[i % 2]);                                          */
/* 			ep2_rand(q[i % 2]);                                         */
/* 			pp_map_oatep_k12(e1, p[i % 2], q[i % 2]);                   */
/* 			ep_rand(p[1 - (i % 2)]);                                    */
/* 			ep2_set_infty(q[1 - (i % 2)]);                              */
/* 			pp_map_sim_oatep_k12(e2, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep_set_infty(p[1 - (i % 2)]);                               */
/* 			ep2_rand(q[1 - (i % 2)]);                                   */
/* 			pp_map_sim_oatep_k12(e2, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep2_set_infty(q[i % 2]);                                    */
/* 			pp_map_sim_oatep_k12(e2, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp_dig(e2, 1) == RLC_EQ, end);            */
/* 			ep_rand(p[0]);                                              */
/* 			ep2_rand(q[0]);                                             */
/* 			pp_map_oatep_k12(e1, p[0], q[0]);                           */
/* 			ep_rand(p[1]);                                              */
/* 			ep2_rand(q[1]);                                             */
/* 			pp_map_oatep_k12(e2, p[1], q[1]);                           */
/* 			fp12_mul(e1, e1, e2);                                       */
/* 			pp_map_sim_oatep_k12(e2, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp(e1, e2) == RLC_EQ, end);               */
/* 			ep_neg(p[1], p[0]);                                         */
/* 			ep2_copy(q[1], q[0]);                                       */
/* 			pp_map_sim_oatep_k12(e1, p, q, 2);                          */
/* 			TEST_ASSERT(fp12_cmp_dig(e1, 1) == RLC_EQ, end);            */
/* 		} TEST_END;                                                  */
// #endif
	/* }                              */
	/* RLC_CATCH_ANY {                */
	/* 	util_print("FATAL ERROR!\n"); */
	/* 	RLC_ERROR(end);               */
	/* }                              */
	code = RLC_OK;
  end:
	bn_free(n);
	bn_free(k);
	fp12_free(e1);
	fp12_free(e2);
	ep2_free(r);

	for (j = 0; j < 2; j++) {
		ep_free(p[j]);
		ep2_free(q[j]);
	}
	return code;
}

int main(void) {
	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}

	util_banner("Tests for the PP module", 0);

	if (ep_param_set_any_pairf() == RLC_ERR) {
		RLC_THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}

	ep_param_print();

	util_banner("Arithmetic", 1);

	if (ep_param_embed() == 12) {
		printf("-----> embed 12\n");

		/* if (doubling12() != RLC_OK) { */
		/* 	core_clean();                */
		/* 	return 1;                    */
		/* }                             */

		/* if (addition12() != RLC_OK) { */
		/* 	core_clean();                */
		/* 	return 1;                    */
		/* }                             */

		if (pairing12() != RLC_OK) {
			core_clean();
			return 1;
		}
	}

	util_banner("All tests have passed.\n", 0);

	core_clean();
	return 0;
}
