/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_ENCODER_PASS2_STRATEGY_H_
#define AOM_AV1_ENCODER_PASS2_STRATEGY_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "common/tools_common.h"

struct AV1_COMP;
struct EncodeFrameParams;
// structure of accumulated stats and features in a gf group
typedef struct {
  double gf_group_err;
  double gf_group_raw_error;
  double gf_group_skip_pct;
  double gf_group_inactive_zone_rows;

  double mv_ratio_accumulator;
  double decay_accumulator;
  double zero_motion_accumulator;
  double loop_decay_rate;
  double last_loop_decay_rate;
  double this_frame_mv_in_out;
  double mv_in_out_accumulator;
  double abs_mv_in_out_accumulator;

  double avg_sr_coded_error;
  double avg_tr_coded_error;
  double avg_pcnt_second_ref;
  double avg_pcnt_third_ref;
  double avg_pcnt_third_ref_nolast;
  double avg_new_mv_count;
  double avg_wavelet_energy;
  double avg_raw_err_stdev;
  int non_zero_stdev_count;
} GF_GROUP_STATS;

static AOM_INLINE void cal_arf_ratio_MLE(const GF_GROUP_STATS *stats, const int num_mbs, RATE_CONTROL *rc) {
  const double MLEweight[] = { 
  0.00292498,-0.01291512, 0.20006554, 0.11155664, 0.00063053, 0.10770296,
  0.18069474, 0.06173544, 0.06092842,-0.18044636,-0.10560494,-0.22200811,
  0.02643944, 0.02235282, 0.19308184,-0.1407296,  0.00356003,-0.01174851,
  -0.00225225 };
  const double MLEbias = 0.00088391;

  fprintf(stderr, "\n calculating arf_ratio");
  double x[19];

  x[0]  = (double)(stats->gf_group_err);
  x[1]  = (double)(stats->gf_group_raw_error    / num_mbs / 16.0);
  x[2]  = (double)(stats->gf_group_skip_pct);
  x[3]  = (double)(stats->gf_group_inactive_zone_rows);
  x[4]  = (double)(stats->mv_ratio_accumulator);
  x[5]  = (double)(stats->decay_accumulator);
  x[6]  = (double)(stats->zero_motion_accumulator);
  x[7]  = (double)(stats->loop_decay_rate);
  x[8]  = (double)(stats->last_loop_decay_rate);
  x[9]  = (double)(stats->this_frame_mv_in_out);
  x[10] = (double)(stats->mv_in_out_accumulator);
  x[11] = (double)(stats->abs_mv_in_out_accumulator);
  x[12] = (double)(stats->avg_sr_coded_error    / num_mbs / 16.0);
  x[13] = (double)(stats->avg_tr_coded_error    / num_mbs / 16.0);
  x[14] = (double)(stats->avg_pcnt_second_ref);
  x[15] = (double)(stats->avg_new_mv_count      / num_mbs / 16.0);
  x[16] = (double)(stats->avg_wavelet_energy    / num_mbs / 16.0);
  x[17] = (double)(stats->avg_raw_err_stdev     / num_mbs / 16.0);
  x[18] = (double)(stats->non_zero_stdev_count);

  double delta_s = 0;
  fprintf(stderr,"\n");
  for(int i = 0; i < 19; i++){
    delta_s += x[i] * MLEweight[i];
    fprintf(stderr," x[%02d]=%6.3lf||ds=%6.3lf", i, x[i], delta_s);
  }
  delta_s += MLEbias;
  // use last_arf_ratio as s
  rc->last_arf_ratio += delta_s;
  // use this_arf_ratio as y(or r)
  rc->this_arf_ratio = 1/(1 + exp(-rc->last_arf_ratio));
  rc->this_arf_ratio = fclamp(rc->this_arf_ratio, 0.0, 1.0);
  fprintf(stderr,"\n delta_s=%6.3lf s=%6.3lf y=%6.3lf", 
    delta_s, rc->last_arf_ratio, rc->this_arf_ratio);
}

// calculate arf_ratio for this gf_group using *normalised* gf_group_stat and last arf_ratio
static AOM_INLINE void cal_arf_ratio(const GF_GROUP_STATS *stats, const int num_mbs, RATE_CONTROL *const rc) {
    fprintf(stderr, "\n calculating arf_ratio");
    double zero_motion_accumulator = stats->zero_motion_accumulator;
    double gf_group_raw_error_N    = stats->gf_group_raw_error / num_mbs;
    double avg_sr_coded_error_N    = stats->avg_sr_coded_error / num_mbs;
    double avg_raw_err_stdev_N     = stats->avg_raw_err_stdev  / num_mbs;
    double decay_accumulator       = stats->decay_accumulator;
    double avg_new_mv_count_N      = stats->avg_new_mv_count   / num_mbs;
    fprintf(stderr, "\n normalized gfstat: [6] %12.4lf [1] %12.4lf [12] %12.4lf [17] %12.4lf [5] %12.4lf [15] %12.4lf",
    zero_motion_accumulator, gf_group_raw_error_N, avg_sr_coded_error_N, avg_raw_err_stdev_N, decay_accumulator, avg_new_mv_count_N);
    
    double arf_ratio_delta = 0;
    if (zero_motion_accumulator <= 0.387) {
      if (gf_group_raw_error_N <= 65.431) {
        if (avg_raw_err_stdev_N > 17.649) arf_ratio_delta = -0.1;
      } else {
        if (decay_accumulator <= 0.19)    arf_ratio_delta = -0.1;
      }
    } else {
      if (avg_sr_coded_error_N <= 4.545) {
        if (avg_new_mv_count_N > 0.01)    arf_ratio_delta =  0.1; 
      } else {
        if (avg_raw_err_stdev_N <= 5.536) arf_ratio_delta =  0.1;
      }
    }

    rc->this_arf_ratio = fclamp(rc->last_arf_ratio + arf_ratio_delta, 0.0, 1.0);
    fprintf(stderr,"\n last_arf_ratio=%4.2lf this_arf_ratio=%4.2lf", rc->last_arf_ratio, rc->this_arf_ratio);
    rc->last_arf_ratio = rc->this_arf_ratio;
}

static AOM_INLINE void output_gf_gf_stats(const GF_GROUP_STATS *stats) {
  {
    fprintf(stderr, "\n GF_STATS: ");
    fprintf(stderr,
            " %12.4lf"
            " %12.4lf"
            " %12.4lf"
            " %12.4lf"

            " %12.4lf"
            " %12.4lf"
            " %12.4lf"
            " %12.4lf"
            " %12.4lf"
            " %12.4lf"
            " %12.4lf"
            " %12.4lf"

            " %12.4lf"
            " %12.4lf"
            " %12.4lf"
            " %12.4lf"
            " %12.4lf"
            " %12.4lf"
            " %12.4lf"
            " %12.4lf"

            " %12d",

            stats->gf_group_err,
            stats->gf_group_raw_error,
            stats->gf_group_skip_pct,
            stats->gf_group_inactive_zone_rows,

            stats->mv_ratio_accumulator,
            stats->decay_accumulator,
            stats->zero_motion_accumulator,
            stats->loop_decay_rate,
            stats->last_loop_decay_rate,
            stats->this_frame_mv_in_out,
            stats->mv_in_out_accumulator,
            stats->abs_mv_in_out_accumulator,

            stats->avg_sr_coded_error,
            stats->avg_tr_coded_error,
            stats->avg_pcnt_second_ref,
            stats->avg_pcnt_third_ref,
            stats->avg_pcnt_third_ref_nolast,
            stats->avg_new_mv_count,
            stats->avg_wavelet_energy,
            stats->avg_raw_err_stdev,

            stats->non_zero_stdev_count);
  }
}

typedef struct {
  double frame_err;
  double frame_coded_error;
  double frame_sr_coded_error;
  double frame_tr_coded_error;
} GF_FRAME_STATS;

void av1_init_second_pass(struct AV1_COMP *cpi);

void av1_init_single_pass_lap(AV1_COMP *cpi);

void av1_get_second_pass_params(struct AV1_COMP *cpi,
                                struct EncodeFrameParams *const frame_params,
                                const EncodeFrameInput *const frame_input,
                                unsigned int frame_flags);

void av1_twopass_postencode_update(struct AV1_COMP *cpi);

void av1_gop_bit_allocation(const AV1_COMP *cpi, RATE_CONTROL *const rc,
                            GF_GROUP *gf_group, int is_key_frame, int use_arf,
                            int64_t gf_group_bits);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_PASS2_STRATEGY_H_
