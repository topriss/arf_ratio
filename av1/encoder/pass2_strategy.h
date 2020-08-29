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
