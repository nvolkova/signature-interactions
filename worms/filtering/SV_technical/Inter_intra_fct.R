# Inter-intra p-value
inter_vs_intra_p_value = function(index_chr, index_bkpt_l, index_bkpt_h, score_cutoff, chr_sizes) {
  # Index rearrangement is intra-chromosomal, other rearrangement is inter-chromosomal
  chr_size = chr_sizes[[index_chr]]
  
  max_dist = inv_p(score_cutoff)
  min_pos_l = pmax(1, index_bkpt_l - max_dist)
  max_pos_l = pmin(chr_size, index_bkpt_l + max_dist)
  min_pos_h = pmax(1, index_bkpt_h - max_dist)
  max_pos_h = pmin(chr_size, index_bkpt_h + max_dist)
  
  prob_one_end_in_index_chr = 2 * (chr_size/genome_size)
  if (min_pos_l <= max_pos_h && min_pos_h <= max_pos_l) {
    out = prob_one_end_in_index_chr * (max(max_pos_l, max_pos_h) - min(min_pos_l, min_pos_h) + 1) / chr_size
  }
  else {
    out = prob_one_end_in_index_chr * (max_pos_l-min_pos_l+1 + max_pos_h-min_pos_h+1) / chr_size
  }
  
  if (out > 1 || out <= 0) {
    stop(sprintf("Failure at inter_vs_intra_p_value() with index_chr=%s, index_bkpt_l=%s, index_bkpt_h=%s, score_cutoff=%s, chr_sizes=%s, OUT=%s", index_chr, index_bkpt_l, index_bkpt_h, score_cutoff, "chr_sizes", out))
  }
  return(out)
}