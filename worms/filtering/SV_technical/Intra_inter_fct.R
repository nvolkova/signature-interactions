# Intra-inter

# This is helper function for intra_vs_inter_p_value()
intra_vs_single_bkpt_p_value = function (index_bkpt, rg_sizes, chr_size, score_cutoff) {
  # Range of positions giving sufficient proximity
  max_dist = inv_p(score_cutoff)
  
  if (all(rg_sizes > chr_size)) {
    rg_sizes = chr_size
  }
  
  out = mean(sapply(
    rg_sizes[rg_sizes <= chr_size],
    function(size) {
      min_pos_l = pmin(pmax(1, index_bkpt - max_dist), chr_size - size)
      max_pos_l = pmin(chr_size - size, index_bkpt + max_dist)
      min_pos_h = pmax(1, index_bkpt - size - max_dist)
      max_pos_h = pmax(1, pmin(chr_size - size, index_bkpt - size + max_dist))
      # cat(paste(size, min_pos_l, max_pos_l, min_pos_h, max_pos_h, "\n"))
      
      # Sanity check
      if (min_pos_l > max_pos_l || min_pos_h > max_pos_h) {
        stop()
      }
      
      if (
        min_pos_l <= max_pos_h &&  # Low end range start vs. high end range end
        max_pos_l >= min_pos_h     # Low end range end vs. high end range start
      ) {
        (pmax(max_pos_l, max_pos_h) - pmin(min_pos_l, min_pos_l) + 1)/chr_size
      }
      else {
        (max_pos_l-min_pos_l+1 + max_pos_h-min_pos_h+1) / chr_size
      }
    }
  ))
  
  if (out > 1 || out < 0) {
    stop(sprintf("Failed with index_bkpt=%s; chr_size=%s; score_cutoff=%s; rg_sizes=%s", index_bkpt, chr_size, score_cutoff, paste(rg_sizes, collapse = ",")))
  }
  out
}

intra_vs_inter_p_value = function(index_chr_l, index_bkpt_l, index_chr_h, index_bkpt_h, score_cutoff, chr_sizes, rg_sizes) {
  # Index rearrangement is inter-chromosomal, other rearrangement is intra-chromosomal
  out = chr_sizes[[index_chr_l]]/genome_size * intra_vs_single_bkpt_p_value(index_bkpt_l, rg_sizes, chr_sizes[[index_chr_l]], score_cutoff) +
    chr_sizes[[index_chr_h]]/genome_size * intra_vs_single_bkpt_p_value(index_bkpt_h, rg_sizes, chr_sizes[[index_chr_h]], score_cutoff)
  
  if (out > 1 || out < 0) {
    stop()
  }
  out
}