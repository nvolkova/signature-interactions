# Inter-inter
inter_vs_inter_p_value = function(index_chr_l, index_bkpt_l, index_chr_h, index_bkpt_h, score_cutoff, chr_sizes) {
  p_match_at_low_chr   = 2 * (chr_sizes[[index_chr_l]]/genome_size) * (1 - chr_sizes[[index_chr_h]]/(genome_size-chr_sizes[[index_chr_l]])) * bkpt_proximity_p_value(index_bkpt_l, chr_sizes[[index_chr_l]], score_cutoff)
  p_match_at_high_chr  = 2 * (chr_sizes[[index_chr_h]]/genome_size) * (1 - chr_sizes[[index_chr_l]]/(genome_size-chr_sizes[[index_chr_h]])) * bkpt_proximity_p_value(index_bkpt_h, chr_sizes[[index_chr_h]], score_cutoff)
  
  # Calculating the area where two inter-chromosomal events produce a p smaller than threshold. 
  # 1: d1 > 0, d2 > 0. Calculate where d1 is when d2 = chr_len - index_bkpt_h
  max_d2 = chr_sizes[[index_chr_h]] - index_bkpt_h
  a = score_cutoff / max_d2 * L^2
  b = chr_sizes[[index_chr_l]] - index_bkpt_l
  area_sum = pmin(b, a) * max_d2
  if (a < b) { 
    area_sum = area_sum + score_cutoff*L*L*(log(b)-log(a))
  }
  
  # 2: d1 > 0, d2 < 0
  max_d2 = index_bkpt_h - 1
  a = score_cutoff / max_d2 * L^2
  b = chr_sizes[[index_chr_l]] - index_bkpt_l
  area_sum = area_sum + pmin(b, a) * max_d2
  if (a < b) {
    area_sum = area_sum + score_cutoff*L*L*(log(b)-log(a))
  }
  
  # 3: d1 < 0, d2 > 0
  max_d2 = chr_sizes[[index_chr_h]] - index_bkpt_h
  a = score_cutoff / max_d2 * L^2
  b = index_bkpt_l - 1
  area_sum = area_sum + pmin(a, b) * max_d2
  if (a < b) { 
    area_sum = area_sum + score_cutoff*L*L*(log(b)-log(a))
  }
  
  # 2: d1 < 0, d2 < 0
  max_d2 = index_bkpt_h - 1
  a = score_cutoff / max_d2 * L^2
  b = index_bkpt_l - 1
  area_sum = area_sum + pmin(a, b) * max_d2
  if (a < b) {
    area_sum = area_sum + score_cutoff*L*L*(log(b)-log(a))
  }
  
  p_match_at_both_chrs = 2 * (chr_sizes[[index_chr_l]]/genome_size) * (chr_sizes[[index_chr_h]]/(genome_size-chr_sizes[[index_chr_l]])) * (area_sum/chr_sizes[[index_chr_l]]/chr_sizes[[index_chr_h]])
  
  out = p_match_at_low_chr + p_match_at_high_chr + p_match_at_both_chrs
  
  if (out > 1 || out <= 0) {
    stop(sprintf("Failed at inter_vs_inter_p_value() with index_chr_l=%s index_bkpt_l=%s index_chr_h=%s index_bkpt_h=%s score_cutoff=%s chr_sizes=c(%s)", index_chr_l, index_bkpt_l, index_chr_h, index_bkpt_h, score_cutoff, paste(chr_sizes, collapse = ",")))
  }
  return(out)
}