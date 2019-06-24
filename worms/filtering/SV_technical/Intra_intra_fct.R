# intra-intra

# This is a helper function for a single secondary breakpoint
intra_vs_single_intra_p_value = function(il, ih, score_cutoff, s, chr_size) {
  il = as.numeric(il)
  ih = as.numeric(ih)
  score_cutoff = as.numeric(score_cutoff)
  s = pmax(as.numeric(s), 1)
  chr_size = as.numeric(chr_size)
  
  # Is there any chance this size produces good enough score?
  if (p(0) * p(abs(ih-il-s)) == score_cutoff) {
    return(1/(chr_size - s + 1))
  }
  if (p(0) * p(abs(ih-il-s)) > score_cutoff) {
    return(0)
  }
  
  good_regions_start = numeric()
  good_regions_end = numeric()
  
  # 1. ll and hh pairing
  if (il <= ih - s) {  # This is if low ends meet first
    # 1) pos <= il. Solution for (il-pos) * (ih-pos-size) - score_cutoff*L^2 = 0
    roots = quadratic_root(1, -il-ih+s, -il*s + il*ih - score_cutoff*L^2)
    roots = round(roots[roots-1 <= il])  # -1 for some rounding problem
    if (length(roots) == 0) {
      # stop("At ll, il <= ih - s, 1): no roots found???")
      warning("At ll, il <= ih - s, 1): no roots found???")
      roots = il
    }
    else if (length(roots) == 2) {
      stop("At ll, il <= ih - s, 1): found two roots???")
    }
    else {
      good_regions_start = c(good_regions_start, pmax(roots, 1))
      good_regions_end   = c(good_regions_end,   max(c(il, roots, 1)))
    }
    
    # 2) il < pos <= ih-s. Solution for (pos-il) * (ih-pos-size) - score_cutoff*L^2 = 0
    roots = quadratic_root(-1, ih-s+il, -il*ih + il*s - score_cutoff*L^2)
    if (length(roots) <= 1) {
      good_regions_start = c(good_regions_start, il)
      good_regions_end   = c(good_regions_end,   ih-s)
    }
    else {
      roots = round(roots)
      good_regions_start = c(good_regions_start, il,       roots[2])
      good_regions_end   = c(good_regions_end,   roots[1], ih-s)
    }
    
    # 3) ih-s < pos. Solution for (pos-il) * (pos+size-ih) - score_cutoff*L^2 = 0
    roots = quadratic_root(1, -ih+s-il, il*ih - il*s - score_cutoff*L^2)
    roots = round(roots[roots+1 >= ih-s])  # +1 for some rounding problems
    if (length(roots) == 0) {
      # stop("At hh, il <= ih - s, 3): no roots found???")
      warning("At hh, il <= ih - s, 3): no roots found???")
      roots = ih - s
    }
    else if (length(roots) == 2) {
      stop("At hh, il <= ih - s, 3): found two roots???")
    }
    else {
      good_regions_start = c(good_regions_start, min(c(ih-s, roots, chr_size)))
      good_regions_end   = c(good_regions_end,   pmin(roots, chr_size))
    }
  }
  else {  # This is if high ends meet first
    # 1) pos+size <= ih. Solution for (il-pos) * (ih-pos-size) - score_cutoff*L^2 = 0
    roots = quadratic_root(1, -il-ih+s, -il*s + il*ih - score_cutoff*L^2)
    roots = round(roots[roots-1 <= ih-s])  # -1 for some rounding problem
    if (length(roots) == 0) {
      # stop("At ll, il > ih - s, 1): no roots found???")
      warning("At ll, il > ih - s, 1): no roots found???")
      roots = ih - s
    }
    else if (length(roots) == 2) {
      stop("At ll, il > ih - s, 1): found two roots???")
    }
    else {
      good_regions_start = c(good_regions_start, pmax(1, roots))
      good_regions_end   = c(good_regions_end,   max(c(ih-s, 1, roots)))
    }
    
    # 2) ih-size < pos <= il. Solution for (il-pos) * (pos+size-ih) - score_cutoff*L^2 = 0
    roots = quadratic_root(-1, ih-s+il, -il*ih + il*s - score_cutoff*L^2)
    if (length(roots) <= 1) {
      good_regions_start = c(good_regions_start, ih-s)
      good_regions_end   = c(good_regions_end,   il)
    }
    else {
      round(roots)
      good_regions_start = c(good_regions_start, ih-s,   roots[2])
      good_regions_end   = c(good_regions_end,   roots[1], il)
    }
    
    # 3) il < pos. Solution for (pos-il) * (pos+size-ih) - score_cutoff*L^2 = 0
    roots = quadratic_root(1, -ih+s-il, il*ih - il*s - score_cutoff*L^2)
    roots = round(roots[roots+1 >= il])  # +1 for some rounding problem
    if (length(roots) == 0) {
      # stop("At hh, il > ih - s, 3): no roots found???")
      warning("At hh, il > ih - s, 3): no roots found???")
      roots = il
    }
    else if (length(roots) == 2) {
      stop("At hh, il > ih - s, 3): found two roots???")
    }
    else {
      good_regions_start = c(good_regions_start, min(c(il+1, roots, chr_size)))
      good_regions_end   = c(good_regions_end,   pmin(roots, chr_size))
    }
  }
  
  # Now lh and hl pairing. Index low end will always meet with second rg high end first.
  # Need to entere here only if there's any chance some of the achieved scores
  # are more extreme than score_cutoff. 
  if (p(0) * p(abs(ih-il+s)) <= score_cutoff) {
    # 1) pos+size <= il: Solution for (il-pos-size) * (ih-pos) - score_cutoff*L^2 = 0
    roots = quadratic_root(1, -ih+s-il, il*ih - ih*s - score_cutoff*L^2)
    roots = round(roots[roots-1 <= il-s])  # -1 for some rounding problem
    if (length(roots) == 0) {
      # stop("At lh, 1): no roots found???")
      warning("At lh, 1): no roots found???")
      roots = il - s
    }
    else if (length(roots) == 2) {
      stop("At lh, 1): found two roots???")
    }
    else {
      good_regions_start = c(good_regions_start, pmax(roots, 1))
      good_regions_end   = c(good_regions_end,   max(c(il-s, roots, 1)))
    }
    
    # 2) pos <= ih: Solution for (pos+size-il) * (ih-pos) - score_cutoff*L^2 = 0
    roots = quadratic_root(-1, ih-s+il, -il*ih + ih*s - score_cutoff*L^2)
    roots = roots[roots+s+1 >= il & roots-1 <= ih]  # +1 and -1 for som erounding problems
    if (length(roots) <= 1) {
      good_regions_start = c(good_regions_start, pmin(il-s+1, ih))
      good_regions_end   = c(good_regions_end,   ih)
    }
    else {
      roots = round(roots)
      good_regions_start = c(good_regions_start, il-s+1,   roots[2])
      good_regions_end   = c(good_regions_end,   roots[1], ih)
    }
    
    # 3) pos > ih: Solution for (pos+size-il) * (pos-ih) - score_cutoff*L^2 = 0
    roots = quadratic_root(1, -ih+s-il, il*ih - ih*s - score_cutoff*L^2)
    roots = round(roots[roots+1 >= ih])  # +1 for some rounding problems
    if (length(roots) == 0) {
      # stop("At hl, 3): no roots found???")
      warning("At hl, 3): no roots found???")
      roots = ih
    }
    else if (length(roots) == 2) {
      stop("At hl, 3): found two roots???")
    }
    else {
      good_regions_start = c(good_regions_start, min(c(roots, ih+1, chr_size)))
      good_regions_end   = c(good_regions_end,   pmin(roots, chr_size))
    }
  }
  
  # Remove regions that are out of bounds
  bad_idx = good_regions_start > chr_size | good_regions_end < 1
  good_regions_start = good_regions_start[!bad_idx]
  good_regions_end   = good_regions_end[!bad_idx]
  if (any(good_regions_end < good_regions_start)) {
    print(good_regions_start)
    print(good_regions_end)
    print(which(good_regions_end < good_regions_start))
    stop(sprintf("Failed after removing out of bounds regions: il=%s; ih=%s; score_cutoff=%s; s=%s; chr_size=%s", il, ih, score_cutoff, s, chr_size))
  }
  
  # Merge regions
  o = order(good_regions_start, good_regions_end)
  good_regions_start = good_regions_start[o]
  good_regions_end   = good_regions_end[o]
  i = 1
  while(i < length(good_regions_start)) {
    if (
      good_regions_start[i]   <= good_regions_end[i+1] &&
      good_regions_start[i+1] <= good_regions_end[i]
    ) {
      good_regions_end[i] = good_regions_end[i+1]
      good_regions_start = good_regions_start[-(i+1)]
      good_regions_end = good_regions_end[-(i+1)]
    }
    else {
      i = i + 1
    }
  }
  
  out = sum(pmin(good_regions_end, chr_size) - pmax(1, good_regions_start) + 1) / chr_size
  if (out > 1 || out < 0) {
    print(good_regions_start)
    print(good_regions_end)
    print(good_regions_end - good_regions_start + 1)
    print(sum(good_regions_end - good_regions_start + 1))
    print(chr_size)
    print(sum(good_regions_end - good_regions_start + 1) > chr_size)
    stop(sprintf("Failed after summing regions: out=%s; il=%s; ih=%s; score_cutoff=%s; s=%s; chr_size=%s", out, il, ih, score_cutoff, s, chr_size))
  }
  
  return(out)
}

intra_vs_intra_p_value = function(index_chr, index_bkpt_l, index_bkpt_h, score_cutoff, chr_size, rg_sizes) {
  if (all(rg_sizes > chr_size)) {
    # return(0)
    rg_sizes = chr_size
  }
  
  res = sapply(
    rg_sizes[rg_sizes <= chr_size],
    function(size) {
      # print(size)
      intra_vs_single_intra_p_value(index_bkpt_l, index_bkpt_h, score_cutoff, size, chr_size)
    }
  )
  
  out = (chr_size/genome_size) * mean(res)
  if (out > 1 || out < 0) {
    stop()
  }
  return(out)
}