# BRASS output classification NEW

library(gtools)
library(VariantAnnotation)
library(IRanges)
library(seqinr)
library(Biostrings)

# get chromosome sizes
WBcel235 <- readDNAStringSet("Caenorhabditis_elegans.WBcel235.dna.toplevel.fa") # get worm reference genome
chr_sizes <- width(WBcel235)
names(chr_sizes) <- c("I","II","III","IV","MtDNA","V","X")
genome_size = sum(as.numeric(chr_sizes))

# Helper functions
L = 1.5e9
p = function(dist) {
  if (dist < 0) {
    stop("Negative input for p()")
  }
  return(pmax(dist, 1)/L)
}
inv_p = function(p) {
  return(pmax(p*L, 1))
}
p_bkpts = function(bkpt_1, bkpt_2) {
  p(abs(bkpt_1 - bkpt_2))
}
bkpt_proximity_p_value = function(index_bkpt, chr_size, score_cutoff) {
  # Returns the fraction of the coordinates in the current chromosome
  # the produces a distance likelihood smaller than score_cutoff.
  
  distance_cutoff = inv_p(score_cutoff)
  min_coord = max(1, index_bkpt - distance_cutoff)
  max_coord = min(chr_size, index_bkpt + distance_cutoff)
  if (min_coord >= max_coord) {
    print(min_coord)
    print(max_coord)
    stop(sprintf("Failed with index_bkpt=%s; chr_size=%s; score_cutoff=%s", index_bkpt, chr_size, score_cutoff))
  } else {
    return((max_coord-min_coord)/chr_size)
  }
}
rg_paired_dist = function(chr_l_1, pos_l_1, chr_h_1, pos_h_1, chr_l_2, pos_l_2, chr_h_2, pos_h_2) {
  p_l = ifelse(
    chr_l_1 == chr_l_2,
    p_bkpts(pos_l_1, pos_l_2),
    1
  )
  p_h = ifelse(
    chr_h_1 == chr_h_2,
    p_bkpts(pos_h_1, pos_h_2),
    1
  )
  p_l * p_h
}
rg_dist = function(chr_l_1, pos_l_1, chr_h_1, pos_h_1, chr_l_2, pos_l_2, chr_h_2, pos_h_2) {
  pmin(
    rg_paired_dist(chr_l_1, pos_l_1, chr_h_1, pos_h_1, chr_l_2, pos_l_2, chr_h_2, pos_h_2),
    rg_paired_dist(chr_l_1, pos_l_1, chr_h_1, pos_h_1, chr_h_2, pos_h_2, chr_l_2, pos_l_2)
  )
}
quadratic_root = function(a, b, c) {
  D = b^2 - 4*a*c
  if (D < 0) {
    return(c())
  }
  else if (D == 0) {
    return(-b/2*a)
  }
  else {
    return(sort(c(-b + sqrt(D), -b - sqrt(D)) / (2 * a)))
  }
}

# 4 types of p-values: inter-inter, inter-intra, intra-inter, intra-intra
source("SV_technical/Inter_inter_fct.R")
source("SV_technical/Inter_intra_fct.R")
source("SV_technical/Intra_inter_fct.R")
source("SV_technical/Intra_intra_fct.R")
# Returns a matrix, where each entry is the expected number of FPs if the
# observed proximity P-value of each pair of SVs was used a the P-value threshold
# for the given pair. FP values are calculated by SV pair category, i.e.
#   del-del
#   del-inter
#   del-inv
#   del-td
#   inter-inter
#   inter-inv
#   inter-td
#   inv-inv
#   inv-td
#   td-td

# cluster matrix????

get_clust_mat_new = function(d,is_already_clustered = rep(F, nrow(d))) {
  n_inter = sum(d[!is_already_clustered, 1] != d[!is_already_clustered, 3])
  is_del = d[,6] == 'DEL'
  is_td  = d[,6] == 'DUP'
  is_inv = d[,6] == 'INV'
  rg_type = ifelse(is_del, "del",
                   ifelse(is_td,  "td",
                          ifelse(is_inv, "inv", "inter")))
  bkpt_l = d[,2]
  bkpt_h = d[,4]
  if (any((is_del | is_td) & bkpt_l >= bkpt_h)) {
    warning(sprintf("is_td & bkpt_l >= bkpt_h at lines c(%s)", paste(which((is_del | is_td) & bkpt_l >= bkpt_h), collapse = ", ")))
    idx = (is_del | is_td) & bkpt_l >= bkpt_h
    temp = bkpt_l[idx]
    bkpt_l[idx] = bkpt_h[idx]
    bkpt_h[idx] = temp
  }
  temp_m = pmin(bkpt_l, bkpt_h)
  temp_M = pmax(bkpt_l, bkpt_h)
  bkpt_l = ifelse(d[,1] == d[,3], temp_m, bkpt_l)
  bkpt_h = ifelse(d[,1] == d[,3], temp_M, bkpt_h)
  
  rg_sizes = list(
    "del" = setNames(bkpt_h-bkpt_l, d[,6])[is_del & !is_already_clustered], #???
    "td"  = setNames(bkpt_h-bkpt_l, d[,6])[is_td & !is_already_clustered], #???
    "inv" = setNames(bkpt_h-bkpt_l, d[,6])[is_inv & !is_already_clustered] #???
  )
  
  m = sapply(
    1:nrow(d),
    function(i) sapply(
      1:nrow(d),
      function(j) {
        if (j == 1) cat(paste("[", i, "] ", d[i,4], "\n", sep = "")) #???
        score_cutoff = rg_dist( d[i,1], bkpt_l[i], d[i,3], bkpt_h[i],
                                d[j,1], bkpt_l[j], d[j,3], bkpt_h[j])
        
        if (i == j) {
          out = NA
        }
        else if (rg_type[i] == "inter" && rg_type[j] == "inter") {
          temp_n = n_inter + sum(is_already_clustered[c(i,j)])
          out = temp_n * (temp_n-1) / 2
          if (score_cutoff < 1) {
            out = out * inter_vs_inter_p_value(d[i,1], bkpt_l[i], d[i,3], bkpt_h[i], score_cutoff, chr_sizes)
          }
        }
        else if (rg_type[i] != "inter" && rg_type[j] == "inter") {
          out = (length(rg_sizes[[rg_type[i]]]) + is_already_clustered[i]) * (n_inter + is_already_clustered[j])
          if (score_cutoff < 1) {
            out = out * inter_vs_intra_p_value(d[i,1], bkpt_l[i], bkpt_h[i], score_cutoff, chr_sizes)
          }
        }
        else if (rg_type[i] == "inter" && rg_type[j] != "inter") {
          out = (n_inter+is_already_clustered[i]) * (length(rg_sizes[[rg_type[j]]]) + is_already_clustered[j])
          cur_rg_sizes = rg_sizes[[rg_type[j]]]
          if (is_already_clustered[j]) {
            cur_rg_sizes = c(cur_rg_sizes, bkpt_h[j] - bkpt_l[j])
          }
          if (score_cutoff < 1) {
            out = out * intra_vs_inter_p_value(d[i,1], bkpt_l[i], d[i,3], bkpt_h[i], score_cutoff, chr_sizes, cur_rg_sizes)
          }
        }
        else if (rg_type[i] != "inter" && rg_type[j] != "inter") {
          if (rg_type[i] == rg_type[j]) {
            temp_n = length(rg_sizes[[rg_type[i]]]) + is_already_clustered[i] + is_already_clustered[j]
            out = temp_n * (temp_n-1) / 2
            cur_rg_sizes = rg_sizes[[rg_type[j]]]
            if (!is_already_clustered[i]) cur_rg_sizes = cur_rg_sizes[-(which(cur_rg_sizes == bkpt_h[i]-bkpt_l[i])[1])]
            if (is_already_clustered[j]) cur_rg_sizes = c(cur_rg_sizes, bkpt_h[j] - bkpt_l[j])
          } else {
            out = (length(rg_sizes[[rg_type[i]]])+is_already_clustered[i]) * (length(rg_sizes[[rg_type[j]]])+is_already_clustered[j])
            cur_rg_sizes = rg_sizes[[rg_type[j]]]
            if (is_already_clustered[j]) cur_rg_sizes = c(cur_rg_sizes, bkpt_h[j] - bkpt_l[j])
          }
          if (score_cutoff < 1) {
            out = out * intra_vs_intra_p_value(d[i,1], bkpt_l[i], bkpt_h[i], score_cutoff, chr_sizes[[d[i,1]]], cur_rg_sizes)
          }
        }
        else {
          stop()
        }
        
        out
      }
    )
  )
  
  m
   
}

get_clust_mat = function(d, is_already_clustered = rep(F, nrow(d))) {
  n_inter = sum(d[!is_already_clustered, 1] != d[!is_already_clustered, 5])
  is_del = d[,1] == d[,5] & d[,2] == "+" & d[,6] == "-"
  is_td  = d[,1] == d[,5] & d[,2] == "-" & d[,6] == "+"
  is_inv = d[,1] == d[,5] & d[,2] == d[,6]
  rg_type = ifelse(is_del, "del",
                   ifelse(is_td,  "td",
                          ifelse(is_inv, "inv", "inter")))
  
  # Fix coordinate columns and perform some sanity checks
  bkpt_l = d[,3]
  bkpt_h = d[,8]
  if (any((is_del | is_td) & bkpt_l >= bkpt_h)) {
    warning(sprintf("is_td & bkpt_l >= bkpt_h at lines c(%s)", paste(which((is_del | is_td) & bkpt_l >= bkpt_h), collapse = ", ")))
    idx = (is_del | is_td) & bkpt_l >= bkpt_h
    temp = bkpt_l[idx]
    bkpt_l[idx] = bkpt_h[idx]
    bkpt_h[idx] = temp
  }
  temp_m = pmin(bkpt_l, bkpt_h)
  temp_M = pmax(bkpt_l, bkpt_h)
  bkpt_l = ifelse(d[,1] == d[,5], temp_m, bkpt_l)
  bkpt_h = ifelse(d[,1] == d[,5], temp_M, bkpt_h)
  
  rg_sizes = list(
    "del" = setNames(bkpt_h-bkpt_l, d[,10])[is_del & !is_already_clustered], 
    "td"  = setNames(bkpt_h-bkpt_l, d[,10])[is_td & !is_already_clustered], 
    "inv" = setNames(bkpt_h-bkpt_l, d[,10])[is_inv & !is_already_clustered] 
  )
  
  m = sapply(
    1:nrow(d),
    function(i) sapply(
      1:nrow(d),
      function(j) {
        if (j == 1) cat(paste("[", i, "] ", d[i,7], "\n", sep = "")) 
        
        score_cutoff = rg_dist( d[i,1], bkpt_l[i], d[i,5], bkpt_h[i],
                                d[j,1], bkpt_l[j], d[j,5], bkpt_h[j])
        
        if (i == j) {
          out = NA
        }
        else if (rg_type[i] == "inter" && rg_type[j] == "inter") {
          temp_n = n_inter + sum(is_already_clustered[c(i,j)])
          out = temp_n * (temp_n-1) / 2
          if (score_cutoff < 1) {
            out = out * inter_vs_inter_p_value(d[i,1], bkpt_l[i], d[i,5], bkpt_h[i], score_cutoff, chr_sizes)
          }
        }
        else if (rg_type[i] != "inter" && rg_type[j] == "inter") {
          out = (length(rg_sizes[[rg_type[i]]]) + is_already_clustered[i]) * (n_inter + is_already_clustered[j])
          if (score_cutoff < 1) {
            out = out * inter_vs_intra_p_value(d[i,1], bkpt_l[i], bkpt_h[i], score_cutoff, chr_sizes)
          }
        }
        else if (rg_type[i] == "inter" && rg_type[j] != "inter") {
          out = (n_inter+is_already_clustered[i]) * (length(rg_sizes[[rg_type[j]]]) + is_already_clustered[j])
          cur_rg_sizes = rg_sizes[[rg_type[j]]]
          if (is_already_clustered[j]) {
            cur_rg_sizes = c(cur_rg_sizes, bkpt_h[j] - bkpt_l[j])
          }
          if (score_cutoff < 1) {
            out = out * intra_vs_inter_p_value(d[i,1], bkpt_l[i], d[i,5], bkpt_h[i], score_cutoff, chr_sizes, cur_rg_sizes)
          }
        }
        else if (rg_type[i] != "inter" && rg_type[j] != "inter") {
          if (rg_type[i] == rg_type[j]) {
            temp_n = length(rg_sizes[[rg_type[i]]]) + is_already_clustered[i] + is_already_clustered[j]
            out = temp_n * (temp_n-1) / 2
            cur_rg_sizes = rg_sizes[[rg_type[j]]]
            if (!is_already_clustered[i]) cur_rg_sizes = cur_rg_sizes[-(which(cur_rg_sizes == bkpt_h[i]-bkpt_l[i])[1])]
            if (is_already_clustered[j]) cur_rg_sizes = c(cur_rg_sizes, bkpt_h[j] - bkpt_l[j])
          } else {
            out = (length(rg_sizes[[rg_type[i]]])+is_already_clustered[i]) * (length(rg_sizes[[rg_type[j]]])+is_already_clustered[j])
            cur_rg_sizes = rg_sizes[[rg_type[j]]]
            if (is_already_clustered[j]) cur_rg_sizes = c(cur_rg_sizes, bkpt_h[j] - bkpt_l[j])
          }
          if (score_cutoff < 1) {
            out = out * intra_vs_intra_p_value(d[i,1], bkpt_l[i], bkpt_h[i], score_cutoff, chr_sizes[[d[i,1]]], cur_rg_sizes)
          }
        }
        else {
          stop()
        }
        
        out
      }
    )
  )
  
  m
}

get_footprints = function(pos, chr, cutoff = 0.001) {
  # Compute footprints given all breakpoints of a cluster
  
  # Edge case
  if (length(pos) == 2 || length(pos) == length(unique(chr))) {
    return(list(
      footprint_idx = 1:length(pos),
      footprint_coords = paste(pos, pos, sep = "-"),
      footprint_bounds = paste(pos, pos, sep = "-")
    ))
  }
  
  # Calculate distances between breakpoints
  o = order(pos, chr)
  reverse_o = match(1:length(pos), o)
  pos = pos[o]
  chr = chr[o]
  dists = diff(pos)
  dists[chr[-length(chr)] != chr[-1]] = NA
  is_new_footprint = rep(0, length(dists) + 1)
  is_new_footprint[which(is.na(dists)) + 1] = 1
  
  # Just a sanity check...
  if (any(pos != sort(pos))) {
    stop("Input positions are not sorted")
  }

  # Below P-values are computed using likelihood ratio test instead
  while (sum(!is.na(dists)) > 2) {
    is_max = which.max(dists)
    p0 = 1/mean(dists, na.rm = T)
    p1 = 1/mean(dists[-is_max], na.rm = T)
    l0 = sum(dexp(dists, p0, log = T), na.rm = T)
    l1 = sum(dexp(dists[-is_max], p1, log = T), na.rm = T) + dexp(dists[is_max], 1/dists[is_max], log = T)
    D = 2 * (l1 - l0)
    pval = pchisq(D, 1, lower.tail = F) * sum(!is.na(dists))
    if (pval < cutoff) {
      is_new_footprint[is_max+1] = 1
      dists[is_max] = NA
    } else {
      break
    }
  }
  footprint_idx = rep(1, length(pos)) + cumsum(is_new_footprint)
  lambda = 1/mean(dists, na.rm = T)
  
  # Chromosome of each footprint
  chr_of_footprint_idx = c()
  for (k in 1:max(footprint_idx)) {
    temp = unique(chr[footprint_idx == k])
    if (length(temp) != 1) { stop("Sanity check failed") }
    chr_of_footprint_idx[k] = temp
  }
  
  footprint_coords_of_idx = c()
  footprint_bounds_of_idx = c()
  
  footprint_idx = footprint_idx[reverse_o]
  return(list(
    footprint_idx = footprint_idx,
    footprint_coords = footprint_coords_of_idx[footprint_idx],
    footprint_bounds = footprint_bounds_of_idx[footprint_idx]
  ))
}

# what does this function do?
clust_rgs = function(d) {
  # set cutoffs for FDR
  cutoff_1 = 0.001
  cutoff_2 = 0.01
  
  # cluster the table
  m = get_clust_mat(d)
  m2 = m/2 + t(m)/2
  h = hclust(as.dist(m2), method = "single")
  fdr = h$height * 10 / seq_along(h$height)
  fdr_cutoff = cutoff_1
  if (fdr[1] > fdr_cutoff) {
    ct = 1:nrow(d)
  } else if (all(fdr <= fdr_cutoff)) {
    ct = rep(1, nrow(d))
  } else {
    hc = h$height[which.min(fdr <= fdr_cutoff) - 1]
    ct = cutree(h, h = hc)
  }
  is_already_clustered = duplicated(ct) | duplicated(ct, fromLast = T)
  
  m = get_clust_mat(d, is_already_clustered)
  m2 = m/2 + t(m)/2
  h = hclust(as.dist(m2), method = "single")
  fdr = h$height * 10 / seq_along(h$height)
  fdr_cutoff = cutoff_2
  if (fdr[1] > fdr_cutoff) {
    ct = 1:nrow(d)
  } else if (all(fdr <= fdr_cutoff)) {
    ct = rep(1, nrow(d))
  } else {
    hc = h$height[which.min(fdr <= fdr_cutoff) - 1]
    ct = cutree(h, h = hc)
  }
  
  cutree_res = ct
  list(
    hclust = h,
    fdr = fdr,
    cutree = ct,
    m = m
  )
}

clust_rgs_new = function(d) {
  # set cutoffs for FDR
  cutoff_1 = 0.001
  cutoff_2 = 0.01
  
  # cluster the table
  m = get_clust_mat_new(d)
  m2 = m/2 + t(m)/2
  h = hclust(as.dist(m2), method = "single")
  fdr = h$height * 10 / seq_along(h$height)
  fdr_cutoff = cutoff_1
  if (fdr[1] > fdr_cutoff) {
    ct = 1:nrow(d)
  } else if (all(fdr <= fdr_cutoff)) {
    ct = rep(1, nrow(d))
  } else {
    hc = h$height[which.min(fdr <= fdr_cutoff) - 1]
    ct = cutree(h, h = hc)
  }
  is_already_clustered = duplicated(ct) | duplicated(ct, fromLast = T)
  
  m = get_clust_mat_new(d, is_already_clustered)
  m2 = m/2 + t(m)/2
  h = hclust(as.dist(m2), method = "single")
  fdr = h$height * 10 / seq_along(h$height)
  fdr_cutoff = cutoff_2
  if (fdr[1] > fdr_cutoff) {
    ct = 1:nrow(d)
  } else if (all(fdr <= fdr_cutoff)) {
    ct = rep(1, nrow(d))
  } else {
    hc = h$height[which.min(fdr <= fdr_cutoff) - 1]
    ct = cutree(h, h = hc)
  }
  
  cutree_res = ct
  list(
    hclust = h,
    fdr = fdr,
    cutree = ct,
    m = m
  )
}

#########################################################################
