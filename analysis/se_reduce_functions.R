# functions to cluster se based on percentage of overlaps
extractClustersFromSelfHits <- function(hits) {
  #hits <- union(hits, t(hits))
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  cid <- seq_len(queryLength(hits))  # cluster ids
  while (TRUE) {
    h <- Hits(qh, cid[sh],
              queryLength(hits), subjectLength(hits))
    cid2 <- pmin(cid, selectHits(h, "first"))
    if (identical(cid2, cid))
      break
    cid <- cid2
  }
  unname(splitAsList(seq_len(queryLength(hits)), cid))
}

mergeConnectedRanges <- function(x, hits) {
  clusters <- extractClustersFromSelfHits(hits)
  ans <- range(extractList(x, clusters))
  ans <- unlist(ans)
  mcols(ans)$cluster <- clusters
  ans
}

group_se_by_overlap <- function(gr,overlap_p=0.5) {
  
  # selfhits
  hits <- findOverlaps(gr)
  x <- gr[queryHits(hits)]
  y <- gr[subjectHits(hits)]
  
  # overlap percentage
  relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y))
  hits_out <- hits[relative_overlap >= overlap_p]
  
  gr_out <- mergeConnectedRanges(gr, hits_out)
  gr_out
  
}

