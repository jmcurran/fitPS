jsd = function(p, q){
  ## Jensen-Shannon Divergence
  ## Not exported because I don't want to encourage misuse.

  p = p / sum(p)
  q = q / sum(q)
  m = 0.5 * (p + q)

  kl = function(p, q){
    klt = p * log(p / q, base = 2)
    klt[p == 0] = 0
    sum(klt)
  }

  return(0.5 * kl(p, m) + 0.5 * kl(q, m))
}
