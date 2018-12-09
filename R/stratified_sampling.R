#' @title stratified sampling
#' @author Jacques.van-Helden\@univ-amu.fr
#' @description Perform a stratified sampling from a vector of group labels:
#' each initial group is replaced by random labels proportional to the initial
#' group proportions.
#'
#' @param labels a vector with labels to be resampled
#' @param size=length(labels) size of the output (resampled) vector
#'
#' @examples
#'
#' ## Generate a vector of group labels
#' labels <- rep(x = c("A", "B", "C", "D"), times = c(256, 64, 16, 4))
#'
#' ## Non-stratified sampling
#' nstr <- sample(labels, size = length(labels) , replace = FALSE)
#' table(labels, nstr)
#' ## Repeat this several times: the proportions vary between iterations
#'
#' ## Stratified permutation of group labels
#' str1 <- StratifiedSampling(labels)
#' ## str1 has same proportions as labels (except for rounding of group sizes)
#' table(labels, str1)
#'
#' ## Second stratified permutation
#' str2 <- StratifiedSampling(labels)
#' ## str2 has same proportions as labels (except for rounding of group sizes)
#' table(labels, str2)
#'
#' ## str1 and str2 are different and independent
#' table(str1, str2)
#'
#' ## Subsampling
#' str3 <- StratifiedSampling(labels, size = length(labels) * 2/3)
#' ## str3 has same proportions as labels (except for rounding of group sizes)
#' table(str3)
#'
#' @export
StratifiedSampling <- function(labels,
                               size = length(labels)) {


  group.sizes.table <- table(labels)
  groups <- names(group.sizes.table)
  group.sizes <- as.vector(group.sizes.table)
  names(group.sizes) <- groups
  group.proportions <- group.sizes / sum(group.sizes)
  resampled <- vector()
  for (group in groups) {
    group.size <- group.proportions[group] * size
    target.sizes <- as.vector(round(group.proportions * group.size))
    group.sample <- sample(rep(x = groups, times = target.sizes), replace = FALSE)

    ## Fix possible differences due to rounding of group sizes
    diff <- group.size - length(group.sample)
    if (diff > 0) {
      ## Add random labels to the too small group
      group.sample <- append(group.sample, sample(x = labels, size = diff, replace = FALSE))
    } else if (diff < 0) {
      ## Select a subset of the too big group
      group.sample <- sample(group.sample, size = group.size, replace = FALSE)
    }
    resampled <- append(resampled, group.sample)
  }


  return(resampled)
}
