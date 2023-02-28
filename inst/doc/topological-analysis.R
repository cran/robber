## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  fig.width = 6,
  fig.height = 4
)

## ----setup--------------------------------------------------------------------
library(robber)
library(ggplot2)
library(tidyr)
library(dplyr)
library(forcats)
library(purrr)
library(igraph)
theme_set(theme_bw(base_size = 15))

## -----------------------------------------------------------------------------
data("hostparasite", package = "robber")
data("seeddispersal", package = "robber")
data("pollination", package = "robber")

## -----------------------------------------------------------------------------
print(paste0("Host-Parasite: ", nrow(hostparasite), " ", 
           ncol(hostparasite), " ", mean(hostparasite)), sep = ",")
print(paste0("Seed Dispersal: ", nrow(seeddispersal), " ", 
           ncol(seeddispersal), " ", mean(seeddispersal)), sep = ",")
print(paste0("Pollination: ", nrow(pollination), " ", 
           ncol(pollination), " ", mean(pollination)), sep = ",")

## -----------------------------------------------------------------------------
(rob_emp_unif <- robustness_emp(hostparasite))
plot(rob_emp_unif)

## -----------------------------------------------------------------------------
(lbm_param <- get_lbm_param(hostparasite, ncores = 1L))

## -----------------------------------------------------------------------------
(rob_lbm_unif <- do.call(robustness_lbm, lbm_param))

## -----------------------------------------------------------------------------
rob_emp_unif$auc
rob_lbm_unif$auc

## -----------------------------------------------------------------------------
plot(rob_emp_unif) + 
  plot(rob_lbm_unif, add = TRUE)

## -----------------------------------------------------------------------------
rob_fun <- lapply(
  X = seq_len(10),
  FUN = function(i) {
    A <- simulate_lbm(con = lbm_param$con, 
                      pi = lbm_param$pi, 
                      rho = lbm_param$rho, 
                      nr = lbm_param$nr, 
                      nc = lbm_param$nc)$A
    return(list(unif = robustness_emp(A),
           dec = robustness_emp(A, ext_seq = "decreasing"),
           inc = robustness_emp(A, ext_seq = "increasing")))
  }
)

## -----------------------------------------------------------------------------
pu <- plot(rob_emp_unif) + 
  plot(rob_lbm_unif, add = TRUE) +
  plot(rob_fun[[1]]$unif, add = TRUE, lty = "dotted", col = "red") +
  plot(rob_fun[[2]]$unif, add = TRUE, lty = "dotted", col = "red") +
  plot(rob_fun[[3]]$unif, add = TRUE, lty = "dotted", col = "red") +
  plot(rob_fun[[4]]$unif, add = TRUE, lty = "dotted", col = "red") +
  plot(rob_fun[[5]]$unif, add = TRUE, lty = "dotted", col = "red") +
  ggtitle("Uniform")
pu

## -----------------------------------------------------------------------------
rob_emp_dec <- robustness_emp(hostparasite, ext_seq = "decreasing")
rob_emp_inc <- robustness_emp(hostparasite, ext_seq = "increasing")
rob_emp_dec$auc
rob_emp_inc$auc

## -----------------------------------------------------------------------------
rob_emp_dec_lin <- robustness_emp(hostparasite, ext_seq = "decreasing", 
                                  method = "linear")
rob_emp_inc_lin <- robustness_emp(hostparasite, ext_seq = "increasing", 
                                  method = "linear")
rob_emp_dec_lin$auc
rob_emp_inc_lin$auc

## -----------------------------------------------------------------------------
rob_lbm_dec <- do.call(robustness_lbm, c(lbm_param, ext_seq = "decreasing"))
rob_lbm_inc <- do.call(robustness_lbm, c(lbm_param, ext_seq = "increasing"))
rob_lbm_dec$auc
rob_lbm_inc$auc

## -----------------------------------------------------------------------------
plot(rob_lbm_unif, col = "black", lty = 1) + 
  plot(rob_emp_unif, add = TRUE, col = "blue", lty = 2) +
  plot(rob_emp_dec, add = TRUE, lty = 3, col = "blue") +
  plot(rob_emp_dec_lin, add = TRUE, lty = 4, col = "blue") +
  plot(rob_emp_inc, add = TRUE, lty = 3, col = "blue") +
  plot(rob_emp_inc_lin, add = TRUE, lty = 4, col = "blue") +
  plot(rob_lbm_dec, add = TRUE, lty = 1, col = "black") +
  plot(rob_lbm_inc, add = TRUE, lty = 1, col = "black") +
  ggtitle("Robustness of hostparasite")

## -----------------------------------------------------------------------------
pl_param <- get_lbm_param(pollination, ncores = 1L)
sd_param <- get_lbm_param(seeddispersal, ncores = 1L)

## -----------------------------------------------------------------------------
pl_param

## -----------------------------------------------------------------------------
sd_param

## -----------------------------------------------------------------------------
pl_lbm_unif <- do.call(robustness_lbm, pl_param)
pl_lbm_dec <- do.call(robustness_lbm, c(pl_param, ext_seq = "decreasing"))
pl_lbm_inc <- do.call(robustness_lbm, c(pl_param, ext_seq = "increasing"))
sd_lbm_unif <- do.call(robustness_lbm, sd_param)
sd_lbm_dec <- do.call(robustness_lbm, c(sd_param, ext_seq = "decreasing"))
sd_lbm_inc <- do.call(robustness_lbm, c(sd_param, ext_seq = "increasing"))

## ----echo=FALSE---------------------------------------------------------------
df <- data.frame(
  Extinction = c("Uniform", "Increasing", "Decreasing"),
  hostparasite = c(rob_lbm_unif$auc, rob_lbm_inc$auc, rob_lbm_dec$auc),
  pollination = c(pl_lbm_unif$auc, pl_lbm_inc$auc, pl_lbm_dec$auc),
  seeddispersal = c(sd_lbm_unif$auc, sd_lbm_inc$auc, sd_lbm_dec$auc))
knitr::kable(df)

## -----------------------------------------------------------------------------
comp_unif <- compare_robustness(list(lbm_param, pl_param, sd_param), 
                   dens = .148, new_nr = 36, new_nc = 51)
comp_inc <- compare_robustness(list(lbm_param, pl_param, sd_param), 
                   dens = .148, new_nr = 36, new_nc = 51, ext_seq = "increasing")
comp_dec <- compare_robustness(list(lbm_param, pl_param, sd_param), 
                   dens = .148, new_nr = 36, new_nc = 51, ext_seq = "decreasing")

## ----echo = FALSE-------------------------------------------------------------
dfc <- data.frame(
  Extinction = c("Uniform", "Increasing", "Decreasing"),
  hostparasite = c(comp_unif[[1]], comp_inc[[1]], comp_dec[[1]]),
  pollination = c(comp_unif[[2]], comp_inc[[2]], comp_dec[[2]]),
  seeddispersal = c(comp_unif[[3]], comp_inc[[3]], comp_dec[[3]]))
knitr::kable(dfc)

## ----set_constant-------------------------------------------------------------
dens <- .0156
nr <-  100
nc <- 100
rob_er <- auc_robustness_lbm(matrix(dens, 1,1), 1, 1, nr, nc)
rob_er

## ----echo=FALSE---------------------------------------------------------------
mod <- matrix(c("a", 1, 1, "a"), 2, 2)
cp  <- matrix(c("a", "a", "a", 1), 2, 2)

## ----echo=FALSE---------------------------------------------------------------
     knitr::kable(mod)

## ----echo=FALSE---------------------------------------------------------------
     knitr::kable(cp)

## ----parameters---------------------------------------------------------------
pi <- c(1/4, 3/4)

## ----analyse_topo, echo = FALSE-----------------------------------------------
if (! file.exists("res_topo.rds")) {
  robust_topology <- tibble::tibble()
  eps <- c(1/seq(8, 1.5, by = -.5), seq(1, 8, by = .5))
  for(i in seq(19)) {
    rho <- c(i*.05, (20-i)*.05)
    list_con_mod <- lapply( eps,
                            function(j) {
                              list(con = matrix(c(j, 1,
                                                  1, j), 2, 2),
                                   pi = pi,
                                   rho = rho)})
    rob_mod <- purrr::map_dfc(
      .x = purrr::set_names(c("uniform", "increasing", "decreasing")), 
      .f = function(x) 
        unlist(compare_robustness(list_param = list_con_mod, 
                                  dens = dens, 
                                  new_nr = nr, 
                                  new_nc = nc, 
                                  ext_seq = x))) %>%  
      dplyr::mutate(Topology = "Modular", 
                    rho = rho[1], 
                    j = seq_along(eps))
    list_con_nest <- lapply( eps,
                             function(j) {
                               list(con = matrix(c(j, j,
                                                   j, 1), 2, 2),
                                    pi = pi,
                                    rho = rho)})
    rob_nest <- purrr::map_dfc(.x = purrr::set_names(c("uniform", "increasing", "decreasing")), 
                               .f = function(x) 
                                 unlist(compare_robustness(list_param = list_con_nest, 
                                                           dens = dens, 
                                                           new_nr = nr, 
                                                           new_nc = nc, 
                                                           ext_seq = x))) %>%  
      dplyr::mutate(Topology = "Core-Periphery", 
                    rho = rho[1], 
                    j = seq_along(eps))
    robust_topology <- dplyr::bind_rows(robust_topology, rob_mod, rob_nest)
    saveRDS(robust_topology, "res_topo.rds")
  }
} else {
  robust_topology <- readRDS("res_topo.rds")
}


## ----eval=FALSE---------------------------------------------------------------
#    robust_topology <- tibble::tibble()
#    eps <- c(1/seq(8, 1.5, by = -.5), seq(1, 8, by = .5))
#    for(i in seq(19)) {
#      rho <- c(i*.05, (20-i)*.05)
#      list_con_mod <- lapply( eps,
#                              function(j) {
#                                list(con = matrix(c(j, 1,
#                                                    1, j), 2, 2),
#                                     pi = pi,
#                                     rho = rho)})
#      rob_mod <- purrr::map_dfc(
#        .x = purrr::set_names(c("uniform", "increasing", "decreasing")),
#        .f = function(x)
#          unlist(compare_robustness(list_param = list_con_mod,
#                                    dens = dens,
#                                    new_nr = nr,
#                                    new_nc = nc,
#                                    ext_seq = x))) %>%
#        dplyr::mutate(Topology = "Modular",
#                      rho = rho[1],
#                      j = seq_along(eps))
#      list_con_nest <- lapply( eps,
#                               function(j) {
#                                 list(con = matrix(c(j, j,
#                                                     j, 1), 2, 2),
#                                      pi = pi,
#                                      rho = rho)})
#      rob_nest <- purrr::map_dfc(.x = purrr::set_names(c("uniform", "increasing", "decreasing")),
#                                 .f = function(x)
#                                   unlist(compare_robustness(list_param = list_con_nest,
#                                                             dens = dens,
#                                                             new_nr = nr,
#                                                             new_nc = nc,
#                                                             ext_seq = x))) %>%
#        dplyr::mutate(Topology = "Core-Periphery",
#                      rho = rho[1],
#                      j = seq_along(eps))
#      robust_topology <- dplyr::bind_rows(robust_topology, rob_mod, rob_nest)
#    }

## ----print_topo---------------------------------------------------------------
prt <- robust_topology %>%
#  rename(Uniform = uniform, Increasing = increasing, Decreasing = decreasing) +
  tidyr::pivot_longer(cols = c("decreasing","uniform", "increasing"), 
                      names_to = "Extinction", 
                      values_to = "Robustness") %>%
  mutate(Extinction = forcats::as_factor(
    case_when(Extinction == "decreasing" ~ "Decreasing",
              Extinction == "uniform" ~ "Uniform",
              Extinction == "increasing" ~ "Increasing"))) %>% 
  ggplot(ggplot2::aes(x = j, y = rho)) +
  geom_tile(ggplot2::aes(fill = Robustness)) +
  facet_grid(Topology  ~ Extinction) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                midpoint = rob_er, 
                                guide_colorbar(title = "R")) +
  annotate(x = 15, y = .5, geom = "point", col = "black", size = 3) +
  annotate(x = 15, y = .55, geom = "text", label = "ER",
                    col = "black", size = 3) +
  scale_x_continuous(breaks = c(0, 7, 15, 22, 30), 
                              labels = c(1/8, 1/4, 1, 4, 8)) +
  coord_fixed(ratio = 30) +
  theme_bw(base_size = 15, base_line_size = 0)

# df_shape <- data.frame(
#   x = c(1, 1, 29, 29, 1, 14),
#   y = c(.05, .05, .05, .5, .95, .05),
#   button = c("a", "d", "b","e", "c", "f"),
#   Topology = rep(c("Core-Periphery", "Modular"), 3),
#   Extinction = c(rep("Decreasing", 2), rep("Increasing", 2), rep("Uniform", 2)))
df_shape <- data.frame(
  x = c(1  , 29 , 29,   29,  1, 1),
  y = c(.85, .05, .35, .95, .6, .4),
  button = c(15L, 17L, 18L, 0L, 2L, 5L),#,
  Topology = rep(c("Core-Periphery", "Modular"), each = 3))
#  Extinction = c(rep("Decreasing", 2), rep("Increasing", 2), rep("Uniform", 2)))
prt + geom_point(data = df_shape, mapping = aes(x = x, y = y, shape = button), 
                 size = 4, show.legend = FALSE)+
  scale_shape_identity() 

## ----compute_gnp--------------------------------------------------------------
if (! file.exists("res_gnp.rds")) {
  n_iter <- 300
  pi <- c(.1,.9)
  rho <- c(.2,.8)
  con <- matrix(c(.8, .4, .4, .1), 2, 2)
  con <- 0.1*con/as.vector(pi%*%con%*%rho)
  nr <- 100
  nc <- 100
  res <- pbmcapply::pbmclapply(
    X = seq(500),
    FUN = function(i) {
      G <- simulate_lbm(nr = nr, nc = nc, pi = pi, rho = rho, con = con)
      auc_mc <- robustness_emp(G$A, nb_iter =  n_iter)$auc
      par <- get_lbm_param(G$A)
      auc_lbm <- robustness_lbm(par$con, par$pi, par$rho, nr, nc)$auc
      return(tibble(type = "cp", mc = auc_mc, lbm = auc_lbm ))
    }, mc.cores = 8
  )
  res_cp <- bind_rows(res)
  
  n_iter <- 300
  pi <- 1
  rho <- 1
  con <- matrix(0.1, 1, 1)
  nr <- 100
  nc <- 100
  res <- pbmcapply::pbmclapply(
    X = seq(500),
    FUN = function(i) {
      G <- simulate_lbm(nr = nr, nc = nc, pi = pi, rho = rho, con = con)
      auc_mc <- robustness_emp(G$A, nb_iter =  n_iter)$auc
      par <- get_lbm_param(G$A)
      auc_lbm <- robustness_lbm(par$con, par$pi, par$rho, nr, nc)$auc
      return(tibble(type = "er", mc = auc_mc, lbm = auc_lbm ))
    }, mc.cores = 8
  )
  res_er <- bind_rows(res)
  
  pi  <- c(.2, .8)
  rho <- c(.9, .1)
  con <- matrix(c(.6, .1, .2, .2), 2, 2)
  con <- 0.1*con/as.vector(pi%*%con%*%rho)
  res <- pbmcapply::pbmclapply(
    X = seq(500),
    FUN = function(i) {
      G <- simulate_lbm(nr = nr, nc = nc, pi = pi, rho = rho, con = con)
      auc_mc <- robustness_emp(G$A, nb_iter =  n_iter)$auc
      par <- get_lbm_param(G$A)
      auc_lbm <- robustness_lbm(par$con, par$pi, par$rho, nr, nc)$auc
      return(tibble(type = "mod", mc = auc_mc, lbm = auc_lbm ))
    }, mc.cores = 8
  )
  res_mod <- bind_rows(res)
  res_gnp <- bind_rows(res_cp, res_mod, res_er)
  saveRDS(res_gnp, "res_gnp.rds")
} else {
  res_gnp <- readRDS("res_gnp.rds")
}

## ----compute_gnm--------------------------------------------------------------
if (! file.exists("res_gnm.rds")) {
  n_iter <- 300
  nr <- 100
  nc <- 100
  
  pi <- c(.1,.9)
  rho <- c(.2,.8)
  con <- matrix(c(.8, .4, .4, .1), 2, 2)
  con <- 0.1*con/as.vector(pi%*%con%*%rho)
  pi <- nr * pi
  rho <- nc * rho
  con <- round(outer(pi, rho) * con)
  
  res <- pbmcapply::pbmclapply(
    X = seq(500),
    FUN = function(i) {
      G <- simulate_lbm(nr = nr, nc = nc, pi = pi, rho = rho, con = con, method = "gnm")
      auc_mc <- robustness_emp(G$A, nb_iter =  n_iter)$auc
      par <- get_lbm_param(G$A)
      auc_lbm <- robustness_lbm(par$con, par$pi, par$rho, nr, nc)$auc
      return(tibble(type = "cp", mc = auc_mc, lbm = auc_lbm ))
    }, mc.cores = 8
  )
  res_cp_gnm <- bind_rows(res)
  
  pi <- 1
  rho <- 1
  con <- matrix(0.1, 1, 1)
  pi <- nr * pi
  rho <- nc * rho
  con <- round(outer(pi, rho) * con)
  res <- pbmcapply::pbmclapply(
    X = seq(500),
    FUN = function(i) {
      G <- simulate_lbm(nr = nr, nc = nc, pi = pi, rho = rho, con = con, method = "gnm")
      auc_mc <- robustness_emp(G$A, nb_iter =  n_iter)$auc
      par <- get_lbm_param(G$A)
      auc_lbm <- robustness_lbm(par$con, par$pi, par$rho, nr, nc)$auc
      return(tibble(type = "er", mc = auc_mc, lbm = auc_lbm ))
    }, mc.cores = 8
  )
  res_er_gnm <- bind_rows(res)
  
  pi  <- c(.2, .8)
  rho <- c(.9, .1)
  con <- matrix(c(.6, .1, .2, .2), 2, 2)
  con <- 0.1*con/as.vector(pi%*%con%*%rho)
  pi <- nr * pi
  rho <- nc * rho
  con <- round(outer(pi, rho) * con)
  res <- pbmcapply::pbmclapply(
    X = seq(500),
    FUN = function(i) {
      G <- simulate_lbm(nr = nr, nc = nc, pi = pi, rho = rho, con = con, method = "gnm")
      auc_mc <- robustness_emp(G$A, nb_iter =  n_iter)$auc
      par <- get_lbm_param(G$A)
      auc_lbm <- robustness_lbm(par$con, par$pi, par$rho, nr, nc)$auc
      return(tibble(type = "mod", mc = auc_mc, lbm = auc_lbm ))
    }, mc.cores = 8
  )
  res_mod_gnm <- bind_rows(res)  
  res_gnm <- bind_rows(res_cp_gnm, res_mod_gnm, res_er_gnm)
  saveRDS(res_gnm, "res_gnm.rds")
} else {
  res_gnm <- readRDS("res_gnm.rds")
}


## ----compute_fixed_block------------------------------------------------------
if (! file.exists("res_block.rds")) {
  n_iter <- 300
  nr <- 100
  nc <- 100
  
  pi <- c(.1,.9)
  rho <- c(.2,.8)
  con <- matrix(c(.8, .4, .4, .1), 2, 2)
  con <- 0.1*con/as.vector(pi%*%con%*%rho)
  simz <- function() {
    A <- diag(0, 100)
     for (k in seq_along(pi)) {
      for (q in seq_along(rho)) {
        A[Z == k, W == q] <- igraph::as_incidence_matrix(
          igraph::sample_bipartite(sum(Z == k), sum(W == q),
                                   type = "gnp", p = con[k,q]))
      }
     }
    return(A)
  }
  
  Z <-  rep(seq(2), times = 100*pi)
  W <-  rep(seq(2), times = 100*rho)
  res <- pbmcapply::pbmclapply(
    X = seq(500),
    FUN = function(i) {
      G <- simz()
      auc_mc <- robustness_emp(G, nb_iter =  n_iter)$auc
      par <- get_lbm_param(G)
      auc_lbm <- robustness_lbm(par$con, par$pi, par$rho, nr, nc)$auc
      return(tibble(type = "cp", mc = auc_mc, lbm = auc_lbm ))
    }, mc.cores = 8
  )
  res_cp_block <- bind_rows(res)
  
  pi <- 1
  rho <- 1
  con <- matrix(0.1, 1, 1)
  Z <-  rep(seq(1), times = 100*pi)
  W <-  rep(seq(1), times = 100*rho)
  res <- pbmcapply::pbmclapply(
    X = seq(500),
    FUN = function(i) {
      G <- simz()
      auc_mc <- robustness_emp(G, nb_iter =  n_iter)$auc
      par <- get_lbm_param(G)
      auc_lbm <- robustness_lbm(par$con, par$pi, par$rho, nr, nc)$auc
      return(tibble(type = "er", mc = auc_mc, lbm = auc_lbm ))
    }, mc.cores = 8
  )
  res_er_block <- bind_rows(res)
  
  pi  <- c(.2, .8)
  rho <- c(.9, .1)
  con <- matrix(c(.6, .1, .2, .2), 2, 2)
  con <- 0.1*con/as.vector(pi%*%con%*%rho)
  W <-  rep(seq(2), times = 100*rho)
  Z <-  rep(seq(2), times = 100*pi)
  res <- pbmcapply::pbmclapply(
    X = seq(500),
    FUN = function(i) {
      G <- simz()
      auc_mc <- robustness_emp(G, nb_iter =  n_iter)$auc
      par <- get_lbm_param(G)
      auc_lbm <- robustness_lbm(par$con, par$pi, par$rho, nr, nc)$auc
      return(tibble(type = "mod", mc = auc_mc, lbm = auc_lbm ))
    }, mc.cores = 8
  )
  res_mod_block <- bind_rows(res)  
  res_block <- bind_rows(res_cp_block, res_mod_block, res_er_block)
  saveRDS(res_block, "res_block.rds")
} else {
  res_block <- readRDS("res_block.rds")
}


## -----------------------------------------------------------------------------
res_gnp %>% 
  mutate(Model = "Canonical") %>%  
  pivot_longer(cols = - c(type, Model), names_to = "Method", values_to = "Robustness") %>% 
#  rename(Topology = type) %>% 
  mutate(
    Topology = case_when(
      type == "cp" ~ "Core-Periphery",
      type == "er" ~ "ER",
      type == "mod" ~ "Modular"
    )) %>% 
  filter(Method == "mc") %>% 
  ggplot(aes(x = Robustness, fill = Topology, lty = Topology)) +
  geom_density(alpha = .3) +
#  ggridges::geom_density_ridges(aes(y = Model), alpha = .3) +
#  geom_density(alpha = .2) +
  ylab(label = "") +
#  facet_wrap(~ Model) +
  theme_minimal(base_size = 15) 

