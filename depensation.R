# Based on:
# @article{perala2022,
#   title = {Allee Effects and the {{Allee-effect}} Zone in Northwest {{Atlantic}} Cod},
#   author = {Per{\"a}l{\"a}, Tommi and Hutchings, Jeffrey A. and Kuparinen, Anna},
#   journal = {Biology Letters},
#   volume = {18},
#   number = {2},
#   pages = {20210439},
#   doi = {10.1098/rsbl.2021.0439},
# }

bh <- function(St, Rinf, S50) {
  Rinf / (1 + (S50/St))
}

bhc <- function(St, Rinf, S50, .c) {
  Rinf / (1 + (S50/St)^.c)
}

St <- seq(1, 200, length.out = 100)
bh(St, 500, 200)

Rt <- bh(St, 500, 200)

plot(St, Rt)
plot(St, log(Rt/St))
plot(St, Rt/St)

St <- seq(1, 200, length.out = 100)
bhc(St, 500, 200, .c = 1.5)

out <- purrr::map_dfr(c(0.5, 1.2, 1, 1.5, 2), function(.x) {
  Rtc <- bhc(St, 500, 200, .c = .x)
  data.frame(St = St, Rt = Rtc, .c = .x)
})

library(ggplot2)

ggplot(out, aes(St, Rt/St, colour = as.factor(.c))) +
  geom_line()

ggplot(out, aes(St, Rt/St, colour = as.factor(.c))) +
  geom_line() +
  scale_y_log10()

ggplot(out, aes(St, Rt, colour = as.factor(.c))) +
  geom_line()

###

#
# bhc <- function(St, alpha, K, .c) {
#   (alpha * St) / (1 + (S50/St)^.c)
# }
# # (S50/St)^.c
#
# bh_ram <- function(St, alpha, K, .c) {
#   (alpha * St^.c) / (1 + ((St^.c)/K))
# }
#
# bhc(100, 0.1, 1000, 1.2)
# bh_ram(100, 0.1, 1000, 1.2)
