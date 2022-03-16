# Fletcher, Q&D:
pt <- function(B, m, Binf, n) {
  gamma <- n^(n / (n - 1)) / (n - 1)
  gamma * m * (B / Binf) - gamma * m * (B / Binf)^n
}

# H&W 8.4.4:
pt_hw <- function(B, r, k, m) {
  r * B - (r / k) * B^m
}

# https://www.sciencedirect.com/science/article/abs/pii/S0304380015002732
# https://tinyurl.com/yuvhtvjt
pt2 <- function(B, r, n, K) {
  (r / n) * B * (1 - (B / K)^n)
}

# https://slideplayer.com/slide/16152688/95/images/5/Pella-Tomlison+Surplus+Production+Model.jpg
pt3 <- function(B, r, n, K) {
  (r / (n - 1)) * B * (1 - (B / K)^(n - 1))
}

plot(1:100, sapply(1:100, function(x) pt(B = x, 60, Binf = 100, n = 2)), type = "l", col = "black")
lines(1:100, sapply(1:100, function(x) pt(B = x, 60, Binf = 100, n = 0.5)), col = "red")
lines(1:100, sapply(1:100, function(x) pt(B = x, 60, Binf = 100, n = 4)), col = "blue")
# :)

plot(1:100, sapply(1:100, function(x) pt_hw(B = x, r = 0.5, k = 100, m = 2)), type = "l", col = "black")
B <- 1:100
lines(B, sapply(B, function(x) pt_hw(B = x, r = 0.2, k = 100, m = 1.90)), type = "l", col = "red")
# :(

plot(1:100, sapply(1:100, function(x) pt2(B = x, r = 0.5, K = 100, n = 2)), type = "l", col = "black")
lines(1:100, sapply(1:100, function(x) pt2(B = x, r = 0.32, K = 100, n = 0.5)), type = "l", col = "black")
# :(
# (correct but different definition of 'n')
# i.e., need to do this:
# note the need to change 'r' compared to Fletcher

plot(1:100, sapply(1:100, function(x) pt2(B = x, r = 0.5, K = 100, n = 1)), type = "l", col = "black")
lines(1:100, sapply(1:100, function(x) pt2(B = x, r = 0.25, K = 100, n = -0.5)), type = "l", col = "black")

plot(1:100, sapply(1:100, function(x) pt3(B = x, r = 0.5, K = 100, n = 2)), type = "l", col = "black")
lines(1:100, sapply(1:100, function(x) pt3(B = x, r = 0.25, K = 100, n = 0.5)), type = "l", col = "black")
# :)
# note the need to change 'r' compared to Fletcher
