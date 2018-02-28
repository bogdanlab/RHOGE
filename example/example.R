library(RHOGE)

# Estimate \eqn{\rho_ge} genome-wide for BMI and Triglyerides
res <- rhoge.gw(bmi, tg, 14000, 5000)
res2 <- rhoge.bd(bmi, tg, 14000, 5000, p1 = 0.05 / nrow(bmi), p2 = 0.05 / nrow(tg))
