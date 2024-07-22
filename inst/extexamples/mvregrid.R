data(AnimalVegetation)
AV <- as.data.frame(AnimalVegetation)

AVmod <- lm(ilr(cbind(disc,spick,din,spin)) ~ regA, data = AV)

AVRG <- ref_grid(AVmod)   |> suppressWarnings()
AVrg

mvregrid(AVRG, newname = "comp", newlevs = c("disc","spick","din","spin")) |>
    confint(by = "regA")
     