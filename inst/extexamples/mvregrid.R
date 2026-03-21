require("compositions")
data(AnimalVegetation)
AV <- as.data.frame(AnimalVegetation)
AV$region <- factor(ifelse(AV$regA == 1, "A", "B"))

AVmod <- lm(ilr(cbind(disc,spick,din,spin)) ~ region, data = AV)

AVRG <- ref_grid(AVmod)   |> suppressWarnings()
AVRG

lvls <- c("disc","spick","din","spin")
mvregrid(AVRG, newname = "comp", newlevs = lvls) |>
    confint(by = "region")

# The ref grid we'd have gotten if the 'clr' transform had been used...
mvregrid(AVRG, transform = "clr", newname = "comp", newlevs = lvls)
