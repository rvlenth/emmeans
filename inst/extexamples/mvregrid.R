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

### Use of fcn argument...
arr.lm <- lm(as.matrix(USArrests) ~ state.region)
arr.emm <- emmeans(arr.lm, ~ crime | state.region, mult.name = "crime")

PC <- princomp(USArrests, cor = TRUE)
std.and.rot <- function(x, pc, ...) {
    z <- sweep(sweep(x, 2, pc$center, "-"), 2, pc$scale, "/")
    z %*% pc$loadings
}

# EMMs obtained if model had been fitted with PC$scores as the response ...
mvregrid(arr.emm, newname = "prin.comp", fcn = std.and.rot, pc = PC)

