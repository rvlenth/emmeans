wine.clm <- clm(rating ~ temp * contact, data = wine)

ref_grid(wine.clm)

# verify that we get the same thing via qdrg:
qdrg(object = wine.clm, ordinal.dim = 5)
