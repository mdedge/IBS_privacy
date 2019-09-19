#Take in a genetic map and linearly interpolate base pairs corresponding to given cM

map <- read.table("../IBS_tiling/plink.GRCh37.map/plink.chr19.GRCh37.map")

diffs <- (1:60)/20
interps <- approx(map[,3], map[,4], xout = sort(c((71.4788 - diffs), (71.4788 + diffs) ) ) )

interp.dat <- cbind(round(interps$x - 71.4788, 2), round(interps$y, 2))

write.table(interp.dat, "cM_basepair_interpolation.txt", row.names = FALSE, col.names = FALSE)


