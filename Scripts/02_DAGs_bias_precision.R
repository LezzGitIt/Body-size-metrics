## Conceptual figure for the paper 'Bergmann's Rule - A practical guide for selecting body size metrics'

# Libraries ---------------------------------------------------------------

library(dagitty)
library(png)

# Define DAGs -------------------------------------------------------------

# Increase precision
model_a <- dagitty('
dag {

  bb = "-.5,-.5,.5,.5"

  Latitude  [exposure, pos="-0.25, 0"]
  Age_sex   [pos="0.25,0"]
  Size      [outcome, pos="0,-0.45"]

  Latitude -> Size
  Age_sex -> Size
}
')
plot(model_a)

# Decrease precision
model_b <- dagitty('
dag {

  bb = "-.5,-.5,.5,.5"

  Latitude  [exposure, pos="-0.25, 0"]
  Age_sex   [pos="0.25,0"]
  Size      [outcome, pos="0,-0.45"]
  " "       [pos="0,0.1"]

  Latitude <-> Age_sex
  Latitude -> Size
}
')
plot(model_b)

# Age_sex are confounders - Differential migration
model_c <- dagitty('
dag {

  bb = "-.5,-.5,.5,.5"

  Latitude  [exposure, pos="-0.25, 0"]
  Age_sex   [pos="0.25,0"]
  Size      [outcome, pos="0,-0.45"]

  Age_sex -> Latitude
  Latitude -> Size
  Age_sex -> Size
}
')
plot(model_c)

# Sampling causes non-random latitudinal gradient in age and sex distribution
model_d <- dagitty('
dag {

  bb = "-.5,-.5,.5,.5"

  Sampling  [pos="0, 0.45"]
  Latitude  [exposure, pos="-0.25, 0"]
  Age_sex   [pos="0.25,0"]
  Size      [outcome, pos="0,-0.45"]

  Sampling -> Latitude
  Sampling -> Age_sex
  Latitude -> Size
  Age_sex -> Size
}
')
plot(model_d)

# Plot --------------------------------------------------------------------
# Open graphics device and set layout
png("Figures/dags_precision_bias.png",
    width = 1400,
    height = 1400,
    res = 200)

par(
  mfrow = c(2, 2),
  mar = c(2, 2, 3, 2)
)

# Common plotting arguments
dag_args <- list(
  node.names = c("Age_sex" = "Age / sex"),
  cex = 1.6
)

# Panel a
plot(model_a,
     node.names = dag_args$node.names,
     cex = dag_args$cex)
mtext("a", side = 3, adj = 0, font = 2, cex = 1.5)

# Panel b
plot(model_b,
     node.names = dag_args$node.names,
     cex = dag_args$cex)
mtext("b", side = 3, adj = 0, font = 2, cex = 1.5)

# Panel c
plot(model_c,
     node.names = dag_args$node.names,
     cex = dag_args$cex)
mtext("c", side = 3, adj = 0, font = 2, cex = 1.5)

# Panel d
plot(model_d,
     node.names = dag_args$node.names,
     cex = dag_args$cex)
mtext("d", side = 3, adj = 0, font = 2, cex = 1.5)

# Close graphics device
dev.off()

# Return graphing parameters to normal
par(mfrow = c(1,1))

# Supplementary -----------------------------------------------------------
# >Temperature-dependent sex determination --------------------------------

# Example with temperature-dependent sex determination in reptiles
# Open graphics device and set layout
png("Figures/Dag_tsd.png")
model_tsd <- dagitty('
dag {

  bb = "-.5,-.5,.5,.5"

  Temperature  [exposure, pos="-0.25, 0"]
  Sex          [pos="0.25,0"]
  Size         [outcome, pos="0,-0.45"]

  Temperature -> Sex
  Temperature -> Size
  Sex -> Size
}
')
plot(model_tsd, cex = dag_args$cex)
#mtext("a", side = 3, adj = 0, font = 2, cex = 1.5)
dev.off()


# >Collider ---------------------------------------------------------------


# Fitting the model Size ~ Latitude + Resources would lead to a biased estimate of Latitude
model_collider <- dagitty('
dag {

  bb = "-.5,-.5,.5,.5"

  Latitude  [exposure, pos="-0.25, 0"]
  Resources   [pos="0.25,0"]
  Size      [outcome, pos="0,-0.45"]

  Latitude -> Resources
  Size -> Resources
}
')
plot(model_collider, cex = dag_args$cex)
mtext("b", side = 3, adj = 0, font = 2, cex = 1.5)

par(mfrow = c(1,1))
