

PlotMeanStddev <- function(pred){
  ncolors <- 200
  greencols.fn <- colorRampPalette(brewer.pal(9, "Greens"))
  greencols <- greencols.fn(ncolors)
  bluecols.fn <- colorRampPalette(brewer.pal(9, "Blues"))
  bluecols <- bluecols.fn(ncolors)
  map.mean <- mapview::mapview(pred, zcol = c("mean"), legend = TRUE,
                               col.regions = greencols)
  map.stddev <- mapview::mapview(pred, zcol = c("stddev"), legend = TRUE, alpha = 0.3, 
                                 col.regions = bluecols)
  leafsync::sync(map.mean, map.stddev)
}


PlotSpatialFields <- function(fitmod, biasfield = TRUE){
  random_sp_effects <- data.frame(cbind(fitmod$model$summary.random$unstr_field$mean,
                                        fitmod$model$summary.random$unstr_field$sd,
                                        Mesh$mesh$loc[,1:2]))
  
  if(biasfield){
    bias_stuff <- data.frame(cbind(fitmod$model$summary.random$bias_field$mean,
                                  fitmod$model$summary.random$bias_field$sd))
    random_sp_effects <- cbind(bias_stuff, random_sp_effects)
  }
  
  column_names <- c("shared_field", "sf_sd", "longitude", "latitude")
  if(biasfield){
    column_names <- c("bias_field", "bf_sd", column_names)
  }
  
  colnames(random_sp_effects) <- column_names
  
  p <- ggplot(random_sp_effects) +
    geom_map(data = norway, map = norway, aes(long, lat, map_id=region), 
             color="#2b2b2b", fill = "white") + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank())
  
  p1 <- p + geom_point(aes(x = longitude, y = latitude, color = shared_field), alpha = 0.6) +
    scale_color_viridis()
  p1.1 <- p + geom_point(aes(x = longitude, y = latitude, color = sf_sd), alpha = 0.6) +
    scale_color_viridis()
  
  if(biasfield){
    p2 <- p + geom_point(aes(x = longitude, y = latitude, color = bias_field), alpha = 0.6) +
      scale_color_viridis()
    p2.1 <- p + geom_point(aes(x = longitude, y = latitude, color = bf_sd), alpha = 0.6) +
      scale_color_viridis()
  }
  
  if(biasfield){
    figure <- ggarrange(p1, p1.1, p2, p2.1, ncol = 2, nrow = 2)
    return(figure)
  }
  
  figure <- ggarrange(p1, p1.1, ncol = 2, nrow = 1)
  return(figure)
}

plot_obs <- function(dataset){
  nosefi <- map_data("world", region = c("Norway(?!:Svalbard)", 
                                         "Sweden", "Finland")) 
  p <- ggplot(dataset) +
    geom_polygon(data = nosefi, aes(long, lat, group = group), 
                 color="#2b2b2b", fill = "white") +
    geom_point(aes(x = decimalLongitude, y = decimalLatitude), 
               color = 'hotpink4', alpha = 0.6, size = 0.5) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    ggtitle(deparse(substitute(dataset)))
  return(p)
}

proj_random_field <- function(model, sp_polygon, mesh){
  Projection <- CRS("+proj=longlat +ellps=WGS84")
  
  r0 <- diff(range(bbox(sp_polygon)[1,]))/diff(range(bbox(sp_polygon)[2,]))
  prj <- inla.mesh.projector(mesh, 
                             xlim=bbox(sp_polygon)[1,],
                             ylim=bbox(sp_polygon)[2,],
                             dims=c(200*r0, 200))
  
  m.bias <- inla.mesh.project(prj, model$summary.ran$bias_field$mean)
  sd.bias <- inla.mesh.project(prj, model$summary.ran$bias_field$sd)
  m.shared <- inla.mesh.project(prj, model$summary.ran$shared_field$mean)
  sd.shared <- inla.mesh.project(prj, model$summary.ran$shared_field$sd)
  proj.points <- SpatialPoints(prj$lattice$loc)
  proj4string(proj.points) <- Projection 
  proj4string(norway.poly) <- Projection
  ov <- sp::over(proj.points, sp_polygon)
  sd.bias[is.na(ov)] <- m.bias[is.na(ov)] <- NA
  sd.shared[is.na(ov)] <- m.shared[is.na(ov)] <- NA
  
  spat_coords <- expand.grid(x = prj$x, y = prj$y)
  bias_df <- data.frame(decimalLongitude = spat_coords$x, decimalLatitude = spat_coords$y,
                        mean = as.vector(m.bias), sd = as.vector(sd.bias))
  shared_df <- data.frame(decimalLongitude = spat_coords$x, decimalLatitude = spat_coords$y,
                          mean = as.vector(m.shared), sd = as.vector(sd.shared))
  tidyr::drop_na(bind_rows(bias = bias_df, shared = shared_df, .id = "field"))
}

make_random_field_df <- function(Model, sp_polygon, mesh){
  spat_fields_df <- proj_random_field(Model$model, sp_polygon = sp_polygon, 
                                      mesh = mesh)
  spat_fields <- spat_fields_df %>% tidyr::gather(key = statistic, value = value, mean:sd)
  return(spat_fields)
}

prediction_df <- function(stk.pred, model){
  pred_df <- data.frame(decimalLongitude = stk.pred$predcoords[,1], 
                        decimalLatitude = stk.pred$predcoords[,2],
                        model$predictions)
  #pred_df$precision <- pred_df$stddev^-2
  Pred <- pred_df %>% tidyr::gather(key = statistic, value = value, mean:stddev)
  return(Pred)
}
  


dots_whiskers_inla <- function(model){
  num_rows <- nrow(model$model$summary.fixed)
  coefs <- model$model$summary.fixed[3:num_rows, 1:5]
  coefs$coefficient <- rownames(coefs)
  colnames(coefs) <- c("mean", "sd", "quant0.025", "quant0.5", "quant0.975", "coefficient")

  p <- ggplot(coefs, aes(x = mean, y = coefficient)) +
    geom_vline(xintercept = 0, linetype = "dotted", size = 0.8, color = "grey") +
    geom_segment(aes(x = quant0.025, xend = quant0.975, y = coefficient, yend = coefficient), 
                 size = 0.8) +
    geom_point(shape = 21, size = 4, fill = "orange") +
    #scale_y_discrete(labels = var_labs) +
    theme_minimal() +
    theme(axis.title = element_blank())
  
  return(p)
}
