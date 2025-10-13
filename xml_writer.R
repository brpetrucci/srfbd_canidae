###########################################
## XML writing script for SRFBD and FBDS ##
## Chapter 2 - Bruno do Rosario Petrucci ##
###########################################

###
# packages

# ape
library(ape)

###
# little auxiliary functions

# add spaces easily
add_spaces <- function(n) {
  paste(rep(' ', n), collapse = '')
}

# draw from a dirichlet distribution
rdirichlet <- function (n, alpha) {
  # length of parameter vector
  l <- length(alpha)
  
  # make the draws from the gamma
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  
  # normalize the draws
  sm <- x %*% rep(1, l)
  
  # return normalized vector
  x / as.vector(sm)
}

# get line if a name is present in list, and nothing if not
line_if_there <- function(list, name) {
  ifelse(name %in% names(list),
         paste0(' ', name, '="', list[name], '"'),
         '')
}

# function to switch 0/1 to {01}
rewrite_polymorphisms <- function(morpho) {
  # create new morpho list
  new_morpho <- vector("list", length(morpho))
  
  # iterate through morpho
  for (i in seq_along(morpho)) {
    # modify it
    new_morpho[[i]] <- ifelse(grepl("^(\\d+(?:/\\d+)+)$", 
                                    morpho[[i]], perl = TRUE),
                              paste0("{", gsub("/", "", morpho[[i]]), "}"),
                              morpho[[i]])
  }
  
  # name it
  names(new_morpho) <- names(morpho)
  
  # return new morpho
  return(new_morpho)
}

# make vector of numbers into character list
range_string <- function(num_list) {
  # choose only unique values of cha
  num_list <- sort(unique(num_list))
  
  # start result vector
  result <- c()
  
  # start current range
  start <- num_list[1]
  end <- num_list[1]
  
  # iterate through rest of list
  for (i in 2:length(num_list)) {
    # if this value is just one more than the previous end, expand end
    if (num_list[i] == end + 1) {
      end <- num_list[i]
    } else {
      # if not, add previous range
      if (start == end) {
        # if just one number, just add it
        result <- c(result, as.character(start))
      } else {
        # if not, add it as a range
        result <- c(result, paste0(start, "-", end))
      }
      
      # start new range from this value
      start <- end <- num_list[i]
    }
  }
  
  # add final range
  if (start == end) {
    result <- c(result, as.character(start))
  } else {
    result <- c(result, paste0(start, "-", end))
  }
  
  # combine into one string
  paste(result, collapse = ",")
}

# make list of morpho partitions and parameters
morpho_partitions <- function(morpho, ascertain) {
  # list to hold character number for each number of states
  char_nums <- list()
  
  # make morpho a data frame
  morpho <- as.data.frame(morpho)
  
  # iterate through columns of morpho
  for (i in 1:nrow(morpho)) {
    # this character data
    char_data <- unlist(morpho[i, ])
    
    # states in this data (not including ? or {})
    char_states <- unique(unlist(strsplit(char_data, "")))
    char_states <- char_states[!(char_states %in% c("?", "{", "}"))]
    
    # number of states for this character
    n_states <- length(char_states)
    
    # name of this character partition
    char_part <- as.character(n_states)
    
    # check if morphoN is already in
    if (char_part %in% names(char_nums)) {
      # add this character to that element
      char_nums[[char_part]] <- c(char_nums[[char_part]], i)
    } else {
      # if not, add it
      char_nums[char_part] <- i
    }
  }
  
  # make results list
  res <- list()
  
  # if ascertained, get ascertainment state numbers
  if (ascertain) {
    # get maximum number of states
    max_states <- max(as.numeric(names(char_nums)))
    
    # ascertain states
    n_ascertain <- 
      (max(unlist(char_nums)) + 1):(max(unlist(char_nums)) + max_states)
  }
  
  # iterate through char_nums to fill results list
  for (i in 1:length(char_nums)) {
    # check number of states
    n_states <- as.numeric(names(char_nums)[i])
    
    # start exclude_from and _to to NULL
    exclude_from <- exclude_to <- NULL
    
    # if ascertained, get exclude to and from, and add n_ascertain
    if (ascertain) {
      # exclude to and from
      exclude_from <- length(char_nums[[i]])
      exclude_to <- exclude_from + n_states
      
      # add n_ascertain to list
      char_nums[[i]] <- c(char_nums[[i]], n_ascertain[1:n_states])
    } 
    
    # add to result
    res[[paste0("morpho", n_states)]] <- list(
      ascertain = ifelse(ascertain, "true", "false"),
      excludefrom = exclude_from,
      excludeto = exclude_to,
      filter = range_string(char_nums[[i]])
    )
  }
  
  # return res
  return(res)
}

# recursively add to names of lists
recursive_names <- function(l, add, mark) {
  # iterate through list
  for (i in 1:length(l)) {
    # check if there are names at all
    if (length(names(l)) > 0) {
      # check if name has the mark
      if (grepl(mark, names(l)[i])) {
        # if so, add to it
        names(l)[i] <- paste0(names(l)[i], add)
      }
      
      # check if element is a list
      if (typeof(l[[i]]) == "list") {
        # if so, recursively call the function
        l[[i]] <- recursive_names(l[[i]], add, mark)
      }
    }
  }
  
  # return list
  return(l)
}

# make operators for a given partition
part_operators <- function(part, partitions) {
  # if not morpho, start list with gamma shape operators
  if (part != "morpho") {
    res <- list(
      `gammaShapeScaler.s:` = list(
        spec = "AdaptableOperatorSampler",
        weight = 0.05,
        parameter = paste0("@gammaShape.s:", part),
        operator = list(
          `AMVNOperator.s:` = list(
            spec = "kernel.AdaptableVarianceMultivariateNormalOperator",
            allowNonsense = "true",
            beta = 0.05,
            burnin = 400,
            initial = 800,
            weight = 0.1,
            transformations = list(
              `AVMNSumTransform.s:` = list(
                spec = "operator.kernel.Transform$LogConstrainedSumTransform",
                f = paste0("@freqParameter.s:", part)
              ),
              `AVMNLogTransform.s:` = list(
                spec = "operator.kernel.Transform$LogTransform",
                f = setNames(as.list(paste0(c("gammaShape.s:", "kappa.s:"), part)),
                             rep("f", 2))
              ),
              `AVMNNoTransform.s:` = list(
                spec = "operator.kernel.Transform$NoTransform",
                f = "@Tree.t:tree"
              )
            )
          ),
          `gammaShapeScalerX.s:` = list(
            spec = "kernel.BactrianScaleOperator",
            parameter = paste0("@gammaShape.s:", part),
            scaleFactor = 0.5,
            upper = 10,
            weight = 0.1
          )
        )
      ))
  } else {
    res <- list()
  }
  
  res[["ucldStdevScaler.c:"]] <- list(
             spec = "ScaleOperator",
             parameter = paste0("@ucldStdev.c:", part),
             scaleFactor = 0.5,
             weight = 3
           )
  res[["CategoriesRandomWalk.c:"]] <- list(
             spec = "operator.IntRandomWalkOperator",
             parameter = paste0("@rateCategories.c:", part),
             weight = 10,
             windowSize = 1
           )
  res[["CategoriesSwapOperator.c:"]] <- list(
             spec = "operator.SwapOperator",
             intparameter = paste0("@rateCategories.c:", part),
             weight = 10
           )
  res[["CategoriesUniform.c:"]] <- list(
             spec = "operator.UniformOperator",
             parameter = paste0("@rateCategories.c:", part),
             weight = 10
           )
  
  # if part is not the first in partitions, create mean operators
  if (part != names(partitions)[1]) {
    res[["ucldMeanScaler.c:"]] <- list(
      spec = "ScaleOperator",
      parameter = paste0("@ucldMean.c:", part),
      scaleFactor = 0.5,
      weight = 1
    )
    
    res[["relaxedUpDownOperator.c:"]] <- list(
      spec = "operator.UpDownOperator",
      scaleFactor = 0.75,
      weight = 3,
      up = paste0("@ucldMean.c:", part),
      down = "@Tree.t:tree"
    )
  }
  
  # if part is not morpho, do kappa and frequencies
  if (part != "morpho") {
    res[["KappaScaler.s:"]] <- list(
      spec = "AdaptableOperatorSampler",
      weight = 0.05,
      parameter = paste0("@kappa.s:", part),
      operator = paste0("@AMVNOperator.s:", part),
      operator = list(
        `KappaScalerX.s:` = list(
          spec = "kernel.BactrianScaleOperator",
          parameter = paste0("@kappa.s:", part),
          scaleFactor = 0.1,
          upper = 10,
          weight = 0.1
        )
      )
    )
    
    res[["FrequenciesExchanger.s:"]] <- list(
      spec = "AdaptableOperatorSampler",
      weight = 0.05, 
      parameter = paste0("@freqParameter.s:", part),
      operator = paste0("@AMVNOperator.s:", part),
      operator = list(
        `FrequenciesExchangerX.s:` = list(
          spec = "operator.kernel.BactrianDeltaExchangeOperator",
          parameter = paste0("@freqParameter.s:", part),
          delta = "0.01",
          weight = "0.1"
        )
      )
    )
  }
  
  # name everything
  res <- recursive_names(res, part, ":")
  
  # return res
  return(res)
}

# make logs based on partitions
make_logs <- function(idref, parts) {
  # create result
  res <- list()
  
  # iterate through desired parts
  for (i in 1:length(parts)) {
    # if it's rate, it's a bit more complicated
    if (idref == "rate.c:") {
      res[[paste0(idref, parts[i])]] <- list(
        spec = "beast.base.evolution.RateStatistic",
        branchratemodel = paste0("@RelaxedClock.c:", parts[i]),
        tree = "@Tree.t:tree"
      )
    } else {
      # otherwise, just make it a normal one
      res[[paste0(idref, parts[i])]] <- paste0(idref, parts[i])
    }
  }
  
  return(res)
}

###
# functions to write scripts

# break up data into first and last occurrences
break_ranges <- function(ranges, mol, morpho, 
                         attach, one_sp = FALSE) {
  # variables to hold results
  new_ranges <- data.frame(matrix(nrow = 0, ncol = 4))
  new_mol <- new_morpho <- list()
  
  # make ages vector
  ages <- c()
  
  # empty mol and morpho vectors for non-attached data
  mol_empty <- paste(rep("?", length(mol[[1]])))
  morpho_empty <- paste(rep("?", length(morpho[[1]])))
  
  # iterate through ranges
  for (i in 1:nrow(ranges)) {
    # this range
    range <- ranges[i, ]
    
    # get taxon
    tx <- range$taxon
    
    if (one_sp) {
      # specimen will be just from fa_max to la_min
      new_ranges <- rbind(new_ranges, c(tx, tx, range$fa_max, range$la_min))
      
      # draw age and add it to vector
      ages <- c(ages, runif(1, range$la_min, range$fa_max))
      
      # add mol and morpho
      new_mol[[tx]] <- mol[[tx]]
      new_morpho[[tx]] <- morpho[[tx]]
    } else {
      # if fa_min is lower than la_min, change it to be la_min
      if (range$fa_min < range$la_min) {
        range$fa_min <- range$la_min
      }
      
      # add first specimen to new_ranges
      new_ranges <- rbind(new_ranges, c(paste0(tx, "_first"), tx,
                                        range$fa_max, range$fa_min))
      
      # get first age
      first_age <- runif(1, range$fa_min, range$fa_max)
      
      # draw age and add it to vector
      ages <- c(ages, first_age)
      
      # add mol and morpho
      new_mol[[paste0(tx, "_first")]] <- mol_empty
      new_morpho[[paste0(tx, "_first")]] <- morpho[[tx]]
      
      # if it is not a singleton
      if (!(range$fa_max == range$la_max && range$fa_min == range$la_min)) {
        # add last specimen to new_ranges
        new_ranges <- rbind(new_ranges, c(paste0(tx, "_last"), tx,
                                          range$la_max, range$la_min))
        
        # draw age and add it to vector
        ages <- c(ages, runif(1, range$la_min, 
                              min(range$la_max, first_age)))
        
        # add mol and morpho
        new_mol[[paste0(tx, "_last")]] <- mol[[tx]]
        new_morpho[[paste0(tx, "_last")]] <- ifelse(rep(attach == "both", 
                                                        length(morpho[[1]])),
                                                    morpho[[tx]],
                                                    morpho_empty)
      }
    }
  }
  
  # name everything
  colnames(new_ranges) <- c("specimen", "taxon", "max_age", "min_age")
  names(new_mol) <- names(new_morpho) <- names(ages) <- new_ranges$specimen
  
  # return
  return(list(RANGES = new_ranges, AGES = ages, 
              MOL = new_mol, MORPHO = new_morpho))
}

# write data to script
write_data <- function(data, type) {
  # start data block
  data_block <- paste0(add_spaces(4), '<data')
  
  # add id and datatype if necessary
  data_block <- c(data_block,
                  paste0('id="', type, '"'),
                  'spec="Alignment"')
  if (type == "morpho") data_block <- c(data_block,
                                        'dataType="standard"')
  data_block[length(data_block)] <- paste0(data_block[length(data_block)],
                                           ">")
  # get number of states
  n_states <- sum(!unique(unlist(strsplit(unique(unlist(data)), ""))) %in% 
                    strsplit("-?nyrswkm{}", "")[[1]])

  # iterate through data
  for (i in 1:length(data)) {
    # start line
    data_line <- paste0(add_spaces(8),
                        '<sequence id="',
                        type,
                        '_seq_',
                        names(data)[i],
                        '" spec="Sequence" taxon="',
                        names(data)[i],
                        '" totalcount="',
                        n_states,
                        '" value="',
                        ifelse(type == "morpho",
                        paste(c(data[[i]], 0:(n_states - 1)), collapse = ""),
                        paste(data[[i]], collapse = "")),
                        '"/>')
    
    # add to data block
    data_block <- c(data_block, data_line)
  }
  
  # if type is morpho, add user data type
  if (type == "morpho") {
    data_block <- c(data_block,
                    paste0(add_spaces(8),
                           '<userDataType id="StandardData.0" ',
                           'spec="beast.base.evolution.datatype.StandardData" ',
                           'ambiguities="" nrOfStates="5"/>'))
  }
  
  # add final data call
  data_block <- c(data_block, 
                  paste0(add_spaces(4),
                         '</data>'))
  
  # return data block
  return(c(data_block, ""))
}

# write prior blocks
prior_blocks <- function(priors) {
  # start text blocks
  prior_block <- c()
  param_block <- c()
  init_block <- c()
  
  # iterate through priors list
  for (i in 1:length(priors)) {
    # if dist doesn't exist, make param line and skip
    if (!("dist" %in% names(priors[[i]]))) {
      # param_line
      param_block <- c(param_block,
                       paste0(add_spaces(12),
                              '<stateNode id="',
                              names(priors)[i],
                              ifelse("spec" %in% names(priors[[i]]),
                                     paste0('" spec="', priors[[i]]$spec, '"'),
                                     '" spec="parameter.RealParameter"'),
                              line_if_there(priors[[i]], "dimension"),
                              '>', priors[[i]]$value, '</stateNode>'))
      next
    }
    
    # get parameter for which we're starting a prior
    param <- names(priors)[i]
    
    # get type of distribution
    dist <- priors[[i]]$dist
    
    # get lower and upper
    if (!is.null(priors[[i]]$lower)) {
      low <- priors[[i]]$lower
    } else {
      low <- -Inf
    }
    if (!is.null(priors[[i]]$upper)) {
      upp <- priors[[i]]$upper
    } else {
      upp <- Inf
    }
    
    # get distribution parameter lines and initial values
    switch(dist,
           "Exponential" = {
             # make parameter text
             param_text <- paste0(add_spaces(20),
                                  '<Exponential id="prior_',
                                  names(priors)[i],
                                  '_Exponential" name="distr" mean="', 
                                  priors[[i]]$mean, '"/>')
             
             # initialize to Inf
             init_val <- Inf
             
             # while init_val is not within the bounds
             while (!(init_val > low && init_val < upp)) {
               init_val <- rexp(1, 1 / priors[[i]]$mean)
             }
           },
           "Uniform" = {
             # make parameter text
             param_text <- paste0(add_spaces(20),
                                  '<Uniform id="prior_',
                                  names(priors)[i],
                                  '_Uniform" name="distr" lower="',
                                  priors[[i]]$lower,
                                  '" upper="',
                                  priors[[i]]$upper, '"/>')
             
             # get initial value
             init_val <- runif(1, priors[[i]]$lower, priors[[i]]$upper)
           },
           "LogNormal" = {
             # make parameter text
             param_text <- paste0(add_spaces(20),
                                  '<LogNormal id="prior_',
                                  names(priors)[i],
                                  '_LogNormal" name="distr" M="',
                                  priors[[i]]$M,
                                  '" S="',
                                  priors[[i]]$S, '"/>')
             # initialize to Inf
             init_val <- Inf
             
             # while init_val is not within the bounds
             while (!(init_val > low && init_val < upp)) {
               # get initial value
               init_val <- rlnorm(1, priors[[i]]$M, priors[[i]]$S)
             }
           },
           "Gamma" = {
             # make parameter text
             param_text <- paste0(add_spaces(20),
                                  '<Gamma id="prior_',
                                  names(priors)[i],
                                  '_Gamma" name="distr" alpha="',
                                  priors[[i]]$alpha,
                                  '" beta="',
                                  priors[[i]]$beta, '"/>')
             
             # initialize to Inf
             init_val <- Inf
             
             # while init_val is not within the bounds
             while (!(init_val > low && init_val < upp)) {
               # get initial value
               init_val <- rgamma(1, priors[[i]]$alpha, 
                                  scale = priors[[i]]$beta)
             }
           },
           "Beta" = {
             # make parameter text
             param_text <- paste0(add_spaces(20),
                                  '<Beta id="prior_',
                                  names(priors)[i],
                                  '_Beta" name="distr" alpha="',
                                  priors[[i]]$alpha,
                                  '" beta="',
                                  priors[[i]]$beta, '"/>')
             
             # initialize to Inf
             init_val <- Inf
             
             # while init_val is not within the bounds
             while (!(init_val > low && init_val < upp)) {
               # get initial value
               init_val <- rbeta(1, priors[[i]]$alpha, priors[[i]]$beta)
             }
           },
           "Dirichlet" = {
             # make parameter text
             param_text <- paste0(add_spaces(20),
                                  '<distr id="prior_',
                                  names(priors)[i],
                                  '_Dirichlet" spec="distribution.Dirichlet"',
                                  ' alpha="',
                                  paste(priors[[i]]$alpha, collapse = " "),
                                  '"/>')
             
             # initialize to Inf
             init_val <- rep(Inf, length(priors[[i]]$alpha))
             
             # while init_val is not within the bounds
             while (!(all(init_val > low & init_val < upp))) {
               # get initial value
               init_val <- rdirichlet(1, priors[[i]]$alpha)
             }
           })
    
    # add to prior block
    prior_block <- c(prior_block,
                     paste0(add_spaces(16),
                            '<prior id="prior_',
                            param, 
                            '" name = "distribution" x="@',
                            param, '">'),
                     param_text,
                     paste0(add_spaces(16), '</prior>'))
    
    # if origin, the init will just be a little lower than upper
    if (names(priors)[i] == "originFBD.t:tree") {
      init_val <- priors[[i]]$upper - 0.1
    }
    
    # add to parameter block
    param_block <- c(param_block,
                     paste0(add_spaces(12),
                            '<stateNode id="',
                            names(priors)[i],
                            '" spec="parameter.RealParameter"',
                            line_if_there(priors[[i]], "dimension"),
                            line_if_there(priors[[i]], "lower"),
                            line_if_there(priors[[i]], "upper"),
                            '>', 
                            paste(init_val, collapse = " "), '</stateNode>'))
    
    # if in the list of feast parameters, add to init block
    if (param %in% c("diversificationRateFBD.t:tree", "turnoverFBD.t:tree",
                     "samplingProportionFBD.t:tree")) {
      init_block <- c(init_block,
                      paste0(add_spaces(8),
                             '<init spec="feast.parameter.RandomRealParameter" ',
                             'initial="@', param, '">'),
                      paste0(add_spaces(12),
                             '<distr idref="prior_', param, '_', dist, '"/>'),
                      paste0(add_spaces(8), '</init>'))
    }
  }
  
  # return blocks
  return(list(PRIORS = prior_block, PARAMS = param_block, INITS = init_block))
}

# write likelihood blocks
likelihood_blocks <- function(partitions) {
  # start block
  block <- c()
  
  # iterate through partitions
  for (i in 1:length(partitions)) {
    # get partition name
    part_name <- names(partitions)[i]
    
    # get info for this partition
    part <- partitions[[i]]
    
    # start tree likelihood
    block <- c(block,
               paste0(add_spaces(16),
                      '<distribution id="treeLikelihood.',
                      part_name,
                      '" spec="',
                      ifelse(partitions[[i]]$threaded,
                             'Threaded', ''),
                      'TreeLikelihood" tree="@Tree.t:tree"',
                      ifelse(partitions[[i]]$data == "morpho" && part_name != "morpho2",
                             ' branchRateModel="@RelaxedClock.c:morpho"', ''),
                      ifelse(i == 1,
                             paste0(' data="@', part_name, '">'), 
                             '>')))
    
    # if not the first partition, need to add data
    if (i != 1) {
      # start data line
      data_line <- paste0(add_spaces(20),
                          '<data id="', part_name, '" ',
                          'spec="FilteredAlignment" data="@',
                          part$data, '" ',
                          ifelse(part$ascertain == "true",
                                 paste0('ascertained="true" excludefrom="',
                                        part$excludefrom,
                                        '" excludeto="',
                                        part$excludeto, '" '), ''),
                          'filter="', part$filter, '"')
      
      # if it is morpho data, need to add data type statement
      if (part$data == "morpho") {
        # finish data line
        data_line <- paste0(data_line, '>')
        
        # add data line and userdatatype
        block <- c(block, data_line,
                   paste0(add_spaces(24),
                          '<userDataType id="morphDataType.',
                          part_name,
                          '" spec="beast.base.evolution.datatype.',
                          part$dataType$type,
                          '" ambiguities="',
                          part$dataType$ambiguities,
                          '" nrOfStates="',
                          part$dataType$nrOfStates, '"/>'),
                   paste0(add_spaces(20), '</data>'))
      } else {
        # just finish data line
        data_line <- paste0(data_line, '/>')
        block <- c(block, data_line)
      }
    }
    
    # finish data block and start site model
    block <- c(block,
               paste0(add_spaces(20),
                      '<siteModel id="',
                      part$siteModel$id,
                      '" spec="SiteModel" ',
                      'mutationRate="',
                      part$siteModel$mutationRate,
                      '" gammaCategoryCount="',
                      part$siteModel$gammaCategoryCount,
                      '" shape="',
                      part$siteModel$shape,
                      '" proportionInvariant="',
                      part$siteModel$proportionInvariant,
                      '">'))
    
    # get substitution model line
    substmodel_line <- switch(part$siteModel$substModel$model,
                              "hky" = {
                                c(paste0(add_spaces(24),
                                         '<substModel id="hky.s:',
                                         part_name,
                                         '" spec="HKY" kappa="',
                                         part$siteModel$substModel$kappa,
                                         '">'),
                                  paste0(add_spaces(28),
                                         '<frequencies id="estimatedFreqs.s:',
                                         part_name,
                                         '" spec="Frequencies" frequencies="',
                                         part$siteModel$substModel$frequencies,
                                         '"/>'),
                                  paste0(add_spaces(24), '</substModel>'))
                              },
                              "LewisMK" = {
                                c(paste0(add_spaces(24),
                                         '<substModel id="LewisMK.s:',
                                         part_name,
                                         '" spec="morphmodels.evolution.substitutionmodel.LewisMK',
                                         '" datatype="@morphDataType.',
                                         part_name, '"/>'))
                              })
    
    # add substitution model line and finish site model
    block <- c(block, substmodel_line,
               paste0(add_spaces(20),
                      '</siteModel>'))
    
    # get branch rate model line
    if (!is.null(part$branchRateModel$model)) {
      brmodel_line <- switch(part$branchRateModel$model,
                             "ucld" = {
                               c(paste0(add_spaces(20),
                                        '<branchRateModel id="RelaxedClock.c:',
                                        ifelse(partitions[[i]]$data == "morpho",
                                               "morpho", part_name),
                                        '" spec="beast.base.evolution.branchratemodel',
                                        '.UCRelaxedClockModel" ',
                                        'clock.rate="',
                                        part$branchRateModel$clock.rate,
                                        '" rateCategories="',
                                        part$branchRateModel$rateCategories,
                                        '" tree="@Tree.t:tree">'),
                                 paste0(add_spaces(24),
                                        '<LogNormal id="LogNormalDistributionModel.c:',
                                        part_name,
                                        '" M="',
                                        part$branchRateModel$M,
                                        '" S="',
                                        part$branchRateModel$S,
                                        '" meanInRealSpace="true" name="distr"/>'),
                                 paste0(add_spaces(20), '</branchRateModel>'))
                             })
    } else {
      brmodel_line <- character(0)
    }
    
    
    # add to block and finish distribution
    block <- c(block, brmodel_line,
               paste0(paste0(add_spaces(16), '</distribution>')))
  }
  
  # return likelihood block
  return(block)
}

# write operator block
operator_block <- function(operators, op_name = "operator", spaces) {
  # create operator block
  block <- c()
  
  # iterate through operators list
  for (i in 1:length(operators)) {
    # logical on whether the operator is in one line
    oneline <- TRUE
    
    # logical on whether first line is ended
    line_ended <- FALSE
    
    # get this particular operator and its name
    op <- operators[[i]]
    decl <- op_name
    
    # begin this line
    line <- paste0(add_spaces(spaces),
                   '<', decl, ' ')
    
    # if operator has length one, just an idref line
    if (length(op) == 1) {
      line <- paste0(line, 'idref="',
                     op[[1]], '"/>')
      block <- c(block, line)
    } else {
      # add id
      line <- paste0(line, 'id="',
                     names(operators)[i], '" ')
      
      # iterate through members of op
      for (j in 1:length(op)) {
        # check if this element is a list
        if (typeof(op[[j]]) == "list") {
          # check if line is not already ended
          if (!line_ended) {
            # remove the last space
            line <- substr(line, 1, nchar(line) - 1)
            
            # add the line ending
            line <- paste0(line, '>')
            
            # add to block
            block <- c(block, line)
            
            # set line_ended to true
            line_ended <- TRUE
          }
          
          # set oneline to false
          oneline <- FALSE
          
          # get the next line
          block <- c(block,
                     operator_block(op[[j]], names(op)[j], spaces + 4))
        } else {
          # add to line
          line <- paste0(line, names(op)[j],
                         '="', op[[j]], '" ')
        }
      }
      
      # add final line, or finish line if oneline is FALSE
      if (oneline) {
        # remove the last space
        line <- substr(line, 1, nchar(line) - 1)
        
        # finish line and add it to block
        line <- paste0(line, '/>')
        block <- c(block, line)
      } else {
        block <- c(block,
                   paste0(add_spaces(spaces), '</', decl, '>'))
      }
    }
  }
  
  # return block
  return(block)
}

# write one whole script
write_script <- function(model, template, maps, 
                         priors, operators, logs,
                         extra_params,
                         ranges, mol, morpho,
                         partitions,
                         attach, one_sp,
                         gens, store_every, log_every,
                         out_file) {
  # variable to hold full script, started with just first line of template
  script <- c(template[1], "")
  
  # break ranges and data up into first and last occurrences
  broken_ranges <- break_ranges(ranges, mol, morpho, attach, one_sp)
  
  # separate that into each object
  new_ranges <- broken_ranges$RANGES
  ages <- broken_ranges$AGES
  new_mol <- broken_ranges$MOL
  new_morpho <- broken_ranges$MORPHO
  
  # get list of taxa from ranges
  taxa <- unique(new_ranges$taxon)
  
  # add molecular and morphological data
  script <- c(script, 
              write_data(new_mol, "mol"),
              write_data(new_morpho, "morpho"))
  
  # if model is srfbd, add taxonset certain
  if (model == "srfbd") {
    # add taxonset declaration
    script <- c(script,
                paste0(add_spaces(4),
                       '<taxonset id="TaxonSet.certain" spec="TaxonSet">'))
    
    # iterate through new ranges
    for (i in 1:nrow(new_ranges)) {
      # if taxon is certain, add to taxonset
      if (new_ranges$max_age[i] == new_ranges$min_age[i]) {
        script <- c(script, paste0(add_spaces(8),
                                   '<taxon id="',
                                   new_ranges$specimen[i],
                                   '" spec="Taxon"/>'))
      }
    }
    
    # finish taxonset
    script <- c(script, paste0(add_spaces(4),
                               '</taxonset>'), "")
  }
  
  # add maps
  script <- c(script, maps, "")
  
  # start MCMC block
  script <- c(script,
              paste0(add_spaces(4),
                     '<run id="mcmc" spec="MCMC" chainLength="',
                     gens, 
                     '" numInitializationAttempts="100000">'),
              paste0(add_spaces(8),
                     '<state id="state" spec="State" storeEvery="',
                     store_every, '">'),
              paste0(add_spaces(12),
                     '<tree id="Tree.t:tree" spec="',
                     ifelse(model == "srfbd",
                            'sr.evolution.tree.SRTree" nodetype="sr.evolution.tree.SRNode" ',
                            'beast.base.evolution.tree.Tree" '),
                     'name="stateNode">'),
              paste0(add_spaces(16),
                     '<trait id="dateTrait.t:tree" ',
                     'spec="beast.base.evolution.tree.TraitSet" ',
                     'traitname="date-backward" value="'))
  
  # iterate through ages
  for (i in 1:length(ages)) {
    # add age
    script <- c(script,
                paste0(add_spaces(16),
                       names(ages)[i],
                       ' = ', ages[i],
                       ifelse(i == length(ages), '">', ',')))
  }
  
  # add taxonset for first DNA partition
  script <- c(script,
              paste0(add_spaces(20),
                     '<taxa id="TaxonSet.',
                     names(partitions)[1],
                     '" spec="TaxonSet">'),
              paste0(add_spaces(24),
                     '<alignment id="',
                     names(partitions)[1],
                     '" spec="FilteredAlignment" filter="',
                     partitions[[1]]$filter, '">'),
              paste0(add_spaces(28),
                     '<data idref="mol"/>'),
              paste0(add_spaces(24), '</alignment>'),
              paste0(add_spaces(20), '</taxa>'),
              paste0(add_spaces(16), '</trait>'),
              paste0(add_spaces(16), 
                     '<taxonset idref="TaxonSet.',
                     names(partitions)[1], 
                     ifelse(model == "srfbd",
                            '">',
                            '"/>')))
  
  # if model is srfbd, need to add sranges to tree
  if (model == "srfbd") {
    # create stratigrahic range count
    sr_count <- 1
    
    # iterate through taxa
    for (i in 1:length(taxa)) {
      # taxon
      tx <- taxa[i]
      
      # if there is no last occurrence, pass
      #if (!any(grepl(paste0(tx, '_last'), new_ranges$specimen))) {
      #  next
      #}
      
      # first and last occurrences
      first_occ <- paste0('@', tx, '_first')
      last_occ <- ifelse(any(grepl(paste0(tx, '_last'), new_ranges$specimen)),
                         paste0('@', tx, '_last'),
                         paste0('@', tx, '_first'))
      
      # add to script
      script <- c(script,
                  paste0(add_spaces(16),
                         '<stratigraphicRange id="r', sr_count,
                         '" spec="StratigraphicRange" ',
                         'firstOccurrence="',
                         first_occ,
                         '" lastOccurrence="',
                         last_occ, '"/>'))
      
      # add to sr_count
      sr_count <- sr_count + 1
    }
    
    # add final taxonset call
    script <- c(script,
                paste0(add_spaces(16),
                       '</taxonset>'))
  }
  
  # add final tree call
  script <- c(script,
              paste0(add_spaces(12),
                     '</tree>'))
  
  # make prior text blocks
  prior_blocks_all <- prior_blocks(priors)
  
  # get each type of block
  prior_block <- prior_blocks_all$PRIORS
  param_block <- prior_blocks_all$PARAMS
  init_block <- prior_blocks_all$INITS
  
  # add parameter block
  script <- c(script,
              param_block,
              paste0(add_spaces(8), '</state>'))
  
  # start init
  script <- c(script,
              paste0(add_spaces(8),
                     '<init id="RandomTree.t:tree" initial="@Tree.t:tree" ',
                     'estimate="false" taxonset="@TaxonSet.',
                     names(partitions)[1],
                     '" ',
                     ifelse(model == "srfbd",
                            'spec="RandomSRangeTree" nodetype="sr.evolution.tree.SRNode"',
                            'spec="RandomTree"'),
                     '>'),
              paste0(add_spaces(12),
                     '<populationModel id="ConstantPopulation0.t:tree" ',
                     'spec="ConstantPopulation">'),
              paste0(add_spaces(16),
                     '<parameter id="randomPopSize.t:tree" name="popSize">1.0</parameter>'),
              paste0(add_spaces(12), '</populationModel>'))
  
  # if srfbd, add stratigraphic range refs here
  if (model == "srfbd") {
    # iterate through SRs
    for (i in 1:(sr_count - 1)) {
      # add to script
      script <- c(script,
                  paste0(add_spaces(12),
                         '<stratigraphicRange idref="r',
                         i,
                         '"/>'))
    }
  }
  
  # end init for tree
  script <- c(script,
              paste0(add_spaces(8), '</init>'),
              init_block)
  
  # start distribution block
  script <- c(script,
              paste0(add_spaces(8),
                     '<distribution id="posterior" spec="CompoundDistribution">'),
              paste0(add_spaces(12),
                     '<distribution id="prior" spec="CompoundDistribution">'),
              paste0(add_spaces(16),
                     '<distribution ',
                     ifelse(model == "srfbd",
                            'id="srfbd" spec="sr.speciation.SRangesBirthDeathModel" ',
                            'id="FBD" spec="sa.evolution.speciation.SABirthDeathModel" '),
                     'tree="@Tree.t:tree" conditionOnRhoSampling="true" ',
                     'origin="@originFBD.t:tree" ',
                     'diversificationRate="@diversificationRateFBD.t:tree" ',
                     'samplingProportion="@samplingProportionFBD.t:tree" ',
                     'turnover="@turnoverFBD.t:tree" ',
                     'rho="', extra_params["rho"], '" ',
                     'removalProbability="', extra_params["removalProbability"], '"/>'))
  
  # add priors
  script <- c(script, prior_block)
  
  # if model is fbds, need to also add uniform for fossil ages
  if (model == "fbds") {
    # number of specimens
    specimen_count <- 1
    
    # iterate through ranges
    for (i in 1:nrow(new_ranges)) {
      # get specimen
      specimen <- new_ranges$specimen[i]
      
      # if both are the same, skip
      if (new_ranges$min_age[i] == new_ranges$max_age[i]) {
        next
      }
      
      # add fossil age uniform
      script <- c(script, 
                  paste0(add_spaces(16),
                         '<distribution id="',
                         specimen, 'Set.prior" ',
                         'spec="sa.math.distributions.SAMRCAPrior" ',
                         'tipsonly="true" tree="@Tree.t:tree">'),
                  paste0(add_spaces(20),
                         '<taxonset id="',
                         specimen, 'Set" spec="TaxonSet">'),
                  paste0(add_spaces(24),
                         '<taxon id="',
                         specimen, '" spec="Taxon"/>'),
                  paste0(add_spaces(20), '</taxonset>'),
                  paste0(add_spaces(20),
                         '<Uniform id="FossilAges.', specimen_count,
                         '" lower="', new_ranges$min_age[i],
                         '" upper="', new_ranges$max_age[i],
                         '" name="distr"/>'),
                  paste0(add_spaces(16),
                         '</distribution>'))
      
      # add to specimen count
      specimen_count <- specimen_count + 1
    }
  }
  
  # close out prior
  script <- c(script,
              paste0(add_spaces(12),
                     '</distribution>'))
  
  # start likelihood
  script <- c(script,
              paste0(add_spaces(12),
                     '<distribution id="likelihood" spec="CompoundDistribution"',
                     ' useThreads="true">'))
  
  # make likelihood blocks for partitions
  lik_blocks <- likelihood_blocks(partitions)
  
  # add to script, and finish likelihood and posterior distributions
  script <- c(script, lik_blocks,
              paste0(add_spaces(12), '</distribution>'),
              paste0(add_spaces(8), '</distribution>'))
  
  # make base operators block
  op_block <- operator_block(operators, spaces = 8)
  
  # add dates sampler for either model
  if (model == "srfbd") {
    # start node date sampler
    op_block <- c(op_block,
                  paste0(add_spaces(8), '<operator ',
                         'id="SRNodeDateSampler" ',
                         'spec="sa.evolution.operators.SampledNodeDateRandomWalker" ',
                         'windowSize="1" tree="@Tree.t:tree" weight="10">'),
                  paste0(add_spaces(12), '<taxonset spec="TaxonSet">'))
    
    # start sampling date lines
    sampling_dates <- c()
    
    # iterate through ranges
    for (i in 1:nrow(new_ranges)) {
      # if it's certain, skip
      if (new_ranges$max_age[i] == new_ranges$min_age[i]) {
        next
      }
      # add taxa to block
      op_block <- c(op_block,
                    paste0(add_spaces(16),
                           '<taxon id="', new_ranges$specimen[i], 
                           '" spec="Taxon"/>'))
      
      # add to sampling dates
      sampling_dates <- c(sampling_dates,
                          paste0(add_spaces(12),
                                 '<samplingDates id="samplingDate_', i,
                                 '" spec="sa.evolution.tree.SamplingDate" ',
                                 'taxon="@',
                                 new_ranges$specimen[i],
                                 '" lower="',
                                 new_ranges$min_age[i],
                                 '" upper="',
                                 new_ranges$max_age[i],
                                 '"/>'))
    }
    
    # finish taxonset and add sampling dates
    op_block <- c(op_block,
                  paste0(add_spaces(12), '</taxonset>'),
                  sampling_dates,
                  paste0(add_spaces(8), '</operator>'))
  } else if (model == "fbds") {
    # iterate through ranges
    for (i in 1:nrow(new_ranges)) {
      # check if it is not a certain taxon
      if (new_ranges$max_age[i] != new_ranges$min_age[i]) {
        # add tip dates sampler
        op_block <- c(op_block,
                      paste0(add_spaces(8),
                             '<operator id="tipDatesSampler.', i, '" ',
                             'spec="sa.evolution.operators.SampledNodeDateRandomWalker" ',
                             'taxonset="@',
                             new_ranges$specimen[i],
                             'Set" tree="@Tree.t:tree" weight="1" ',
                             'windowSize="1"/>'))
      }
    }
  } else {
    stop("Model must be srfbd or fbds")
  }
  
  # add operator block to script
  script <- c(script, op_block)
  
  # start logger
  log_block <- c(paste0(add_spaces(8),
                        '<logger id="tracelog" spec="Logger" ',
                        'fileName="output/$(filebase).log" ',
                        'logEvery="', log_every, '" ',
                        'model="@posterior" sanitiseHeaders="true" sort="smart">'))
  
  # iterate through logs
  for (i in 1:length(logs)) {
    # get this log
    log <- logs[[i]]
    
    # if length is 1, just add with idref
    if (length(log) == 1) {
      log_block <- c(log_block,
                     paste0(add_spaces(12),
                            '<log idref="', log, '"/>'))
    } else {
      # if not, start this log line
      log_line <- paste0(add_spaces(12), '<log id="', names(logs)[i], '" ')
      
      # iterate through log
      for (j in 1:length(log)) {
        log_line <- paste0(log_line,
                           names(log)[j], '="',
                           log[[j]], '" ')
      }
      
      # delete last space
      substr(log_line, 1, nchar(log_line) - 1)
    
      # add end of line
      log_line <- paste0(log_line, '/>')
      
      # add to log block
      log_block <- c(log_block, log_line)
    }
  }
  
  # add ending to log_block
  log_block <- c(log_block,
                 paste0(add_spaces(8), '</logger>'))
  
  # add log_block to script, and add screen and tree log
  script <- c(script, log_block,
              paste0(add_spaces(8),
                     '<logger id="screenlog" spec="Logger" logEvery="',
                     log_every, '">'),
              paste0(add_spaces(12),
                     '<log idref="posterior"/>'),
              paste0(add_spaces(12),
                     '<log idref="likelihood"/>'),
              paste0(add_spaces(12),
                     '<log idref="prior"/>'),
              paste0(add_spaces(8), '</logger>'),
              paste0(add_spaces(8),
                     '<logger id="treelog.t:tree" spec="Logger" ',
                     'fileName="output/$(filebase).trees" ',
                     'logEvery="', log_every, '" mode="tree">'),
              paste0(add_spaces(12),
                     '<log id="TreeWithMetaDataLogger.t:tree" ',
                     'spec="',
                     ifelse(model == "srfbd",
                            'sr.evolution.tree.TreeWithMetadataLogger" ',
                            'beast.base.evolution.TreeWithMetaDataLogger" '),
                     'branchratemodel="@',
                     switch(partitions$morpho2$branchRateModel$model,
                       "ucld" = {
                         "RelaxedClock.c:morpho"
                       }
                     ),
                     '" tree="@Tree.t:tree"/>'),
              paste0(add_spaces(8), '</logger>'))
  
  # add final lines
  script <- c(script,
              paste0(add_spaces(8),
                     '<operatorschedule id="OperatorSchedule" spec="OperatorSchedule"/>'),
              paste0(add_spaces(4), '</run>'),
              '</beast>')
  
  # write script
  writeLines(script, out_file)
}

# write all scripts
write_n_scripts <- function(n_scripts, model, template, maps, 
                            priors, operators, logs,
                            extra_params,
                            ranges, mol, morpho,
                            partitions,
                            attach, one_sp,
                            gens, store_every, log_every,
                            out_dir, out_filename) {
  # iterate through number of scripts
  for (i in 1:n_scripts) {
    # out file
    out_file_n <- paste0(out_dir, out_filename, "_", i, ".xml")
    
    # write script
    write_script(model, template, maps, 
                 priors, operators, logs,
                 extra_params,
                 ranges, mol, morpho,
                 partitions,
                 attach, one_sp,
                 gens, store_every, log_every,
                 out_file_n)
  }
}