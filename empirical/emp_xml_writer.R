###########################################
##   Write scripts for Canidae analyses  ##
## Chapter 2 - Bruno do Rosario Petrucci ##
###########################################

###
# source xml writer

# get base directory
base_dir <- "/Users/petrucci/Documents/research/srfbd_canidae/"
emp_dir <- paste0(base_dir, "empirical/")

# source xml_writer.R
source(paste0(base_dir, "xml_writer.R"))

###
# write scripts

# templates
maps <- readLines(paste0(base_dir, "maps.xml"))

# function to write scripts for a given setup
scripts_setup <- function(n_scripts, model, one_sp, attach) {
  # template
  template <- readLines(paste0(base_dir, "template_", model, ".xml"))
  
  # data
  ranges <- read.delim(paste0(emp_dir, "data/ranges.tsv"))
  mol <- read.nexus.data(paste0(emp_dir, "data/mol.nex"))
  morpho <- read.nexus.data(paste0(emp_dir, "data/morpho.nex"))
  
  # rewrite polymorphisms
  morpho <- rewrite_polymorphisms(morpho)
  
  # get morpho partitions data
  morpho_partitions <- morpho_partitions(morpho, ascertain = TRUE)
  
  # start partitions vector with molecular partitions
  partitions <- list(nDNA = list(
    data = "mol",
    filter = "1-14742",
    ascertained = "false",
    dataType = NULL,
    siteModel = list(
      id = "SiteModel.s:nDNA",
      gammaCategoryCount = 4,
      shape = "@gammaShape.s:nDNA",
      mutationRate = 1,
      proportionInvariant = 0,
      substModel = list(
        model = "hky",
        kappa = "@kappa.s:nDNA",
        frequencies = "@freqParameter.s:nDNA"
      )),
    branchRateModel = list(
      id = "RelaxedClock.c:nDNA",
      model = "ucld",
      clock.rate = 1,
      rateCategories = "@rateCategories.c:nDNA",
      dist = "LogNormal",
      M = 1,
      S = "@ucldStdev.c:nDNA"
    ),
    threaded = TRUE),
    mt12DNA = list(
      data = "mol",
      filter = paste0("14743-15426\\3,15427-16570\\3,",
                      "14744-15426\\3,15428-16570\\3"),
      ascertained = "false",
      dataType = NULL,
      siteModel = list(
        id = "SiteModel.s:mt12DNA",
        gammaCategoryCount = 4,
        shape = "@gammaShape.s:mt12DNA",
        mutationRate = 1,
        proportionInvariant = 0,
        substModel = list(
          model = "hky",
          kappa = "@kappa.s:mt12DNA",
          frequencies = "@freqParameter.s:mt12DNA"
        )),
      branchRateModel = list(
        id = "RelaxedClock.c:mt12DNA",
        model = "ucld",
        clock.rate = "@ucldMean.c:mt12DNA",
        rateCategories = "@rateCategories.c:mt12DNA",
        dist = "LogNormal",
        M = 1,
        S = "@ucldStdev.c:mt12DNA"
      ),
      threaded = TRUE),
    mt3DNA = list(
      data = "mol",
      filter = "14745-15426\\3,15429-16570\\3",
      ascertained = "false",
      dataType = NULL,
      siteModel = list(
        id = "SiteModel.s:mt3DNA",
        gammaCategoryCount = 4,
        shape = "@gammaShape.s:mt3DNA",
        mutationRate = 1,
        proportionInvariant = 0,
        substModel = list(
          model = "hky",
          kappa = "@kappa.s:mt3DNA",
          frequencies = "@freqParameter.s:mt3DNA"
        )),
      branchRateModel = list(
        id = "RelaxedClock.c:mt3DNA",
        model = "ucld",
        clock.rate = "@ucldMean.c:mt3DNA",
        rateCategories = "@rateCategories.c:mt3DNA",
        dist = "LogNormal",
        M = 1,
        S = "@ucldStdev.c:mt3DNA"
      ),
      threaded = TRUE))
  
  # add partitions for morpho
  for (i in 1:length(morpho_partitions)) {
    # get this partition
    part <- morpho_partitions[[i]]
    
    # get the name
    part_name <- names(morpho_partitions)[i]
    
    # add partition
    partitions[[part_name]] <- list(
      data = "morpho",
      filter = part$filter,
      ascertained = part$ascertain,
      dataType = list(
        type = "StandardData",
        nrOfStates = as.numeric(sub('.*([0-9]+)', '\\1', part_name))
      ),
      excludefrom = part$excludefrom,
      excludeto = part$excludeto,
      siteModel = list(
        id = paste0("SiteModel.s:", part_name),
        gammaCategoryCount = 0,
        shape = 1,
        mutationRate = paste0("@mutationRate.s:", part_name),
        proportionInvariant = 0,
        substModel = list(
          model = "LewisMK"
        )
      ),
      threaded = FALSE)
    
    # branch rate model only for the first partition
    if (i == 1) {
      partitions[[part_name]]$branchRateModel <- list(
        id = "RelaxedClock.c:morpho",
        model = "ucld",
        clock.rate = "@ucldMean.c:morpho",
        rateCategories = "@rateCategories.c:morpho",
        dist = "LogNormal",
        M = 1, 
        S = "@ucldStdev.c:morpho"
      )
    }
  }
  
  # parameter list (with priors)
  priors = list(
    `mutationRate.s:morpho2` = list(
      spec = "parameter.RealParameter",
      value = 1
    ),
    `mutationRate.s:morpho3` = list(
      spec = "parameter.RealParameter",
      value = 1
    ),
    `mutationRate.s:morpho4` = list(
      spec = "parameter.RealParameter",
      value = 1
    ),
    `mutationRate.s:morpho5` = list(
      spec = "parameter.RealParameter",
      value = 1
    ),
    `diversificationRateFBD.t:tree` = list(
      dist = "Exponential",
      lower = 0,
      mean = 0.0508
    ),
    `turnoverFBD.t:tree` = list(
      dist = "Beta",
      lower = 0,
      upper = 1,
      alpha = 2,
      beta = 1
    ),
    `samplingProportionFBD.t:tree` = list(
      dist = "Beta",
      lower = 0,
      upper = 1,
      alpha = 3.1,
      beta = 6.9
    ),
    `originFBD.t:tree` = list(
      dist = "Uniform",
      lower = 37,
      upper = 50
    ),
    `gammaShape.s:nDNA` = list(
      dist = "Exponential",
      lower = 0.1,
      mean = 1
    ),
    `gammaShape.s:mt12DNA` = list(
      dist = "Exponential",
      lower = 0.1,
      mean = 1
    ),
    `gammaShape.s:mt3DNA` = list(
      dist = "Exponential",
      lower = 0.1,
      mean = 1
    ),
    `ucldMean.c:mt12DNA` = list(
      dist = "LogNormal",
      lower = 0,
      M = -3.5,
      S = 1
    ),
    `ucldMean.c:mt3DNA` = list(
      dist = "LogNormal",
      lower = 0,
      M = -3.5,
      S = 1
    ),
    `ucldMean.c:morpho` = list(
      dist = "LogNormal",
      lower = 0,
      M = -3.5,
      S = 1
    ),
    `ucldStdev.c:nDNA` = list(
      dist = "Gamma",
      lower = 0,
      alpha = 0.5396,
      beta = 0.3819
    ),
    `ucldStdev.c:mt12DNA` = list(
      dist = "Gamma",
      lower = 0,
      alpha = 0.5396,
      beta = 0.3819
    ),
    `ucldStdev.c:mt3DNA` = list(
      dist = "Gamma",
      lower = 0,
      alpha = 0.5396,
      beta = 0.3819
    ),
    `ucldStdev.c:morpho` = list(
      dist = "Gamma",
      lower = 0,
      alpha = 0.5396,
      beta = 0.3819
    ),
    `kappa.s:nDNA` = list(
      dist = "LogNormal",
      lower = 0,
      M = 1,
      S = 1.25
    ),
    `kappa.s:mt12DNA` = list(
      dist = "LogNormal",
      lower = 0,
      M = 1,
      S = 1.25
    ),
    `kappa.s:mt3DNA` = list(
      dist = "LogNormal",
      lower = 0,
      M = 1,
      S = 1.25
    ),
    `freqParameter.s:nDNA` = list(
      dist = "Dirichlet",
      dimension = 4,
      lower = 0,
      upper = 1,
      alpha = c(4, 4, 4, 4)
    ),
    `freqParameter.s:mt12DNA` = list(
      dist = "Dirichlet",
      dimension = 4,
      lower = 0,
      upper = 1,
      alpha = c(4, 4, 4, 4)
    ),
    `freqParameter.s:mt3DNA` = list(
      dist = "Dirichlet",
      dimension = 4,
      lower = 0,
      upper = 1,
      alpha = c(4, 4, 4, 4)
    )
  )
  
  # if srfbd or fbds with one sample, just use number of taxa
  if (model == "srfbd" || (model == "fbds" && one_sp)) {
    priors <- append(priors, list(
      `rateCategories.c:nDNA` = list(
        spec = "parameter.IntegerParameter",
        dimension = 2 * length(unique(ranges$taxon)) - 2,
        value = 1
      ),
      `rateCategories.c:mt12DNA` = list(
        spec = "parameter.IntegerParameter",
        dimension = 2 * length(unique(ranges$taxon)) - 2,
        value = 1
      ),
      `rateCategories.c:mt3DNA` = list(
        spec = "parameter.IntegerParameter",
        dimension = 2 * length(unique(ranges$taxon)) - 2,
        value = 1
      ),
      `rateCategories.c:morpho` = list(
        spec = "parameter.IntegerParameter",
        dimension = 2 * length(unique(ranges$taxon)) - 2,
        value = 1
      )
    ))
  } else {
    # otherwise, use number of samples
    priors <- append(priors, list(
      `rateCategories.c:nDNA` = list(
        spec = "parameter.IntegerParameter",
        dimension = 2 * length(unique(ranges$taxon)) - 
          sum(ranges$fa_max == ranges$la_max & ranges$fa_min == ranges$la_min) - 2,
        value = 1
      ),
      `rateCategories.c:mt12DNA` = list(
        spec = "parameter.IntegerParameter",
        dimension = 2 * length(unique(ranges$taxon)) - 
          sum(ranges$fa_max == ranges$la_max & ranges$fa_min == ranges$la_min) - 2,
        value = 1
      ),
      `rateCategories.c:mt3DNA` = list(
        spec = "parameter.IntegerParameter",
        dimension = 2 * length(unique(ranges$taxon)) - 
          sum(ranges$fa_max == ranges$la_max & ranges$fa_min == ranges$la_min) - 2,
        value = 1
      ),
      `rateCategories.c:morpho` = list(
        spec = "parameter.IntegerParameter",
        dimension = 2 * length(unique(ranges$taxon)) - 
          sum(ranges$fa_max == ranges$la_max & ranges$fa_min == ranges$la_min) - 2,
        value = 1
      )
    ))
  }
  
  # extra parameters we might need around the script
  extra_params <- c(rho = 34/36, removalProbability = 0)
  
  # operators
  operators <- list(
    FixMeanMutationRatesOperator = list(
      spec = "operator.kernel.BactrianDeltaExchangeOperator",
      delta = 0.75,
      weight = 2,
      parameter = setNames(as.list(paste0("mutationRate.s:morpho", 2:5)),
                           rep("parameter", 4)),
      weightvector = list(
        `weightparameter` = list(
          spec = "parameter.IntegerParameter",
          dimension = 4,
          estimate = "false",
          lower = 0,
          upper = 0,
          value = "71 35 8 4")
      )
    ))
  operators <- c(operators,
                 part_operators("nDNA", partitions),
                 part_operators("mt12DNA", partitions),
                 part_operators("mt3DNA", partitions),
                 part_operators("morpho", partitions))
  operators <- append(operators, list(
    `originScalerFBD.t:tree` = list(
      spec = "ScaleOperator",
      parameter = "@originFBD.t:tree",
      weight = 3
    ),
    `divRateScalerFBD.t:tree` = list(
      spec = "ScaleOperator",
      parameter = "@diversificationRateFBD.t:tree",
      weight = 10
    ),
    `turnoverScalerFBD.t:tree` = list(
      spec = "ScaleOperator",
      parameter = "@turnoverFBD.t:tree",
      weight = 10
    ),
    `samplingPScalerFBD.t:tree` = list(
      spec = "ScaleOperator",
      parameter = "@samplingProportionFBD.t:tree",
      weight = 10
    )))
  
  
  # logs
  logs <- c(list(
    `posterior` = "posterior",
    `likelihood` = "likelihood",
    `prior` = "prior",
    `diversificationRateFBD.t:tree` = "diversificationRateFBD.t:tree",
    `turnoverFBD.t:tree` = "turnoverFBD.t:tree",
    `samplingProportionFBD.t:tree` = "samplingProportionFBD.t:tree",
    `originFBD.t:tree` = "originFBD.t:tree"),
    make_logs("treeLikelihood.", names(partitions)),
    make_logs("mutationRate.s:", paste0("morpho", 2:5)),
    make_logs("gammaShape.s:", c("nDNA", "mt12DNA", "mt3DNA")),
    make_logs("ucldMean.c:", c("mt12DNA", "mt3DNA", "morpho")),
    make_logs("ucldStdev.c:", c("nDNA", "mt12DNA", "mt3DNA", "morpho")),
    make_logs("rate.c:", c("nDNA", "mt12DNA", "mt3DNA", "morpho")),
    make_logs("kappa.s:", c("nDNA", "mt12DNA", "mt3DNA")),
    make_logs("freqParameter.s:", c("nDNA", "mt12DNA", "mt3DNA"))
  )
  
  # check model
  if (model == "srfbd") {
    # srfbd operators
    operators <- append(operators, list(
      `SRTreeRootScaler` = list(
        spec = "SAScaleOperator",
        rootOnly = "true",
        scaleFactor = 0.95,
        tree = "@Tree.t:tree",
        weight = 1
      ),
      `SRWilsonBalding` = list(
        spec = "SRWilsonBalding",
        tree = "@Tree.t:tree",
        weight = 20
      ),
      `LeftRightChildSwap` = list(
        spec = "LeftRightChildSwap",
        tree = "@Tree.t:tree",
        weight = 3
      ),
      `LeafToSampledAncestorJump` = list(
        spec = "SRLeafToSampledAncestorJump",
        tree = "@Tree.t:tree",
        weight = 10
      ),
      `SRUniformOperator` = list(
        spec = "SAUniform",
        tree = "@Tree.t:tree",
        weight = 20
      ),
      `SRTreeScaler` = list(
        spec = "SAScaleOperator",
        scaleFactor = 0.95,
        tree = "@Tree.t:tree",
        weight = 3
      )
    ))
    
    # srfbd logs
    logs <- c(logs, list(
      `srfbd` = "srfbd",
      `SACountFBD.t:tree` = list(
        spec = "sr.evolution.tree.SampledAncestorLogger",
        tree = "@Tree.t:tree"
      ),
      `treeHeight.t:tree` = list(
        spec = "beast.base.evolution.tree.TreeHeightLogger",
        tree = "@Tree.t:tree"
      )
    ))
  } else {
    # fbds operators
    operators <- append(operators, list(
      `SATreeRootScaler` = list(
        spec = "sa.evolution.operators.SAScaleOperator",
        rootOnly = "true",
        scaleFactor = 0.95,
        tree = "@Tree.t:tree",
        weight = 1
      ),
      `SAWilsonBalding` = list(
        spec = "sa.evolution.operators.SAWilsonBalding",
        tree = "@Tree.t:tree",
        weight = 10
      ),
      `LeafToSampledAncestorJump` = list(
        spec = "sa.evolution.operators.LeafToSampledAncestorJump",
        tree = "@Tree.t:tree",
        weight = 10
      ),
      `SAUniformOperator` = list(
        spec = "sa.evolution.operators.SAUniform",
        tree = "@Tree.t:tree",
        weight = 20
      ),
      `SATreeScaler` = list(
        spec = "sa.evolution.operators.SAScaleOperator",
        scaleFactor = 0.95,
        tree = "@Tree.t:tree",
        weight = 3
      ),
      `SAWide` = list(
        spec = "sa.evolution.operators.SAExchange",
        isNarrow = "false",
        tree = "@Tree.t:tree",
        weight = 10
      ),
      `SANarrow` = list(
        spec = "sa.evolution.operators.SAExchange",
        tree = "@Tree.t:tree",
        weight = 10
      )
    ))
    
    # fbds logs
    logs <- c(logs, list(
      `FBD` = "FBD",
      `SACountFBD.t:tree` = list(
        spec = "sa.evolution.tree.SampledAncestorLogger",
        tree = "@Tree.t:tree"
      ),
      `treeHeight.t:tree` = list(
        spec = "beast.base.evolution.tree.TreeStatLogger",
        tree = "@Tree.t:tree"
      )
    ))
  }
  
  # get script name
  script_name <- paste0(model, "_",
                        ifelse(model == "srfbd", attach, ""),
                        ifelse(model == "fbds", 
                               ifelse(one_sp, "one", "two"), ""))
  
  # write scripts
  write_n_scripts(n_scripts, model, template, maps,
                  priors, operators, logs,
                  extra_params, ranges, mol, morpho,
                  partitions, attach, one_sp,
                  1000000000, 1000000, 1000000,
                  paste0(emp_dir, "scripts/"),
                  script_name)
  
}

# set script number
n_scripts <- 1

# write scripts for all combinations
scripts_setup(n_scripts, "srfbd", FALSE, "both")
scripts_setup(n_scripts, "srfbd", FALSE, "first")
scripts_setup(n_scripts, "fbds", FALSE, "both")
scripts_setup(n_scripts, "fbds", TRUE, "both")