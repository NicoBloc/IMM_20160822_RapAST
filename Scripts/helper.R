#############
# Load libs #
#############

library(gdata) # read.xls()
library(ggplot2)
library(reshape2)  # melt()
library(scales)  # pretty_breaks()
library(plyr)  # revalue()
library(dplyr)
library(RColorBrewer)
library(printr)
library(knitr)
library(cowplot)
library(pander)
library(tidyr)
library(pROC)

##########################
# Configure ggplot-theme #
##########################

scale_colour_discrete <- function(...) scale_colour_brewer(..., palette="Set1")
scale_fill_discrete <- function(...) scale_fill_brewer(..., palette="Set1")

theme_set(theme_bw(base_size=10) + theme(panel.grid.major=element_line(size=0.25), panel.grid.minor=element_blank(), plot.title=element_text(size=10, face="bold")))

#########################
# Define some functions #
#########################

myQ <- function(x) {
  r <- quantile(x, probs = c(0, 0.05, 0.5, 0.95, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

addBoxplot <- function(plt) {
  plt + stat_summary(fun.data=myQ, geom='boxplot', position='dodge', aes(width=0.5))
}

masterplot <- function(d, wt, allNonWt) {
  # Duplicate wt-data, we use the "phenotypes" for facetting and "isWt" for filling:
  dTmpDupl <- d %>% select(phenotypes, hDiscr, diameter) %>% filter(phenotypes != wt) %>% mutate(isWt="non-wild-type")
  dTmpWt <- d %>% select(phenotypes, hDiscr, diameter) %>% filter(phenotypes == wt) %>% select(-phenotypes) %>% mutate(isWt="wild-type")
  for (pt in allNonWt) {
    dTmpDupl <- rbind(dTmpDupl, data.frame(dTmpWt, phenotypes=pt))
  }
  
  addBoxplot(ggplot(dTmpDupl, aes(x=hDiscr, y=diameter, fill=isWt)) + facet_wrap(~ phenotypes) + scale_fill_discrete("", drop=FALSE)) %>% polishPlot
}

lineplotHombachPanel <- function(d) {
  polishPanel(ggplot(d, aes(x=h, y=diameter, group=checkin, colour=isWt)) + geom_point() + geom_line() + scale_x_continuous(breaks=c(6, 8, 12, 18), minor=NULL))
}

lineplotHombachMaster <- function(d) {
  pltWt <- lineplotHombachPanel(dTmp[dTmp$phenotypes == wt, ])
  
  pltNonWt <- list()
  pltWtAndNonWt <- list()
  for (i in 1:(nlevels(dTmp$phenotypes) - 1)) {
    nonWt <- allNonWt[i]
    pltNonWt[[i]] <- lineplotHombachPanel(dTmp[dTmp$phenotypes == nonWt, ])
    pltWtAndNonWt[[i]] <- lineplotHombachPanel(dTmp[dTmp$phenotypes %in% c(wt, nonWt), ])
  }
  
  allPlots <- list()
  for (i in 1:(3 * (nlevels(dTmp$phenotypes) - 1))) {
    if (i %% 3 == 1) {
      allPlots[[i]] = pltWt
    } else if (i %% 3 == 2) {
      allPlots[[i]] = pltNonWt[[i %/% 3 + 1]]
    } else if (i %% 3 == 0) {
      allPlots[[i]] = pltWtAndNonWt[[i %/% 3]]
    }
  }
  
  labels <- as.vector(t(array(c(paste0('(', letters[1:nlevels(dTmp$phenotypes)], ') ', c(allNonWt, 'all phenotypes combined')), rep('', 2 * nlevels(dTmp$phenotypes))), c(nlevels(dTmp$phenotypes), 3))))
  
  do.call("plot_grid", c(allPlots, ncol=3, labels=list(labels), label_size=6, hjust=0, vjust=0.8))
}

lineplotWFacets <- function(d) {
  polishPlot(ggplot(dTmp, aes(x=h, y=diameter, group=checkin, colour=isWt)) + geom_line() + geom_point() + facet_wrap(~ phenotypes) + scale_x_continuous(breaks=c(6, 8, 12, 18), minor=NULL) + scale_colour_discrete(""))
}

myBoxplot <- function(d) {
  ggplot(d, aes(x=hDiscr, y=diameter, fill=phenotypes)) %>% addBoxplot %>% polishPanel
}

polishPanel <- function(plt) {
  polishPlot(plt + guides(colour=FALSE) + scale_colour_discrete(drop=FALSE) + scale_fill_discrete(name=""))
}

polishPlot <- function(plt) {
  plt + xlab('time / h') + ylab('diameter / mm') + coord_cartesian(ylim=c(6, 40)) + theme(legend.position="top", strip.text.x=element_text(face='bold'))
}

plotOverviewAntibiotic <- function(d, type="boxplot", groupRarePhenotypes=FALSE) {
  stopifnot(type %in% c("boxplot", "lineplot"))
  importantSpecies <- c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae", "Staphylococcus epidermidis", "Staphylococcus aureus")
  
  dSummary <- d[d$organism %in% importantSpecies, ]
  dSummary$organism <- factor(dSummary$organism, levels=importantSpecies)
  
  if (groupRarePhenotypes) {
    # Group rare phenotypes together:
    dSummary$phenotypes <- factor(dSummary$phenotypes)
    rarePhenotypes <- names(which(table(dSummary$phenotypes) < 20 * 4))
    if (length(rarePhenotypes) > 1) {
      repl = rep('others', length(rarePhenotypes))
      names(repl) = rarePhenotypes
      dSummary$phenotypes <- revalue(dSummary$phenotypes, replace=repl)
    }
  }
  
  plt <- if (type == "boxplot") {
    addBoxplot(ggplot(dSummary, aes(x=hDiscr, y=diameter, fill=organism)) + scale_fill_discrete(drop=FALSE))
  } else {
    ggplot(dSummary, aes(x=h, y=diameter, group=checkin, colour=organism)) + geom_point() + geom_line() + scale_colour_discrete(drop=FALSE) + scale_x_continuous(breaks=c(6, 8, 12, 18))
  }
  plt + facet_wrap(~ phenotypes) + theme(legend.position="top") + xlab('time/h') + ylab('diameter/mm')
}


##############
# Masterloop #
##############

masterLoop <- function(manuscript) {
  stopifnot(manuscript %in% 2:3)
  
  showAllPlots <- FALSE  # not tested whether script works if this is set to TRUE
  groupRarePhenotypes <- FALSE  # not tested whether script works if this is set to TRUE
  
  cat('\n\n\\pagebreak\n')
  
  for (ab in levels(dRaw$antibiotic)) {
    pandoc.header(paste(sep=', ', ab), 1)
    
    d <- filter(dRaw, antibiotic == ab) %>% select(-antibiotic)
    
    if (showAllPlots) {
      print(plotOverviewAntibiotic(d, groupRarePhenotypes=groupRarePhenotypes))  # Summary figure for individual antibiotics
      cat('\n\n\\pagebreak\n')
    }
    
    # Group Enterobacteriaceae with few data together:
    smallEnterobacteriaceae <- setdiff(names(which(table(d$organism) < 20 * 4)), c("Staphylococcus aureus", "Staphylococcus epidermidis"))
    repl <- rep('other Enterobacteriaceae', length(smallEnterobacteriaceae))
    names(repl) <- smallEnterobacteriaceae
    d$organism <- revalue(d$organism, replace=repl)
    
    for (sp in levels(d$organism)) {
      dTmp <- filter(d, organism == sp) %>% select(-organism)
      if (nrow(dTmp) == 0) {
        next
      }

      pandoc.header(paste0(ab, ", *", sp, "*"), 2)

      dTmp$phenotypes <- factor(dTmp$phenotypes)
      wt <- levels(dTmp$phenotypes)[c(grep(" wild-type", levels(dTmp$phenotypes)), which("Wild-type" == levels(dTmp$phenotypes)), which("wild-type" == levels(dTmp$phenotypes)), which("Wild type" == levels(dTmp$phenotypes)), which("wild type" == levels(dTmp$phenotypes)))]
      if (length(wt) == 0 && manuscript == 2) {
        cat('\n\n', paste(sep=', ', ab, sp), ': No data for wild type available.\n\n')
        next
      }

      if (groupRarePhenotypes && manuscript == 2) {
        stopifnot(length(wt) == 1)
        # Group rare phenotypes together:
        rarePhenotypes <- setdiff(names(which(table(dTmp$phenotypes) < 5 * 4)), wt)
        if (length(rarePhenotypes) > 1) {
          repl <- rep('others', length(rarePhenotypes))
          names(repl) <- rarePhenotypes
          dTmp$phenotypes <- revalue(dTmp$phenotypes, replace=repl)
        }
      }

      allNonWt <- setdiff(levels(dTmp$phenotypes), wt)

      isWtLvls <- c("non-wild-type", "wild-type")
      if (length(wt) == 1) {
        dTmp$isWt <- factor(isWtLvls[(dTmp$phenotypes == wt) + 1], levels=isWtLvls)
      } else if (length(wt) == 0) {
        dTmp$isWt <- factor(isWtLvls[1], levels=isWtLvls)
      }

      #######################################
      # Number of data points per phenotype #
      #######################################

      cat('\n\n\n\n')
      dTmp %>% group_by(phenotypes) %>% summarise(n=n_distinct(checkin)) %>% rename(phenotype=phenotypes) %>% kable %>% print
      cat('\n\n\n\n')

      ##############
      # Masterplot #
      ##############

      if (manuscript == 2) {
        if (length(allNonWt) == 0) {
          addBoxplot(ggplot(dTmp, aes(x=hDiscr, y=diameter, fill=isWt)) + facet_wrap(~ phenotypes) + scale_fill_discrete("", drop=FALSE)) %>% polishPlot %>% print
          cat(sep='', '\n\n\\emph{Caption or footnote for ', paste(sep=' / ', ab, sp), ', to be written}.\n\n')
          cat('\n\n\\pagebreak\n')
          next
        }
        
        masterplot(dTmp, wt=wt, allNonWt=allNonWt) %>% print
      } else if (manuscript == 3) {
        myBoxplot(dTmp) %>% print
      }

      ######################
      # Caption / footnote #
      ######################

      cat(sep='', '\n\n\\emph{Caption or footnote for ', paste(sep=' / ', ab, sp), ', to be written}. Boxes range from the 5th to the 95th percentiles, whiskers cover the whole range of the data, and thick lines indicate medians.\n\n')

      ####################
      # Additional plots #
      ####################

      if (showAllPlots) {
        cat('\n\n\\pagebreak\n')
        print(lineplotHombachMaster(dTmp))  # lineplot requested by M. Hombach and M. Jetter
        cat('\n\n\\pagebreak\n')
        print(lineplotWFacets(dTmp))  # simple lineplot

        # More accessible figure based on boxplots:
        bxplt <- list()
        for (i in 1:(nlevels(dTmp$phenotypes) - 1)) {
          nonWt <- allNonWt[i]
          bxplt[[i]] <- myBoxplot(dTmp[dTmp$phenotypes %in% c(wt, nonWt), ]) + ggtitle(nonWt)
        }
        cat('\n\n\\pagebreak\n')
        print(do.call("plot_grid", c(bxplt)))
      }

      cat('\n\n\\pagebreak\n')
    }
  }
}

##########
# Tab. 2 #
##########

tab2 <- function() {
  cutoff <- 1  # adjust manually in figure legend below!!
  discrDeltaDiam <- function(d) {
    factor(ifelse(d < -cutoff, 1, ifelse(d > cutoff, 3, 2)), levels=1:3)
  }
  
  dTab2 <- dRaw %>% filter(phenotypes == "wild-type") %>%
    group_by(organism, antibiotic, h) %>% summarise(q05Diam=quantile(diameter, 0.05)) %>%
    mutate(deltaDiam=(q05Diam[h == 18] - q05Diam[h == 6])) %>% mutate(deltaDiscr=discrDeltaDiam(deltaDiam))
  
  print(ggplot(dTab2, aes(x=h, y=q05Diam, colour=deltaDiscr)) +
          geom_line() + geom_point() +
          facet_grid(antibiotic ~ organism) +
          theme(strip.text.y=element_text(angle=0), legend.position="top") +
          scale_x_continuous(breaks=c(6, 8, 12, 18)) + scale_y_continuous(breaks=pretty_breaks(n=3)) +
          xlab("time / h") + ylab("5th percentile of diameter for wild-type strains / mm") +
          scale_colour_manual(name=expression(Delta==d[18*h]-d[6*h]), values=c("orangered3", "skyblue4", "orange"),
                              labels=c(expression("decrease "*(Delta < -1*mm)),
                                       expression("stable "*(-1*mm<=paste(Delta<=1*mm))),
                                       expression("increase "*(Delta>1*mm))), drop=FALSE))
  
  dTab2 %>% select(organism, antibiotic, deltaDiscr) %>% unique %>% spread(organism, deltaDiscr) %>% write.table("tab2.csv", sep=";", row.names=FALSE)
  ggsave(paste(sep="", "Figs/fig3.png"), width=18, height=22, units="cm", dpi=1200, limitsize=TRUE)
}