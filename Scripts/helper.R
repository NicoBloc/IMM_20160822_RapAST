#############
# Load libs #
#############

library(multcomp)  # glht
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
library(gtable)
library(gridExtra)
library(grid)
library(prevalence)  # propCI()

##########################
# Configure ggplot-theme #
##########################

scale_colour_discrete <- function(...) scale_colour_brewer(..., palette="Set1")
scale_fill_discrete <- function(...) scale_fill_brewer(..., palette="Set1")

theme_set(theme_bw(base_size=10) + theme(panel.grid.major=element_line(size=0.25), panel.grid.minor=element_blank()))

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
  
  addBoxplot(ggplot(dTmpDupl, aes(x=hDiscr, y=diameter, fill=isWt)) + facet_grid(~ phenotypes) + scale_fill_discrete("", drop=FALSE) + theme(aspect.ratio=0.45 + 0.5 * length(allNonWt) * 0.9 / 7)) %>% polishPlot
}

aucBreaks <- function(dAUC) {
  # breaks for y-axis of aucPlot()
  a <- min(0.95, 0.1 * ceiling(10 * min(dAUC$AUC, dAUC$lower, na.rm=TRUE)))
  if (a == 0.95) {
    c(0.95, 1)
  } else if (a >= 0.6) {
    seq(from=a, to=1, by=0.1)
  } else {
    c(a, (a + 1) / 2, 1)
  }
}

aucPlot <- function(d, manuscript, wt=NULL, allNonWt=NULL) {
  dAUC <- data.frame()
  for (hh in unique(d$h)) {
    if (manuscript == 2) {
      for (pt in allNonWt) {
        idx <- d$h == hh & d$phenotypes %in% c(pt, wt)
        auc <- auc(d$isWt[idx], d$diameter[idx])
        ci <- ci.auc(auc)
        dAUC <- rbind(dAUC, data.frame(phenotypes=pt, h=hh, AUC=auc, lower=ci[1], upper=ci[3]))
      }
    } else if (manuscript == 3) {
      idx <- d$h == hh
      auc <- auc(d$isWt[idx], d$diameter[idx])
      ci <- ci.auc(auc)
      dAUC <- rbind(dAUC, data.frame(h=hh, AUC=auc, lower=ci[1], upper=ci[3]))
    }
  }
  
  c1 <- "black"
  c2 <- "grey60"
  
  plt <- ggplot(dAUC, aes(x=h, y=AUC, ymin=lower, ymax=upper)) + geom_line(colour=c2) + geom_errorbar(width=1.5, colour=c1) + geom_point(colour=c1) + xlab('time / h') + ylab('AUC') + scale_x_continuous(breaks=c(6, 8, 12, 18)) + scale_y_continuous(breaks=aucBreaks(dAUC))  + coord_cartesian(ylim=c(min(0.95, dAUC$AUC, dAUC$lower, na.rm=TRUE), 1)) # + scale_y_continuous(breaks=c(0.5, 0.75, 1))
  
  if (manuscript == 2) {
    plt <- plt + facet_grid(~ phenotypes) + theme(strip.text.x=element_text(face='bold'), aspect.ratio=(0.23 + 0.5 * length(allNonWt) * 0.46 / 7))  # , plot.background=element_rect(fill="grey90")
  } else if (manuscript == 3) {
    plt <- plt + theme(aspect.ratio=(0.5))
  }
  return(plt)
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
  polishPanel((ggplot(d, aes(x=hDiscr, y=diameter, fill=phenotypes)) %>% addBoxplot) + theme(aspect.ratio=(0.7)))
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

alignPlots <- function(p1, p2) {
  # Align two plots
  # http://stackoverflow.com/questions/26159495/align-multiple-ggplot-graphs-with-and-without-legends
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  grob_combined <- gtable:::rbind_gtable(g1, g2, "last")
  grid.newpage()
  grid.draw(grob_combined)
}

##############
# Masterloop #
##############

masterLoop <- function(manuscript) {
  stopifnot(manuscript %in% 2:3)
  
  showAllPlots <- FALSE  # not tested whether script works if this is set to TRUE
  groupRarePhenotypes <- FALSE  # not tested whether script works if this is set to TRUE
  
  cat('\n\n\\pagebreak\n')
  
  for (ab in levels(dRaw$antibiotic)) {  # [c(1, 4, 7, 24)]
    pandoc.header(paste(sep=', ', ab), 1)
    
    d <- filter(dRaw, antibiotic == ab) %>% select(-antibiotic)
    # Some modifications requested by Böttger (02.12.2016):
    if (manuscript == 2) {
      # d <- d %>% filter(!(organism == "S. epidermidis" & h == 6))
    } else if (manuscript == 3) {
      d <- d %>% filter(!(organism == "P. aeruginosa" & h == 6))
    }
    
    if (showAllPlots) {
      print(plotOverviewAntibiotic(d, groupRarePhenotypes=groupRarePhenotypes))  # Summary figure for individual antibiotics
      cat('\n\n\\pagebreak\n')
    }
    
    # # Group Enterobacteriaceae with few data together:
    # smallEnterobacteriaceae <- setdiff(names(which(table(d$organism) < 20 * 4)), c("Staphylococcus aureus", "Staphylococcus epidermidis"))
    # repl <- rep('other Enterobacteriaceae', length(smallEnterobacteriaceae))
    # names(repl) <- smallEnterobacteriaceae
    # d$organism <- revalue(d$organism, replace=repl)
    
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
      # dTmp %>% group_by(phenotypes) %>% summarise(n=n_distinct(checkin)) %>% rename(phenotype=phenotypes) %>% kable %>% print
      # dTmp %>% group_by(phenotypes, hDiscr) %>% summarise(n=n_distinct(checkin), readable=(100 * sum(!is.na(diameter)) / n)) %>% mutate(tmp=factor(paste0("r", hDiscr), levels=paste0("r", levels(dTmp$hDiscr)))) %>% select(-hDiscr) %>% spread(tmp, readable) %>% kable(digits=0) %>% print
      cat('\n\n\n\n')
      
      tbl <- dTmp %>% group_by(phenotypes, hDiscr) %>% summarise(n=n_distinct(checkin), readable=round(100 * sum(!is.na(diameter)) / n)) %>% mutate(tmp=factor(paste0("r", hDiscr), levels=paste0("r", levels(dTmp$hDiscr)))) %>% select(-hDiscr) %>% spread(tmp, readable)
      
      cat("\\begin{center}\n")
      cat("\n\nSample sizes and readabilities for different phenotypes.\n\n")
      # cat("\n\n\\smallskip\n\n")
      cat("\\begin{tabular}[c]{@{}l", paste(collapse="", rep("r", ncol(tbl) - 1)), "@{}}\n", sep="")
      cat("\\toprule\n")
      cat(" & & \\multicolumn{", ncol(tbl) - 2, "}{c}{readability / \\%}\\\\\n", sep="")
      cat("\\cmidrule{3-", ncol(tbl), "}\n", sep="")
      cat("phenotype & n &", tbl %>% ungroup() %>% select(-1, -2) %>% names %>% substring(2) %>% paste(collapse=" h & "), "h\\tabularnewline\n")
      cat("\\midrule\n")
      tbl %>% apply(1, function(x) {cat(sep="", paste(x, collapse="&"), "\\tabularnewline\n")})
      cat("\\bottomrule\n")
      cat("\\end{tabular}\n")
      cat("\\end{center}\n")
      cat("\n\n")
      
      ##############
      # Masterplot #
      ##############
      
      dTmp <- dTmp[!is.na(dTmp$diameter), ]

      if (manuscript == 2) {
        if (length(allNonWt) == 0) {
          addBoxplot(ggplot(dTmp, aes(x=hDiscr, y=diameter, fill=isWt)) + facet_wrap(~ phenotypes) + scale_fill_discrete("", drop=FALSE) + theme(aspect.ratio=1)) %>% polishPlot %>% print
        } else {
          p1 <- masterplot(dTmp, wt=wt, allNonWt=allNonWt)
          p2 <- aucPlot(dTmp, manuscript=manuscript, wt=wt, allNonWt=allNonWt)
          alignPlots(p1, p2)
        }
      } else if (manuscript == 3) {
        p1 <- myBoxplot(dTmp)
        p2 <- aucPlot(dTmp, manuscript=manuscript)
        alignPlots(p1, p2)
        
        # plot_grid(p1, p2, ncol=1, rel_heights=c(2, 1)) %>% print  # , labels=c("a", "b")
      }

      ######################
      # Caption / footnote #
      ######################
      
      if (length(allNonWt) == 0) {
        cat(sep='', '\n\nBoxes represent diameter ranges from the 5th to the 95th percentile, whiskers represent the full range of diameters, and bold lines indicate median diameter values.\n\n')
      } else {
        cat(sep='', '\n\n\\textbf{(Top)} Boxes represent diameter ranges from the 5th to the 95th percentile, whiskers represent the full range of diameters, and bold lines indicate median diameter values. \\textbf{(Bottom)} The area under the receiver operating characteristic curve (AUC) quantifies how well a phenotype is separated from the wild type. An AUC of 1 indicates perfect separation while an AUC of 0.5 corresponds to complete overlap of the two populations. Error bars display 95\\% confidence intervals.\n\n')
      }

      # ####################
      # # Additional plots #
      # ####################
      # 
      # if (showAllPlots) {
      #   cat('\n\n\\pagebreak\n')
      #   print(lineplotHombachMaster(dTmp))  # lineplot requested by M. Hombach and M. Jetter
      #   cat('\n\n\\pagebreak\n')
      #   print(lineplotWFacets(dTmp))  # simple lineplot
      # 
      #   # More accessible figure based on boxplots:
      #   bxplt <- list()
      #   for (i in 1:(nlevels(dTmp$phenotypes) - 1)) {
      #     nonWt <- allNonWt[i]
      #     bxplt[[i]] <- myBoxplot(dTmp[dTmp$phenotypes %in% c(wt, nonWt), ]) + ggtitle(nonWt)
      #   }
      #   cat('\n\n\\pagebreak\n')
      #   print(do.call("plot_grid", c(bxplt)))
      # }

      cat('\n\n\\pagebreak\n')
    }
  }
}

##########
# Tab. 2 #
##########

tab2 <- function() {
  cutoff <- 2  # adjust manually in figure legend below!!
  discrDeltaDiam <- function(d) {
    factor(ifelse(d < -cutoff, 1, ifelse(d > cutoff, 3, 2)), levels=1:3)
  }
  
  dTab2 <- dRaw %>% filter(phenotypes == "wild-type") %>%
    group_by(organism, antibiotic, h) %>% summarise(q05Diam=quantile(diameter, 0.05)) %>%
    mutate(deltaDiam=(q05Diam[h == 18] - q05Diam[h == 6])) %>% mutate(deltaDiscr=discrDeltaDiam(deltaDiam))
  
  ggplot(dTab2, aes(x=h, y=q05Diam, colour=deltaDiscr)) +
          geom_line() + geom_point() +
          facet_grid(antibiotic ~ organism) +
          theme(strip.text.y=element_text(angle=0), legend.position="top") +
          scale_x_continuous(breaks=c(6, 8, 12, 18)) + scale_y_continuous(breaks=pretty_breaks(n=3)) +
          xlab("time / h") + ylab("5th percentile of diameter for wild-type strains / mm") +
          scale_colour_manual(name=expression(Delta==d[18*h]-d[6*h]), values=c("orangered3", brewer.pal(3, 'Set1')[3], "orange"),
                              labels=c(expression("decrease "*(Delta < -2*mm)),
                                       expression("stable "*(-2*mm<=paste(Delta<=2*mm))),
                                       expression("increase "*(Delta>2*mm))), drop=FALSE)  # "skyblue4"
  
  ggsave(paste(sep="", "Figs/fig3.png"), width=18, height=22, units="cm", dpi=1200, limitsize=TRUE)
  dTab2 %>% select(organism, antibiotic, deltaDiscr) %>% unique %>% spread(organism, deltaDiscr) %>% write.table("tab2.csv", sep=";", row.names=FALSE)
}