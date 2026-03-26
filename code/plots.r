# Copyright 2025 OTH - Laboratory for Digitalisation (LfD)
# Written by Lukas Schmidbauer
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

library(tidyverse)
library(patchwork)
library(tikzDevice)
library(scales)
library(ggh4x)

INCH.PER.CM <- 1 / 2.54
WIDTH <- 18.1 * INCH.PER.CM

COL.WIDTH <- 6.95 * INCH.PER.CM
PLOT.HEIGHT.S6 <- 1.1*COL.WIDTH

BASE.SIZE <- 9
FORMAT <- "tex"
theme_paper_base <- function() {
    return(theme_bw(base_size = BASE.SIZE) +
        theme(
            axis.title.x = element_text(size = BASE.SIZE),
            axis.title.y = element_text(size = BASE.SIZE),
            legend.title = element_text(size = BASE.SIZE, margin=margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")),
            legend.position = "top",
            legend.box = "vertical",
            legend.spacing.y = unit(-0.2, "cm"),
            plot.margin = unit(c(0, 0.1, 0, 0.1), "cm")
        ))
}

guide_paper_base <- function() {
    return(guides(
        shape = guide_legend(order = 1, nrow = 1, byrow = TRUE),
        col = guide_legend(order = 2, nrow = 1, byrow = TRUE, reverse = TRUE)
    ))
}

guide_paper_base_section6 <- function() {
    return(guides(
        shape = guide_legend(order = 1, nrow = 1, byrow = TRUE),
        col = guide_legend(order = 2, nrow = 2, byrow = TRUE, reverse = FALSE)
    ))
}

options(
    tikzDocumentDeclaration = "\\documentclass[conference]{IEEEtran}",
    tikzLatexPackages = c(
        getOption("tikzLatexPackages"),
        "\\usepackage{amsmath}"
    ),
    tikzSanitizeCharacters = c("%", "#", "_"),
    tikzReplacementCharacters = c("\\%", "\\#", "\\_")
)
do.save.tikz <- function(g, out.name, .width, .height, .sanitize=FALSE) {
    tikz(str_c(OUT.PATH, out.name, ".tex"),
        width = .width, height = .height,
        sanitize = .sanitize
    )
    print(g)
    dev.off()
}

extract_monomial_density <- function(text) {
    matches <- gregexpr("\\d.\\d", text, perl= TRUE)
    extracted <- regmatches(text, matches)
    return(as.numeric(extracted))
}


OUT.PATH <- "./RPlots/"
ifelse(!dir.exists(file.path(OUT.PATH)),
        dir.create(file.path(OUT.PATH)),
        "")
POINT.ALPHA <- 0.6
POINT.SIZE <- 0.5

BOXPLOT.WIDTH <- 1
BOXPLOT.LW <- 0.5 
BOXPLOT.Med <- 0.7
BOXPLOT.OUTSH <- 20
BOXPLOT.OUTAL <- 0.7

LFD.COLOURS <- c("black", "#E69F00", "#999999", "#beaed4", "#009371", "#ed665a", "#1f78b4", "#009371")
LFD.COLOURS.S6 <- c("black", "#E69F00", "#999999", "#009371", "#beaed4", "#ed665a", "#1f78b4", "#009371")
LFD.SHAPES <- c(15, 16, 17, 4, 5, 8, 9, 20)

###################################################################
####### Industrial k-SAT Input data
###################################################################


read_kSAT_csvs <- function(pattern= "doublepower", type="c", folder = "Experiments_Industrial") {
  df <- list.files(path = folder, pattern = paste("^", pattern, ".*", type, "\\.csv$", sep=""), full.names = TRUE) %>%
        lapply(read_csv, col_types=cols(bitstr = col_character(), SA_bitstr = col_character())) %>%
        bind_rows()

  return(df)
}

relable_k <- function(orig) {
    return(paste0("$k=", orig, "$"))
}

relable_path <- function(x) {
    case_when(
        grepl("npVPubo",x) ~ "direct PUBO",
        grepl("npVQubo",x) ~ "direct QUBO",
        grepl("pVPubo",x) ~ "optimised PUBO",
        grepl("pVQubo",x) ~ "optimised QUBO",
        TRUE ~ "NA"
    )
}

relable_clauses <- function(x) {
    return(paste0("$|C|=",x,"$"))
}

relable_clauses_reverse <- function(x) {
    return(gsub("\\D+","",x))
}

relable_SA_steps <-function(x){
    case_when(
        grepl("linear",x) ~ "Linear # steps",
        grepl("quadratic",x) ~ "Quadratic # steps",
        TRUE ~ "NA"
    )
}

#######################################
# Params
#######################################
SAT_gen_variant = "power" # "doublepower", "power", "uniform"
SAT_clauses= c(13, 53, 107, 163, 263)


if(file.exists("Experiments_Industrial/dfMetrics.csv")){
    dfMetrics <- read_csv("Experiments_Industrial/dfMetrics.csv")
    
}else{
    dfMetrics <- read_kSAT_csvs(SAT_gen_variant) %>%
        filter(SAT_c %in% SAT_clauses) %>%
        mutate(SAT_avgk = relable_k(SAT_avgk), SAT_c=relable_clauses(SAT_c))
    write.csv(dfMetrics, "Experiments_Industrial/dfMetrics.csv")
}
dfMetrics$SAT_avgk = factor(dfMetrics$SAT_avgk, levels=relable_k(c(3,5,7,11,13))) 
dfMetrics$SAT_c = factor(dfMetrics$SAT_c, levels=relable_clauses(SAT_clauses))

if(file.exists("Experiments_Industrial/dfSA.csv")){
    dfSA <- read_csv("Experiments_Industrial/dfSA.csv")
    
}else{
    dfSA <- read_kSAT_csvs(SAT_gen_variant, type="c_SA")
    dfSA <- dfSA %>%
        left_join(dfMetrics, by="SAT_filename")
    write.csv(dfSA, "Experiments_Industrial/dfSA.csv")
}

dfSA$SAT_avgk = factor(dfSA$SAT_avgk, levels=relable_k(c(3,5,7,11,13))) 


if(file.exists("Experiments_Industrial/dfQAOA.csv")){
    dfQAOA <- read_csv("Experiments_Industrial/dfQAOA.csv")
    
}else{
    dfQAOA <- read_kSAT_csvs(SAT_gen_variant, type="c_QAOA") %>%
        filter(QA_status == "ok")
    dfQAOA <- dfQAOA %>%
    left_join(dfMetrics, by="SAT_filename")
    write.csv(dfQAOA, "Experiments_Industrial/dfQAOA.csv")
}

if(file.exists("Experiments_Industrial_paperWithoutSA/dfMetrics.csv")){
    print("Detected scaling experiments. Continuing with different dataframe for those except for SA and QAOA.")
    dfMetrics <- read_csv("Experiments_Industrial_paperWithoutSA/dfMetrics.csv")
}
dfMetrics$SAT_avgk = factor(dfMetrics$SAT_avgk, levels=relable_k(c(3,5,7,11,13))) 
dfMetrics$SAT_c = factor(dfMetrics$SAT_c, levels=relable_clauses(SAT_clauses))


#Approx # Monomials for optimised SAT
create.grid <- function(x.grid, f, t) {
    dat <- tibble(x = x.grid, y = f(x, t), t=t)
    return(dat)
}

tab.ApproxMonomials <- bind_rows(lapply(c(2,4,8,16), function(t) {
    x.grid <- seq(0, min(t,6), 0.1)
    create.grid(x.grid, \(x, t) 2**x + 3*(t-x), t)
})) |> mutate(t=as.factor(t))

plt.ApproxMonomials <- ggplot(tab.ApproxMonomials, aes(x=x, y=y, colour=t, linetype=t)) + 
    geom_line() +
    scale_colour_manual(values=LFD.COLOURS.S6) +
    xlab("Remaining positive literals $\\alpha$") +
    ylab("Approx. # monomials") +
    theme_paper_base() 

ggsave(plot = plt.ApproxMonomials, filename = "SATApproxMonomials.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.ApproxMonomials, "genSATApproxMonomials", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)

#Power Law distributed
tab.PowerLawPhase <- bind_rows(lapply(c(1), function(t) {
    x.grid <- seq(3, 30, 0.1)
    create.grid(x.grid, \(x, t) (2*x-1)/(x-1), t)
})) |> mutate(t=as.factor(t))

plt.PowerLawPhase <- ggplot(tab.PowerLawPhase, aes(x=x, y=y)) + 
    geom_line() +
    xlab("$k$") +
    ylab("$\\frac{2k-1}{k-1}$") +
    theme_paper_base() 

ggsave(plot = plt.PowerLawPhase, filename = "SATPowerLawPhases.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.PowerLawPhase, "genSATPowerLawPhase", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)


#SAT
tab.SAT <- dfMetrics

plt.SATvars <- ggplot(tab.SAT, aes(x=SAT_v, y=pVSAT_v, colour=SAT_c)) +
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line() +
    facet_grid(rows = vars(SAT_avgk)) +
    theme_paper_base() +
    scale_colour_manual("# Clauses", values=LFD.COLOURS.S6) +
    xlab("# Variables in $k$-SAT") +
    ylab("# Variables in optimised SAT")

plt.SATpLits <- ggplot(tab.SAT, aes(x=SAT_pLit, y=pVSAT_pLit, colour=SAT_c)) +
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line() +
    facet_grid(rows = vars(SAT_avgk)) +
    theme_paper_base() +
    scale_colour_manual("# Clauses", values=LFD.COLOURS.S6) +
    xlab("# Positive literals in $k$-SAT") +
    ylab("# Positive literals in optimised SAT")

plt.SATnLits <- ggplot(tab.SAT, aes(x=SAT_nLit, y=pVSAT_nLit, colour=SAT_c)) +
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line() +
    facet_grid(rows = vars(SAT_avgk)) +
    theme_paper_base() +
    scale_colour_manual("# Clauses", values=LFD.COLOURS.S6) +
    xlab("# Negative literals in $k$-SAT") +
    ylab("# Negative literals in optimised SAT")

tab.SATClauses <- tab.SAT %>% 
    mutate(SAT_c = as.numeric(relable_clauses_reverse(SAT_c)))

design.SATcombined <- "
111
234
"
plt.SATcombined <- guide_area() + plt.SATvars + plt.SATpLits + plt.SATnLits + plot_layout(design=design.SATcombined, guides="collect", heights = c(0.5,3) )
ggsave(plot = plt.SATcombined, filename = "SATMetrics.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.SATcombined, "genSATVarsLits", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)

#Number of variables PBF
tab.NumVarsPBF <- dfMetrics %>% 
    pivot_longer(cols = c(npVPubo_num_vars, pVPubo_num_vars, npVQubo_num_vars, pVQubo_num_vars), names_to="Variant", values_to="Num_Vars") %>%
    mutate(Variant = relable_path(Variant))
    
plt.NumVarsPBF <- ggplot(tab.NumVarsPBF, aes(x=SAT_v, y=Num_Vars, colour=Variant)) + #colour=as.factor(nCVariant), linetype=as.factor(nLVariant))) +
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line() +
    facet_grid(cols = vars(SAT_c), rows = vars(SAT_avgk)) +
    theme_paper_base() +
    scale_y_log10(breaks = c(1e1, 1e2, 1e3), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    scale_colour_manual("Transformation path", values=LFD.COLOURS.S6) +
    guide_paper_base_section6() +
    xlab("# Variables in $k$-SAT") +
    ylab("# Variables in PBF [log]")
ggsave(plot = plt.NumVarsPBF, filename = "NumVarsPBF.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.NumVarsPBF, "genSATNumVarsPBF", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)


#Number of Monomials PBF
tab.NumMonomialsPBF <- dfMetrics %>% 
    pivot_longer(cols = c(npVPubo_num_monomials, pVPubo_num_monomials, npVQubo_num_monomials, pVQubo_num_monomials), names_to="Variant", values_to="Num_monomials") %>%
    mutate(Variant = relable_path(Variant))

plt.NumMonomialsPBF <- ggplot(tab.NumMonomialsPBF, aes(x=SAT_v, y=Num_monomials, colour=Variant)) +
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line() +
    facet_grid(cols = vars(SAT_c), rows = vars(SAT_avgk)) +
    theme_paper_base() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    scale_colour_manual("Transformation path", values=LFD.COLOURS.S6) +
    guide_paper_base_section6() +
    xlab("# Variables in $k$-SAT") +
    ylab("# Monomials in PBF [log]")
ggsave(plot = plt.NumMonomialsPBF, filename = "NumMonomialsPBF.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.NumMonomialsPBF, "genSATNumMonomialsPBF", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)

#Number of variables Ising
tab.NumVarsIsing <- dfMetrics %>% 
    pivot_longer(cols = c(npVPuboIsing_num_vars, pVPuboIsing_num_vars, npVQuboIsing_num_vars, pVQuboIsing_num_vars), names_to="Variant", values_to="Num_Vars") %>%
    mutate(Variant = relable_path(Variant))

plt.NumVarsIsing <- ggplot(tab.NumVarsIsing, aes(x=SAT_v, y=Num_Vars, colour=Variant)) +
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line() +
    facet_grid(cols = vars(SAT_c), rows = vars(SAT_avgk)) +
    theme_paper_base() +
    scale_y_log10(breaks = c(1e1, 1e2, 1e3), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    scale_colour_manual("Transformation path", values=LFD.COLOURS.S6) +
    guide_paper_base_section6() +
    xlab("# Variables in $k$-SAT") +
    ylab("# Variables in Ising [log]")
ggsave(plot = plt.NumVarsIsing, filename = "NumVarsIsing.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.NumVarsIsing, "genSATNumVarsIsing", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)

#Number of Monomials Ising
tab.NumMonomialsIsing <- dfMetrics %>% 
    pivot_longer(cols = c(npVPuboIsing_num_monomials, pVPuboIsing_num_monomials, npVQuboIsing_num_monomials, pVQuboIsing_num_monomials), names_to="Variant", values_to="Num_monomials") %>%
    mutate(Variant = relable_path(Variant))

plt.NumMonomialsIsing <- ggplot(tab.NumMonomialsIsing, aes(x=SAT_v, y=Num_monomials, colour=Variant)) +
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line() +
    facet_grid(cols = vars(SAT_c), rows = vars(SAT_avgk)) +
    theme_paper_base() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    scale_colour_manual("Transformation path", values=LFD.COLOURS.S6) +
    guide_paper_base_section6() +
    xlab("# Variables in $k$-SAT") +
    ylab("# Monomials in Ising [log]")
ggsave(plot = plt.NumMonomialsIsing, filename = "NumMonomialsIsing.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.NumMonomialsIsing, "genSATNumMonomialsIsing", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)

#Number of qubits QAOA
tab.NumQubitsQAOA <- dfMetrics %>% 
    pivot_longer(cols = c(npVPuboIsingLQaoa_num_qubits, pVPuboIsingLQaoa_num_qubits, npVQuboIsingLQaoa_num_qubits, pVQuboIsingLQaoa_num_qubits), names_to="Variant", values_to="Num_qubits") %>%
    mutate(Variant = relable_path(Variant))

plt.NumQubitsQAOA <- ggplot(tab.NumQubitsQAOA, aes(x=SAT_v, y=Num_qubits, colour=Variant)) +
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line() +
    facet_grid(cols = vars(SAT_c), rows = vars(SAT_avgk)) +
    theme_paper_base() +
    scale_y_log10(breaks = c(1e1,1e2,1e3), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    scale_colour_manual("Transformation path", values=LFD.COLOURS.S6) +
    guide_paper_base_section6() +
    xlab("# Variables in $k$-SAT") +
    ylab("# Qubits QAOA [log]")
ggsave(plot = plt.NumQubitsQAOA, filename = "NumQubitsQAOA.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.NumQubitsQAOA, "genSATNumQubitsQAOA", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)

#Number of 1 qubit gates QAOA
tab.NumOneQubitGates <- dfMetrics %>% 
    pivot_longer(cols = c(npVPuboIsingLQaoa_num_local_gates, pVPuboIsingLQaoa_num_local_gates, npVQuboIsingLQaoa_num_local_gates, pVQuboIsingLQaoa_num_local_gates), names_to="Variant", values_to="Num_OneQubitGates") %>%
    mutate(Variant = relable_path(Variant))

plt.NumOneQubitGates <- ggplot(tab.NumOneQubitGates, aes(x=SAT_v, y=Num_OneQubitGates, colour=Variant)) +
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line() +
    facet_grid(cols = vars(SAT_c), rows = vars(SAT_avgk)) +
    theme_paper_base() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    scale_colour_manual("Transformation path", values=LFD.COLOURS.S6) +
    guide_paper_base_section6() +
    xlab("# Variables in $k$-SAT") +
    ylab("# 1-Qubit gates in QAOA ($p=1$) [log]")
ggsave(plot = plt.NumOneQubitGates, filename = "NumOneQubitGatesQAOA.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.NumOneQubitGates, "genSATNumOneQubitGatesQAOA", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)

#Number of 2 qubit gates QAOA
tab.NumTwoQubitGates <- dfMetrics %>% 
    pivot_longer(cols = c(npVPuboIsingLQaoa_num_non_local_gates, pVPuboIsingLQaoa_num_non_local_gates, npVQuboIsingLQaoa_num_non_local_gates, pVQuboIsingLQaoa_num_non_local_gates), names_to="Variant", values_to="Num_TwoQubitGates") %>%
    mutate(Variant = relable_path(Variant))

plt.NumTwoQubitGates <- ggplot(tab.NumTwoQubitGates, aes(x=SAT_v, y=Num_TwoQubitGates, colour=Variant)) +
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line() +
    facet_grid(cols = vars(SAT_c), rows = vars(SAT_avgk)) +
    theme_paper_base() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    scale_colour_manual("Transformation path", values=LFD.COLOURS.S6) +
    guide_paper_base_section6() +
    xlab("# Variables in $k$-SAT") +
    ylab("# 2-Qubit gates in QAOA ($p=1$) [log]")
ggsave(plot = plt.NumTwoQubitGates, filename = "NumTwoQubitGatesQAOA.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.NumTwoQubitGates, "genSATNumTwoQubitGatesQAOA", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)

#Circuit depth QAOA
tab.CircuitDepth <- dfMetrics %>% 
    pivot_longer(cols = c(npVPuboIsingLQaoa_depth, pVPuboIsingLQaoa_depth, npVQuboIsingLQaoa_depth, pVQuboIsingLQaoa_depth), names_to="Variant", values_to="CircuitDepth") %>%
    mutate(Variant = relable_path(Variant))

plt.CircuitDepth <- ggplot(tab.CircuitDepth, aes(x=SAT_v, y=CircuitDepth, colour=Variant)) +
    geom_point(size=POINT.SIZE, alpha=POINT.ALPHA) +
    geom_line() +
    facet_grid(cols = vars(SAT_c), rows = vars(SAT_avgk)) +
    theme_paper_base() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "l") +
    scale_colour_manual("Transformation path", values=LFD.COLOURS.S6) +
    guide_paper_base_section6() +
    xlab("# Variables in $k$-SAT") +
    ylab("Circuit depth for single QAOA layer [log]")
ggsave(plot = plt.CircuitDepth, filename = "CircuitDepthQAOA.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.CircuitDepth, "genSATCircuitDepthQAOA", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)

tab.SABox <- dfSA %>%
    filter(SAT_c %in% c(relable_clauses(53)) & SAT_v %in% c(13,23,35,37,41,46,47)) %>%
    mutate(SA_pbf_type = relable_path(SA_pbf_type)) %>%
    group_by(SAT_filename, SA_pbf_type, SA_steps_variant, SA_init_temp) %>%
    mutate(min_group_energy = min(SA_energy_original_pbf)) %>%
    group_by(SAT_filename, SA_pbf_type, SA_steps_variant) %>%
    filter(min(min_group_energy) == min_group_energy) %>%
    mutate(SA_steps_variant = relable_SA_steps(SA_steps_variant))

plt.SABox <- ggplot(tab.SABox, aes(x=as.factor(SAT_v), y=SA_energy_original_pbf, colour=SA_pbf_type)) +
    geom_boxplot(linewidth = BOXPLOT.LW, fatten = BOXPLOT.Med, outlier.shape = BOXPLOT.OUTSH, outlier.alpha = BOXPLOT.OUTAL) +
    facet_grid2(cols = vars(SA_steps_variant), rows = vars(SAT_avgk), scales = "free_y") +
    theme_paper_base() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_colour_manual("Transformation path", values=LFD.COLOURS.S6) +
    facetted_pos_scales(y=list(TRUE ~ scale_y_continuous(n.breaks=4)))+
    guide_paper_base_section6() +
    xlab("# Variables in $k$-SAT") +
    ylab("Energy in PUBO ($|C| = 53$)")

ggsave(plot = plt.SABox, filename = "SABox.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.SABox, "genSATSABox", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)

tab.QAOAbox <- dfQAOA %>%
    filter(SAT_c %in% relable_clauses(c(13, 53)) & SAT_avgk %in% relable_k(c(3, 5))) %>%
    mutate(QA_pbf_type = relable_path(QA_pbf_type))

plt.QAOAbox <- ggplot(tab.QAOAbox, aes(x=as.factor(SAT_v), y=energy_original_pbf, colour=QA_pbf_type)) +
    geom_boxplot(linewidth = BOXPLOT.LW, fatten = BOXPLOT.Med, outlier.shape = BOXPLOT.OUTSH, outlier.alpha = BOXPLOT.OUTAL) +
    facet_grid(cols = vars(SAT_c), rows = vars(SAT_avgk), scales="free") + 
    theme_paper_base() +
    scale_colour_manual("Transformation path", values=LFD.COLOURS.S6) +
    guide_paper_base_section6() +
    xlab("# Variables in $k$-SAT") +
    ylab("Pubo energy")
ggsave(plot = plt.QAOAbox, filename = "QAOAbox.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = PLOT.HEIGHT.S6)
do.save.tikz(plt.QAOAbox, "QAOAbox", 2* COL.WIDTH, PLOT.HEIGHT.S6, TRUE)

###################################################################
####### Input data
###################################################################

StypeGB <- c("0.5", "0.8", "1.0")
StypeMB <- c("simple", "medium", "better")
data <- read.csv("Performance_results.csv") %>% 
    mutate(Function_density = extract_monomial_density(gen_type)) %>% 
    filter(Selection_type %in%  StypeGB | Selection_type %in% StypeMB)

GBdata <- data %>% filter(monomial_or_graph_based == "graph_based") 
MBdata <- data %>% filter(monomial_or_graph_based == "monomial_based")

###################################################################
####### Plot variables vs runtime
###################################################################

tab.runtime.data <- filter(data, Function_density %in% c(0.2, 0.4, 0.6, 0.8, 1.0)) 

plt.runtimeMBvsGB <- ggplot(tab.runtime.data, aes(x=num_variables_before, y=time, colour=Selection_type)) +
    geom_point(size = POINT.SIZE, alpha = POINT.ALPHA) +
    facet_grid(cols = vars(Function_density), rows= vars(monomial_or_graph_based)) + 
    theme_paper_base() +
    guide_paper_base() +
    scale_colour_manual("Selection type", values = LFD.COLOURS) +
    xlab("# Variables prior to reduction") +
    ylab("Quadratisation time [s]")

ggsave(plot = plt.runtimeMBvsGB, filename = "RuntimeMBvsGB.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = COL.WIDTH)
do.save.tikz(plt.runtimeMBvsGB, "genRuntimeMBvsGB", 2* COL.WIDTH, COL.WIDTH, TRUE)

plt.RuntimeMBvsGBlog <- ggplot(tab.runtime.data, aes(x=num_variables_before, y=time, colour=Selection_type)) +
    geom_point(size = POINT.SIZE, alpha = POINT.ALPHA) +
    facet_grid(cols = vars(Function_density), rows= vars(monomial_or_graph_based)) + 
    theme_paper_base() +
    guide_paper_base() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "lr") + 
    scale_colour_manual("Selection type", values = LFD.COLOURS) +
    xlab("# Variables prior to reduction") +
    ylab("Quadratisation time [s, log]")

ggsave(plot = plt.RuntimeMBvsGBlog, filename = "RuntimeMBvsGBlog.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = COL.WIDTH)
do.save.tikz(plt.RuntimeMBvsGBlog, "genRuntimeMBvsGBlog", 2* COL.WIDTH, COL.WIDTH, TRUE)

###################################################################
####### Plot terms vs runtime
###################################################################

tab.runtime.data$poly_len = tab.runtime.data$poly_len / 1000000.0

plt.termsRuntimeMBvsGB <- ggplot(tab.runtime.data, aes(x=poly_len, y=time, colour=Selection_type)) +
    geom_point(size = POINT.SIZE, alpha = POINT.ALPHA) +
    facet_grid(cols = vars(Function_density), rows= vars(monomial_or_graph_based)) + 
    theme_paper_base() +
    guide_paper_base() +
    scale_colour_manual("Selection type", values = LFD.COLOURS) +
    scale_x_continuous(labels = function(x) format(x, big.mark = "'", scientific = FALSE)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) + 
    xlab("# Terms in $f$ prior to reduction in Mio.") +
    ylab("Quadratisation time [s]")

ggsave(plot = plt.termsRuntimeMBvsGB, filename = "TermsRuntimeMBvsGB.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = COL.WIDTH)
do.save.tikz(plt.termsRuntimeMBvsGB, "genTermsRuntimeMBvsGB", 2* COL.WIDTH, COL.WIDTH, TRUE)

plt.termsRuntimeMBvsGBlog <- ggplot(tab.runtime.data, aes(x=poly_len, y=time, colour=Selection_type)) +
    geom_point(size = POINT.SIZE, alpha = POINT.ALPHA) +
    facet_grid(cols = vars(Function_density), rows= vars(monomial_or_graph_based)) + 
    theme_paper_base() +
    guide_paper_base() +
    scale_colour_manual("Selection type", values = LFD.COLOURS) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "r") + 
    scale_x_continuous(labels = function(x) format(x, big.mark = "'", scientific = FALSE)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) + 
    xlab("# Terms in $f$ prior to reduction in Mio.") +
    ylab("Quadratisation time [s, log]")

ggsave(plot = plt.termsRuntimeMBvsGBlog, filename = "TermsRuntimeMBvsGBlog.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = COL.WIDTH)
do.save.tikz(plt.termsRuntimeMBvsGBlog, "genTermsRuntimeMBvsGBlog", 2* COL.WIDTH, COL.WIDTH, TRUE)

###################################################################
####### Plot introduced variables
###################################################################

tab.vars.data <- filter(data, Function_density %in% c(0.2, 0.6, 1.0))

guide_Vars <- function(.reverse=TRUE) {
    return(guides(
        shape = guide_legend(order = 1, nrow = 1, byrow = TRUE),
        col = guide_legend(order = 2, nrow = 3, byrow = TRUE, reverse = .reverse)
    ))
}
y.lim <- c(0,max(tab.vars.data$num_variables_after))

plt.VarsBeforeAfterMB <- ggplot(filter(tab.vars.data, monomial_or_graph_based == "monomial_based"), aes(x=num_variables_before, y=num_variables_after, colour = Selection_type)) +
    geom_point(size = POINT.SIZE, alpha = POINT.ALPHA) +
    facet_grid(rows=vars(Function_density)) + 
    theme_paper_base() + 
    guide_Vars(FALSE) +
    scale_colour_manual("Selection type\n(monomial_based)  ", values = LFD.COLOURS[c(4,5,6)]) +
    scale_y_continuous(limits=y.lim) +
    xlab("# Variables prior to reduction") +
    ylab("# Variables after reduction")

plt.VarsBeforeAfterGB <- ggplot(filter(tab.vars.data, monomial_or_graph_based == "graph_based"), aes(x=num_variables_before, y=num_variables_after, colour = Selection_type)) +
    geom_point(size = POINT.SIZE, alpha = POINT.ALPHA) +
    facet_grid(rows=vars(Function_density)) + 
    theme_paper_base() + 
    guide_Vars() +
    scale_colour_manual("Selection type\n(graph_based)  ", values = LFD.COLOURS) +
    scale_y_continuous(limits=y.lim) +
    xlab("# Variables prior to reduction") +
    ylab("# Variables after reduction")

plt.VarsBeforeAfterMBvsGB <- (plt.VarsBeforeAfterGB + plt.VarsBeforeAfterMB + plot_layout(axes = "collect") ) 

ggsave(plot = plt.VarsBeforeAfterMBvsGB, filename = "VarsBeforeAfterMBvsGB.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height=COL.WIDTH)
do.save.tikz(plt.VarsBeforeAfterMBvsGB, "genVarsBeforeAfterMBvsGB", 2*COL.WIDTH, COL.WIDTH, TRUE)

###################################################################
####### Plot densities
###################################################################

tab.densities.data <- data %>% 
                    filter(Function_density %in% c(0.1, 0.3, 0.5, 0.8, 1.0)) %>% 
                    mutate(monomial_or_graph_based = case_when( monomial_or_graph_based == "graph_based" ~ "Selection type (graph\\_based)",
                                                                monomial_or_graph_based == "monomial_based" ~ "Selection type (monomial\\_based)"))

plt.densitiesMBvsGB <- ggplot(tab.densities.data, aes(x=num_variables_before, y=density_deg2r, colour=as.factor(Function_density))) +
    geom_point(size = POINT.SIZE, alpha=POINT.ALPHA) +
    facet_nested(~ monomial_or_graph_based + Selection_type) +
    theme_paper_base() +
    guide_paper_base() +
    scale_colour_manual("Function density", values = LFD.COLOURS) +
    scale_shape_manual("Algorithm", values = LFD.SHAPES) +
    xlab("\\# Variables prior to reduction") +
    ylab("Density $d_2$ after reduction")

ggsave(plot = plt.densitiesMBvsGB, filename = "DensityMBvsGB.pdf", path = OUT.PATH, width = 2*COL.WIDTH, height = COL.WIDTH)
do.save.tikz(plt.densitiesMBvsGB, "genDensityMBvsGB", 2*COL.WIDTH, COL.WIDTH, FALSE)


###################################################################
####### Plot variables vs terms
###################################################################

tab.DegreeSweepVariables <- read.csv("Performance_results_degreeSweep20.csv") %>%
    mutate(Function_density = extract_monomial_density(gen_type)) 

tab.DegreeSweepVariables$Function_degree = as.numeric(str_extract(tab.DegreeSweepVariables$gen_type, "\\s\\d+"))

tab.DegreeSweepVariables = tab.DegreeSweepVariables %>%
     filter(Function_degree %in% c(1,2,4,8,12,16))

plt.DegreeSweepVariables <- ggplot(tab.DegreeSweepVariables, aes(x=num_variables_before, y=poly_len, colour = as.factor(Function_degree))) +
    geom_point(size = POINT.SIZE, alpha = POINT.ALPHA) +
    geom_line() +
    theme_paper_base() +
    guide_paper_base() +
    scale_colour_manual("Function degree", values = LFD.COLOURS) +
    scale_y_continuous(trans="log2", breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
    annotation_logticks(base=2 ,sides = "lr") + 
    xlab("\\# Variables prior to reduction") +
    ylab("\\# Terms in $f$ prior to reduction [log]")

ggsave(plot = plt.DegreeSweepVariables, filename = "DegreeSweepVariables.pdf", path=OUT.PATH, width = 2*COL.WIDTH, height = COL.WIDTH)
do.save.tikz(plt.DegreeSweepVariables, "genDegreeSweepVariables", 2*COL.WIDTH, COL.WIDTH, FALSE)