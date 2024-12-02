library(tidyverse)
library(patchwork)
library(tikzDevice)
library(scales)
library(ggh4x)
library(ggpmisc)


INCH.PER.CM <- 1 / 2.54
WIDTH <- 18.1 * INCH.PER.CM

COL.WIDTH <- 8.85 * INCH.PER.CM
BASE.SIZE <- 9
FORMAT <- "tex"
theme_paper_base <- function() {
    return(theme_bw(base_size = BASE.SIZE) +
        theme(
            axis.title.x = element_text(size = BASE.SIZE),
            axis.title.y = element_text(size = BASE.SIZE),
            legend.title = element_text(size = BASE.SIZE),
            legend.position = "top",
            legend.box = "vertical",
            legend.spacing.y = unit(-0.2, "cm"),
            plot.margin = unit(c(0, 0.2, 0, 0), "cm")
        ))
}

guide_paper_base <- function() {
    return(guides(
        shape = guide_legend(order = 1, nrow = 1, byrow = TRUE),
        col = guide_legend(order = 2, nrow = 1, byrow = TRUE, reverse = TRUE)
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

LFD.COLOURS <- c("black", "#E69F00", "#999999", "#beaed4", "#009371", "#ed665a", "#1f78b4", "#009371")
LFD.SHAPES <- c(15, 16, 17, 4, 5, 8, 9, 20)

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
    scale_colour_manual("Selection type\n(monomial_based)", values = LFD.COLOURS[c(4,5,6)]) +
    scale_y_continuous(limits=y.lim) +
    xlab("# Variables prior to reduction") +
    ylab("# Variables after reduction")

plt.VarsBeforeAfterGB <- ggplot(filter(tab.vars.data, monomial_or_graph_based == "graph_based"), aes(x=num_variables_before, y=num_variables_after, colour = Selection_type)) +
    geom_point(size = POINT.SIZE, alpha = POINT.ALPHA) +
    facet_grid(rows=vars(Function_density)) + 
    theme_paper_base() + 
    guide_Vars() +
    scale_colour_manual("Selection type\n(graph_based)", values = LFD.COLOURS) +
    scale_y_continuous(limits=y.lim) +
    xlab("# Variables prior to reduction") +
    ylab("# Variables after reduction")

plt.VarsBeforeAfterMBvsGB <- (plt.VarsBeforeAfterGB + plt.VarsBeforeAfterMB + plot_layout(axes = "collect") ) 

ggsave(plot = plt.VarsBeforeAfterMBvsGB, filename = "VarsBeforeAfterMBvsGB.pdf", path = OUT.PATH, width = COL.WIDTH, height=COL.WIDTH)
do.save.tikz(plt.VarsBeforeAfterMBvsGB, "genVarsBeforeAfterMBvsGB", COL.WIDTH, COL.WIDTH, TRUE)

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
    scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
    annotation_logticks(sides = "lr") + 
    xlab("\\# Variables prior to reduction") +
    ylab("\\# Terms in $f$ prior to reduction [log]")

ggsave(plot = plt.DegreeSweepVariables, filename = "DegreeSweepVariables.pdf", path=OUT.PATH, width = COL.WIDTH, height = COL.WIDTH)
do.save.tikz(plt.DegreeSweepVariables, "genDegreeSweepVariables", COL.WIDTH, COL.WIDTH, FALSE)