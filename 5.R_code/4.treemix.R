# --------------------------- #
#
# Maerl Genomics Project
#
# Description:
# TreeMix analysis
#
# https://github.com/Tom-Jenkins/maerl_genomics
#
# --------------------------- #

# Load plotting functions from TreeMix
source("TreeMix/plotting_funcs.R")

# Plot variance explained
# Obtained by running the variance_explained.sh script
migration = c(0:7)
variance = c(93.13,94.17,98.38,99.00,99.99,100.00,100.00,100.00)
plot(migration, variance, xlab = "Migration event", ylab = "Variance explained (%)")

# Plot tree with no migration parameters
# treemix -seed 123 -i maerl_SNPs_treemix.gz -o results -root Fal
plot_tree("TreeMix/Results/results", cex = 1.5, lwd = 1.5, ybar = 0.2, disp = 0.0005)
# plot_resid("TreeMix/Results/results", pop_order = "TreeMix/maerl.list")

# Plot tree with one migration event
# treemix -seed 123 -i maerl_SNPs_treemix.gz -o resultsM1 -m 1 -root Fal
plot_tree("TreeMix/Results/resultsM1", cex = 1.5, lwd = 1.5, ybar = 0.2, disp = 0.0005)
# plot_resid("TreeMix/Results/resultsM1", pop_order = "TreeMix/maerl.list")

# Plot tree with two migration events
# treemix -seed 123 -i maerl_SNPs_treemix.gz -o resultsM2 -m 2 -root Fal
plot_tree("TreeMix/Results/resultsM2", cex = 1.5, lwd = 1.5, ybar = 0.2, disp = 0.0005)
plot_resid("TreeMix/Results/resultsM2", pop_order = "TreeMix/maerl.list")

# Plot tree with three migration events
# treemix -seed 123 -i maerl_SNPs_treemix.gz -o resultsM3 -m 3 -root Fal
plot_tree("TreeMix/Results/resultsM3", cex = 1.5, lwd = 1.5, ybar = 0.2, disp = 0.0005)
# plot_resid("TreeMix/Results/resultsM3", pop_order = "TreeMix/maerl.list")

# Plot tree with four migration events
# treemix -seed 123 -i maerl_SNPs_treemix.gz -o resultsM4 -m 4 -root Fal
plot_tree("TreeMix/Results/resultsM4", cex = 1.5, lwd = 1.5, ybar = 0.2, disp = 0.0005)
# plot_resid("TreeMix/Results/resultsM4", pop_order = "TreeMix/maerl.list")


# ----------------- #
#
# Generate figures
#
# ----------------- #

# Plot tree with two migration events
png("Figures/Figure3.png", width = 15, height = 8, units = "in", res = 600)
# pdf("Figures/Figure3.pdf", width = 15, height = 8)
par(mfrow = c(1, 2))
plot_tree("TreeMix/Results/resultsM2", cex = 1.5, lwd = 1.5, ybar = 0.4, disp = 0.0005)
mtext("(a)", side = 3, line = 2, at = par("usr")[1], cex = 1.5, font = 2)
plot_resid("TreeMix/Results/resultsM2", pop_order = "TreeMix/maerl.list", cex = 1.2)
mtext("(b)", side = 3, line = 2, at = par("usr")[1], cex = 1.5, font = 2)
dev.off()

# Plot varianced explained and all trees in composite figure
png("Figures/FigureS5.pdf", width = 10, height = 12, units = "in", res = 600)
par(mfrow = c(3, 2))
plot(migration, variance, xlab = "Migration event", ylab = "Variance explained (%)")
plot_tree("TreeMix/Results/results", cex = 1.5, lwd = 1.5, ybar = 0.2, disp = 0.0005)
plot_tree("TreeMix/Results/resultsM1", cex = 1.5, lwd = 1.5, ybar = 0.2, disp = 0.0005)
plot_tree("TreeMix/Results/resultsM2", cex = 1.5, lwd = 1.5, ybar = 0.2, disp = 0.0005)
plot_tree("TreeMix/Results/resultsM3", cex = 1.5, lwd = 1.5, ybar = 0.2, disp = 0.0005)
plot_tree("TreeMix/Results/resultsM4", cex = 1.5, lwd = 1.5, ybar = 0.2, disp = 0.0005)
dev.off()

