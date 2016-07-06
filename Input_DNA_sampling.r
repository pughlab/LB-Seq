# |nput DNA sampling power curve
# Prepared to support manuscript by Kis et al.
# "Circulating tumor DNA sequence analysis as an alternative to bone marrow aspirates
# for clinical management of multiple myeloma"

# Plot the probability that a tumour fragment at a given concentration (target.conc) is
# present within an aliquot of total DNA extracted from plasma (ng.dna), given the mass of
# a haploid genome (ng.genome).
target.conc <- 0.0005
ng.dna <- seq(1,100, by=1)
ng.genome <- 0.0035
ng.probs <- pbinom(1, as.integer(ng.dna/ng.genome), prob=target.conc, lower.tail=FALSE)
plot(ng.probs,
	 xlab="Aliquot ng of ctDNA",
	 ylab=paste("Probability aliquot contains tumor fragment at ",
	            signif(target.conc*100, 2),
	            "% of total DNA",
	            sep=""),
	 pch=16,
	 col="black")

#Highlight in red aliquot quantities at desired power levels
powers <- c(0.999, 0.9999)
for(power in powers) {
    ng.at.power <- min(which(ng.probs > power))
    points(ng.at.power, ng.probs[ng.at.power], col="red", pch=16)
    print(paste(ng.at.power,
                " ng has ",
                power*100,
                "% chance of containing a tumour fragment at ",
                signif(target.conc*100, 2),
                "% concentration in plasma",
                sep=""))
}
