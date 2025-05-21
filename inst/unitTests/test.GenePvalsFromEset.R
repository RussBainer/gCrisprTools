library(RUnit)
library(gCrisprTools)

test.ct.applyAlpha <- function() {
    p = seq(0, 1, length.out=20)
    fc = seq(-3, 3, length.out=20)
    fc[2] = NA
    fc[3] = -20
    stats = data.frame(
        Depletion.P=p,
        Enrichment.P=rev(p),
        fc=fc
    )

    checkIdentical(
        data.frame(
            Depletion.P=p,
            Enrichment.P=rev(p),
            fc=fc,
            scores.deplete=c(0.10,rep(1.0,19)),
            scores.enrich=c(rep(1.0,18),0.10,0.05)
        ),
        ct.applyAlpha(stats,scoring="combined")
    )

    checkIdentical(
        data.frame(
            Depletion.P=p,
            Enrichment.P=rev(p),
            fc=fc,
            scores.deplete=c(0.05,0.10,rep(1.0,18)),
            scores.enrich=c(rep(1.0,18),0.10,0.05)
        ),
        ct.applyAlpha(stats,scoring="pvalue")
    )
    
    checkIdentical(
        data.frame(
            Depletion.P=p,
            Enrichment.P=rev(p),
            fc=fc,
            scores.deplete=c(0.10,1.0,0.05,rep(1.0,17)),
            scores.enrich=c(rep(1.0,18),0.10,0.05)
        ),
        ct.applyAlpha(stats,scoring="fc")
    )

}
