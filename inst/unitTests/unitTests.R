test_flowCut <- function() {

    res <- flowCut(f=NA)
    checkEquals(NA, res$frame)
    checkEquals(NULL, res$worstChan)
    checkEquals(NULL, res$ind)

    data(flowCutData)

    res <- flowCut(f=flowCutData[[1]], Segment = round(0.6*nrow(flowCutData[[1]])))
    checkEquals(flowCutData[[1]], res$frame)
    checkEquals(NULL, res$worstChan)
    checkEquals(NULL, res$ind)

    res <- flowCut(flowCutData[[2]], MaxPercCut = 0.5)
    checkEquals(7, res$worstChan)
    checkEquals(3238, length(res$ind))
    checkEquals(10189, nrow(res$frame))
}
