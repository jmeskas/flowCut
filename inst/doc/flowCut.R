## ---- label=Load_Library, warning=FALSE, message=FALSE, echo=TRUE, fig=FALSE, results='hide'----
library("flowCut")

## ----label=Load_Data, fig=FALSE, results='hide'--------------------------
data(flowCutData)

## ---- label=SimpleExample_Image, fig=FALSE, echo=TRUE, height=22, width=22, keep.source=TRUE----
res_flowCut <- flowCut(flowCutData[[1]], FileID = 1, Plot = "All")

## ---- label=SimpleExample_Data, fig=FALSE, echo=TRUE---------------------
res_flowCut$data

## ---- label=MaxValleyHgt0p01_Image, fig=FALSE, echo=TRUE, height=22, width=22, keep.source=TRUE----
res_flowCut <- flowCut(flowCutData[[2]], FileID=2, Plot ="All", MaxValleyHgt = 0.01)

## ---- label=MaxValleyHgt0p1_Image, fig=FALSE, echo=TRUE, height=22, width=22, keep.source=TRUE----
res_flowCut <- flowCut(flowCutData[[2]], FileID=3, Plot ="All", MaxValleyHgt = 0.1)

## ---- label=MaxPercCut0p15_Image, fig=FALSE, echo=TRUE, height=22, width=22, keep.source=TRUE----
res_flowCut <- flowCut(flowCutData[[3]], FileID=4, Plot ="All", MaxPercCut = 0.15)

## ---- label=MaxPercCut0p3_Image, fig=FALSE, echo=TRUE, keep.source=TRUE----
res_flowCut <- flowCut(flowCutData[[3]], FileID=5, Plot ="All", MaxPercCut = 0.3)

## ---- label=EventsRemoved_1----------------------------------------------
res_flowCut$data["% of events removed",]

## ---- label=GateLineForce_44_Image, fig=FALSE, echo=TRUE, height=22, width=22, keep.source=TRUE----
res_flowCut <- flowCut(flowCutData[[3]], FileID=6, Plot ="All", GateLineForce=44)

## ---- label=EventsRemoved_2----------------------------------------------
res_flowCut$data["% of events removed",]

## ---- label=AllowFlaggedRerun_UseOnlyWorstChannels_Image_1, fig=FALSE, echo=TRUE, height=22, width=22, keep.source=TRUE----
res_flowCut <- flowCut(flowCutData[[4]], FileID=7, Plot ="All")

## ---- label=AllowFlaggedRerun_UseOnlyWorstChannels_Image_2, fig=FALSE, echo=TRUE, height=22, width=22, keep.source=TRUE----
res_flowCut <- flowCut(flowCutData[[4]], FileID=8, AllowFlaggedRerun = TRUE,
                    UseOnlyWorstChannels = TRUE, Plot="All")

