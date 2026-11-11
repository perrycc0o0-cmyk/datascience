# Render manuscript CellChat panels from calculated joint-model tables only.
# Usage: Rscript render_cellchat.R RUN_ROOT OUTPUT_DIR
suppressPackageStartupMessages({library(CellChat); library(grid)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript render_cellchat.R RUN_ROOT OUTPUT_DIR")
tables <- file.path(args[1], "06_CellChat", "manuscript_tables")
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
files <- c(P = "outgoing_selected_comparisons", A = "outgoing_by_receiver_scope",
           B = "joint_differential_edges", C = "outgoing_pathways_increased",
           D = "outgoing_pathways_decreased", E = "outgoing_pathways_heatmap")
dat <- lapply(file.path(tables, paste0(files, ".csv")), read.csv, fileEncoding = "UTF-8")
names(dat) <- names(files)
stopifnot(nrow(dat$P) == 5L, nrow(dat$A) == 10L, nrow(dat$E) == 65L)
A_display <- dat$A

.cm_hash <- function(f) digest::digest(f,algo="sha256",file=TRUE,serialize=FALSE)
.cm_box <- list(A=c(34,0,1202,333),B=c(34,336,595,308),C=c(657,332,610,320),
                D=c(34,650,595,347),E=c(650,650,617,347),P=c(739,1445,687,306))
.cm_labels <- c(GLU_ADGM="GLU-ADGM",GLU_L2_L3="GLU-L2/3",GLU_L4_5="GLU-L4/5",
                GLU_L5="GLU-L5",GLU_L5_6="GLU-L5/6",GLU_L6_CLA="GLU-L6 CLA")
.cm_blue <- "#4F50FF"; .cm_red <- "#FF5367"
.cm_color <- function(d) ifelse(!is.finite(d),"#D9D9D9",ifelse(d>1e-12,.cm_red,ifelse(d< -1e-12,.cm_blue,"#9A9A9A")))
.cm_canvas <- function(box) {
  X <- function(x) (x-box[1])/box[3]
  Y <- function(y) 1-(y-box[2])/box[4]
  point <- function(t,x,y,size=8.5,just="centre",col="black",rot=0,face="plain")
    grid::grid.text(t,X(x),Y(y),just=just,rot=rot,gp=grid::gpar(fontfamily="Arial",fontsize=size,col=col,fontface=face,lineheight=1.05))
  line <- function(x1,y1,x2,y2,col="black",lwd=.65)
    grid::grid.lines(c(X(x1),X(x2)),c(Y(y1),Y(y2)),gp=grid::gpar(col=col,lwd=lwd))
  rect <- function(x1,y1,x2,y2,fill,col=NA,lwd=.5)
    grid::grid.rect(mean(c(X(x1),X(x2))),mean(c(Y(y1),Y(y2))),abs(X(x2)-X(x1)),abs(Y(y2)-Y(y1)),gp=grid::gpar(fill=fill,col=col,lwd=lwd))
  list(X=X,Y=Y,text=point,line=line,rect=rect)
}
.cm_range <- function(values,annotate=FALSE) {
  z <- c(0,values[is.finite(values)]); span <- diff(range(z))
  if(span==0) return(c(0,1))
  c(if(min(z)<0) min(z)-.05*span else 0,max(z)+(if(annotate) .20 else .05)*span)
}
.cm_axis <- function(cv,range,left,right,top,bottom,digits,n=4,zero=FALSE,tick_size=7.2,zero_col="#777777",zero_lwd=.55) {
  map <- function(v) left+(v-range[1])/diff(range)*(right-left)
  ticks <- pretty(range,n=n);ticks <- ticks[ticks>=range[1] & ticks<=range[2]]
  cv$line(left,top,left,bottom);cv$line(left,bottom,right,bottom)
  for(v in ticks) {x<-map(v);cv$line(x,bottom,x,bottom+4,lwd=.55)
    cv$text(if(digits==0) format(round(v),trim=TRUE) else formatC(v,format="f",digits=digits),x,bottom+15,size=tick_size)}
  if(zero && range[1]<0 && range[2]>0) cv$line(map(0),top,map(0),bottom,col=zero_col,lwd=zero_lwd)
  map
}
state <- new.env(parent = emptyenv())
canvas <- function(box, panel) {
  cv <- .cm_canvas(box)
  original_text <- cv$text
  original_line <- cv$line
  cv$text <- function(t, x, y, size = 8.5, just = "centre", col = "black", rot = 0, face = "plain") {
    t <- gsub("<U+2192>", intToUtf8(8594), t, fixed = TRUE)
    t[t == "GLU-L5 | All other subtypes"] <- "GLU-L5 | All other\nsubtypes"
    original_text(t, x, y, size, just, col, rot, if (panel == "P") "plain" else "bold")
  }
  cv$line <- function(x1, y1, x2, y2, col = "black", lwd = .65) {
    original_line(x1, y1, x2, y2, col, max(lwd * 1.6, 2))
  }
  cv
}

draw_A <- function(cv) {
  state$panel <- "A"; cv$text("A", 45, 35, size = 14)
  cv$rect(1040, 11, 1058, 22, .cm_blue); cv$text("Ext", 1071, 17, size = 8, just = "left")
  cv$rect(1113, 11, 1131, 22, .cm_red); cv$text("Ret-Ext", 1144, 17, size = 8, just = "left")
  for (scope in c("Excitatory", "All")) {
    d <- A_display[A_display$target_scope == scope, , drop = FALSE]; d <- d[order(d$display_order), ]
    left <- if (scope == "Excitatory") 184 else 794; right <- if (scope == "Excitatory") 589 else 1200
    map <- .cm_axis(cv, .cm_range(c(d$outgoing_Ext, d$outgoing_RetExt), TRUE), left, right, 44, 269, 2, 4)
    ys <- seq(64, 250, length.out = 5L)
    for (i in seq_len(5L)) {
      cv$text(d$display_label[[i]], left - 7, ys[[i]], size = 7.3, just = "right")
      for (condition in c("Ext", "Ret-Ext")) {
        v <- if (condition == "Ext") d$outgoing_Ext[[i]] else d$outgoing_RetExt[[i]]
        y <- ys[[i]] + if (condition == "Ext") -6.5 else 6.5
        if (!is.finite(v)) next
        cv$rect(map(0), y - 5.7, map(v), y + 5.7, if (condition == "Ext") .cm_blue else .cm_red)
        cv$text(sprintf("%.3f", v), map(v) + 4, y, size = 6.8, just = "left")
      }
    }
    label <- paste0("Outgoing interaction strength\n", if (scope == "Excitatory")
      "(other excitatory subtypes)" else "(all other included neuronal subtypes)")
    cv$text(label, mean(c(left, right)), 316, size = 8.3)
  }
}

nodes <- c("GLU_L2_L3", "GLU_L4_5", "GLU_L5", "GLU_L5_6", "GLU_L6_CLA")
senders <- c("GLU_L5", "GLU_L5_6", "GLU_L6_CLA")
edge <- dat$B[dat$B$selected_for_display, , drop = FALSE]
Ext <- RetExt <- matrix(0, 5, 5, dimnames = list(unname(.cm_labels[nodes]), unname(.cm_labels[nodes])))
Ext[cbind(match(edge$source, nodes), match(edge$target, nodes))] <- edge$weight_Ext
RetExt[cbind(match(edge$source, nodes), match(edge$target, nodes))] <- edge$weight_RetExt
signed <- RetExt - Ext
container <- methods::new("CellChat")
container@net <- list(Ext = list(weight = Ext), RetExt = list(weight = RetExt))
graph <- igraph::graph_from_adjacency_matrix(signed, mode = "directed", weighted = TRUE)
native_layout <- igraph::in_circle(order = c(2, 1, 5, 4, 3))
maxabs <- max(abs(signed))
if (maxabs == 0) stop("No estimable nonzero edges for the differential network.")
lty <- rep(1L, igraph::ecount(graph))
node_delta <- vapply(nodes, function(n) dat$A$delta_RetExt_minus_Ext[
  dat$A$source == n & dat$A$target_scope == "Excitatory"], numeric(1))
palette <- ifelse(is.finite(node_delta) & node_delta > 0, "#FFD9DE",
                  ifelse(is.finite(node_delta) & node_delta < 0, "#E3E6FF", "#F7F7F7"))

draw_native_B <- function(cv) {
  oldpar <- par(no.readonly = TRUE); oldopt <- igraph::igraph_options(edge.lty = lty)
  on.exit({do.call(igraph::igraph_options, oldopt); par(oldpar)}, add = TRUE)
  par(fig = c(cv$X(128), cv$X(448), cv$Y(640), cv$Y(340)), new = FALSE,
    mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), ps = 12, cex = 1, font = 2, xpd = NA, family = "Arial")
  z <- CellChat::netVisual_diffInteraction(object = container, comparison = c(1, 2), measure = "weight",
    sources.use = unname(.cm_labels[senders]), targets.use = NULL, color.use = palette,
    color.edge = c("#FF5367", "#4F50FF"), title.name = "", remove.isolate = FALSE, top = 1,
    weight.scale = TRUE, vertex.weight = rep(1, 5), vertex.weight.max = 1, vertex.size.max = 5,
    vertex.label.cex = 2/3, vertex.label.color = "transparent", edge.weight.max = maxabs,
    edge.width.max = 4, alpha.edge = .6, label.edge = FALSE, edge.curved = .2,
    layout = native_layout, margin = .3, arrow.width = 1, arrow.size = .15)
  stopifnot(inherits(z, "recordedplot"))
  xy <- igraph::norm_coords(scale(igraph::layout_(graph, native_layout)), xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  xx <- (grconvertX(xy[, 1], from = "user", to = "ndc") - cv$X(0)) / (cv$X(1) - cv$X(0))
  yy <- (grconvertY(xy[, 2], from = "user", to = "ndc") - cv$Y(0)) / (cv$Y(1) - cv$Y(0))
  state$nodes <- data.frame(node = nodes, x = xx + ifelse(xy[, 1] >= 0, 18, -18),
    y = yy + ifelse(abs(xy[, 2]) > .3, -sign(xy[, 2]) * 12, 0), just = ifelse(xy[, 1] >= 0, "left", "right"))
}
draw_B_text <- function(cv) {
  state$panel <- "B"; cv$text("B", 45, 346, size = 14, face = "bold")
  for (i in seq_len(nrow(state$nodes))) {
    p <- state$nodes[i, ]; cv$text(unname(.cm_labels[p$node]), p$x, p$y, size = 8.5, just = p$just, face = "bold")
  }
}
draw_CD <- function(panel, cv) {
  state$panel <- panel; d <- dat[[panel]]; d <- d[order(d$display_order), ]
  if (panel == "C") {left <- 794; right <- 1200; top <- 355; bottom <- 590; ys <- seq(370, 577, length.out = 9); h <- 18; labely <- 634}
  else {left <- 184; right <- 589; top <- 671; bottom <- 922; ys <- seq(682, 913, length.out = 12); h <- 14; labely <- 964}
  cv$text(panel, if (panel == "C") 667 else 45, if (panel == "C") 346 else 663, size = 14)
  map <- .cm_axis(cv, .cm_range(d$plot_value), left, right, top, bottom, if (panel == "C") 4 else 3, 4, TRUE)
  for (i in seq_len(nrow(d))) {
    cv$text(d$display_label[[i]], left - 7, ys[[i]], size = 8, just = "right")
    cv$rect(map(0), ys[[i]] - h/2, map(d$plot_value[[i]]), ys[[i]] + h/2, .cm_color(d$delta_RetExt_minus_Ext[[i]]))
  }
  cv$text(paste0("Outgoing pathway strength difference\n(",
    if (panel == "C") "Ret-Ext - Ext" else "Ext - Ret-Ext", ")"), mean(c(left, right)), labely, size = 8.3)
}

draw_E <- function(cv){
  d<-dat$E;lim<-max(abs(c(1, d$display_value)), na.rm = TRUE);pal<-grDevices::colorRampPalette(c("#0000FF","#FFFFFF","#FF0000"))(2001);zlim<-asinh(lim/5);T<-function(v)asinh(v/5)
  color<-function(v)if(!is.finite(v))"#D9D9D9" else pal[round((T(v)/zlim+1)*1000)+1L]
  cv$text("E",667,663,size=14);left<-815;right<-1178;top<-685;bottom<-896;cw<-(right-left)/13;ch<-(bottom-top)/5
  for(i in seq_len(nrow(d))){x<-left+(d$column_order[i]-.5)*cw;y<-top+(d$row_order[i]-.5)*ch;v<-d$display_value[i];cv$rect(x-cw/2,y-ch/2,x+cw/2,y+ch/2,color(v));if(is.finite(v))cv$text(sprintf("%.2f",v),x,y,size=4.4,col=if(abs(T(v))/zlim>.65)"white"else"black")}
  for(i in 1:5)cv$text(d$row_label[d$row_order==i][1],809,top+(i-.5)*ch,size=6.6,just="right");for(j in 1:13)cv$text(d$pathway[d$column_order==j][1],left+(j-.5)*cw,915,size=7.2,rot=55,just="right")
  xx<-1206;y1<-711;y2<-875;ww<-13;for(i in 1:200)cv$rect(xx-ww/2,y2-(i-1)*(y2-y1)/200,xx+ww/2,y2-i*(y2-y1)/200,pal[round((i-.5)*2000/200)+1L]);cv$line(xx-ww/2,y1,xx-ww/2,y2,lwd=.65);cv$line(xx+ww/2,y1,xx+ww/2,y2,lwd=.65)
  grid::grid.polygon(cv$X(c(xx-ww/2,xx,xx+ww/2)),cv$Y(c(y1,y1-8,y1)),gp=grid::gpar(fill="#FF0000",col="black",lwd=2));grid::grid.polygon(cv$X(c(xx-ww/2,xx,xx+ww/2)),cv$Y(c(y2,y2+8,y2)),gp=grid::gpar(fill="#0000FF",col="black",lwd=2))
  for(t in c(-50,-10,-5,0,5,10,50)[abs(c(-50,-10,-5,0,5,10,50)) <= lim])cv$text(format(t,trim=TRUE),1219,mean(c(y1,y2))-T(t)/zlim*(y2-y1)/2,size=6.2,just="left");cv$text("Difference (\u00D710\u00B3)",1258,793,size=7.2,rot=90);cv$text("asinh scale",1246,793,size=5.5,rot=90)
}

P_draw<-function(cv){
  d<-dat$P[order(dat$P$display_order),];cv$text('P',748,1474,size=10.5,face='bold')
  vals<-d$percent_change;rg<-c(min(c(0,vals))-6,max(c(0,vals))+14)
  boldcv<-cv;boldcv$line<-function(x1,y1,x2,y2,col='black',lwd=.65)cv$line(x1,y1,x2,y2,col=col,lwd=lwd*1.6)
  boldcv$text<-function(t,x,y,size=8.5,just='centre',col='black',rot=0,face='plain')cv$text(t,x,y,size=size,just=just,col=col,rot=rot,face='bold')
  map<-.cm_axis(boldcv,rg,1007,1376,1465,1686,0,5,TRUE,tick_size=7,zero_col='black',zero_lwd=.92)
  ys<-seq(1493,1661,length.out=5)
  # The descriptive vertical heading is omitted to keep bold row labels unobstructed.
  for(i in 1:5){
    cv$text(d$display_label[i],999,ys[i],size=7.1,just='right',face='bold');v<-d$percent_change[i]
    cv$rect(map(0),ys[i]-16,map(v),ys[i]+16,.cm_color(v))
    if(v<0&&abs(v)<10)cv$text(sprintf('%+.1f%%',v),map(v)-5,ys[i],size=7,just='right',face='bold') else
      if(v<0)cv$text(sprintf('%+.1f%%',v),map(v)+2,ys[i],size=7,just='left',col='white',face='bold') else
        cv$text(sprintf('%+.1f%%',v),map(v)+5,ys[i],size=7,just='left',face='bold')
  }
  cv$text('Percent change: Ret-Ext vs Ext',1194,1732,size=8,face='bold')
}

paint <- function(panel, box) {
  cv <- canvas(box, panel)
  if (panel %in% c("B", "combined")) draw_native_B(cv) else grid.newpage()
  pushViewport(viewport(clip = "on"))
  on.exit(popViewport(), add = TRUE)
  switch(panel, P = P_draw(cv), A = draw_A(cv), B = draw_B_text(cv),
         C = draw_CD("C", cv), D = draw_CD("D", cv), E = draw_E(cv),
         combined = {
           draw_A(cv); draw_B_text(cv); draw_CD("C", cv); draw_CD("D", cv); draw_E(cv)
         })
}

stems <- c(P = "outgoing_selected_comparisons", A = "outgoing_by_receiver_scope",
           B = "differential_communication_network", C = "outgoing_pathways_increased",
           D = "outgoing_pathways_decreased", E = "outgoing_pathways_heatmap",
           combined = "supplementary_communication_panels")
if (any(file.exists(outer(file.path(out, stems), c(".pdf", ".png", ".svg"), paste0)))) {
  stop("An output figure already exists. Choose a new output directory.")
}
for (panel in names(stems)) {
  box <- if (panel == "combined") c(0, 0, 1267, 997) else .cm_box[[panel]]
  px <- if (panel == "P") 1426 else 1267
  width <- box[[3L]] * 240 / px / 25.4
  height <- box[[4L]] * 240 / px / 25.4
  stem <- file.path(out, stems[[panel]])
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width, height, family = "Arial", bg = "white")
  tryCatch(paint(panel, box), finally = dev.off())
  ragg::agg_png(paste0(stem, ".png"), width, height, units = "in", res = 600, background = "white")
  tryCatch(paint(panel, box), finally = dev.off())
  svglite::svglite(paste0(stem, ".svg"), width, height, bg = "white")
  tryCatch(paint(panel, box), finally = dev.off())
}
