#' Trajplot for fitted model
#'
#' plot the longitudinal trajectory of features by local and global clusterings
#'
#' @param fit an objective output from BCC.multi() function
#' @param feature.ind a numeric value indicating which feature to plot.
#'                     The number indicates the order of the feature specified
#'                     in mydat argument of the BCC.multi() function
#' @param which.cluster a character value: "global" or "local", indicating
#'                      whether to plot the trajectory by global cluster or
#'                      local cluster indices
#' @param title Title for the trace plot
#' @param xlab Label for x axis
#' @param ylab Label for y axis
#' @param color Color for the trajplot
#' @return A plot object
#' @examples
#' # get data from the package
#' data(epil1)
#' fit.BCC <- epil1
#' # for local cluster
#' trajplot(fit=fit.BCC,feature.ind=1, which.cluster = "local.cluster",
#'          title= "Local Clustering",xlab="time (months)",
#'          ylab="anxiety",color=c("#00BA38", "#619CFF"))
#'
#' # for global cluster
#' trajplot(fit=fit.BCC,feature.ind=1,
#'          which.cluster = "global.cluster",
#'          title="Global Clustering",xlab="time (months)",
#'          ylab="anxiety",color=c("#00BA38", "#619CFF"))
#'
#' @export
#' @import ggplot2
#' @importFrom graphics par
#' @useDynLib BCClong, .registration=TRUE

trajplot <- function(fit,feature.ind=1,which.cluster = "global.cluster",
                     title = NULL, ylab = NULL,xlab = NULL, color=NULL){
  time.org <- y <- plot.cluster <- id <- NULL
  dat <- fit$dat
  if (which.cluster == "local.cluster"){
    number.cluster <- length(unique(dat[[feature.ind]]$cluster.local))
    per <- round(100*table(fit$cluster.local[[feature.ind]])/fit$N,1)
    dat[[feature.ind]]$plot.cluster <- factor(dat[[feature.ind]]$cluster.local,
            labels=paste("Cluster ",1:number.cluster," (",per,"%",")",sep=""))
  }
  if (which.cluster == "global.cluster"){
    number.cluster <- length(unique(dat[[feature.ind]]$cluster.global))
    per <- round(100*table(fit$cluster.global)/fit$N,1)
    dat[[feature.ind]]$plot.cluster <- factor(dat[[feature.ind]]$cluster.global,
            labels=paste("Cluster ",1:number.cluster," (",per,"%",")",sep=""))}

  if(is.null(color)==FALSE) {
    gp <- ggplot(data = dat[[feature.ind]], aes(x =time.org, y =y,
                                                color=plot.cluster,
                                                linetype=plot.cluster,
                                                fill=plot.cluster))+
          geom_point(size=2,alpha=0.2) +
          geom_line(aes(x = time.org, y = y,group=id,color=plot.cluster),
                    linewidth=1.5,alpha=0.2)+
          geom_smooth(method = "loess", linewidth = 3,se = FALSE,span=2) +
          ggtitle(title) +
          theme_bw() +
          ylab(ylab) + xlab(xlab)+
          theme(legend.position ="bottom",
                legend.title=element_blank(),
                plot.title = element_text(size = 16, face = "bold"),
                axis.text=element_text(size=16),
                axis.title=element_text(size=16),
                legend.text=element_text(size=12),
                axis.text.x = element_text(angle = 0 ),
                strip.text.x = element_text(size = 16, angle = 0),
                strip.text.y = element_text(size = 16,face="bold")) +
          guides(color=guide_legend(nrow=1,byrow=FALSE),
                 linetype=guide_legend(nrow=1,byrow=FALSE),
                 fill=guide_legend(nrow=1,byrow=FALSE)) +
          scale_color_manual(values=color)+  scale_fill_manual(values=color)
  }
  else{
    gp <- ggplot(data = dat[[feature.ind]],
                aes(x =time.org, y =y, color=plot.cluster,
                    linetype=plot.cluster,fill=plot.cluster))+
          geom_point(size=2,alpha=0.2) +
          geom_line(aes(x = time.org, y = y,group=id,color=plot.cluster),
                    linewidth=1.5,alpha=0.2)+
          geom_smooth(method = "loess", linewidth = 3,se = FALSE,span=2) +
          ggtitle(title) +
          theme_bw() +
          ylab(ylab) + xlab(xlab)+
          theme(legend.position ="bottom",
                legend.title=element_blank(),
                plot.title = element_text(size = 16, face = "bold"),
                axis.text=element_text(size=16),
                axis.title=element_text(size=16),
                legend.text=element_text(size=12),
                axis.text.x = element_text(angle = 0 ),
                strip.text.x = element_text(size = 16, angle = 0),
                strip.text.y = element_text(size = 16,face="bold")) +
          guides(color=guide_legend(nrow=1,byrow=FALSE),
                 linetype=guide_legend(nrow=1,byrow=FALSE),
                 fill=guide_legend(nrow=1,byrow=FALSE))
  }
  return(gp)
}

# [END]
