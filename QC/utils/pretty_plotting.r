library(ggplot2)
library(dplyr)
library(ggsci)
library(ggrastr)
library(ggrepel)
library(data.table)
library(latex2exp)

chr_lengths_38 <- c(248956422, 242193529, 198295559,190214555,181538259,170805979,159345973,
    145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
    83257441,80373285,58617616,64444167,46709983,50818468,156040895)

chr_lengths_37 <- c(249250621, 243199373,198022430,191154276,180915260,171115067,159138663,
    146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,
    81195210,78077248,59128983,63025520,48129895,51304566,155270560)

create_pretty_forest <- function(df, title, save_figure=FALSE, file='file_out',
    width=160, height=90, scaling=1, y_label='', pass_order=NULL, print_p=TRUE,
    horizontal=TRUE, hline_at=0, ggplot_theme=theme_classic)
{
    if (!is.null(pass_order)) {
        df$label <- factor(df$label, levels = df$label[pass_order])
    }

    p <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
        geom_pointrange() + 
        geom_hline(yintercept=hline_at, lty=2)

    if (horizontal==TRUE) p <- p + coord_flip()

    p <- p + labs(title=title, y=y_label, x='') +
        ggplot_theme()

    if (horizontal == TRUE) {
        p <- p + geom_text(aes(y=mean, label=p_vals), vjust=-1)
    } else {
        p <- p + geom_text(aes(y=mean, label=p_vals), hjust=-0.2)
    }

    if (print_p)  print(p)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm', dpi=500)
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
    }
    
    return(p)
}

create_pretty_hist <- function(df, aest, x_label, threshold=NULL, threshold_max=NULL,
    file='file_out', title='', binwidth=0.002, width=160, height=90, scaling=1,
    save_figure=FALSE, xlim=NULL, key_label='', print_p=TRUE, title.hjust=0.5, ggplot_theme=theme_classic)
{
    p <- ggplot(df, aest)
    geom_histogram(binwidth=binwidth, fill='#aec7e8', color='#1f77b4')
    ylim <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range
    p <- ggplot(df, aest)

    if (is.null(aest$fill)) {
        p <- p + geom_histogram(binwidth=binwidth, fill='#aec7e8', color='#1f77b4')
    } else {
        p <- p + geom_histogram(binwidth=binwidth, color='grey50')
    }

    if (!is.null(xlim)) {
        p <- p + coord_cartesian(xlim=xlim)
    } 

    p <- p + labs(x=x_label, y='Count', title=title, fill=key_label) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
    scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
    ggplot_theme() +
    theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
    axis.title.y = element_text(margin=ggplot2::margin(r=10)),
    plot.title = element_text(hjust=title.hjust))

    if (!is.null(threshold)) p <- p + geom_vline(xintercept=threshold, linetype='dashed')
    if (!is.null(threshold_max)) p <- p + geom_vline(xintercept=threshold_max, linetype='dashed')

    
    if (print_p) print(p)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm', dpi=500)
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
    }
    
    return(p)
}

create_pretty_density <- function(df, aest, x_label, threshold=NULL, threshold_max=NULL,
    file='file_out', title='', binwidth=0.002, width=160, height=90, scaling=1,
    save_figure=FALSE, xlim=NULL, key_label='', print_p=FALSE, title.hjust=0.5,
    ggplot_theme=theme_classic)
{
    p <- ggplot(df, aest)
    ylim <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range
    p <- ggplot(df, aest)

    if (is.null(aest$fill)) {
        p <- p + geom_density(fill='#aec7e8', color='#1f77b4')
    } else {
        p <- p + geom_density(color='grey50')
    }

    if (!is.null(xlim)) {
        p <- p + coord_cartesian(xlim=xlim)
    } 

    p <- p + labs(x=x_label, y='Count', title=title, fill=key_label) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
    scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
    ggplot_theme() +
    theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
    axis.title.y = element_text(margin=ggplot2::margin(r=10)),
    plot.title = element_text(hjust=title.hjust))

    if (!is.null(threshold)) {
        p <- p + geom_vline(xintercept=threshold, linetype='dashed') #+
        # annotate('text', label=paste0(x_label, ' = ', threshold), x=threshold, y=0.8*ylim[2], angle=90, vjust=-1)
    }
    if (!is.null(threshold_max)) {
        p <- p + geom_vline(xintercept=threshold_max, linetype='dashed') #+
        # annotate('text', label=paste0(x_label, ' = ', threshold_max), x=threshold_max, y=0.8*ylim[2], angle=90, vjust=2)
    }

    if (print_p) print(p)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm', dpi=500)
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
    }
    
    return(p)
}

create_pretty_boxplots <- function(df, aes, aes_col, threshold=NULL,
    threshold_max=NULL, file='file_out', title='', x_label='', y_label='',
    key_label='', xlim=NULL, legend=FALSE, save_figure=FALSE, 
    width=160, height=90, scaling=1, facet=FALSE, facet_grid=NULL, jitter_size=0.5,
    outlier.shape=NA, n_ticks=10, print_p=FALSE, alpha=0.6, title.hjust=0.5, ggplot_theme=theme_classic)
{
    p = ggplot(df, aes) +
        geom_boxplot(outlier.shape=outlier.shape, coef=0, color='grey50', fill='grey95', show.legend=FALSE) + 
        geom_jitter_rast(width=0.2, height=0, size=jitter_size, aes_col, show.legend=legend, alpha=alpha, stroke=0.05, raster.dpi=500) + 
        coord_flip(ylim=xlim) +
        labs(title=title, x=y_label, y=x_label, color=key_label) + 
        scale_color_d3('category20') +
        scale_y_continuous(breaks = scales::pretty_breaks(n=n_ticks)) +
        guides(color = guide_legend(override.aes = list(size=2))) +
        ggplot_theme() +
        theme(axis.title.x = element_text(margin = ggplot2::margin(t=10)),
              plot.title = element_text(hjust=title.hjust))

    if (!is.null(threshold)) {
        p <- p + geom_hline(yintercept=threshold, linetype='dashed')
    }
    if (!is.null(threshold_max)) {
        p <- p + geom_hline(yintercept=threshold_max, linetype='dashed')
    }
    if (facet){
        p <- p + facet_grid
    }

    if (print_p) print(p)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm', dpi=500)
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
    }
    return(p)
}

create_pretty_cumulative <- function(df, aes, x_label, threshold, threshold_max=NULL,
    file='file_out', title='', width=160, height=90, scaling=1, save_figure=FALSE,
    xlim=c(0,1), key_label='', print_p=FALSE, title.hjust=0.5, ggplot_theme=theme_classic)
{
    # These next ones are cdfs.
    p = ggplot(df, aes) + 
        stat_ecdf(geom='line', pad=FALSE) +
        geom_vline(xintercept=threshold, linetype='dashed') +
        coord_cartesian(xlim=xlim) +
        labs(x=x_label, y='Cumulative Proportion', title=title, color=key_label) +
        scale_color_d3('category10') +
        scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
        ggplot_theme() +
        theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
              axis.title.y = element_text(margin=ggplot2::margin(r=10)),
              plot.title = element_text(hjust=title.hjust))

    # p <- p + annotate('text', label=paste0(x_label, ' = ', threshold), x=threshold, y=mean(ylim), angle=90, vjust=-1)
    # ylim <- ggplot_build(p)$panel$ranges[[1]]$y.range
    # print(ylim)

    if (!is.null(threshold_max)) {
        p <- p + geom_vline(xintercept=threshold_max, linetype='dashed')# +
        # annotate('text', label=paste0(x_label, ' = ', threshold_max), x=threshold_max, y=mean(ylim), angle=90, vjust=2)
    }
    if (print_p) print(p)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm', dpi=500)
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
    }
    return(p)
}

create_pretty_scatter <- function(dt, aes, file='file_out', save_figure=FALSE, 
    key_label='', title='', limits=NULL, width=160, height=90, presentation=FALSE,
    add_final_layer=FALSE, final_layer=NULL, n_x_ticks=10, n_y_ticks=10, x_label=NULL,
    y_label=NULL, print_p=FALSE, gradient=FALSE, midpoint=0, gradient_title="", title.hjust=0.5,
    legend_max_min=NULL, n_legend_step=4, xlim=NULL, ylim=NULL, ggplot_theme=theme_classic, alpha=0.6, size=1,
    include_qq_ribbon=FALSE, aes_ribbon=NULL, ribbon_p=0.95)
{
    p <- ggplot(dt, aes)
    print(dt)
    if (include_qq_ribbon) {
        p <- p + geom_ribbon(aes_ribbon, fill="grey80", color="grey80")
    }
    
    if ("size" %in% names(aes)) {
        p <- p + geom_point_rast(alpha=alpha, raster.dpi=500)
    } else {
        p <- p + geom_point_rast(alpha=alpha, size=size, raster.dpi=500)
    }
    p <- p +
        scale_color_d3('category20', limits=limits) +
        labs(title=title, color=paste0(key_label)) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=n_x_ticks)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=n_y_ticks)) +
        ggplot_theme() +
        theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
              axis.title.y = element_text(margin=ggplot2::margin(r=10)),
              plot.title = element_text(hjust=title.hjust))

    if (add_final_layer) {
        cat("Adding final layer...\n")
        p <- p + guides(fill=FALSE) + geom_point_rast(mapping=aes, data=final_layer, raster.dpi=500) + scale_color_d3('category20')
    }

    if (gradient) {
        print(gradient_title)
        if(!is.null(legend_max_min)) {
            p <- p + scale_color_gradient2(
                low="blue", high="red", mid='grey50', midpoint=midpoint,
                name=gradient_title,
                breaks=seq(legend_max_min[1], legend_max_min[2], (legend_max_min[2]-legend_max_min[1])/n_legend_step))
        } else {
            p <- p + scale_color_gradient2(
                low="blue", high="red", mid='grey50', midpoint=midpoint,
                name=gradient_title)
        }
    }
    
    cat("Adding axis labels...\n")
    if (!is.null(x_label)) p <- p + labs(x=x_label)
    if (!is.null(y_label)) p <- p + labs(y=y_label)
    if (!is.null(xlim)) p <- p + xlim(xlim[1], xlim[2])
    if (!is.null(ylim)) p <- p + ylim(ylim[1], ylim[2])

    cat("Saving figure...\n")
    if (print_p) print(p)

    if (save_figure) {
        
        if (presentation == TRUE) {
            width <- 160
            height <- 90
            ggsave(paste0(file, '.jpg'), p, width=width, height=height, units='mm', dpi=500)
        } else {
            ggsave(paste0(file, '.jpg'), p, width=width, height=height, units='mm', dpi=500)
            ggsave(paste0(file, '.pdf'), p, width=width, height=height, units='mm')
        }
    
    }

    return(p)
}

create_pretty_qq_plot <- function(dt, aes, file='file_out', save_figure=FALSE,
    plot_title='', limits=NULL, width=110, height=100, n_x_ticks=10, n_y_ticks=10,
    x_label=TeX("$-\\log_{10}(\\mathit{p}_{permutation})$"), 
    y_label=TeX("$-\\log_{10}(\\mathit{p}_{observed})$"),
    n_to_include=NULL, cex_labels=1, print_p=TRUE, gradient=FALSE,
    gradient_title="", title.hjust=0.5, legend_max_min=NULL, n_legend_step=4,
    xlim=NULL, ylim=NULL, y_x_col="grey40", ggplot_theme=theme_classic, alpha=0.6,
    include_qq_ribbon=TRUE, aes_ribbon=NULL, ribbon_p=0.95, key_cols=c("pval"))
{
    cat("Creating scatter-plot...\n")
    dt <- data.table(dt)
    setkeyv(dt, cols=key_cols)
    print(dt)
    p <- create_pretty_scatter(dt, aes, file=file, save_figure=FALSE,
        title=plot_title, limits=limits,
        width=width, height=height, n_x_ticks=n_x_ticks, n_y_ticks=n_y_ticks,
        x_label=x_label, y_label=y_label, gradient=gradient, gradient_title=gradient_title,
        title.hjust=title.hjust, legend_max_min=legend_max_min, n_legend_step=n_legend_step,
        xlim=xlim, ylim=ylim, ggplot_theme=ggplot_theme, alpha=alpha,
        include_qq_ribbon=include_qq_ribbon, aes_ribbon=aes_ribbon, ribbon_p=ribbon_p)

    cat("Adding y=x line...\n")
    p <- p + geom_abline(intercept=0, slope=1, color=y_x_col) #+ coord_fixed()

    if (!is.null(n_to_include)) {
        cat("Adding labels...\n")
        p <- p + geom_label_repel(data=dt[(nrow(dt)-n_to_include+1):nrow(dt), ],
            aes(label=labels), box.padding = 0.4, label.padding=0.1, point.padding = 0.2,
            color = 'grey30', segment.color = 'grey50', max.overlaps=Inf,
            size=cex_labels, segment.size=0.1, show.legend = FALSE)
    }

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width, height=height, units='mm', dpi=500)
        ggsave(paste0(file, '.pdf'), p, width=width, height=height, units='mm')
    }
    
    cat("Created scatter-plot...\n")
    if (print_p) print(p)

    return(p)
}

create_pretty_pointrange <- function(dt, y, ymin, ymax, colors="phenotype",
    break_on="phenotype", x_labels="phenotype", threshold=NULL, file='file_out',
    title='', width=165, height=110, save_figure=FALSE, xlim=NULL, key_label='',
    print_p=FALSE, title.hjust=0.5, split_on=NULL, spacing=7, x_label="", y_label="",
    title_size=NULL, manual_colours=NULL, colour_levels=NULL, break_on_levels=NULL,
    ggplot_theme=theme_classic, ylim=NULL)
{   
    if(!is.null(colour_levels)) {
        dt[[colors]] <- factor(dt[[colors]], levels=colour_levels)
    } else {
        dt[[colors]] <- factor(dt[[colors]])
    }

    if(!is.null(break_on_levels)) {
        dt[[break_on]] <- factor(dt[[break_on]], levels=break_on_levels)
    } else {
        dt[[break_on]] <- factor(dt[[break_on]])
    }

    dt <- dt %>% arrange(dt[[break_on]], dt[[colors]])

    dt1 <- cbind(dt, aux=rep(1,length(dt[,1]))) 
    dt1 <- within(dt1, {aux = unlist(by(aux, dt1[[break_on]], cumsum))})
    dt1$aux <- dt1$aux + as.numeric(factor(dt1[[break_on]])) * (length(unique(factor(dt1[[break_on]]))) + spacing)

    print(head(dt1))

    define_breaks <- function(dt, split_on) {
        breaks <- c()
        for (split in unique(dt[[split_on]])) {
            breaks <- c(breaks, mean(unlist(dt %>% filter(dt[[split_on]] == split) %>% select(aux))))
        }
        return(breaks)
    }

    p <- ggplot(dt1, aes(x=dt1$aux, y = dt1[[y]], ymin = dt1[[ymin]], ymax = dt1[[ymax]])) +
    geom_pointrange(aes(color = dt1[[colors]]), size = 0.5) + 
    scale_x_continuous("", breaks=define_breaks(dt1, break_on), labels=unique(dt1[[break_on]])) +
    labs(x=x_label, y=y_label, title=title, color=key_label)

    if (!is.null(manual_colours)) {
        p <- p + scale_color_manual(labels=colour_levels, values=manual_colours)
    } else {
        p <- p + scale_color_d3('category10')
    }
    p <- p +
    ggplot_theme() +
    theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
          axis.title.y = element_text(margin=ggplot2::margin(r=10)),
          plot.title = element_text(hjust=title.hjust))

    if(!is.null(threshold)) {
        p <- p + geom_hline(yintercept=threshold, linetype='dashed', col="grey40")
    }

    if (!is.null(title_size)) 
        p <- p + theme(plot.title = element_text(size=title_size))

    if (!is.null(ylim)) p <- p + ylim(ylim[1], ylim[2])

    if(print_p) print(p)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width, height=height, units='mm', dpi=500)
        ggsave(paste0(file, '.pdf'), p, width=width, height=height, units='mm')
    }
    return(p)
}

get_start_and_end <- function(chr_lengths) {
    start <- rep(0, length(chr_lengths))
    start[1] <- 1
    end <- rep(0, length(chr_lengths))
    end[1] <- chr_lengths[1]
    for(chr in 2:length(chr_lengths)) {
        start[chr] <- start[chr-1] + chr_lengths[chr-1]
        end[chr] <- end[chr-1] + chr_lengths[chr]
    }
    return(list(start=start, end=end))
}

make_manhattan_plot = function(contigs, positions, pvals, log_OR=NULL, labels=NULL, size_by_p=FALSE, buffer=100000000, title='', threshold=5, chr_lengths=chr_lengths_38,
  colour_aes=NULL, log_p_vals=FALSE, significance_T=5e-8, ggplot_theme=theme_bw, two_tone=TRUE, by_OR=FALSE, colour_1='#2b59a1', colour_2='#5fb756',
  save_figure=FALSE, file='file_out', scaling=1, width=230, height=100, title_size=NULL, minus_log_p_max=NULL, print_p=FALSE)
{  
    contigs_ <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
    start_end <- get_start_and_end(chr_lengths)
    dt_contigs <- data.frame(contig=contigs_, start=start_end$start, end=start_end$end) %>%
    mutate(middle = floor(start + (end-start)/2),
           length = (end-start)) %>%
    mutate(shifted_position=middle + (contig - 1) * buffer)

    dt_plot <- data.frame(contig=contigs, position=as.integer(positions), pval=as.numeric(pvals), labels=labels) %>%
    mutate(x = dt_contigs[gsub('X', '23', contig), 'start'] + position + (as.integer(gsub('X', '23', contig))-1)*buffer)

    # Include two tone chromosome plotting
    if (two_tone) {
        dt_plot <- dt_plot %>% mutate(colour=ifelse((as.integer(gsub('X', '23', contig)) %% 2) == 0, colour_1, colour_2))
    } else {
        dt_plot <- dt_plot %>% mutate(colour=colour_1) 
    }

    if (by_OR & !is.null(log_OR)) {
        dt_plot$log_OR <- log_OR
        dt_plot <- dt_plot %>% mutate(colour = case_when(
            dt_plot$log_OR < -1 ~ "blue3",
            ((dt_plot$log_OR > -1) & (dt_plot$log_OR <= 0)) ~ "cornflowerblue",
            ((dt_plot$log_OR > 0) & (dt_plot$log_OR <= 1)) ~ "indianred3",
            dt_plot$log_OR > 1 ~ "red",
            TRUE ~ "grey40")
        )
        dt_plot$colour <- factor(dt_plot$colour, levels = c("red", "indianred3", "cornflowerblue", "blue3"))
    }

    # Were log p-values passed?
    if(!log_p_vals) {
        dt_plot <- dt_plot %>% mutate(y = -log10(pval)) %>% select(x, y, colour, labels)
    } else {
        dt_plot <- dt_plot %>% mutate(y= pval) %>% select(x, y, colour, labels)
    }

    if (size_by_p) {
        dt_plot <- dt_plot %>% mutate(size = case_when(
            dt_plot$y > 3 ~ 20,
            ((dt_plot$y > -log10(0.05)) & (dt_plot$y <= 3)) ~ 5,
            TRUE ~ 0.5)
        )
    }

    if (size_by_p) {
        p <- ggplot(dt_plot, aes(x=x,y=y,label=labels, col=colour, size=size)) + geom_point_rast()
    } else {
        p <- ggplot(dt_plot, aes(x=x,y=y,label=labels, col=colour)) + geom_point_rast(size=0.5)
    }

    p <- p + geom_hline(yintercept=-log10(significance_T), color='#E15759', linetype='dashed') +
        scale_x_continuous(breaks=dt_contigs$shifted_position, labels=gsub(23, 'X', dt_contigs$contig)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
        labs(x='Chromosome', y=TeX("$-\\log_{10}(\\mathit{p})$"), title=title) + ggplot_theme()

    if ((by_OR == FALSE) & (size_by_p==FALSE)) {
        p <- p + theme(legend.position = "none")
    }

    if (two_tone) { p <- p + scale_color_manual(values=c(colour_1,colour_2)) }
    if (by_OR) { 
        p <- p + scale_colour_manual(name="Odds ratio",
            values = levels(dt_plot$colour), labels=c("> 10", "1 - 10", "0.1 - 1", "< 0.1"))
    }

    if (!is.null(labels))
        p <- p + geom_label_repel(
            data=subset(dt_plot, y > threshold), size = 5,
            aes(label=labels), color='grey30', box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50')

    if (!is.null(title_size)) 
        p <- p + theme(plot.title = element_text(size=title_size))

    if (!is.null(minus_log_p_max))
        p <- p + ylim(0, minus_log_p_max)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
    }

    if (print_p) { print(p) }

    return(list(p=p, dt=dt_plot))
}