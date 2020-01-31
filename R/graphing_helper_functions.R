#' @importFrom grid unit
#' @importFrom grid grid.text
#' @importFrom grid gpar

## add gridlines to line graph or scatterplot
.plotGridlines <- function(
    ylim, yaxis_increment, show_gridlines, gridline_colour="#aaaaaa")
{
    if (show_gridlines == TRUE) {
        if (yaxis_increment >= ylim[2] | yaxis_increment <= 0) {
            stop(
                "Invalid yaxis_increment: should be more than zero but less
                than the upper limit of ylim")
        }

        ## draw grid lines
        for (h in seq(ylim[1], ylim[2], by=yaxis_increment)) {
            circos.lines(
                CELL_META$cell.xlim, c(h, h), lty=3, col=gridline_colour)
        }
    }
}

## add y axis to line graph or scatterplot
.plotYaxis <- function(
    show_yaxis, ylim, yaxis_increment, yaxis_label_size=0.6,
    yaxis_tick_size=0.5, yaxis_side="left",
    yaxis_location=CELL_META$sector.index, yaxis_colour="black")
{
    if (show_yaxis == TRUE) {
        if (yaxis_increment >= ylim[2] | yaxis_increment <= 0) {
            stop(
                "Invalid yaxis_increment: should be more than zero but less
                than the upper limit of ylim")
        }

        circos.yaxis(
            side=yaxis_side, sector.index=yaxis_location, col=yaxis_colour,
            at=round(seq(ylim[1], ylim[2], by=yaxis_increment), digits=2),
            labels.cex=yaxis_label_size, labels.col=yaxis_colour,
            tick.length=convert_x(
                yaxis_tick_size, "mm", sector.index=yaxis_location))
    }
}

## add title to the graph
.plotTitle <- function(
    title=NULL, title_x_position=0.5, title_y_position=0.9,
    title_font_face="bold", title_size=1.1, title_colour="black")
{
    if (!is.null(title)) {
        x <- grid::unit(title_x_position, "npc")
        y <- grid::unit(title_y_position, "npc")
        grid::grid.text(
            label=title, x=x, y=y,
            gp=grid::gpar(
                col=title_colour, cex=title_size, font=title_font_face))
    }
}

## helper function to reverse the direction of a sector
## from: https://github.com/jokergoo/circlize/issues/124
.reverseXaxis <- function(x)
{
    x_range <- CELL_META$xlim
    x_range[2] - x + x_range[1]
}


## helper function to plot text on the next track inwards
.plotNewTrackText <- function(
    label, x_start, x_end, this_track, sector_index, feature_label_size,
    label_track_height, max_track)
{

    ## labels are centred in the middle of the track. Thus, if we have a mix of
    ## long and short labels the short ones will be far from their feature and
    ## the long ones overlapping. To fix this, make an adjustment based on the
    ## length of the label
    if (this_track == max_track) {
        label_adjust <- (nchar(label)^0.7) * feature_label_size * 0.035
    } else {
        label_adjust <- (nchar(label)^0.7) * feature_label_size * 0.035 * 6
    }

    ## try to plot text on the next track in. if this fails (the track doesn't)
    ## exist, create the new track & retry
    tryCatch(expr={
        circos.text(
            x=((x_start + x_end) / 2), y=(CELL_META$ylim[2] - label_adjust),
            track.index=this_track + 1, facing="clockwise", niceFacing=TRUE,
            sector.index=sector_index, labels=as.character(label),
            cex=feature_label_size)
    },

    error = function(e){
        if (this_track == max_track) {
            label_track_height <- CELL_META$yplot[1] - 0.1
        }
        circos.track(
            ylim=c(0, 1), bg.border=NA, bg.col=NA,
            track.height=label_track_height)
        circos.text(
            x=((x_start + x_end) / 2), y=(CELL_META$ylim[2] - label_adjust),
            labels=as.character(label), facing="clockwise", niceFacing=TRUE,
            cex=feature_label_size, track.index=this_track + 1,
            sector.index=sector_index)
    })
}

.plotIdeogram <- function(
    current_sectors, coverage_data=NULL, coverage_rectangle=NULL,
    current_track, colour, border_colour, sector_labels=TRUE,
    sector_label_size, sector_label_colour)
{

    for (i in seq_along(current_sectors)) {
        ## get the correct sector
        circos.update(
            sector.index=current_sectors[i], bg.border=NA, bg.col=NA,
            track.index=current_track)

        ## coverage rectangle
        if (as.character(current_sectors[i]) %in%
            as.character(coverage_rectangle)) {
            this_sector_coverage <- .getThisSectorCoverage(
                coverage_data, current_sectors[i])
            circos.lines(
                x=this_sector_coverage$start, col=colour[i], area=TRUE,
                y=this_sector_coverage$coverage, border=border_colour[i],
                sector.index=current_sectors[i], track.index=current_track)

        } else {
            ## regular rectangle
            circos.rect(
                xleft=CELL_META$xlim[1], xright=CELL_META$xlim[2],
                col=colour[i], ybottom=CELL_META$ylim[1],
                border=border_colour[i], ytop=CELL_META$ylim[2] / 2,
                sector.index=current_sectors[i], track.index=current_track)
        }
        ## add label if needed
        if (sector_labels == TRUE) {
            circos.text(
                x=mean(CELL_META$xlim), sector.index=current_sectors[i],
                cex=sector_label_size,
                y=mean(c(CELL_META$ylim[1], CELL_META$ylim[2] / 2)),
                labels=as.character(CELL_META$sector.index),
                col=sector_label_colour, facing="bending.inside",
                niceFacing=TRUE, track.index=current_track)
        }
    }
}

## helper function to make the ylim and add a new track
.plotNewTrack <- function(
    current_sectors, custom_ylim, coverage_rectangle,
    coverage_data, track_height=0.15)
{
    ## if it's a coverage rectangle the ylim needs to match the coverage
    ## max/min otherwise it should just be c(0,1)
    if (any(as.character(current_sectors) %in% as.character(
        coverage_rectangle)) & is.null(custom_ylim)) {
        this_sector_coverage <- .getThisSectorCoverage(
            coverage_data, current_sectors)
        ylim <- c(0, max(this_sector_coverage$coverage))

    } else if (!is.null(coverage_rectangle) & !is.null(custom_ylim)) {
        ylim <- custom_ylim

    } else {
        ylim <- c(0, 1)
    }
    circos.track(
        ylim=ylim, track.height=track_height, bg.border=NA, bg.col=NA,
        track.margin=c(0, 0.01))
}

## helper function to plot an appropriate x axis
.plotXaxis <- function(
    xaxis = TRUE, xaxis_label_size=0.35, xaxis_colour="#747577",
    xaxis_spacing=10, flipped_sector=FALSE, sectors=get.all.sector.index(),
    track_index=get.cell.meta.data("track.index"))
{
    if (xaxis == FALSE){ # if not plotting an xaxis stop here
        return(xaxis)
    }

    ## if we're using the "start_end" xaxis_spacing, we can just make the
    ## space between each tick really massive (360) degrees so that we only
    ## get the first and last ticks (which always plot irrespective of the
    ## sector's size)
    if (xaxis_spacing == "start_end") {
        spacing_angle <- 360
    } else {
        spacing_angle <- xaxis_spacing
    }

    units <- c("bp" = 0, "kb" = 1, "mb" = 2, "gb" = 3)
    ## go sector by sector to account for different scales
    for (i in seq_along(sectors)) {
        ## heavily based on circos.axis and .default.major.by functions
        cell_xlim <- get.cell.meta.data("xlim", sectors[i])
        cell_start_degree <- get.cell.meta.data(
            "cell.start.degree", sectors[i], track_index)
        ## minus for normal (clockwise) sectors, because angles increase as you
        ## go anticlockwise around the unit circle (opposite for reversed
        ## sectors). The ticks are 'spacing_angle' degrees apart
        if (sectors[i] %in% flipped_sector) {
            first_tick_degree <- cell_start_degree + spacing_angle
        } else {
            first_tick_degree <- cell_start_degree - spacing_angle
        }

        ## idk why this works but it does- otherwise the x values come out of
        ## the reverse.circlize() function are really off
        if (cell_start_degree <= circos.par("start.degree")) {
            first_tick_degree <- 360 + first_tick_degree
        }

        ## find the x position of the first tick
        first_tick_x <- reverse.circlize(
            x=first_tick_degree,
            y=get.cell.meta.data(
                "cell.bottom.radius", sectors[i], track_index),
            sector.index=sectors[i], track.index=track_index)[1]

        ## find all of the ticks using the first tick as a guide for how far
        ## apart in bp they need to be
        if (cell_xlim[1] < 0) {
            ## if the start of the sector is negative (common for
            ## insertionDiagram) then correct the first tick to begin at 0
            tick_start <- 0
        } else {
            tick_start <- cell_xlim[1]
        }
        increment_by <- abs(cell_xlim[1] - first_tick_x)
        ticks <- seq(tick_start, cell_xlim[2], by=increment_by)
        ticks <- c(ticks, ticks[length(ticks)] + increment_by)

        ## make the last tick at the very end of the sector:
        ticks[length(ticks)] <- cell_xlim[2]

        ## but if the last two ticks are now too close delete one
        last_tick_degree <- circlize(
            ticks[length(ticks)], 0, sectors[i], CELL_META$track.index)[1]
        second_last_tick_degree <- circlize(
            ticks[length(ticks) - 1], 0, sectors[i], CELL_META$track.index)[1]

        if (abs(last_tick_degree - second_last_tick_degree) < 6) {
            ticks <- ticks[-(length(ticks) - 1)]
        }

        ## need to flip axes for flipped sectors
        if (sectors[i] %in% flipped_sector) {
            ticks <- .reverseXaxis(ticks)
        }

        original_ticks <- round(ticks, 2)
        ## if the number for ticks is more than 4 digits, divide it until its
        ## 4 or less so we dont have massive numbers crowding the whole diagram
        for (j in seq_along(ticks)) {
            names(ticks)[j] <- 0
            while (nchar(ticks[j]) >= 4 & ticks[j] > 1000) {
                ticks[j] <- ticks[j] / 1000
                ## keep track of how many times we divded it so we know what
                ## units to put on
                names(ticks)[j] <- as.numeric(names(ticks)[j]) + 1
            }
            for (k in seq_along(units)) {
                if (names(ticks)[j] == units[k]) {
                    names(ticks)[j] <- names(units)[k]
                    break
                }
            }
        }

        ticks <- signif(ticks, 4)
        labels <- paste0(ticks, names(ticks))

        ## plot the axis
        circos.axis(
            h="top", labels.col=xaxis_colour, labels.facing="clockwise",
            labels.cex=xaxis_label_size, major.at=original_ticks,
            labels=labels, sector.index=sectors[i])

    }
}
