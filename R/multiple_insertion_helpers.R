## helpers for the multipleInsertionDiagram function

#' @importFrom grid bezierGrob
#' @importFrom grid bezierPoints
#' @importFrom grid convertX
#' @importFrom grid convertY

## the plotting area for the MID is made up of 9 slots:
## 1 2 3
## 4 5 6
## 7 8 9
## with slot 5 always being the central genome circle
## we divide the plotting region into 4 'areas' (quarters) each of which
## contains 3 slots (1 corner slot and 2 others). For example:
## area 1 = slots 2, 1, 4
## area 2 = slots 4, 7, 8
## area 3 = slots 8, 9, 6
## area 4 = slots 5, 3, 2
## generally, it is preferred to place plots into the corner slot (2nd slot)

## the .equivSlotFinder function finds the 'equivalent slot' in either the
## next area or the previous area (for example, equivalent slot for slot 7 is
## slot 9 in the next area or slot 1 in the previous area, because they're both
## corner slots). This is important when the allocation of plots/events to
## slots is dependent on the adjacent areas (because the areas have overlap)
.equivSlotFinder <- function(current_slots, which_side){
    ## this is the order of the slots as given above
    all_slots <- c(2, 1, 4, 7, 8, 9, 6, 3)
    new_slot_indexes <- vector()
    for (i in seq_along(current_slots)){
        ## for each slot we need to find an equivalent for
        current_slot_index <- grep(current_slots[i], all_slots)
        if (which_side == "previous") {
            new_slot_index <- current_slot_index - 2
            if (new_slot_index < 1) {
                ## if we've gone negative loop back to end
                new_slot_index <- length(all_slots) + new_slot_index
            }
        } else {
            new_slot_index <- current_slot_index + 2
            if (new_slot_index > length(all_slots)) {
                ## if we've gone off the end come back to the start
                new_slot_index <- new_slot_index - length(all_slots)
            }
        }
        new_slot_indexes[i] <- new_slot_index
    }
    return(all_slots[new_slot_indexes])
}

## .firstOrLastTwo takes the plot/event names you give it and determines
## whether they are in the first two slots of their area (e.g. 2 & 1 or 4 & 7)
## or the last two slots (e.g. 2 & 3 or 7 & 8)
.firstOrLastTwo <- function(event_names, assigned_slots){
    all_slots <- c(2, 1, 4, 7, 8, 9, 6, 3)
    first_event <- event_names[1]

    ## all of the corner slots in the all_slots vector have an even index. So
    ## we can find whether a plot/event is in the corner or not by seeing if it
    ## has an even index:
    first_event_slot <- assigned_slots[first_event]
    first_event_index <- grep(first_event_slot, all_slots)
    if (first_event_index %% 2 == 0){
        ## if the first event is in the corner, the second event must be after
        ## it so these events are in the last two slots
        return("last_two")
    } else {
        ## otherwise the second event must be in the corner (because we always
        ## allocate one of the plots/events to a corner) and so the first event
        ## is in the first slot (so the first two slots are occupied)
        return("first_two")
    }
}

## the .chooseTwoSlots function is used to assign two plots/events in the
## same area into the 3 possible slots based on their position in the genome
.chooseTwoSlots <- function(
    event_angles, area_slots, corner_slot_angle, assigned_slots){
    ## the plot that is closest to the middle of the corner slot goes into the
    ## corner and the other plot goes either before or after it depending on
    ## whether it is before or after
    if (abs(event_angles[1] - corner_slot_angle) <
        abs(event_angles[2] - corner_slot_angle)){
        assigned_slots[names(event_angles)[1]] <- area_slots[2]
        assigned_slots[names(event_angles)[2]] <- area_slots[3]
    } else {
        assigned_slots[names(event_angles)[2]] <- area_slots[2]
        assigned_slots[names(event_angles)[1]] <- area_slots[1]
    }
    return(assigned_slots)
}


## .slotAssigner does what you'd expect. Given the contents of each of the four
## areas, it assigns all of the plots/events into a slot.
.slotAssigner <- function(area_one, area_two, area_three, area_four) {
    ## these basic pieces of information are always the same. They are:
    ## all_area_to_assign: the events we're assigning
    all_area_to_assign <- list(
        area_one, area_two, area_three, area_four)

    ## all_area_neighbours: the neighbouring areas for each area
    all_area_neighbours <- list(
        list(area_four, area_two), list(area_one, area_three),
        list(area_two, area_four), list(area_three, area_one))

    ## all_area_slots: the slots contained within each area
    all_area_slots <- list(
        c(2, 1, 4), c(4, 7, 8), c(8, 9, 6), c(6, 3, 2))

    ## all_corner_slot_angles: the angles for the middle of the corner slot
    ## in each of the four areas
    all_corner_slot_angles <- c(135, 225, 315, 45)

    ## all_area_opposites: the non-neighbouring (opposite) area for each of
    ## these four areas
    all_area_opposites <- list(
        area_three, area_four, area_one, area_two)

    ## assign slots one area at a time
    assigned_slots <- vector()
    for (i in seq(1,4)) {
        ## set up the key information for this area so we don't have to have
        ## so many [[i]] everywhere
        area_to_assign <- all_area_to_assign[[i]]
        area_neighbours <- all_area_neighbours[[i]]
        area_slots <- all_area_slots[[i]]
        corner_slot_angle <- all_corner_slot_angles[[i]]
        area_opposite <- all_area_opposites[[i]]

        ## the way we assign slots depends on how many plots/events are in that
        ## area
        if (length(area_to_assign) == 1) {
            ## if there's only one plot/event, just chuck it in the corner as
            ## this makes the best use of the space
            assigned_slots[names(area_to_assign)[1]] <- area_slots[2]
        } else if (length(area_to_assign) == 3) {
            ## three plots/events are the maximum we can fit in one area so
            ## we don't really have a choice about how we assign them. Just put
            ## the one that comes first in the first slot and so on
            assigned_slots[names(area_to_assign)[1]] <- area_slots[1]
            assigned_slots[names(area_to_assign)[2]] <- area_slots[2]
            assigned_slots[names(area_to_assign)[3]] <- area_slots[3]
        } else if (length(area_to_assign) == 2){
            ## having two plots/events in one area is a little more tricky. In
            ## this case, where we can/can't allocate our plots depends on the
            ## plots/events in neighbouring areas.
            neighbour_events <- c(
                length(area_neighbours[[1]]), length(area_neighbours[[2]]))
            if (neighbour_events[1] == 3 & neighbour_events[2] == 3) {
                ## if both neighbours have three events then 2/3 of the slots
                ## in this area are already full due to the overlap. Thus we
                ## can't plot two events here.
                stop("Sorry, there are too many plots/events close together.")
            } else if (neighbour_events[1] == 3 | neighbour_events[2] == 3) {
                ## again, because of the overlap if we have a single neighbour
                ## with three plots/events then 1/3 of our slots are already
                ## full. Thus, we don't have a choice about where to put the
                ## plots in this area.
                if (neighbour_events[1] == 3) {
                    ## the exact two slots we can use depend on whether its
                    ## the previous or next neighbour who has 3 plots
                    assigned_slots[names(area_to_assign)[1]] <- area_slots[2]
                    assigned_slots[names(area_to_assign)[2]] <- area_slots[3]
                } else if (neighbour_events[2] == 3){
                    assigned_slots[names(area_to_assign)[1]] <- area_slots[1]
                    assigned_slots[names(area_to_assign)[2]] <- area_slots[2]
            }
            } else {
                ## if neither of our neighbours have three plots, then we can
                ## make a choice about which 2 of the 3 slots we want to use
                ## (keeping in mind we'd like to use the second (corner) slot)
                if (neighbour_events[1] <= 1 & neighbour_events[2] <= 1) {
                    ## if both neighbours have 1 or 0 events in them we can
                    ## use whichever slots we like as the overlap won't be a
                    ## problem (neighbours with 0 events won't use any slots,
                    ## neighbours with 1 event will always use the corner
                    ## slot). choose_two_slots chooses the best two slots to
                    ## suit these particular events
                    assigned_slots <- .chooseTwoSlots(
                        area_to_assign, area_slots, corner_slot_angle,
                        assigned_slots)
                } else {
                    ## the only remaining case is where one or both of the
                    ## neighbours have 2 events in them. In this case, we need
                    ## to consider whether the other (non-neighbouring/
                    ## opposite) area has three events, as this will influence
                    ## where the neighbouring areas can place their plots
                    if (length(area_opposite) == 3){
                        if (neighbour_events[1] == 2 &
                            neighbour_events[2] == 2) {
                            ## the figure can accept a maximum of 8 plots total
                            ## and we already have 7 (so can't plot 2 more)
                            stop(
                                "Sorry, multipleInsertionDiagram can only
                                accept a maximum of 8 plots/events")
                        } else {
                            ## otherwise, (only one neighbour has 2 events) we
                            ## need to know which one. This tells us which of
                            ## the slots in the current area we can use
                            if (neighbour_events[1] == 2) {
                                ## if its the area before, the first slot is
                                ## full already so use the 2nd and 3rd slots
                                assigned_slots[names(area_to_assign)[1]] <-
                                    area_slots[2]
                                assigned_slots[names(area_to_assign)[2]] <-
                                    area_slots[3]
                            } else {
                                ## if its the area after, the last slot is
                                ## full already so use the 1st and 2nd slots
                                assigned_slots[names(area_to_assign)[1]] <-
                                    area_slots[1]
                                assigned_slots[names(area_to_assign)[2]] <-
                                    area_slots[2]
                            }
                        }
                    } else {
                        ## otherwise (opposite area has < 3 plots/events and
                        ## one or both neighbours has 2) there are two possible
                        ## ways we can assign our slots.
                        if (neighbour_events[1] == 2 &
                            neighbour_events[2] == 2){
                            ## when both neighbours have 2 events, we'll start
                            ## assigning slots from the area before the current
                            ## one, then the current area, and finally the next
                            ## area
                            previous_neighbour <- area_neighbours[[1]]
                            next_neighbour <- area_neighbours[[2]]

                            ## assign previous_neighbour
                            assigned_slots <- .chooseTwoSlots(
                                previous_neighbour,
                                area_slots = .equivSlotFinder(
                                    area_slots, "previous"),
                                corner_slot_angle = corner_slot_angle - 90,
                                assigned_slots)

                            ## the current area depends on whether the previous
                            ## area uses the first two or last two slots
                            position <- .firstOrLastTwo(
                                names(previous_neighbour), assigned_slots)
                            if (position == "last_two") {
                                ## if the previous area uses the last two
                                ## slots, we must also use the last two slots
                                ## for this area and the area after
                                assigned_slots[names(area_to_assign)[1]] <-
                                    area_slots[2]
                                assigned_slots[names(area_to_assign)[2]] <-
                                    area_slots[3]
                                assigned_slots[names(next_neighbour)[1]] <-
                                    .equivSlotFinder(area_slots[2], "next")
                                assigned_slots[names(next_neighbour)[2]] <-
                                    .equivSlotFinder(area_slots[3], "next")
                            } else {
                                ## otherwise can choose either of the two slots
                                ## for this current area
                                assigned_slots <- .chooseTwoSlots(
                                    area_to_assign, area_slots,
                                    corner_slot_angle, assigned_slots)

                                position <- .firstOrLastTwo(
                                    names(area_to_assign), assigned_slots)

                                ## finally allocate the next area
                                if (position == "last_two") {
                                    ## again, if we use the last two slots on
                                    ## this area, we must use the last two
                                    ## slots on the next area
                                    assigned_slots[names(next_neighbour)[1]] <-
                                        .equivSlotFinder(area_slots[2], "next")
                                    assigned_slots[names(next_neighbour)[2]] <-
                                        .equivSlotFinder(area_slots[3], "next")
                                } else {
                                    ## otherwise choose slots for the next area
                                    assigned_slots <- .chooseTwoSlots(
                                        next_neighbour,
                                        area_slots = .equivSlotFinder(
                                            area_slots, "next"),
                                        corner_slot_angle =
                                            (corner_slot_angle + 90),
                                        assigned_slots)
                                }
                            }
                        } else if (neighbour_events[1] == 2){
                            ## as above, we'll assign the neighbouring areas
                            ## first then do the current one based on that
                            previous_neighbour <- area_neighbours[[1]]
                            assigned_slots <- .chooseTwoSlots(
                                previous_neighbour,
                                area_slots = .equivSlotFinder(
                                    area_slots, "previous"),
                                corner_slot_angle = corner_slot_angle - 90,
                                assigned_slots)

                            position <- .firstOrLastTwo(
                                names(previous_neighbour), assigned_slots)
                            if (position == "last_two") {
                                assigned_slots[names(area_to_assign)[1]] <-
                                    area_slots[2]
                                assigned_slots[names(area_to_assign)[2]] <-
                                    area_slots[3]
                            } else {
                                assigned_slots <- .chooseTwoSlots(
                                    area_to_assign, area_slots,
                                    corner_slot_angle, assigned_slots)
                            }

                        } else {
                            ## exactly the same as above except in this case
                            ## it's the next neighbour not the previous
                            ## neighbour with 2 events
                            next_neighbour <- area_neighbours[[2]]
                            assigned_slots <- .chooseTwoSlots(
                                next_neighbour,
                                area_slots = .equivSlotFinder(
                                    area_slots, "next"),
                                corner_slot_angle = corner_slot_angle + 90,
                                assigned_slots)

                            position <- .firstOrLastTwo(
                                names(next_neighbour), assigned_slots)

                            if (position == "first_two") {
                                ## if the next area uses the first two slots
                                ## then we must also because the first slot of
                                ## the next area is the last slot of this area
                                assigned_slots[names(area_to_assign)[1]] <-
                                    area_slots[1]
                                assigned_slots[names(area_to_assign)[2]] <-
                                    area_slots[2]
                            } else {
                                assigned_slots <- .chooseTwoSlots(
                                    area_to_assign, area_slots,
                                    corner_slot_angle, assigned_slots)
                            }
                        }
                    }
                }
            }
        } else if (length(area_to_assign) > 3){
            ## currently there is no support for more than 3 events in an area
            ## because it would just be too crowded.
            stop("Sorry, there are too many plots/events close together.")
        }
    }
    return(assigned_slots)
}

## make_link_points generates the x and y coordinates for the points that make
## up the link (which is just a polygon)
.makeLinkPoints <- function(slot, x_grid, y_grid, mc_x, mc_y){
    ## a link is made up of three parts: 2 bezier curves and then the points
    ## that run along the event circle. We already have the points along the
    ## event circle as we generated those during the process of plotting the
    ## events themselves, so let's make the bezier curves
    ## the control points for the bezier curve depend on the slot
    if (slot == 1) {
        control_x <- c(0.28, 0.36)
        control_y <- c(0.75, 0.65)
    } else if (slot == 2) {
        control_x <- c(0.50, 0.50)
        control_y <- c(0.66, 0.665)
    } else if (slot == 3) {
        control_x <- c(0.72, 0.62)
        control_y <- c(0.75, 0.65)
    } else if (slot == 4) {
        control_x <- c(0.30, 0.36)
        control_y <- c(0.50, 0.50)
    } else if (slot == 6) {
        control_x <- c(0.65, 0.70)
        control_y <- c(0.50, 0.50)
    } else if (slot == 7) {
        control_x <- c(0.28, 0.35)
        control_y <- c(0.28, 0.35)
    } else if (slot == 8) {
        control_x <- c(0.50, 0.5)
        control_y <- c(0.33, 0.33)
    } else { # slot 9
        control_x <- c(0.72, 0.62)
        control_y <- c(0.28, 0.35)
    }

    ## create a bezier grob and then extract the points from it
    r_bezier <- grid::bezierGrob(
        x = c(x_grid[1], control_x, mc_x), y = c(y_grid[1], control_y, mc_y))
    l_bezier <- grid::bezierGrob(
        x = c(x_grid[length(x_grid)], control_x, mc_x),
        y = c(y_grid[length(y_grid)], control_y, mc_y))
    r_bezier_points <- grid::bezierPoints(r_bezier)
    l_bezier_points <- grid::bezierPoints(l_bezier)

    ## because each bezier curve ends at (mc_x, mc_y), when we join them we'll
    ## have the same point twice. Fix by removing it from one of the curves
    l_bezier_points_no_mc_x <- utils::head(l_bezier_points$x, -1)
    l_bezier_points_no_mc_y <- utils::head(l_bezier_points$y, -1)

    ## bezierPoints outputs in inches, we want npc:
    l_bezier_points_no_mc_x <- grid::convertX(l_bezier_points_no_mc_x, "npc")
    l_bezier_points_no_mc_y <- grid::convertY(l_bezier_points_no_mc_y, "npc")
    r_bezier_points$x <- grid::convertX(r_bezier_points$x, "npc")
    r_bezier_points$y <- grid::convertY(r_bezier_points$y, "npc")

    ## now put all the points together in order
    link_x_points <- c(
        r_bezier_points$x, # the first bezier
        rev(l_bezier_points_no_mc_x), # second bezier
        rev(x_grid)) # the points along the circle
    link_y_points <- c( # same as above
        r_bezier_points$y, rev(l_bezier_points_no_mc_y), rev(y_grid))
    return(list(link_x_points, link_y_points))
}
