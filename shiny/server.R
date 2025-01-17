function(input, output, session)
{

  values <- reactiveValues(
    readCounts = emptyReadCounts,
    readCountSummary = "",
    snvs = emptySnvs,
    snvSummary = "",
    selectedSubstitution = NULL,
    selectedId = NULL
  )


  # read counts tab pane

  observe({
    file <- input$readCountFile
    if (is.null(file))
    {
      values$readCounts <- emptyReadCounts
      values$readCountSummary <- ""
      return(NULL)
    }

    readCounts <- read.delim(file$datapath, colClasses = "character", stringsAsFactors = FALSE, check.names = FALSE)

    missingColumns <- setdiff(readCountColumns, colnames(readCounts))
    if (length(missingColumns) > 0)
    {
      values$readCounts <- emptyReadCounts
      values$readCountSummary <- str_c("Error: missing columns: ", str_c(missingColumns, collapse = ", "))
      return(NULL)
    }

    readCounts <- readCounts %>%
      select(one_of(readCountColumns)) %>%
      mutate(Position = as.integer(Position)) %>%
      mutate(`Depth unfiltered` = as.integer(`Depth unfiltered`)) %>%
      mutate(Depth = as.integer(Depth)) %>%
      mutate(`A count` = as.integer(`A count`)) %>%
      mutate(`C count` = as.integer(`C count`)) %>%
      mutate(`G count` = as.integer(`G count`)) %>%
      mutate(`T count` = as.integer(`T count`))

    librarySampleMapping <- readCounts %>%
      select(ID, Sample) %>%
      distinct
    ids <- librarySampleMapping %>%
      select(ID) %>%
      distinct
    if (nrow(ids) != nrow(librarySampleMapping))
    {
      values$readCounts <- emptyReadCounts
      values$readCountSummary <- "Error: multiple samples specified for a single library"
      return(NULL)
    }

    values$readCounts <- readCounts
    values$readCountSummary <- str_c(
      "<b>",
      readCounts %>% nrow(),
      "</b> rows for <b>",
      readCounts %>% select(ID) %>% distinct() %>% nrow(),
      "</b> libraries (<b>",
      readCounts %>% select(Sample) %>% distinct() %>% nrow(),
      "</b> samples) at <b>",
      readCounts %>% select(Amplicon, Chromosome, Position) %>% distinct() %>% nrow(),
      "</b> target locations"
    )
  })

  output$readCountSummary <- renderUI({
    HTML(values$readCountSummary)
  })

  output$readCountTable <- DT::renderDataTable(
    datatable(
      values$readCounts,
      rownames = FALSE,
      selection = "single",
      filter = "bottom",
      extensions = "Buttons",
      options = list(
        dom = 'Bfrtip',
        buttons = I('pageLength'),
        pageLength = 10,
        lengthMenu = list(c(10, 25, 50, 100), c('Show 10 rows', 'Show 25 rows', 'Show 50 rows', 'Show 100 rows')),
        searchHighlight = TRUE,
        columnDefs = list(list(className = 'dt-center', targets = 5))
      )
    ),
    server = TRUE
  )

  observe({
    selectedRow <- input$readCountTable_rows_selected
    if (is.null(selectedRow)) return(NULL)

    selectedReadCounts <- slice(values$readCounts, selectedRow)

    values$selectedId <- selectedReadCounts$ID

    referenceBase <- selectedReadCounts %>%
      select(`Reference base`) %>%
      unlist(use.names = FALSE)
    alternateAllele <- setdiff(bases, referenceBase)[1]

    values$selectedSubstitution <- selectedReadCounts %>%
      select(Amplicon, Chromosome, Position, `Reference base`) %>%
      mutate(`Alternate allele` = alternateAllele)
  })


  # SNVs tab pane

  observe({
    file <- input$snvFile
    if (is.null(file))
    {
      values$snvs <- emptySnvs
      values$snvSummary <- ""
      return(NULL)
    }

    snvs <- read.delim(file$datapath, colClasses = "character", stringsAsFactors = FALSE, check.names = FALSE)

    idColumns <- grep("^ID \\d+$", colnames(snvs), value = TRUE)
    numberOfReplicates <- length(idColumns)
    if (numberOfReplicates > 0)
    {
      expectedIdColumns <- str_c("ID ", 1:numberOfReplicates)
      missingColumns <- setdiff(expectedIdColumns, colnames(snvs))
      if (length(missingColumns) > 0)
      {
        values$snvs <- emptySnvs
        values$snvSummary <- "Error: ID columns should be numbered consecutively, i.e. ID 1, ID 2, etc."
        return(NULL)
      }

      filterColumns <- str_c("Filters ", 1:numberOfReplicates)
      otherColumns <- setdiff(colnames(snvs), c(idColumns, filterColumns))

      expandedSnvs <- NULL

      for (i in 1:numberOfReplicates)
      {
        idColumn <- str_c("ID ", i)
        filterColumn <- str_c("Filters ", i)

        if (filterColumn %in% colnames(snvs)) {
          replicateSnvs <- snvs %>%
            select(ID = matches(idColumn), one_of(otherColumns), Filters = matches(filterColumn)) %>%
            filter(Filters != "no call")
        } else {
          replicateSnvs <- snvs %>%
            select(ID = matches(idColumn), one_of(otherColumns))
        }

        replicateSnvs <- replicateSnvs %>%
          filter(ID != "")

        expandedSnvs <- bind_rows(expandedSnvs, replicateSnvs)
      }

      snvs <- expandedSnvs
    }

    if (!"Reference base" %in% colnames(snvs) && "Ref" %in% colnames(snvs))
      snvs <- rename(snvs, `Reference base` = Ref)

    if (!"Alternate allele" %in% colnames(snvs) && "Alt" %in% colnames(snvs))
      snvs <- rename(snvs, `Alternate allele` = Alt)

    if (!"Filters" %in% colnames(snvs))
      snvs <- mutate(snvs, Filters = "")

    if (!"Confidence" %in% colnames(snvs))
      snvs <- mutate(snvs, Confidence = "")

    missingColumns <- setdiff(snvColumns, colnames(snvs))
    if (length(missingColumns) > 0)
    {
      values$snvs <- emptySnvs
      values$snvSummary <- str_c("Error: missing columns: ", str_c(missingColumns, collapse = ", "))
      return(NULL)
    }

    snvs <- snvs %>%
      select(one_of(snvColumns)) %>%
      mutate(Position = as.integer(Position)) %>%
      distinct %>%
      arrange(ID, Chromosome, Position, `Reference base`, `Alternate allele`, Amplicon)

    values$snvs <- snvs
    values$snvSummary <- str_c(
      "<b>",
      nrow(snvs),
      "</b> single nucleotide variants"
    )
  })

  output$snvSummary <- renderUI({
    HTML(values$snvSummary)
  })

  snvTableData <- reactive({
    snvs <- values$snvs
    readCounts <- values$readCounts

    alleleFractions <- snvs %>%
      select(ID, Amplicon, Chromosome, Position, `Reference base`) %>%
      distinct %>%
      inner_join(readCounts, by = c("ID", "Amplicon", "Chromosome", "Position", "Reference base")) %>%
      mutate(
        A = `A count` / Depth,
        C = `C count` / Depth,
        G = `G count` / Depth,
        T = `T count` / Depth
      ) %>%
      select(-`Depth unfiltered`, -Depth) %>%
      gather(`Alternate allele`, `Allele fraction`, A, C, G, T)

    snvs %>%
      left_join(alleleFractions, by = c("ID", "Amplicon", "Chromosome", "Position", "Reference base", "Alternate allele")) %>%
      select(ID, Sample, Amplicon, Chromosome, Position, `Reference base`, `Alternate allele`, `Allele fraction`, `A count`, `C count`, `G count`, `T count`, Filters, Confidence)
  })

  output$snvTable <- DT::renderDataTable(
    datatable(
      snvTableData(),
      rownames = FALSE,
      selection = "single",
      filter = "bottom",
      extensions = "Buttons",
      options = list(
        dom = 'Bfrtip',
        buttons = I('pageLength'),
        pageLength = 10,
        lengthMenu = list(c(10, 25, 50, -1), c('Show 10 rows', 'Show 25 rows', 'Show 50 rows', 'Show all rows')),
        searchHighlight = TRUE,
        columnDefs = list(list(className = 'dt-center', targets = c(4, 5)))
      )
    ) %>%
      formatRound(c("Allele fraction"), digits = 5),
    server = TRUE
  )

  observe({
    selectedRow <- input$snvTable_rows_selected
    if (is.null(selectedRow)) return(NULL)

    selectedSnv <- slice(snvTableData(), selectedRow)

    values$selectedId <- selectedSnv$ID

    values$selectedSubstitution <- selectedSnv %>%
      select(Amplicon, Chromosome, Position, `Reference base`, `Alternate allele`)
  })


  # locations tab pane

  updatingLocationControls <- FALSE

  locationTableData <- reactive({
    values$readCounts %>%
      select(Chromosome, Position, `Reference base`) %>%
      distinct %>%
      arrange(Chromosome, Position, `Reference base`)
  })

  output$locationSummary <- renderUI({
    HTML(
      str_c(
        "Read counts available for <b>",
        nrow(locationTableData()),
        "</b> distinct target locations"
      )
    )
  })

  # location selection table should only get updated once immediately after
  # a new set of read counts is uploaded
  output$locationTable <- DT::renderDataTable({

    data <- locationTableData()

    selection <- "single"
    isolate(selectedSubstitution <- values$selectedSubstitution)
    if (!is.null(selectedSubstitution))
    {
      selectedRow <- which(data$Chromosome == selectedSubstitution$Chromosome & data$Position == selectedSubstitution$Position)
      if (length(selectedRow) == 1)
        selection <- list(mode = "single", selected = selectedRow)
    }

    datatable(
      data,
      rownames = FALSE,
      selection = selection,
      extensions = 'Scroller',
      options = list(
        dom = 'ft',
        scroller = TRUE,
        deferRender = TRUE,
        scrollY = 300,
        searchHighlight = TRUE,
        columnDefs = list(list(className = 'dt-center', targets = 2))
      )
    )
  },
    server = FALSE
  )

  locationTableProxy <- dataTableProxy('locationTable')

  selectedLocation <- reactive({
    selectedRow <- input$locationTable_rows_selected
    # if (is.null(selectedRow))
    #   NULL
    # else
    #   slice(locationTableData(), selectedRow)
    if (is.null(selectedRow))
    {
      NULL
       message("Location table: no row selected")
    }
    else
    {
      result <- slice(locationTableData(), selectedRow)
      message("Location table: ", unlist(result))
      result
    }
  })

  selectedLocationReferenceBase <- reactive({
    selectedLocation <- selectedLocation()
    if (is.null(selectedLocation))
      NULL
    else
      selectedLocation %>%
      select(`Reference base`) %>%
      unlist(use.names = FALSE)
  })

  locationAlternateAllele <- reactive({
    referenceBase <- selectedLocationReferenceBase()
    selectedAllele <- input$locationAlternateAllele
    if (is.null(referenceBase) || is.null(selectedAllele) || referenceBase == selectedAllele)
      NULL
    else
      selectedAllele
  })

  # updates location substitution selection controls in response to
  # update of values$selectedSubstitution
  observe({
    selectedSubstitution <- values$selectedSubstitution
    if (is.null(selectedSubstitution)) return(NULL)
    message("Updating location selection controls for ", unlist(selectedSubstitution))

    data <- locationTableData()

    selectedRow <- which(data$Chromosome == selectedSubstitution$Chromosome & data$Position == selectedSubstitution$Position)
    if (length(selectedRow) == 0) return(NULL)

    updatingLocationControls <- TRUE

    isolate(previouslySelectedRow <- input$locationTable_rows_selected)

    if (is.null(previouslySelectedRow) || previouslySelectedRow != selectedRow)
      selectRows(locationTableProxy, selectedRow)

    amplicons <- values$readCounts %>%
      filter(Chromosome == selectedSubstitution$Chromosome & Position == selectedSubstitution$Position) %>%
      select(Amplicon) %>%
      distinct %>%
      arrange(Amplicon) %>%
      unlist(use.names = FALSE)

    selectedAmplicon <- selectedSubstitution$Amplicon

    ampliconLabel <- "Amplicon"
    if (length(amplicons) > 1)
      ampliconLabel <- str_c("Amplicons (", length(amplicons), ")")

    updateSelectInput(
      session,
      "locationAmplicon",
      label = ampliconLabel,
      choices = amplicons,
      selected = selectedAmplicon
    )

    referenceBase <- selectedSubstitution$`Reference base`
    alternateAllele <- selectedSubstitution$`Alternate allele`
    alleles <- setdiff(bases, selectedSubstitution$`Reference base`)

    message("Updating locationAlternateAllele radio buttons ", alternateAllele, " [", alleles, "]")
    updateRadioButtons(
      session,
      "locationAlternateAllele",
      choices = alleles,
      selected = alternateAllele,
      inline = TRUE
    )
  })

  # updates amplicon drop-down list in response to change in the selected location
  observe({
    selectedLocation <- selectedLocation()
    if (is.null(selectedLocation)) return(NULL)

    if (updatingLocationControls) return(NULL)

    amplicons <- values$readCounts %>%
      filter(Chromosome == selectedLocation$Chromosome & Position == selectedLocation$Position) %>%
      select(Amplicon) %>%
      distinct %>%
      arrange(Amplicon) %>%
      unlist(use.names = FALSE)

    ampliconLabel <- "Amplicon"
    if (length(amplicons) > 1)
      ampliconLabel <- str_c("Amplicons (", length(amplicons), ")")

    selectedAmplicon <- NULL
    if (length(amplicons) > 0)
    {
      isolate(selectedSubstitution <- values$selectedSubstitution)
      if (!is.null(selectedSubstitution))
      {
        amplicon <- selectedSubstitution$Amplicon
        if (!is.null(amplicon) && is.element(amplicon, amplicons))
          selectedAmplicon <- amplicon
      }

      if (is.null(selectedAmplicon))
      {
        isolate(amplicon <- input$locationAmplicon)
        if (!is.null(amplicon) && is.element(amplicon, amplicons))
          selectedAmplicon <- amplicon
      }

      if (is.null(selectedAmplicon))
        selectedAmplicon <- amplicons[1]
    }

    message("Updating amplicon drop down: ", selectedAmplicon, " (", str_c(amplicons, collapse = ", "), ")")
    updateSelectInput(
      session,
      "locationAmplicon",
      label = ampliconLabel,
      choices = amplicons,
      selected = selectedAmplicon
    )
  })

  # updates alternate allele radio buttons in response to change in the selected location
  observe({
    referenceBase <- selectedLocationReferenceBase()

    if (updatingLocationControls) return(NULL)

    alleles <- bases
    if (!is.null(referenceBase)) alleles <- setdiff(alleles, referenceBase)
    message("Updating alternate allele radio buttons")
    message("Reference base: ", referenceBase)

    selectedAllele <- NULL

    isolate(selectedSubstitution <- values$selectedSubstitution)
    if (!is.null(selectedSubstitution))
    {
      allele <- selectedSubstitution$`Alternate allele`
      message("Alternate allele for selected substitution: ", allele)
      if (is.element(allele, alleles))
        selectedAllele <- allele
      message("Selection allele: ", selectedAllele)
    }

    if (is.null(selectedAllele))
    {
      isolate(allele <- input$locationAlternateAllele)
      if (is.element(allele, alleles))
        selectedAllele <- allele
    }
    message("Selection allele: ", selectedAllele)

    updateRadioButtons(session, "locationAlternateAllele", choices = alleles, selected = selectedAllele, inline = TRUE)
  })

  # updates values$selectedSubstitution in response to selections made in the
  # location substitution controls
  observe({
    location <- selectedLocation()
    amplicon <- input$locationAmplicon
    alternateAllele <- input$locationAlternateAllele

    isolate(selectedSubstitution <- values$selectedSubstitution)

    if (updatingLocationControls)
    {
      if (is.null(selectedSubstitution) ||
          is.null(location) ||
          is.null(amplicon) ||
          is.null(alternateAllele) ||
          location$Chromosome != selectedSubstitution$Chromosome ||
          location$Position != selectedSubstitution$Position ||
          amplicon != selectedSubstitution$Amplicon ||
          alternateAllele != selectedSubstitution$AlternateAllele)
      {
        updatingLocationControls <- FALSE
      }
    }
    else
    {
      if (is.null(location) ||
          is.null(amplicon) ||
          is.null(alternateAllele))
        return(NULL)

      amplicons <- values$readCounts %>%
        filter(Chromosome == location$Chromosome & Position == location$Position) %>%
        select(Amplicon) %>%
        distinct %>%
        arrange(Amplicon) %>%
        unlist(use.names = FALSE)
      if (!is.element(amplicon, amplicons)) return(NULL)

      if (alternateAllele == location$`Reference base`) return(NULL)

      substitution <- location %>%
        select(Chromosome, Position, `Reference base`) %>%
        mutate(`Alternate allele` = alternateAllele) %>%
        mutate(Amplicon = amplicon)

      if (is.null(selectedSubstitution) ||
          selectedSubstitution$Chromosome != substitution$Chromosome ||
          selectedSubstitution$Position != substitution$Position ||
          selectedSubstitution$`Reference base` != substitution$`Reference base` ||
          selectedSubstitution$`Alternate allele` != substitution$`Alternate allele` ||
          selectedSubstitution$Amplicon != substitution$Amplicon)
      {
        message("Updating values$selectedSubstitution (1) ", unlist(selectedSubstitution), " -> ", unlist(substitution))
        values$selectedSubstitution <- substitution
      }
    }
  })

  locationAlternateAlleleData <- reactive({

    selectedLocation <- selectedLocation()
    alternateAllele <- locationAlternateAllele()

    result <- emptyReadCounts

    if (!is.null(selectedLocation) && !is.null(alternateAllele))
    {
      result <- values$readCounts %>%
        inner_join(selectedLocation, by = c("Chromosome", "Position", "Reference base"))
    }

    result <- result %>%
      mutate(
        A = `A count` / Depth,
        C = `C count` / Depth,
        G = `G count` / Depth,
        T = `T count` / Depth
      ) %>%
      gather(`Alternate allele`, `Allele fraction`, A, C, G, T)

    if (!is.null(alternateAllele))
    {
      result <- result %>%
        filter(`Alternate allele` == alternateAllele)
    }

    result <- result %>%
      filter(Depth >= input$locationMinimumDepth)

    snvs <- values$snvs
    result <- result %>%
      left_join(
        mutate(
          snvs,
          Called = Filters != "no call",
          Filtered = ifelse(!Called, NA, !(Filters %in% c("pass", "PASS", "")))
        ),
        by = c("ID", "Amplicon", "Chromosome", "Position", "Reference base", "Alternate allele")
      ) %>%
      mutate(Called = ifelse(is.na(Called), ifelse(nrow(snvs) == 0, NA, FALSE), Called)) %>%
      mutate(Series = ifelse(is.na(Called), "Call status unknown", ifelse(Called, ifelse(Filtered, "Filtered", "Called"), "Not called")))

    result
  })

  locationAmpliconAlternateAlleleData <- reactive({
    result <- locationAlternateAlleleData()
    selectedAmplicon <- input$locationAmplicon
    if (is.null(selectedAmplicon))
      slice(result, 0)
    else
      filter(result, Amplicon == selectedAmplicon)
  })

  output$locationSubstitutionSummary <- renderUI({

    summary <- "No location selected"

    selectedLocation <- selectedLocation()
    if (!is.null(selectedLocation))
    {
      summary <- selectedLocation %>%
        select(Chromosome, Position, `Reference base`) %>%
        unlist(use.names = FALSE)

      summary <- str_c(summary, collapse = " ")

      alternateAllele <- locationAlternateAllele()
      if (!is.null(alternateAllele))
        summary <- str_c(summary, alternateAllele, sep = ">")
      summary <- str_c("<b>", summary, "</b>")

      selectedAmplicon <- input$locationAmplicon
      if (!is.null(selectedAmplicon))
        summary <- str_c(summary, selectedAmplicon, sep = "&nbsp;&nbsp;")

      data <- locationAmpliconAlternateAlleleData()
      if (nrow(data) > 0)
        summary <- str_c(summary, str_c(nrow(data), " libraries"), sep = "&nbsp;&nbsp;")
    }

    HTML(summary)
  })

  locationFittedDistribution <- reactive({
    data <- locationAmpliconAlternateAlleleData()
    alleleFractions <- data$`Allele fraction`
    result <- fitDistribution(
      alleleFractions,
      excludeHighestProportion = input$locationExcludeHighestProportionForFitting,
      maximumAlleleFraction = input$locationMaximumAlleleFractionForFitting,
      includeZeroAlleleFractionValues = input$includeZeroAlleleFractionValues,
      distribution = input$locationDistribution,
      thresholdProbabilities = thresholdProbabilities
    )
    result$alternateAlleleData <- data
    result
  })

  locationMaximumAlleleFraction <- reactive({
    fittedDistribution <- locationFittedDistribution()
    maximumAlleleFraction <- fittedDistribution$maximumAlleleFraction
    if (input$setLocationMaximumAlleleFraction) maximumAlleleFraction <- input$locationMaximumAlleleFraction
    maximumAlleleFraction
  })

  locationThresholdAlleleFraction <- reactive({
    thresholdAlleleFraction <- -1.0
    fittedDistribution <- locationFittedDistribution()
    thresholds <- fittedDistribution$thresholds
    if (!is.null(thresholds))
    {
      thresholdProbability <- input$locationThresholdProbability
      thresholdIndex <- which(thresholdProbabilities == thresholdProbability)
      if (length(thresholdIndex) == 1)
        thresholdAlleleFraction <- thresholds[thresholdIndex]
    }
    thresholdAlleleFraction
  })

  output$locationScatterBoxPlot <- renderHighchart({

    fittedDistribution <- locationFittedDistribution()

    data <- fittedDistribution$alternateAlleleData
    if (!is.null(data))
    {
      snvs <- values$snvs

      selectedId <- values$selectedId
      if (!is.null(selectedId))
      {
        selectedSample <- data %>%
          filter(ID == selectedId) %>%
          select(Sample) %>%
          distinct %>%
          unlist(use.names = FALSE)

        if (length(selectedSample) == 1)
        {
          data <- data %>%
            mutate(Series = ifelse(Sample == selectedSample, "Selected replicate", Series))
        }

        data <- data %>%
          mutate(Series = ifelse(ID == selectedId, "Selected", Series))
      }

      data <- data %>%
        arrange(ID) %>%
        transmute(
          Sample,
          `Allele fraction`,
          Series,
          CallStatus = ifelse(!is.na(Filters) & !(Filters %in% c("pass", "PASS", "")), str_c("<br>Filters: ", Filters),
                       ifelse(Series == "Selected" | Series == "Selected replicate", ifelse(is.na(Called), "<br>Call status unknown", ifelse(Called, "<br>Called", "<br>Not called")), "")),
          tooltip = str_c(
            "ID: ", ID,
            "<br>Sample: ", Sample,
            "<br>AF: ", format(round(`Allele fraction`, digits = 5), scientific = FALSE),
            CallStatus,
            ifelse(Called, str_c("<br>Confidence: ", Confidence), "")
          ),
          ID = str_c(ID, Chromosome, Position, `Reference base`, `Alternate allele`, Amplicon, sep = "\t")
        )
    }

    maximumAlleleFraction <- locationMaximumAlleleFraction()
    thresholdAlleleFraction <- locationThresholdAlleleFraction()

    alleleFractionScatterBoxPlot(
      data,
      series = c("Call status unknown", "Not called", "Filtered", "Called", "Selected", "Selected replicate"),
      colours = c(hex_to_rgba("#7cb5ec", 0.5), hex_to_rgba("#7cb5ec", 0.5), hex_to_rgba("#8085e9", 0.5), "#90ed7d", "#f15c80", hex_to_rgba("#f15c80", 0.55)),
      xlabel = "Library",
      maximumAlleleFraction = maximumAlleleFraction,
      thresholdAlleleFraction = thresholdAlleleFraction,
      thresholdLabel = str_c("p = ", input$locationThresholdProbability, ", AF = ", format(round(thresholdAlleleFraction, digits = 4), scientific = FALSE)),
      clicked = "locationScatterBoxPlot_clicked"
    )
  })

  observe({
    selected <- input$locationScatterBoxPlot_clicked
    if (is.null(selected)) return(NULL)

    selected <- selected %>%
      strsplit("\t") %>%
      unlist

    id <- selected[1]
    chromosome <- selected[2]
    position <- as.integer(selected[3])
    referenceBase <- selected[4]
    alternateAllele <- selected[5]
    amplicon = selected[6]

    isolate(selectedId <- values$selectedId)

    if (is.null(selectedId) || selectedId != id)
      values$selectedId <- id

    substitution <- tibble(
      Amplicon = amplicon,
      Chromosome = chromosome,
      Position = position,
      `Reference base` = referenceBase,
      `Alternate allele` = alternateAllele
    )

    isolate(selectedSubstitution <- values$selectedSubstition)

    if (is.null(selectedSubstitution) ||
        selectedSubstitution$Chromosome != substitution$Chromosome ||
        selectedSubstitution$Position != substitution$Position ||
        selectedSubstitution$`Reference base` != substitution$`Reference base` ||
        selectedSubstitution$`Alternate allele` != substitution$`Alternate allele` ||
        selectedSubstitution$Amplicon != substitution$Amplicon)
    {
      message("Updating values$selectedSubstitution (2) ", unlist(selectedSubstitution), " -> ", unlist(substitution))
      values$selectedSubstitution <- substitution
    }
  })

  observe({
    input$locationCheckSelectedButton
    message(isolate(values$selectedId))
    message(isolate(values$selectedSubstitution) %>% unlist(use.names = FALSE))
  })

  observe({
    input$locationClearSelectedButton
    values$selectedId <- NULL
  })

  locationSelectedSampleData <- reactive({

    data <- locationAlternateAlleleData()

    selectedId <- values$selectedId
    if (is.null(selectedId))
      data <- slice(data, 0)
    else
    {
      selectedSample <- data %>%
        filter(ID == selectedId) %>%
        select(Sample) %>%
        distinct %>%
        unlist(use.names = FALSE)

      if (length(selectedSample) == 1)
        data <- filter(data, Sample == selectedSample)
      else
        data <- slice(data, 0)
    }

    data <- data %>%
      mutate(Called = ifelse(is.na(Called), "", ifelse(Called, "yes", "no"))) %>%
      mutate(Filtered = ifelse(is.na(Filtered), "", ifelse(Filtered, "yes", "no")))

    data
  })

  output$locationSelectedSampleTable <- DT::renderDataTable({

    data <- locationSelectedSampleData() %>%
      select(ID, Sample, Amplicon, Chromosome, Position, `Reference base`, `Alternate allele`, Called, Filtered, `Allele fraction`, `A count`, `C count`, `G count`, `T count`) %>%
      arrange(ID, Sample, Amplicon)

    columnNames <- c("ID", "Sample", "Amplicon", "Chr", "Pos", "Ref", "Alt", "Called", "Filtered", "AF", "A", "C", "G", "T")

    isolate(snvs <- values$snvs)

    if (nrow(snvs) == 0)
    {
      data <- select(data, -Called, -Filtered)
      columnNames <- setdiff(columnNames, c("Called", "Filtered"))
    }

    # TODO
    # selection <- "single"
    # isolate(selectedSubstitution <- values$selectedSubstitution)
    # selectedAmplicon <- NULL
    # if (!is.null(selectedSubstitution)) selectedAmplicon <- selectedSubstitution$Amplicon
    # message("Selected amplicon: ", selectedAmplicon)
    # isolate(selectedId <- values$selectedId)
    # message("Selected ID: ", selectedId)
    # selectedRow <- which(data$ID == selectedId & data$Amplicon == selectedAmplicon)
    # message("Selected rows: ", selectedRow)
    # if (length(selectedRow) == 1)
    #   selection <- list(mode = "single", selected = selectedRow)

    datatable(
      data,
      colnames = columnNames,
      rownames = FALSE,
      selection = "none",
      options = list(
        dom = 'ti',
        pageLength = 1000,
        searchHighlight = TRUE,
        columnDefs = list(list(className = 'dt-center', targets = c(5, 6, 7)))
      )
    ) %>%
      formatRound(c("Allele fraction"), digits = 5)
  },
  server = FALSE
  )

  # Disabled selection in table because of problems caused by having a selected
  # row when selecting a replicate sample in the library details table
  # Note that it is possible to select a replicate by clicking on the point in
  # the scatter plot and it is possible to select an overlapping amplicon through
  # the drop-down menu so being able to select a specific replicate and amplicon
  # in this table, while convenient, is not absolutely necessary
  # observe({
  #   selectedRow <- input$locationSelectedSampleTable_rows_selected
  #   if (is.null(selectedRow)) return(NULL)
  #
  #   data <- locationSelectedSampleData()
  #   if (nrow(data) == 0) return(NULL)
  #
  #   selected <- slice(data, selectedRow)
  #   if (nrow(selected) != 1) return(NULL)
  #
  #   id <- selected$ID
  #   isolate(selectedId <- values$selectedId)
  #   if (is.null(selectedId) || selectedId != id) values$selectedId <- id
  #
  #   substitution <- selected %>%
  #     select(Amplicon, Chromosome, Position, `Reference base`, `Alternate allele`)
  #   isolate(selectedSubstitution <- values$selectedSubstitution)
  #   if (is.null(selectedSubstitution) ||
  #       selectedSubstitution$Chromosome != substitution$Chromosome ||
  #       selectedSubstitution$Position != substitution$Position ||
  #       selectedSubstitution$`Reference base` != substitution$`Reference base` ||
  #       selectedSubstitution$`Alternate allele` != substitution$`Alternate allele` ||
  #       selectedSubstitution$Amplicon != substitution$Amplicon)
  #   {
  #     message("Updating substitution ", unlist(substitution))
  #     values$selectedSubstitution <- substitution
  #   }
  # })

  output$locationDensityPlot <- renderHighchart({

    data <- NULL
    series <- NULL

    fittedDistribution <- locationFittedDistribution()

    density <- fittedDistribution$density
    if (!is.null(density))
    {
      densitySeries <- str_c("Unfiltered (", length(fittedDistribution$alleleFractions), " libraries)")
      data <- mutate(density, Series = densitySeries)
      series <- densitySeries

      filteredDensity <- fittedDistribution$filteredDensity
      if (!is.null(filteredDensity))
      {
        filteredDensitySeries <- densitySeries <- str_c("Filtered (", length(fittedDistribution$filteredAlleleFractions), " libraries)")
        data <- data %>%
          bind_rows(filteredDensity %>%
          mutate(Series = filteredDensitySeries))
        series <- c(series, filteredDensitySeries)
      }

      fitted <- fittedDistribution$fitted
      if (!is.null(fitted))
      {
        fittedDistributionSeries <- "Fitted distribution"
        data <- data %>%
          bind_rows(fitted %>%
          mutate(Series = fittedDistributionSeries))
        series <- c(series, fittedDistributionSeries)
      }
    }

    maximumAlleleFraction <- locationMaximumAlleleFraction()
    thresholdAlleleFraction <- locationThresholdAlleleFraction()

    alleleFractionDensityPlot(
      data,
      series = series,
      type = c("areaspline", "areaspline", "spline"),
      colours = c("#8085e9", "#7cb5ec", "#434348"),
      maximumAlleleFraction = maximumAlleleFraction,
      thresholdAlleleFraction = thresholdAlleleFraction,
      thresholdLabel = str_c("p = ", input$locationThresholdProbability, ", AF = ", format(round(thresholdAlleleFraction, digits = 4), scientific = FALSE)),
      decimalPlaces = 4
    )
  })

  output$locationCullenFreyGraph <- renderPlot({

    data <- locationAmpliconAlternateAlleleData()

    alleleFractions <- data$`Allele fraction`
    nonZeroAlleleFractions <- alleleFractions[alleleFractions != 0]

    if (length(nonZeroAlleleFractions) > 0)
      descdist(nonZeroAlleleFractions)
    else
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  })


  # libraries tab pane

  libraryTableData <- reactive({
    values$readCounts %>%
      select(ID, Sample) %>%
      distinct %>%
      arrange(ID)
  })

  output$librarySummary <- renderUI({
    HTML(
      str_c(
        "Read counts available for <b>",
        nrow(libraryTableData()),
        "</b> libraries"
      )
    )
  })

  # library selection table should only get updated once immediately after
  # a new set of read counts is uploaded
  output$libraryTable <- DT::renderDataTable({

    data <- libraryTableData()

    selection <- "single"
    isolate(selectedId <- values$selectedId)
    if (!is.null(selectedId))
    {
      selectedRow <- which(data$ID == selectedId)
      if (length(selectedRow) == 1)
        selection <- list(mode = "single", selected = selectedRow)
    }

    datatable(
      data,
      rownames = FALSE,
      selection = selection,
      extensions = 'Scroller',
      options = list(
        dom = 'ft',
        scroller = TRUE,
        deferRender = TRUE,
        scrollY = 300,
        searchHighlight = TRUE
      )
    )
  },
  server = FALSE
  )

  libraryTableProxy <- dataTableProxy('libraryTable')

  selectedLibrary <- reactive({
    selectedRow <- input$libraryTable_rows_selected
    if (is.null(selectedRow))
      NULL
    else
      slice(libraryTableData(), selectedRow)
  })

  # updates values$selectedId in response to selection made in the library table
  observe({
    library <- selectedLibrary()
    if (!is.null(library)) values$selectedId <- library$ID
  })

  # update selected row in library table following change to values$selectedId
  observe({
    selectedId <- values$selectedId
    if (is.null(selectedId)) return(NULL)

    data <- libraryTableData()

    selectedRow <- which(data$ID == selectedId)
    if (length(selectedRow) == 0) return(NULL)

    isolate(previouslySelectedRow <- input$libraryTable_rows_selected)

    if (is.null(previouslySelectedRow) || previouslySelectedRow != selectedRow)
    {
      selectRows(libraryTableProxy, selectedRow)
    }
  })

  # update substitution drop-down selection following change to values$selectedSubstitution
  observe({
    selectedSubstitution <- values$selectedSubstitution
    if (is.null(selectedSubstitution)) return(NULL)

    referenceBase <- selectedSubstitution$`Reference base`
    alternateAllele <- selectedSubstitution$`Alternate allele`
    selectedSubstitution <- str_c(referenceBase, alternateAllele, sep = ">")

    isolate(previouslySelectedSubstitution <- input$librarySubstitution)
    message("Substitution ", previouslySelectedSubstitution, " -> ", selectedSubstitution)

    if (is.null(previouslySelectedSubstitution) || previouslySelectedSubstitution != selectedSubstitution)
    {
      message("Updating librarySubstitution selectInput ", selectedSubstitution)
      updateSelectInput(
        session,
        "librarySubstitution",
        label = "Sustitution",
        choices = substitutions,
        selected = selectedSubstitution
      )
    }
  })

  sampleSubstitutionData <- reactive({

    result <- emptyReadCounts

    library <- selectedLibrary()
    if (!is.null(library))
      result <- filter(values$readCounts, Sample == library$Sample)

    result <- result %>%
      mutate(
        A = `A count` / Depth,
        C = `C count` / Depth,
        G = `G count` / Depth,
        T = `T count` / Depth
      ) %>%
      gather(`Alternate allele`, `Allele fraction`, A, C, G, T)

    substitution <- input$librarySubstitution %>%
      strsplit(">") %>%
      unlist
    referenceBase <- substitution[1]
    alternateAllele <- substitution[2]
    result <- result %>%
      filter(`Reference base` == referenceBase) %>%
      filter(`Alternate allele` == alternateAllele)

    result <- filter(result, Depth >= input$libraryMinimumDepth)

    snvs <- values$snvs
    result <- result %>%
      left_join(
        mutate(
          snvs,
          Called = Filters != "no call",
          Filtered = ifelse(!Called, NA, !(Filters %in% c("pass", "PASS", "")))
        ),
        by = c("ID", "Amplicon", "Chromosome", "Position", "Reference base", "Alternate allele")
      ) %>%
      mutate(Called = ifelse(is.na(Called), ifelse(nrow(snvs) == 0, NA, FALSE), Called)) %>%
      mutate(Series = ifelse(is.na(Called), "Call status unknown", ifelse(Called, ifelse(Filtered, "Filtered", "Called"), "Not called")))

    result
  })

  librarySubstitutionData <- reactive({
    result <- sampleSubstitutionData()
    library <- selectedLibrary()
    if (is.null(library))
      slice(result, 0)
    else
      filter(result, ID == library$ID)
  })

  libraryFittedDistribution <- reactive({
    data <- librarySubstitutionData()
    alleleFractions <- data$`Allele fraction`
    result <- fitDistribution(
      alleleFractions,
      excludeHighestProportion = input$libraryExcludeHighestProportionForFitting,
      maximumAlleleFraction = input$libraryMaximumAlleleFractionForFitting,
      distribution = input$libraryDistribution,
      thresholdProbabilities = thresholdProbabilities
    )
    result$substitutionData <- data
    result
  })

  libraryMaximumAlleleFraction <- reactive({
    fittedDistribution <- libraryFittedDistribution()
    maximumAlleleFraction <- fittedDistribution$maximumAlleleFraction
    if (input$setLibraryMaximumAlleleFraction) maximumAlleleFraction <- input$libraryMaximumAlleleFraction
    maximumAlleleFraction
  })

  libraryThresholdAlleleFraction <- reactive({
    thresholdAlleleFraction <- -1.0
    fittedDistribution <- libraryFittedDistribution()
    thresholds <- fittedDistribution$thresholds
    if (!is.null(thresholds))
    {
      thresholdProbability <- input$libraryThresholdProbability
      thresholdIndex <- which(thresholdProbabilities == thresholdProbability)
      if (length(thresholdIndex) == 1)
        thresholdAlleleFraction <- thresholds[thresholdIndex]
    }
    thresholdAlleleFraction
  })

  output$librarySubstitutionSummary <- renderUI({

    summary <- "No library/substitution selected"

    library <- selectedLibrary()
    if (!is.null(library))
    {
      summary <- str_c(
        "<b>",
        library$ID,
        " ",
        library$Sample,
        " ",
        input$librarySubstitution,
        "</b>"
      )

      librarySubstitutionData <- librarySubstitutionData()
      if (nrow(librarySubstitutionData) > 0)
        summary <- str_c(summary, str_c(nrow(librarySubstitutionData), " target locations"), sep = "&nbsp;&nbsp;")
    }

    HTML(summary)
  })

  output$libraryScatterBoxPlot <- renderHighchart({

    fittedDistribution <- libraryFittedDistribution()

    data <- fittedDistribution$substitutionData
    if (!is.null(data))
    {
      selectedSubstitution <- values$selectedSubstitution
      if (!is.null(selectedSubstitution))
      {
        data <- data %>%
          mutate(Series = ifelse(
            Chromosome == selectedSubstitution$Chromosome &
              Position == selectedSubstitution$Position,
            # `Reference base` == selectedSubstitution$`Reference base` &
            # `Alternate allele` == selectedSubstitution$`Alternate allele`,
            "Selected overlapping amplicon", Series)) %>%
          mutate(Series = ifelse(
            Chromosome == selectedSubstitution$Chromosome &
              Position == selectedSubstitution$Position &
              # `Reference base` == selectedSubstitution$`Reference base` &
              # `Alternate allele` == selectedSubstitution$`Alternate allele`,
              Amplicon == selectedSubstitution$Amplicon,
            "Selected", Series))
      }

      data <- data %>%
        arrange(Chromosome, Position) %>%
        transmute(
          Chromosome,
          Position,
          Amplicon,
          `Allele fraction`,
          Series,
          CallStatus = ifelse(!is.na(Filters) & !(Filters %in% c("pass", "PASS", "")), str_c("<br>Filters: ", Filters),
                              ifelse(Series == "Selected" | Series == "Selected replicate", ifelse(is.na(Called), "<br>Call status unknown", ifelse(Called, "<br>Called", "<br>Not called")), "")),
          tooltip = str_c(
            Chromosome, " ", Position, " ", input$librarySubstitution,
            "<br>Amplicon: ", Amplicon,
            "<br>AF: ", format(round(`Allele fraction`, digits = 5), scientific = FALSE),
            CallStatus,
            ifelse(Called, str_c("<br>Confidence: ", Confidence), "")
          ),
          ID = str_c(ID, Chromosome, Position, `Reference base`, `Alternate allele`, Amplicon, sep = "\t")
        )
    }

    maximumAlleleFraction <- libraryMaximumAlleleFraction()
    thresholdAlleleFraction <- libraryThresholdAlleleFraction()

    alleleFractionScatterBoxPlot(
      data,
      series = c("Call status unknown", "Not called", "Filtered", "Called", "Selected", "Selected overlapping amplicon"),
      colours = c(hex_to_rgba("#7cb5ec", 0.5), hex_to_rgba("#7cb5ec", 0.5), hex_to_rgba("#8085e9", 0.5), "#90ed7d", "#f15c80", hex_to_rgba("#f15c80", 0.55)),
      xlabel = "Location",
      maximumAlleleFraction = maximumAlleleFraction,
      thresholdAlleleFraction = thresholdAlleleFraction,
      thresholdLabel = str_c("p = ", input$libraryThresholdProbability, ", AF = ", format(round(thresholdAlleleFraction, digits = 4), scientific = FALSE)),
      clicked = "libraryScatterBoxPlot_clicked"
    )
  })

  observe({
    selected <- input$libraryScatterBoxPlot_clicked
    if (is.null(selected)) return(NULL)

    selected <- selected %>%
      strsplit("\t") %>%
      unlist

    id <- selected[1]
    chromosome <- selected[2]
    position <- as.integer(selected[3])
    referenceBase <- selected[4]
    alternateAllele <- selected[5]
    amplicon = selected[6]

    isolate(selectedId <- values$selectedId)

    if (is.null(selectedId) || selectedId != id)
      values$selectedId <- id

    substitution <- tibble(
      Amplicon = amplicon,
      Chromosome = chromosome,
      Position = position,
      `Reference base` = referenceBase,
      `Alternate allele` = alternateAllele
    )

    isolate(selectedSubstitution <- values$selectedSubstition)

    if (is.null(selectedSubstitution) ||
        selectedSubstitution$Chromosome != substitution$Chromosome ||
        selectedSubstitution$Position != substitution$Position ||
        selectedSubstitution$`Reference base` != substitution$`Reference base` ||
        selectedSubstitution$`Alternate allele` != substitution$`Alternate allele` ||
        selectedSubstitution$Amplicon != substitution$Amplicon)
    {
      message("Updating values$selectedSubstitution (3) ", unlist(selectedSubstitution), " -> ", unlist(substitution))
      values$selectedSubstitution <- substitution
    }
  })

  observe({
    input$libraryCheckSelectedButton
    message(isolate(values$selectedId))
    message(isolate(values$selectedSubstitution) %>% unlist(use.names = FALSE))
  })

  observe({
    input$libraryClearSelectedButton
    values$selectedSubstitution <- NULL
  })

  librarySelectedLocationData <- reactive({

    data <- sampleSubstitutionData()

    selectedSubstitution <- values$selectedSubstitution
    if (is.null(selectedSubstitution))
      data <- slice(data, 0)
    else
    {
      data <- data %>%
        filter(Chromosome == selectedSubstitution$Chromosome & Position == selectedSubstitution$Position)
    }

    data <- data %>%
      mutate(Called = ifelse(is.na(Called), "", ifelse(Called, "yes", "no"))) %>%
      mutate(Filtered = ifelse(is.na(Filtered), "", ifelse(Filtered, "yes", "no")))

    data
  })

  output$librarySelectedLocationTable <- DT::renderDataTable({

    data <- librarySelectedLocationData() %>%
      select(ID, Sample, Amplicon, Chromosome, Position, `Reference base`, `Alternate allele`, Called, Filtered, `Allele fraction`, `A count`, `C count`, `G count`, `T count`) %>%
      arrange(ID, Sample, Amplicon)

    columnNames <- c("ID", "Sample", "Amplicon", "Chr", "Pos", "Ref", "Alt", "Called", "Filtered", "AF", "A", "C", "G", "T")

    isolate(snvs <- values$snvs)

    if (nrow(snvs) == 0)
    {
      data <- select(data, -Called, -Filtered)
      columnNames <- setdiff(columnNames, c("Called", "Filtered"))
    }

    datatable(
      data,
      colnames = columnNames,
      rownames = FALSE,
      selection = "single",
      options = list(
        dom = 'tip',
        pageLength = 50,
        searchHighlight = TRUE,
        columnDefs = list(list(className = 'dt-center', targets = c(3, 4, 5)))
      )
    ) %>%
      formatRound(c("Allele fraction"), digits = 5)
  },
  server = FALSE
  )

  observe({
    selectedRow <- input$librarySelectedLocationTable_rows_selected
    if (is.null(selectedRow)) return(NULL)

    data <- librarySelectedLocationData()
    if (nrow(data) == 0) return(NULL)

    selected <- slice(data, selectedRow)
    if (nrow(selected) != 1) return(NULL)

    id <- selected$ID
    isolate(selectedId <- values$selectedId)
    if (is.null(selectedId) || selectedId != id) values$selectedId <- id

    substitution <- selected %>%
      select(Chromosome, Position, `Reference base`, `Alternate allele`, Amplicon)
    isolate(selectedSubstitution <- values$selectedSubstitution)
    if (is.null(selectedSubstitution) ||
        selectedSubstitution$Chromosome != substitution$Chromosome ||
        selectedSubstitution$Position != substitution$Position ||
        selectedSubstitution$`Reference base` != substitution$`Reference base` ||
        selectedSubstitution$`Alternate allele` != substitution$`Alternate allele` ||
        selectedSubstitution$Amplicon != substitution$Amplicon)
    {
      message("Updating values$selectedSubstitution (4) ", unlist(selectedSubstitution), " -> ", unlist(substitution))
      values$selectedSubstitution <- substitution
    }
  })

  output$libraryDensityPlot <- renderHighchart({

    data <- NULL
    series <- NULL

    fittedDistribution <- libraryFittedDistribution()

    density <- fittedDistribution$density
    if (!is.null(density))
    {
      densitySeries <- str_c("Unfiltered (", length(fittedDistribution$alleleFractions), " locations)")
      data <- mutate(density, Series = densitySeries)
      series <- densitySeries

      filteredDensity <- fittedDistribution$filteredDensity
      if (!is.null(filteredDensity))
      {
        filteredDensitySeries <- densitySeries <- str_c("Filtered (", length(fittedDistribution$filteredAlleleFractions), " locations)")
        data <- data %>%
          bind_rows(filteredDensity %>%
          mutate(Series = filteredDensitySeries))
        series <- c(series, filteredDensitySeries)
      }

      fitted <- fittedDistribution$fitted
      if (!is.null(fitted))
      {
        fittedDistributionSeries <- "Fitted distribution"
        data <- data %>%
          bind_rows(fitted %>%
          mutate(Series = fittedDistributionSeries))
        series <- c(series, fittedDistributionSeries)
      }
    }

    maximumAlleleFraction <- libraryMaximumAlleleFraction()
    thresholdAlleleFraction <- libraryThresholdAlleleFraction()

    alleleFractionDensityPlot(
      data,
      series = series,
      type = c("areaspline", "areaspline", "spline"),
      colours = c("#8085e9", "#7cb5ec", "#434348"),
      maximumAlleleFraction = maximumAlleleFraction,
      thresholdAlleleFraction = thresholdAlleleFraction,
      thresholdLabel = str_c("p = ", input$libraryThresholdProbability, ", AF = ", format(round(thresholdAlleleFraction, digits = 4), scientific = FALSE)),
      decimalPlaces = 4
    )
  })

  output$libraryCullenFreyGraph <- renderPlot({

    data <- librarySubstitutionData()

    alleleFractions <- data$`Allele fraction`
    nonZeroAlleleFractions <- alleleFractions[alleleFractions != 0]

    if (length(nonZeroAlleleFractions) > 0)
      descdist(nonZeroAlleleFractions)
    else
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  })

}

