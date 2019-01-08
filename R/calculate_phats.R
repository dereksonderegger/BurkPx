#' Calculate phat values for an input file of assay values
#'
#' @param data_in The input data set. This could be either a data frame or a character
#'                string pointing to the name of an input file. If it is a file name, then
#'                it should either be the file name in
#'                the current working directory or be the full path to the file.
#'                This file may either be a .csv, .xls, or .xlsx file. See below
#'                for the format of the input file or input data frame.
#' @param file_out The name of the output .csv file that contains all of the phat
#'                 values. If it isn't given, the default is to return the output
#'                 as a data frame.
#' @param Species The species we want to predict for.  Accepted values are
#'                "Human" or "NHP". This defaults to "Human".
#' @param Trained_on The data set which the model was trained on.  Typically users
#'                   will want the model trained on the full dataset, but for cross-validation
#'                   purposes, we might want to specify the model built just using the
#'                   designated training data. Valid options are 'Full' and 'Training' and
#'                   defaults to the full data set.
#' @param Isotype Which isotypes should we report. Valid options are any combination
#'                of 'IgG', 'IgM' or 'IgGM'. By default, the procedure produces all of them.
#'                If the necessary data is not presented, the results are NAs or NaNs.
#' @param Method What model selection method was used. Valid options are 'LASSO', 'Ridge' or
#'               'HCP1' or any combination of these.
#' @param average_replicates If a serum sample has replicate observations, should we
#'                average the replicate values and produce a single phat or should we
#'                treat them seperatately and produce phat values for each replicate. The
#'                default (FALSE) is to not average and produce phat values for each replicate.
#' @param simplify_model_name The column names for the resulting phat values are the model names
#'                of the model that produced it. However, because the model names is often quite
#'                long (e.g. Human_Training_IgGM_LASSO) it is often desirable to reduce the model
#'                name to just the aspects that change so that the resulting column names
#'                are simpler. The defaults to TRUE.  If there is only one model requested,
#'                the simplified column name is just "phat".
#'
#' @details The input file should be a spreadsheet file (either .csv, .xls, .xlsx file) are data
#'          frame with the same column requirements.  The input should have a column
#'          'Serum' which denotes the Serum ID. If the Serum ID values start with IgG or IgM, then we
#'          will use those as the Isotype.  If the Serum ID values do not start with IgG or IgM, then there
#'          must be an 'Isotype' column that contains that information. The remaining columns are the
#'          antibody values and should include some subset of the column names.
#'
#'          The following antibodies can be used:
#'          BPSL2096_AhpC, BPSL1404_ClpX, BPSS0476_GroS, BPSS0477_GroEL2, BPSS0135,
#'          BPSL1743_Arg, BPSL2827_DNAK, BPSL3222_rpIL, MSHR5855.WCL, BPSL1201_IMPS,
#'          BPSL3396_AtpD, BPSS0530, BPSL2522_OmpA, BPSS1850, LPSA, LPSB, CPS,
#'          BPSS1769_NADH, BPSS1652, BPSL2697_GroEL, BPSS1498_HCP1.B
#'
#'          If a particular model does not use an antibody (e.g. the HCP1 models only use HCP1
#'          values), then the input data could be missing all of the other columns and the
#'          function will still work.
#'
#' @export
calculate_phats <- function(data_in,
            file_out=NULL,
            Species = 'Human',
            Trained_on = 'Full',
            Isotype = c('IgG', 'IgM', 'IgGM'),
            Method = 'LASSO',
            average_replicates = FALSE,
            simplify_model_name = TRUE){

  if( is.character(data_in) ){
    # load the data
    if( str_detect(data_in, fixed('.csv')) ){
      data <- read.csv(data_in)
    }else{
      data <- readxl::read_excel(data_in)
    }
  }else{
    data <- data_in
  }


  # The first column should be named 'Serum', but if not,
  # warn the user and change it.
  if( colnames(data)[1] != 'Serum' ){
    colnames(data)[1] <- 'Serum'
    message( 'First column in the data should be named Serum!')
  }

  # If the Isotype is found via the first characters of the Serum
  if( str_sub( data$Serum[1], 1,3) %in% c('IgG', 'IgM') ){
    data <- data %>%
      mutate(Isotype = str_sub(Serum, 1, 3)) %>%
      mutate( Serum = str_sub(Serum, start=5) )
    message( 'Pulling Isotype from first characters of Serum' )
  }

  # If there is no Isotype column, stop
  if( !('Isotype' %in% colnames(data)) ){
    stop('There is no Isotype column!')
  }

  # If there are replicates, add a column for that for that
  data <- data %>%
    group_by(Serum, Isotype) %>%
    mutate(Rep = 1:n() )  %>%
    select(Serum, Isotype, Rep, everything()) %>%
    group_by()

  # average the replicate values, if desired.  This keeps the Rep column
  # even though it is now useless because the rest of the script assumes it is
  # present
  if( average_replicates ){
    data <- data %>%
      group_by(Serum, Isotype ) %>%
      summarize_all( mean )
  }

  # Build the name vector of the set of models the user wants
  models_to_use <- expand.grid(Species=Species, Trained_on=Trained_on,
                               Isotype=Isotype, Method=Method ) %>%
    mutate(model_name = str_c(Species, Trained_on, Isotype, Method, sep='_')) %>%
    mutate(simple_model_name = '')

  if( length(unique(models_to_use$Species)) > 1 )
    models_to_use <- models_to_use %>% mutate( simple_model_name = str_c( simple_model_name, Species, '_' ) )
  if( length(unique(models_to_use$Trained_on)) > 1 )
    models_to_use <- models_to_use %>% mutate( simple_model_name = str_c( simple_model_name, Trained_on, '_' ) )
  if( length(unique(models_to_use$Isotype)) > 1 )
    models_to_use <- models_to_use %>% mutate( simple_model_name = str_c( simple_model_name, Isotype, '_' ) )
  if( length(unique(models_to_use$Method)) > 1 )
    models_to_use <- models_to_use %>% mutate( simple_model_name = str_c( simple_model_name, Method, '_' ) )

  models_to_use <- models_to_use %>%
    mutate( simple_model_name = str_sub(simple_model_name,1, -2) )
  if( nrow(models_to_use) == 1 ){
    models_to_use$simple_model_name = 'p.hat'
  }

  # Split the data into IgG and IgM and make the IgGM combination
  IgG <- data %>% filter(Isotype == 'IgG')
  IgM <- data %>% filter(Isotype == 'IgM')

  if( nrow(IgG) > 0  & nrow(IgM) > 0 ){
    IgGM <- data %>% group_by(Serum, Isotype, Rep) %>%
      gather('Antigen','Value', -Serum, -Isotype, -Rep) %>%
      tidyr::unite('Antigen', Isotype, Antigen) %>% ungroup() %>%
      spread(Antigen, Value) %>% arrange(Serum, Rep) %>%
      group_by(Serum, Rep) %>% summarise_all(max, na.rm=TRUE)
  }else{
    IgGM <- NULL
  }


  # Now apply all the models
  out <- NULL
  for( i in 1:nrow( models_to_use) ){
    df <- get(as.character(models_to_use[i, 'Isotype'] ) )  # get the data set to predict on
    if( !is.null(df) && nrow(df) > 0 ){      # if we actually have the data
      model <- models[[ models_to_use[i, 'model_name'] ]]
      out <-
        df %>%
        select(Serum, Rep) %>%
        mutate(model_name = models_to_use[i, 'model_name']) %>% ungroup() %>%
        mutate( phat = predict( model, newdata=df, type='response' )) %>%
        rbind( out, . )
    }
  }


  if( simplify_model_name ){
    out <- left_join( out, models_to_use[, c('model_name','simple_model_name')], by='model_name' ) %>%
      mutate(model_name = simple_model_name) %>% select(-simple_model_name)
  }

  out <- out %>%
    group_by(Serum, Rep) %>%
    spread(model_name, phat) %>%
    arrange(Serum, Rep)

  if( is.null(file_out) ){
    return(out)
  }else{
    write.csv(out, file=file_out)
    return(invisible(out))
  }

}

