server <- function(input, output) {
  
  observeEvent(input$seqs,{
    if (input$seqs == '') {
      hide('clip')
      hide('downloadgb')
    } else {
      show('clip')
      show('downloadgb')
    }
  })
  
  mods <- suppressWarnings(read_csv('all_IDT_modifications.csv', col_names = TRUE, show_col_types = FALSE) %>% as.data.frame())
  
  processed <- eventReactive(input$seqs, {
    if (input$seqs != '') {
      df <- read.table(text = input$seqs, sep = '\t')
      df_mods <- df[,1:2]
      df_mods$V2 <- str_remove(df_mods$V2, '^/')
      df_mods <- suppressWarnings(separate(df_mods, 'V2', into = as.character(1:max(str_count(df_mods$V2, '/'))), sep = '/'))
      df_mods[is.na(df_mods)] <- ''
      
      seqs <- str_replace_all(df$V2, '/[^/]+/', '')
      
      df_mods[apply(df_mods, 2, function(x) str_detect(x, '^[AGTCUN][AGTCUN]*$'))] <- ''
      df_mods$seqs <- seqs
      # remove empty columns
      df_mods <- df_mods[,as.numeric(which(apply(df_mods,2, function(x) !all(str_detect(x, '^$')))))]
      df_out <- df_mods[,c('V1', 'seqs')]
      df_out$length <- str_count(df_out$seqs)
      seqs_single_char_as_mod <- gsub('/[^/]+/', '€', df$V2)
      if (any(str_detect(seqs_single_char_as_mod, '€'))) {
        pos <- lapply(str_locate_all(seqs_single_char_as_mod, '€'), function(x) as.numeric(x[,1]) - c(0, 1:(length(x[,1]) -1)))
        pos <- lapply(1:length(pos), function(x) c(pos[[x]][which(pos[[x]] != max(pos[[x]]))],pos[[x]][which(pos[[x]] == max(pos[[x]]))] -1))
  
        df_rest <- df_mods[str_detect(colnames(df_mods), '^[0-9]+$')]
        for (i in 1:ncol(df_rest)) {
          df_out[,paste0('mod_', i)] <- df_rest[, i]
          df_out[df_rest[,i] != "", paste0('mod_', i, '_name')] <- mods[mods$code %in% df_rest[,i], 'name']
          df_out[, paste0('mod_', i, '_pos')] <- sapply(pos, function(x) x[i])
        }
      }
      colnames(df_out)[1:2] <- c("name", "sequence")
      df_out
    }
  })
  
  out_text <- eventReactive(processed(), {
    if (input$seqs != '') {
      df <- processed()
      
      bb_text <- paste0(
        'LOCUS       ', df$name, '              ', df$length, ' bp    DNA     linear   ', toupper(format(Sys.Date(), format="%d-%b-%Y")), '
FEATURES             Location/Qualifiers
     primer_bind     1..', df$length, '
                     /Sequence="', df$sequence, '"
                     /label="Binding Region"
     PLACE_HOLDER_FOR_MODS
ORIGIN      
        1 ', df$sequence, '
//')
      
      df$out <- paste0(
'LOCUS       ', df$name, '              ', df$length, ' bp    DNA     linear   ', toupper(format(Sys.Date(), format="%d-%b-%Y")), '
FEATURES             Location/Qualifiers
     primer_bind     1..', df$length, '
                     /Sequence="', df$sequence, '"
                     /label="Binding Region"
PLACE_HOLDER_FOR_MODS
ORIGIN      
        1 ', df$sequence, '
//')
      for (i in 1:nrow(df)) {
        max_mods <- sum(!is.na(df[i, str_detect(colnames(df), '_pos$')]))
        if (max_mods > 0) {
          repl <- lapply(1:max_mods, function(x) paste('     misc_feature   ', df[i, paste0('mod_', x, '_pos')], '\n                     /label="', df[i, paste0('mod_', x)] ,'"')) %>% unlist() %>% paste0(collapse = '\n')
        } else {
          repl <- ''
        }
        df[i, 'out'] <- str_replace(df[i, 'out'], 'PLACE_HOLDER_FOR_MODS', repl)
      }
      df$out
    }
  })
  
  output$outtext <- renderUI({
    out <- out_text()
    out <- str_replace_all(out, '\n', '<br/>')
    out <- str_replace_all(out, '\t', '&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp')
    out <- str_replace_all(out, ' ', '&nbsp')
    HTML(paste0('<pre><tt>', paste(out, collapse = '<br/>'), '</tt></pre>'))
  })
  
  output$outtable <- renderDataTable({
    df <- processed()
  })
  
  output$downloadgb <- downloadHandler(
    filename = function() {
      paste("IDT2GB-", Sys.Date(), ".gb", sep="")
    },
    content = function(file) {
      write_file(paste(out_text(), collapse = '\n'), file)
  })  
  
  output$clip <- renderUI({
    rclipButton("clipbtn", HTML('<small><font color="grey">Copy .gb format to clipboard</font></small>'), paste(out_text(), collapse = '\n'), modal = FALSE, icon("copy"), )
  })
}