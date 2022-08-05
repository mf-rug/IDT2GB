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
  
  processed <- eventReactive(c(input$seqs, input$remr), {
    if (input$seqs != '') {
      df <- read.table(text = input$seqs, sep = '\t')
      df_mods <- df[,1:2]
      df_mods$V2 <- str_remove(df_mods$V2, '^/')
      df_mods$V2 <- str_replace_all(df_mods$V2, '//', '/')
      df_mods <- suppressWarnings(separate(df_mods, 'V2', into = as.character(1:max(str_count(df_mods$V2, '/'))), sep = '/'))
      df_mods[is.na(df_mods)] <- ''
      
      seqs <- str_replace_all(df$V2, '/[^/]+/', '')
      print(paste('seqs found to be', paste(seqs, collapse = ',')))
      
      df_mods[apply(df_mods, 2, function(x) str_detect(x, '^[+*mrAGTCUN][+*mrAGTCUN]*$'))] <- ''
      df_mods <- apply(df_mods, 1, function(x) c(x[x != ""], rep("", length(x[x == ""])))) %>% t() %>% as.data.frame()
      if (input$remr) {
        seqs <- str_remove_all(seqs, 'r')
      }
      if (input$remp) {
        seqs <- str_remove_all(seqs, '[+]')
      }
      if (input$rems) {
        seqs <- str_remove_all(seqs, '[*]')
      }
      if (input$remm) {
        seqs <- str_remove_all(seqs, 'm')
      }
      df_mods$seqs <- seqs
      
      # remove empty columns
      df_mods <- df_mods[,as.numeric(which(apply(df_mods,2, function(x) !all(str_detect(x, '^$')))))]
      df_out <- df_mods[,c('V1', 'seqs')]
      df_out$length <- str_count(df_out$seqs)
      seqs_single_char_as_mod <- gsub('/[^/]+/', '€', df$V2)
      if (any(str_detect(seqs_single_char_as_mod, '€'))) {
        print('some mods found')
        
        #this calculates position of the mod by compensating for the extra char that is the modification
        #the first stays, but all others get subtracted by the number of mods that follow
        pos <- lapply(str_locate_all(seqs_single_char_as_mod, '€'), function(x) as.numeric(x[,1]) - c(0, seq_len(max(0, length(x[,1]) -1))))
        
        #this corrects for the fact there is more positions for modifications than bases
        #usually a mod is assigned to the next base, except the 3' end, which gets assigned to the last base
        pos <- lapply(1:length(pos), function(x) if (!identical(pos[[x]], numeric(0))) {
          c(pos[[x]][which(pos[[x]] != max(pos[[x]]) | pos[[x]] == 1)],
            pos[[x]][which(pos[[x]] == max(pos[[x]]) & pos[[x]] != 1)] -1)
        } else { NULL })
        
        df_rest <- df_mods[which(!colnames(df_mods) %in% c('V1', 'seqs'))]
        for (i in 1:ncol(df_rest)) {
          df_out[,paste0('mod_', i)] <- df_rest[, i]
          lookup_mod <- mods[match(df_rest[,i],mods$code), 'name']
          df_out[df_rest[,i] != "", paste0('mod_', i, '_name')] <- lookup_mod[!is.na(lookup_mod)]
          lookup_pos <- sapply(pos, function(x) x[i]) %>% unlist()
          df_out[df_rest[,i] != "", paste0('mod_', i, '_pos')] <- lookup_pos[!is.na(lookup_pos)]
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
          repl <- lapply(1:max_mods, function(x) paste('     Modified_base   ', 
                                                       df[i, paste0('mod_', x, '_pos')], 
                                                       '\n                     /label="', 
                                                       df[i, paste0('mod_', x)],
                                                       '"\n                     /Modification="', 
                                                       df[i, paste0('mod_', x, '_name')])) %>% unlist() %>% paste0(collapse = '\n')
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