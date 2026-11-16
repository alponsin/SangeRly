# Lists the texts of the program that are not yet in www/i18n_id.json, so that they can be translated.
# Run it from the project folder:  source("list_untranslated.R")   ->  writes untranslated_texts.txt
dict <- names(jsonlite::fromJSON("www/i18n_id.json")$dict)
files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
strs <- character()
for (f in files) {
  pd <- tryCatch(utils::getParseData(parse(f, keep.source = TRUE)), error = function(e) NULL)
  if (is.null(pd)) next
  x <- pd$text[pd$token == "STR_CONST"]
  x <- substr(x, 2, nchar(x) - 1)
  strs <- c(strs, x)
}
strs <- unique(trimws(gsub("\\\\n", " ", gsub('\\\\"', '"', strs))))
ui <- strs[nchar(strs) >= 4 & grepl("[A-Za-z]{3}", strs) & !grepl("[{};<>$]|::|^\\.|\\.(R|csv|png|svg|pdf|xlsx|js|css)$", strs) &
           (grepl(" ", strs) | grepl("^[A-Z]", strs))]
todo <- setdiff(ui, dict)
writeLines(sort(todo), "untranslated_texts.txt")
cat(length(todo), "texts are not translated yet (of", length(ui), "). They are in untranslated_texts.txt\n")
