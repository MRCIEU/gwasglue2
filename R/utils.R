

# "Not in" function
 "%ni%" <- Negate("%in%")


# @title Utility function to display warning messages as they occur (from susieR)
# @param ... warning message
# @param style either "warning" or "hint"
#'@importFrom crayon combine_styles
warning_message = function(..., style=c("warning", "hint")) {
  style = match.arg(style)
  if (style=="warning" && getOption("warn")>=0) {
    alert <- combine_styles("bold", "underline", "red")	
    message(alert("WARNING:"), " ", ...)
  } else {
    alert <- combine_styles("bold", "underline", "magenta")	
    message(alert("HINT:"), " ", ...)
  }
}
