# Rana_kode.R - Kodet i UTF-8 for C% stC8tte norske tegn (C%, C8, C&)

# Installer og last inn nC8dvendige pakker hvis de ikke allerede er tilgjengelige
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("fitdistrplus", quietly = TRUE)) install.packages("fitdistrplus")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("knitr", quietly = TRUE)) install.packages("knitr")  # For tabellformatering

# Last inn verktC8y for C% lese Excel-filer, manipulere data, passe til fordelinger og tegne grafer
library(readxl)
library(dplyr)
library(fitdistrplus)
library(ggplot2)
library(knitr)  # For C% lage pene tabeller

# Input-verdier samlet C8verst for enkel endring (alle valg styres her)
input_parameters <- list(
  file_path = "C:\\Users\\Magnu\\OneDrive\\Dokumenter\\NTNU-Emner\\Master\\Rana_data.xlsx",  # Sti til Excel-filen med data
  areas_of_interest = c("U-stollen, N4.5", "U-stollen, N5", "U-stollen, N5.5", "U-stollen, N6", "U-stollen, N6.5"),  # OmrC%der som skal analyseres
  lithologies_of_interest = c("GlimSkif", "AmfSkif", "GranatGlimSkif", "KarbGlimSkif"),  # Bergarter som skal inkluderes
  combined_lithologies = c("GlimSkif", "AmfSkif", "GranatGlimSkif", "KarbGlimSkif"),  # Bergarter som skal kombineres pC% slutten
  max_depth = 20,  # Maksimal dybde (i meter) for C% begrense analysen
  include_faultzone = FALSE,  # Velg om "FaultZone" skal inkluderes (True/False)
  num_simulations = 10000  # Antall simuleringer for Monte Carlo-analyse
)

# Skriv ut input-verdiene i konsollen for C% sjekke at de er korrekte
cat("Input-parametere:\n")
cat("Filsti:", input_parameters$file_path, "\n")
cat("OmrC%der:", paste(input_parameters$areas_of_interest, collapse = ", "), "\n")
cat("Litologier:", paste(input_parameters$lithologies_of_interest, collapse = ", "), "\n")
cat("Kombinerte litologier:", paste(input_parameters$combined_lithologies, collapse = ", "), "\n")
cat("Maks dybde:", input_parameters$max_depth, "meter\n")
cat("Inkluder FaultZone:", input_parameters$include_faultzone, "\n")
cat("Antall simuleringer:", input_parameters$num_simulations, "\n\n")

# Last inn data fra spesifikke ark i Excel-filen og gjC8r dem tilgjengelige i programmet
required_sheets <- c("Vestmalmen_Headers", "Vestmalmen_Lithology", "Vestmalmen_Rock_Mechanics_old")  # Navn pC% arkene som skal leses
data_list <- setNames(lapply(required_sheets, function(sheet) read_excel(input_parameters$file_path, sheet = sheet)), required_sheets)  # Les arkene og gi dem navn
list2env(data_list, envir = .GlobalEnv)  # GjC8r dataene til globale variabler som kan brukes overalt

# Sjekk at alle nC8dvendige kolonner finnes i dataene, stopp hvis noe mangler
check_columns <- function(df, required_cols) {  # Funksjon for C% sjekke kolonner
  missing_cols <- setdiff(required_cols, names(df))  # Finn kolonner som mangler
  if (length(missing_cols) > 0) stop("Manglende kolonner: ", paste(missing_cols, collapse = ", "))  # Gi feilmelding hvis noe mangler
}

check_columns(Vestmalmen_Headers, c("Hole_number", "Area"))  # Sjekk kolonner i "Headers"-arket
check_columns(Vestmalmen_Lithology, c("Hole_number", "From", "To", "Lithology", "Lithology Sub-class"))  # Sjekk kolonner i "Lithology"-arket
check_columns(Vestmalmen_Rock_Mechanics_old, c("Hole_number", "From", "To", "RQD_(calc)", "Jr", "Ja", "Q-Value_(calc)"))  # Sjekk kolonner i "Rock Mechanics"-arket

# Filtrer ut borehull som ikke er i de valgte omrC%dene (uavhengig av store/smC% bokstaver)
valid_holes <- Vestmalmen_Headers$Hole_number[tolower(Vestmalmen_Headers$Area) %in% tolower(input_parameters$areas_of_interest)]  # Finn gyldige borehull

# Filtrer litologidata basert pC% borehull og dybde, juster sluttpunkt til maks dybde
lithology_filtered <- Vestmalmen_Lithology %>%
  filter(Hole_number %in% valid_holes, (From < input_parameters$max_depth | To <= input_parameters$max_depth)) %>%  # Behold data innen dybde
  mutate(To = pmin(To, input_parameters$max_depth))  # Sett sluttpunkt til maks dybde hvis det er hC8yere

# Filtrer mekanikkdata basert pC% borehull og dybde, juster sluttpunkt til maks dybde
rock_mechanics_filtered <- Vestmalmen_Rock_Mechanics_old %>%
  filter(Hole_number %in% valid_holes, (From < input_parameters$max_depth | To <= input_parameters$max_depth)) %>%  # Behold data innen dybde
  mutate(To = pmin(To, input_parameters$max_depth))  # Sett sluttpunkt til maks dybde hvis det er hC8yere

# Lag en tekstfil for C% lagre resultater og skriv ut overskrift
output_file <- file("Rana_output.txt", "w")  # Opprett en fil for C% skrive resultater
cat("Verdier for RQD, Jr, Ja, og Q per meterintervall (alle intervaller, uansett NA-verdier)\n", file = output_file)  # Skriv overskrift

# GC% gjennom hvert borehull og skriv ut data
for (hole in unique(rock_mechanics_filtered$Hole_number)) {  # For hvert unike borehull
  cat("\nBorehull:", hole, "\n", file = output_file)  # Skriv borehullnummer
  hole_mech <- rock_mechanics_filtered %>% filter(Hole_number == hole)  # Hent mekanikkdata for dette borehullet
  hole_lith <- lithology_filtered %>% filter(Hole_number == hole)  # Hent litologidata for dette borehullet
  
  # Bruk alle intervaller uansett om det er manglende verdier
  hole_mech_all <- hole_mech  # Behold alle data for borehullet
  
  # Bestem hvor mange intervaller som skal vises (alle i dette tilfellet)
  n_to_show <- nrow(hole_mech_all)  # Antall rader/intervaller
  
  if (n_to_show == 0) {  # Hvis ingen data
    cat("Ingen data for dette borehullet b