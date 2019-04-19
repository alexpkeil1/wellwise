#' Well water measurements
#'
#' A dataset containing well water measurements and simulated
#' outcomes for 52440 units across 20 iterations
#'
#' @format A data frame with 52440 rows and 29 variables:
#' \describe{
#'    \item{iter}{iteration number}
#'    \item{y}{outcome}
#'    \item{Acidity}{Acidity measurement}
#'    \item{Acidity_lod}{Acidity at LOD (1/0/NA)}
#'    \item{Aluminum}{Aluminum measurement}
#'    \item{Aluminum_lod}{Aluminum at LOD (1/0/NA)}
#'    \item{Ammonia}{Ammonia measurement}
#'    \item{Ammonia_lod}{Ammonia at LOD (1/0/NA)}
#'    \item{Antimony}{Antimony measurement}
#'    \item{Antimony_lod}{Antimony at LOD (1/0/NA)}
#'    \item{Arsenic}{Arsenic measurement}
#'    \item{Arsenic_lod}{Arsenic at LOD (1/0/NA)}
#'    \item{Barium}{Barium measurement}
#'    \item{Barium_lod}{Barium at LOD (1/0/NA)}
#'    \item{Beryllium}{Beryllium measurement}
#'    \item{Beryllium_lod}{Beryllium at LOD (1/0/NA)}
#'    \item{Boron}{Boron measurement}
#'    \item{Boron_lod}{Boron at LOD (1/0/NA)}
#'    \item{Cadmium}{Cadmium measurement}
#'    \item{Cadmium_lod}{Cadmium at LOD (1/0/NA)}
#'    \item{Calcium}{Calcium measurement}
#'    \item{Calcium_lod}{Calcium at LOD (1/0/NA)}
#'    \item{Chloride}{Chloride measurement}
#'    \item{Chloride_lod}{Chloride at LOD (1/0/NA)}
#'    \item{Chromium}{Chromium measurement}
#'    \item{Chromium_lod}{Chromium at LOD (1/0/NA)}
#'    \item{Cobalt}{Cobalt measurement}
#'    \item{Cobalt_lod}{Cobalt at LOD (1/0/NA)}
#'    \item{Color}{Color measurement}
#'    \item{Color_lod}{Color at LOD (1/0/NA)}
#'    \item{Conductivity}{Conductivity measurement}
#'    \item{Conductivity_lod}{Conductivity at LOD (1/0/NA)}
#'    \item{Copper}{Copper measurement}
#'    \item{Copper_lod}{Copper at LOD (1/0/NA)}
#'    \item{Cyanide}{Cyanide measurement}
#'    \item{Cyanide_lod}{Cyanide at LOD (1/0/NA)}
#'    \item{Fluoride}{Fluoride measurement}
#'    \item{Fluoride_lod}{Fluoride at LOD (1/0/NA)}
#'    \item{Gold}{Gold measurement}
#'    \item{Gold_lod}{Gold at LOD (1/0/NA)}
#'    \item{Hexavalent_Chromium}{Hexavalent_Chromium measurement}
#'    \item{Hexavalent_Chromium_lod}{Hexavalent_Chromium at LOD (1/0/NA)}
#'    \item{Insoluble_Iron}{Insoluble_Iron measurement}
#'    \item{Insoluble_Iron_lod}{Insoluble_Iron at LOD (1/0/NA)}
#'    \item{Insoluble_Manganese}{Insoluble_Manganese measurement}
#'    \item{Insoluble_Manganese_lod}{Insoluble_Manganese at LOD (1/0/NA)}
#'    \item{Iron}{Iron measurement}
#'    \item{Iron_lod}{Iron at LOD (1/0/NA)}
#'    \item{Lead}{Lead measurement}
#'    \item{Lead_lod}{Lead at LOD (1/0/NA)}
#'    \item{Lithium}{Lithium measurement}
#'    \item{Lithium_lod}{Lithium at LOD (1/0/NA)}
#'    \item{MSC}{MSC measurement}
#'    \item{MSC_lod}{MSC at LOD (1/0/NA)}
#'    \item{Magnesium}{Magnesium measurement}
#'    \item{Magnesium_lod}{Magnesium at LOD (1/0/NA)}
#'    \item{Manganese}{Manganese measurement}
#'    \item{Manganese_lod}{Manganese at LOD (1/0/NA)}
#'    \item{Mercury}{Mercury measurement}
#'    \item{Mercury_lod}{Mercury at LOD (1/0/NA)}
#'    \item{Molybdenum}{Molybdenum measurement}
#'    \item{Molybdenum_lod}{Molybdenum at LOD (1/0/NA)}
#'    \item{Nickel}{Nickel measurement}
#'    \item{Nickel_lod}{Nickel at LOD (1/0/NA)}
#'    \item{Nitrate}{Nitrate measurement}
#'    \item{Nitrate_lod}{Nitrate at LOD (1/0/NA)}
#'    \item{Nitrite}{Nitrite measurement}
#'    \item{Nitrite_lod}{Nitrite at LOD (1/0/NA)}
#'    \item{Orthophosphate}{Orthophosphate measurement}
#'    \item{Orthophosphate_lod}{Orthophosphate at LOD (1/0/NA)}
#'    \item{Potassium}{Potassium measurement}
#'    \item{Potassium_lod}{Potassium at LOD (1/0/NA)}
#'    \item{Selenium}{Selenium measurement}
#'    \item{Selenium_lod}{Selenium at LOD (1/0/NA)}
#'    \item{Settleable_Solids}{Settleable_Solids measurement}
#'    \item{Settleable_Solids_lod}{Settleable_Solids at LOD (1/0/NA)}
#'    \item{Silica}{Silica measurement}
#'    \item{Silica_lod}{Silica at LOD (1/0/NA)}
#'    \item{Silver}{Silver measurement}
#'    \item{Silver_lod}{Silver at LOD (1/0/NA)}
#'    \item{Sodium}{Sodium measurement}
#'    \item{Sodium_lod}{Sodium at LOD (1/0/NA)}
#'    \item{Soluble_Iron}{Soluble_Iron measurement}
#'    \item{Soluble_Iron_lod}{Soluble_Iron at LOD (1/0/NA)}
#'    \item{Soluble_Manganese}{Soluble_Manganese measurement}
#'    \item{Soluble_Manganese_lod}{Soluble_Manganese at LOD (1/0/NA)}
#'    \item{Strontium}{Strontium measurement}
#'    \item{Strontium_lod}{Strontium at LOD (1/0/NA)}
#'    \item{Sulfate}{Sulfate measurement}
#'    \item{Sulfate_lod}{Sulfate at LOD (1/0/NA)}
#'    \item{Thallium}{Thallium measurement}
#'    \item{Thallium_lod}{Thallium at LOD (1/0/NA)}
#'    \item{Tin}{Tin measurement}
#'    \item{Tin_lod}{Tin at LOD (1/0/NA)}
#'    \item{Total_Alkalinity}{Total_Alkalinity measurement}
#'    \item{Total_Alkalinity_lod}{Total_Alkalinity at LOD (1/0/NA)}
#'    \item{Total_Dissolved_Solids}{Total_Dissolved_Solids measurement}
#'    \item{Total_Dissolved_Solids_lod}{Total_Dissolved_Solids at LOD (1/0/NA)}
#'    \item{Total_Hardness}{Total_Hardness measurement}
#'    \item{Total_Hardness_lod}{Total_Hardness at LOD (1/0/NA)}
#'    \item{Total_Phosphate}{Total_Phosphate measurement}
#'    \item{Total_Phosphate_lod}{Total_Phosphate at LOD (1/0/NA)}
#'    \item{Total_Suspended_Solids}{Total_Suspended_Solids measurement}
#'    \item{Total_Suspended_Solids_lod}{Total_Suspended_Solids at LOD (1/0/NA)}
#'    \item{Turbidity}{Turbidity measurement}
#'    \item{Turbidity_lod}{Turbidity at LOD (1/0/NA)}
#'    \item{Uranium}{Uranium measurement}
#'    \item{Uranium_lod}{Uranium at LOD (1/0/NA)}
#'    \item{Vanadium}{Vanadium measurement}
#'    \item{Vanadium_lod}{Vanadium at LOD (1/0/NA)}
#'    \item{Zinc}{Zinc measurement}
#'    \item{Zinc_lod}{Zinc at LOD (1/0/NA)}
#'    \item{pH}{pH measurement}
#'    \item{pH_lod}{pH at LOD (1/0/NA)}
#'  }
"welldata"

